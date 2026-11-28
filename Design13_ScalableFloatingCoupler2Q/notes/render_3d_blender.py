# -*- coding: utf-8 -*-
"""Build and render a 3D photorealistic Blender 5.3 scene from Design13 r4 GDS polygons.

Run via:
    blender --background --python Design13_ScalableFloatingCoupler2Q/notes/render_3d_blender.py
"""
import json
import math
from pathlib import Path

import bmesh
import bpy
from mathutils import Vector

HERE = Path(__file__).resolve().parent
POLY_JSON = HERE / 'fig' / 'design13_r4_3d_polys.json'
BLEND_OUT = HERE / 'fig' / 'Design13_r4_3d.blend'
PNG_CENTER = HERE / 'fig' / 'render_3d_center.png'
PNG_MACRO = HERE / 'fig' / 'render_3d_macro_bridge.png'
PNG_JJ = HERE / 'fig' / 'render_3d_manhattan_jj.png'
PNG_HERO = HERE / 'fig' / 'render_3d_hero_all.png'


def clear_scene():
    bpy.ops.wm.read_factory_settings(use_empty=True)


def make_material(name, base_color, metallic=0.9, roughness=0.2, coat=0.0, emission=None, emission_strength=0.0):
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get('Principled BSDF')
    bsdf.inputs['Base Color'].default_value = (*base_color, 1.0)
    bsdf.inputs['Metallic'].default_value = metallic
    bsdf.inputs['Roughness'].default_value = roughness
    if 'Coat Weight' in bsdf.inputs:
        bsdf.inputs['Coat Weight'].default_value = coat
    if emission is not None and 'Emission Color' in bsdf.inputs:
        bsdf.inputs['Emission Color'].default_value = (*emission, 1.0)
        bsdf.inputs['Emission Strength'].default_value = emission_strength
    return mat


def _signed_area(pts):
    a = 0.0
    n = len(pts)
    for i in range(n):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % n]
        a += x0 * y1 - x1 * y0
    return 0.5 * a


def add_extruded_layer(name, polys_mm, z_bot, z_top, mat, bevel=0.0):
    """Build a 100% solid extruded mesh via bmesh with enforced CCW winding and ear-clip triangulation."""
    if not polys_mm:
        return None
    bm = bmesh.new()
    dz = z_top - z_bot
    for pts in polys_mm:
        clean = []
        for x, y in pts:
            if not clean or math.hypot(x - clean[-1][0], y - clean[-1][1]) > 1e-6:
                clean.append((float(x), float(y)))
        if len(clean) >= 2 and math.hypot(clean[0][0] - clean[-1][0], clean[0][1] - clean[-1][1]) <= 1e-6:
            clean.pop()
        if len(clean) < 3:
            continue
        if _signed_area(clean) < 0:
            clean.reverse()
        verts = [bm.verts.new((x, y, z_bot)) for x, y in clean]
        try:
            f = bm.faces.new(verts)
            bmesh.ops.triangulate(bm, faces=[f], quad_method='FIXED', ngon_method='CLIP')
        except ValueError:
            pass

    bm.faces.ensure_lookup_table()
    res = bmesh.ops.extrude_face_region(bm, geom=bm.faces[:])
    extruded_verts = [g for g in res['geom'] if isinstance(g, bmesh.types.BMVert)]
    bmesh.ops.translate(bm, verts=extruded_verts, vec=(0.0, 0.0, dz))
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces[:])

    mesh = bpy.data.meshes.new(f'{name}_mesh')
    bm.to_mesh(mesh)
    bm.free()
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    obj.data.materials.append(mat)
    return obj


def add_realistic_arch_airbridges(name, span_polys_mm, z_base, arch_height, ribbon_thick, mat):
    """Build realistic reflowed-photoresist curved Al arch airbridges (smooth cosine/parabolic arch)."""
    if not span_polys_mm:
        return None
    bm = bmesh.new()
    n_seg = 24

    for pts in span_polys_mm:
        if len(pts) < 4:
            continue
        # Extract oriented bounding box (center, long axis u_hat, short axis v_hat, length L, width W)
        arr = [Vector((float(x), float(y))) for x, y in pts]
        cx = sum(v.x for v in arr) / len(arr)
        cy = sum(v.y for v in arr) / len(arr)
        # Find longest edge -> defines the bridge span direction u_hat
        best_len = -1.0
        u_hat = Vector((1.0, 0.0))
        for i in range(len(arr)):
            d = arr[(i + 1) % len(arr)] - arr[i]
            if d.length > best_len:
                best_len = d.length
                u_hat = d.normalized()
        v_hat = Vector((-u_hat.y, u_hat.x))
        # Project vertices to get exact span L and width W
        us = [(v - Vector((cx, cy))).dot(u_hat) for v in arr]
        vs = [(v - Vector((cx, cy))).dot(v_hat) for v in arr]
        L = max(us) - min(us)
        W = max(vs) - min(vs)
        if L < 1e-4 or W < 1e-4:
            continue

        # Foot pads sit flat on base metal at both ends (18% each side), smooth raised-cosine arch in the middle 64%
        rings = []
        for k in range(n_seg + 1):
            t = k / n_seg  # 0 .. 1 along span
            u = (t - 0.5) * L
            if t <= 0.16 or t >= 0.84:
                z = z_base
            else:
                s = (t - 0.16) / 0.68  # 0 .. 1 across arch
                z = z_base + arch_height * 0.5 * (1.0 - math.cos(2.0 * math.pi * s))
            # Local tangent slope for perpendicular ribbon thickness
            if 0.16 < t < 0.84:
                s = (t - 0.16) / 0.68
                dz_du = (arch_height * math.pi / (0.68 * L)) * math.sin(2.0 * math.pi * s)
            else:
                dz_du = 0.0
            n_len = math.hypot(dz_du, 1.0)
            nu, nz = -dz_du / n_len, 1.0 / n_len

            p_bot_l = Vector((cx, cy, z)) + Vector((u_hat.x * u - v_hat.x * 0.5 * W, u_hat.y * u - v_hat.y * 0.5 * W, 0.0))
            p_bot_r = Vector((cx, cy, z)) + Vector((u_hat.x * u + v_hat.x * 0.5 * W, u_hat.y * u + v_hat.y * 0.5 * W, 0.0))
            off = Vector((u_hat.x * nu * ribbon_thick, u_hat.y * nu * ribbon_thick, nz * ribbon_thick))
            v0 = bm.verts.new(p_bot_l)
            v1 = bm.verts.new(p_bot_r)
            v2 = bm.verts.new(p_bot_r + off)
            v3 = bm.verts.new(p_bot_l + off)
            rings.append((v0, v1, v2, v3))

        # Connect consecutive cross-section rings into a watertight smooth arch tube
        for k in range(n_seg):
            a0, a1, a2, a3 = rings[k]
            b0, b1, b2, b3 = rings[k + 1]
            for quad in ((a0, a1, b1, b0), (a1, a2, b2, b1), (a2, a3, b3, b2), (a3, a0, b0, b3)):
                try:
                    bm.faces.new(quad)
                except ValueError:
                    pass
        # End caps
        for cap in ((rings[0][3], rings[0][2], rings[0][1], rings[0][0]),
                    (rings[-1][0], rings[-1][1], rings[-1][2], rings[-1][3])):
            try:
                bm.faces.new(cap)
            except ValueError:
                pass

    bm.faces.ensure_lookup_table()
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces[:])
    mesh = bpy.data.meshes.new(f'{name}_mesh')
    bm.to_mesh(mesh)
    bm.free()
    # Smooth shading on the curved arch faces!
    for poly in mesh.polygons:
        poly.use_smooth = True
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    obj.data.materials.append(mat)
    return obj


def _add_oriented_box(bm, cx, cy, z0, z1, length, width, angle_rad):
    """Add an oriented 3D rectangular prism into a bmesh."""
    ca, sa = math.cos(angle_rad), math.sin(angle_rad)
    hl, hw = 0.5 * length, 0.5 * width
    corners = [(-hl, -hw), (+hl, -hw), (+hl, +hw), (-hl, +hw)]
    b_verts = [bm.verts.new((cx + lx * ca - ly * sa, cy + lx * sa + ly * ca, z0)) for lx, ly in corners]
    t_verts = [bm.verts.new((cx + lx * ca - ly * sa, cy + lx * sa + ly * ca, z1)) for lx, ly in corners]
    bm.faces.new(b_verts[::-1])
    bm.faces.new(t_verts)
    for i in range(4):
        j = (i + 1) % 4
        bm.faces.new((b_verts[i], b_verts[j], t_verts[j], t_verts[i]))


def add_manhattan_jjs(sites, mat_bot, mat_top, mat_oxide):
    """Build 3D crossed orthogonal Al/AlOx/Al Manhattan JJs + opposing 45-deg bandage pads at every site."""
    if not sites:
        return
    bm_bot, bm_top, bm_ox = bmesh.new(), bmesh.new(), bmesh.new()

    for s in sites:
        cx, cy = float(s['x']), float(s['y'])
        ori = math.radians(float(s['ori']))
        # Local frame (u along ori connects opposing tapered Nb leads at u = -0.008 and u = +0.008):
        # Finger 1 at +45 deg lands on -u lead; Finger 2 at -45 deg lands on +u lead;
        # they cross orthogonally (90 deg) right at (cx, cy)!
        a1 = ori + math.radians(45.0)
        a2 = ori - math.radians(45.0)

        # Bottom Al finger (sits cleanly at top of Nb lead z=0.010, crossing at +45 deg with 2 um stub)
        dx1 = -0.0022 * math.cos(a1)
        dy1 = -0.0022 * math.sin(a1)
        _add_oriented_box(bm_bot, cx + dx1, cy + dy1, 0.0100, 0.0126, 0.0115, 0.0014, a1)
        # Compact bandage pad landing cleanly on top of the -u tapered Nb lead
        pad1_x = cx - 0.0058 * math.cos(a1)
        pad1_y = cy - 0.0058 * math.sin(a1)
        _add_oriented_box(bm_bot, pad1_x, pad1_y, 0.0098, 0.0136, 0.0038, 0.0034, a1)

        # Top Al finger (crosses orthogonally at -45 deg over the bottom finger from the +u lead!)
        dx2 = +0.0022 * math.cos(a2)
        dy2 = +0.0022 * math.sin(a2)
        _add_oriented_box(bm_top, cx + dx2, cy + dy2, 0.0128, 0.0154, 0.0115, 0.0014, a2)
        # Compact bandage pad landing cleanly on top of the +u tapered Nb lead
        pad2_x = cx + 0.0058 * math.cos(a2)
        pad2_y = cy + 0.0058 * math.sin(a2)
        _add_oriented_box(bm_top, pad2_x, pad2_y, 0.0098, 0.0160, 0.0038, 0.0034, a2)

        # Glowing AlOx tunnel barrier right at the orthogonal cross (cx, cy)!
        _add_oriented_box(bm_ox, cx, cy, 0.0122, 0.0136, 0.0018, 0.0018, a1)

    for nm, bm, mat in (
        ('Manhattan_JJ_BottomAl', bm_bot, mat_bot),
        ('Manhattan_JJ_TopAl', bm_top, mat_top),
        ('Manhattan_JJ_AlOx', bm_ox, mat_oxide),
    ):
        bmesh.ops.recalc_face_normals(bm, faces=bm.faces[:])
        mesh = bpy.data.meshes.new(f'{nm}_mesh')
        bm.to_mesh(mesh)
        bm.free()
        obj = bpy.data.objects.new(nm, mesh)
        bpy.context.collection.objects.link(obj)
        obj.data.materials.append(mat)


def add_box(name, cx, cy, cz, sx, sy, sz, mat):
    bpy.ops.mesh.primitive_cube_add(size=1.0, location=(cx, cy, cz))
    obj = bpy.context.active_object
    obj.name = name
    obj.scale = (sx, sy, sz)
    obj.data.materials.append(mat)
    return obj


def add_area_light(name, loc, target, energy, size, color=(1.0, 1.0, 1.0)):
    ldata = bpy.data.lights.new(name=name, type='AREA')
    ldata.energy = energy
    ldata.size = size
    ldata.color = color
    lobj = bpy.data.objects.new(name=name, object_data=ldata)
    lobj.location = loc
    direction = Vector(target) - Vector(loc)
    lobj.rotation_euler = direction.to_track_quat('-Z', 'Y').to_euler()
    bpy.context.collection.objects.link(lobj)
    return lobj


def setup_camera(name, loc, target, lens=65.0, fstop=None):
    cdata = bpy.data.cameras.new(name=name)
    cdata.lens = lens
    cdata.clip_start = 0.01
    cdata.clip_end = 500.0
    if fstop is not None:
        cdata.dof.use_dof = True
        cdata.dof.focus_distance = float((Vector(target) - Vector(loc)).length)
        cdata.dof.aperture_fstop = fstop
    cobj = bpy.data.objects.new(name=name, object_data=cdata)
    cobj.location = loc
    direction = Vector(target) - Vector(loc)
    cobj.rotation_euler = direction.to_track_quat('-Z', 'Y').to_euler()
    bpy.context.collection.objects.link(cobj)
    return cobj


def main():
    data = json.loads(POLY_JSON.read_text(encoding='utf-8'))
    clear_scene()

    # World background (bright clean studio cyclorama for crisp silver metallic reflections)
    world = bpy.data.worlds.new('StudioWorld')
    bpy.context.scene.world = world
    world.use_nodes = True
    bg = world.node_tree.nodes.get('Background')
    bg.inputs['Color'].default_value = (0.38, 0.42, 0.48, 1.0)
    bg.inputs['Strength'].default_value = 1.65

    # Materials (bright polished silver Nb on deep blue-grey Si, gold arch airbridges, Manhattan JJs!)
    mat_si = make_material('Silicon_Substrate', (0.045, 0.065, 0.095), metallic=0.35, roughness=0.18, coat=0.5)
    mat_metal = make_material('BaseMetal_Nb', (0.88, 0.91, 0.96), metallic=0.90, roughness=0.22)
    mat_ab = make_material('Airbridge_GoldAl', (0.98, 0.74, 0.18), metallic=0.96, roughness=0.18)
    mat_jj_bot = make_material('Manhattan_BotAl', (0.45, 0.78, 0.98), metallic=0.95, roughness=0.16)
    mat_jj_top = make_material('Manhattan_TopAl', (0.95, 0.42, 0.25), metallic=0.95, roughness=0.16)
    mat_jj_ox = make_material('Manhattan_AlOx', (1.0, 0.85, 0.40), metallic=0.2, roughness=0.15,
                              emission=(1.0, 0.55, 0.15), emission_strength=4.5)
    mat_floor = make_material('StudioFloor', (0.015, 0.018, 0.025), metallic=0.1, roughness=0.6)

    # 1. Silicon die substrate (16 x 16 x 0.50 mm) + subtle studio floor
    die_w, die_h = data['die_mm']
    add_box('Silicon_Die', 0.0, 0.0, -0.25, die_w, die_h, 0.50, mat_si)
    add_box('Studio_Floor', 0.0, 0.0, -0.55, 80.0, 80.0, 0.10, mat_floor)

    # 2. Base metal + realistic smooth reflowed-resist curved Al/Au ARCH airbridges + 65 crossed Manhattan JJs!
    add_extruded_layer('BaseMetal_L1', data['layer_1'], z_bot=0.000, z_top=0.010, mat=mat_metal)
    add_realistic_arch_airbridges(
        'Airbridges_Arch_L20', data['layer_20'],
        z_base=0.009, arch_height=0.024, ribbon_thick=0.0040, mat=mat_ab,
    )
    add_manhattan_jjs(data.get('manhattan_sites', []), mat_jj_bot, mat_jj_top, mat_jj_ox)

    # 3. Bright clean studio lighting (soft overhead + raking key/fill/rim for crisp 3D depth)
    add_area_light('KeyLight', (8.0, -9.0, 5.5), (0.0, 0.0, 0.0), energy=2200.0, size=7.0, color=(1.0, 0.96, 0.90))
    add_area_light('FillLight', (-10.0, -6.0, 7.0), (0.0, 0.0, 0.0), energy=1200.0, size=10.0, color=(0.86, 0.93, 1.0))
    add_area_light('RimLight', (-2.0, 12.0, 4.5), (0.0, 0.0, 0.0), energy=1600.0, size=8.0, color=(0.92, 0.96, 1.0))
    add_area_light('TopSoftbox', (1.0, 0.5, 9.0), (1.0, 0.5, 0.0), energy=1400.0, size=12.0, color=(1.0, 0.99, 0.98))

    # 4. Configure Cycles CPU (denoised)
    scene = bpy.context.scene
    scene.render.engine = 'CYCLES'
    scene.cycles.device = 'CPU'
    scene.cycles.samples = 56
    scene.cycles.use_denoising = True
    scene.render.resolution_x = 1600
    scene.render.resolution_y = 1050
    scene.render.image_settings.file_format = 'PNG'
    scene.view_settings.view_transform = 'AgX'
    scene.view_settings.exposure = 0.65

    # 5. Low-angle zoomed-in 3D cameras: all-in-one hero (Qubit + Coupler + Arch Airbridges + Manhattan JJs)!
    cam_center = setup_camera('Cam_Center3D', loc=(2.65, -3.55, 1.55), target=(0.15, 0.05, 0.03), lens=62.0, fstop=5.6)
    cam_macro = setup_camera('Cam_MacroBridge', loc=(1.88, -0.58, 0.44), target=(1.18, 0.48, 0.04), lens=72.0, fstop=4.0)
    cam_jj = setup_camera('Cam_ManhattanJJ', loc=(1.42, 0.38, 0.19), target=(1.20, 0.60, 0.006), lens=82.0, fstop=8.0)
    cam_hero = setup_camera('Cam_HeroAll', loc=(0.42, -0.58, 1.02), target=(1.14, 0.52, 0.012), lens=52.0, fstop=11.0)

    scene.camera = cam_hero
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_OUT))

    scene.render.filepath = str(PNG_HERO)
    bpy.ops.render.render(write_still=True)
    print('BLENDER_RENDER_OK:', BLEND_OUT, PNG_HERO)


if __name__ == '__main__':
    main()
