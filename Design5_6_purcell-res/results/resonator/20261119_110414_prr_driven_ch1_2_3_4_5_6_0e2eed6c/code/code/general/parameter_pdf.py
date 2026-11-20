"""Render per-design parameter / simulation-result sheets (A4 PDF) from a JSON spec.

    python -I -B code/general/parameter_pdf.py <spec.json> [--out-dir DIR] [--combined report.pdf]

--combined writes every sheet (each starting on a new page, continuous page numbers) into one PDF.
A section may carry "image": path (relative to the spec file) and "image_height": inches.

Spec: {"sheets": [{"file": "X.pdf", "title": str, "subtitle": str, "header": [[key, value], ...],
                   "sections": [{"title": str, "note": str?, "rows": [[parameter, value, unit, status_and_source], ...]}],
                   "missing": [str, ...]}]}
The renderer only lays out text; every number and status comes from the spec, which cites its sources.
"""
from pathlib import Path
import argparse
import json
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

W, H = 8.27, 11.69                     # A4 portrait (inch)
LEFT, RIGHT, TOP, BOTTOM = 0.55, 7.72, 11.1, 0.6
COLS = [(0.55, 30), (2.95, 26), (5.0, 9), (5.62, 36)]  # (x inch, wrap chars) parameter, value, unit, status/source
LINE = 0.145
FONT = 7.2


class Sheet:
    def __init__(self, pdf, title, counter=None):
        self.pdf, self.title, self.page = pdf, title, 0
        self.counter = counter if counter is not None else {'n': 0}
        self.new_page()

    def new_page(self):
        if self.page:
            self.finish()
        self.page += 1
        self.counter['n'] += 1
        self.fig = plt.figure(figsize=(W, H))
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, W); self.ax.set_ylim(0, H); self.ax.axis('off')
        self.y = TOP
        if self.page > 1:
            self.text(LEFT, self.y, f'{self.title} (continued)', size=8, color='0.4')
            self.y -= 0.3

    def finish(self):
        self.text(RIGHT, 0.35, f'page {self.counter["n"]}', size=6.5, color='0.5', ha='right')
        self.pdf.savefig(self.fig); plt.close(self.fig)

    def text(self, x, y, s, size=FONT, **kw):
        self.ax.text(x, y, s, fontsize=size, va='top', family='DejaVu Sans', **kw)

    def need(self, h):
        if self.y - h < BOTTOM:
            self.new_page()

    def paragraph(self, s, size=FONT, width=None, color='0.15', gap=0.06, **kw):
        width = width or int(140 * FONT / size)      # characters per line at this font size
        lines = textwrap.wrap(s, width) or ['']
        self.need(len(lines) * LINE * size / FONT)
        for ln in lines:
            self.text(LEFT, self.y, ln, size=size, color=color, **kw)
            self.y -= LINE * size / FONT
        self.y -= gap

    def table(self, rows, header=('Parameter', 'Value', 'Unit', 'Status / source')):
        def wrapped(row):
            return [textwrap.wrap(str(c), w) or [''] for c, (_, w) in zip(row, COLS)]
        self.need(LINE * 3)
        self.ax.add_patch(plt.Rectangle((LEFT - .05, self.y - LINE - .03), RIGHT - LEFT + .1, LINE + .06,
                                        color='#dfe4ea', lw=0))
        for (x, _), h in zip(COLS, header):
            self.text(x, self.y, h, weight='bold')
        self.y -= LINE + 0.07
        for k, row in enumerate(rows):
            cells = wrapped(row)
            n = max(len(c) for c in cells)
            if self.y - n * LINE < BOTTOM:
                self.new_page()
                self.text(LEFT, self.y, '(table continued)', size=6.5, color='0.5'); self.y -= LINE
            if k % 2:
                self.ax.add_patch(plt.Rectangle((LEFT - .05, self.y - n * LINE - .02), RIGHT - LEFT + .1,
                                                n * LINE + .04, color='#f4f6f8', lw=0))
            for (x, _), c in zip(COLS, cells):
                for i, ln in enumerate(c):
                    self.text(x, self.y - i * LINE, ln, weight='bold' if c is cells[1] else 'normal')
            self.y -= n * LINE + 0.05
        self.y -= 0.1


    def image(self, path, height):
        import matplotlib.image as mpimg
        img = mpimg.imread(str(path))
        h = min(height, TOP - BOTTOM - 0.3)
        w = min(RIGHT - LEFT, h * img.shape[1] / img.shape[0])
        h = w * img.shape[0] / img.shape[1]
        self.need(h + 0.1)
        ax = self.fig.add_axes([LEFT / W, (self.y - h) / H, w / W, h / H])
        ax.imshow(img); ax.axis('off')
        self.y -= h + 0.15


def render(sheet, out, base=Path('.')):
    with PdfPages(out) as pdf:
        render_into(pdf, sheet, base)


def render_into(pdf, sheet, base=Path('.'), counter=None):
    if True:
        s = Sheet(pdf, sheet['title'], counter)
        s.text(LEFT, s.y, sheet['title'], size=15, weight='bold'); s.y -= 0.32
        s.paragraph(sheet.get('subtitle', ''), size=8.5, color='0.3', gap=0.12)
        for k, v in sheet.get('header', []):
            s.need(LINE)
            s.text(LEFT, s.y, k, weight='bold');
            lines = textwrap.wrap(str(v), 118) or ['']
            for i, ln in enumerate(lines):
                s.text(1.95, s.y - i * LINE, ln)
            s.y -= LINE * len(lines) + 0.03
        s.y -= 0.12
        for sec in sheet['sections']:
            s.need(0.6)
            s.text(LEFT, s.y, sec['title'], size=10.5, weight='bold', color='#1f3b5c'); s.y -= 0.24
            if sec.get('note'):
                s.paragraph(sec['note'], size=6.8, color='0.3')
            if sec.get('rows'):
                s.table(sec['rows'])
            if sec.get('image'):
                s.image(base / sec['image'], sec.get('image_height', 4.0))
            for para in sec.get('paragraphs', []):
                s.paragraph(para, size=7.4)
        if sheet.get('missing'):
            s.need(0.6)
            s.text(LEFT, s.y, 'Missing, unvalidated or to decide', size=10.5, weight='bold', color='#8a2b1d')
            s.y -= 0.24
            for m in sheet['missing']:
                s.paragraph('- ' + m, size=7.2)
        s.finish()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('spec', type=Path)
    ap.add_argument('--out-dir', type=Path)
    ap.add_argument('--combined', help='write all sheets into this one PDF (in --out-dir)')
    a = ap.parse_args()
    spec = json.loads(a.spec.read_text(encoding='utf-8'))
    out_dir = a.out_dir or a.spec.parent
    base = a.spec.parent
    if a.combined:
        out, counter = out_dir / a.combined, {'n': 0}
        with PdfPages(out) as pdf:
            for sheet in spec['sheets']:
                render_into(pdf, sheet, base, counter)
        print(out)
        return
    for sheet in spec['sheets']:
        out = out_dir / sheet['file']
        render(sheet, out, base)
        print(out)


if __name__ == '__main__':
    main()
