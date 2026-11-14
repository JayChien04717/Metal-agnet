"""Plot eigenmode frequency trajectories and spatial electric-field scores."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

folder = Path(__file__).resolve().parent / "results/layout3_full/readout_eigenmode"
fig, axes = plt.subplots(2, 2, figsize=(10, 6), layout="constrained")
for i, side in enumerate(("left", "right")):
    path = folder / f"local_eigen_{side}_fields.csv"
    if not path.exists():
        continue
    fields = pd.read_csv(path)
    freq = pd.read_csv(folder / f"local_eigen_{side}_freqs.csv", index_col=0)
    for mode, column in enumerate(freq, 1):
        axes[i,0].plot(freq.index, freq[column], "-o", ms=4, label=f"Mode {mode}")
    axes[i,0].set(xlabel="Adaptive pass", ylabel="Frequency (GHz)", title=f"{side.capitalize()} fixture")
    axes[i,0].set_xticks(freq.index)
    axes[i,0].legend(frameon=False, fontsize=8)
    ax = axes[i,1]
    base = np.zeros(len(fields))
    for key,label,color in (("tee","Coupling end","#317cb1"),("meander","Meander","#80b1ce"),
                            ("readout_pad","Readout pad","#a3c28a"),("qubit","Fluxonium region","#b77863")):
        values=fields[key].to_numpy()
        ax.bar(fields["mode"],values,bottom=base,label=label,color=color)
        base += values
    ax.set(xticks=fields["mode"], xticklabels=[f"{f:.3f} GHz" for f in fields.frequency_GHz],
           ylim=(0,1), ylabel="Fraction of sampled |E|²")
    if i == 0:
        ax.legend(frameon=False, fontsize=8, ncol=2)
fig.savefig(folder/"readout_eigenmode_fields.png",dpi=180)
