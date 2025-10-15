import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle, Patch

rng = np.random.default_rng(7)
rows = []
for chr_ in range(1, 23):  
    n_blocks = rng.integers(60, 110)
    local_rg = rng.normal(0.0, 0.28, n_blocks)
    cov      = rng.uniform(0.75, 1.00, n_blocks)   
    rows += list(zip([chr_]*n_blocks, range(n_blocks), local_rg, cov))

df = pd.DataFrame(rows, columns=["CHR","block_idx","local_rg","coverage"])\
       .sort_values(["CHR","block_idx"], ignore_index=True)

COV_THR = 0.90
df["lowcov"]   = df["coverage"] < COV_THR
df["excluded"] = False
df.loc[(df["CHR"] == 6) & (df["block_idx"] < 8), "excluded"] = True 

n_chr        = 22
block_counts = df.groupby("CHR")["block_idx"].max().add(1)  
max_blocks   = int(block_counts.max())

vals     = np.full((n_chr, max_blocks), np.nan)
mask_low = np.zeros_like(vals, dtype=bool)
mask_exc = np.zeros_like(vals, dtype=bool)

for chr_ in range(1, 23):
    d   = df[df["CHR"] == chr_]
    idx = d["block_idx"].to_numpy()
    vals[chr_-1, idx]     = d["local_rg"].to_numpy()
    mask_low[chr_-1, idx] = d["lowcov"].to_numpy()
    mask_exc[chr_-1, idx] = d["excluded"].to_numpy()

vals_T, mask_low_T, mask_exc_T = vals.T, mask_low.T, mask_exc.T

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "font.size": 11,
    "axes.edgecolor": "#333", "axes.titleweight": "bold"
})

fig, ax = plt.subplots(figsize=(14, 8), constrained_layout=True)

cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad("#E6E6E6")
norm = TwoSlopeNorm(vmin=np.nanmin(vals_T), vcenter=0.0, vmax=np.nanmax(vals_T))
im = ax.imshow(vals_T, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")

ax.set_xticks(np.arange(n_chr))
ax.set_xticklabels([str(i) for i in range(1, 23)])
ax.set_xlabel("Chromosome")
ax.set_ylabel("LD block index (equal width)")
ax.set_title("LD-block heatmap of local r$_g$")

step = 10
ax.set_yticks(np.arange(0, max_blocks, step))
ax.grid(axis="y", color="#000", alpha=0.08, linewidth=0.6)

for c in range(n_chr):
    n_b = int(block_counts.iloc[c])
    for r in np.where(mask_low_T[:n_b, c])[0]:
        ax.add_patch(Rectangle((c-0.5, r-0.5), 1, 1,
                               facecolor="#000000", alpha=0.18, edgecolor="none"))
    for r in np.where(mask_exc_T[:n_b, c])[0]:
        ax.add_patch(Rectangle((c-0.5, r-0.5), 1, 1,
                               facecolor="none", edgecolor="#666",
                               hatch="////", linewidth=0.0, alpha=0.9))

secax = ax.secondary_xaxis('top')
secax.set_xticks(np.arange(n_chr))
secax.set_xticklabels([str(int(n)) for n in block_counts])
secax.set_xlabel("# LD blocks")

ax.set_xlim(-0.5, n_chr - 0.5)
ax.set_ylim(-0.5, max_blocks - 0.5)

cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
cbar.set_label("local r$_g$ (diverging, centered at 0)")
legend_handles = [
    Patch(facecolor="#000000", alpha=0.18, label=f"Coverage < {COV_THR:.2f}"),
    Patch(facecolor="none", edgecolor="#666", hatch="////", label="Excluded region")
]
ax.legend(handles=legend_handles, loc="upper right", frameon=True, framealpha=0.9)

# plt.savefig("LDblock_local_rg_heatmap.pdf", bbox_inches="tight")
plt.show()
