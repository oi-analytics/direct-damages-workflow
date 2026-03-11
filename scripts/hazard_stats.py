"""Plot the hazard statistics as line plots."""
#%%
import pandas as pd
import rasterio
from pathlib import Path
from glob import glob
from tqdm import tqdm
import matplotlib.pyplot as plt

from direct_damages.utils import read_config

config = read_config("../workflow/config.yaml")
hazdir = Path(config["input"]) / "hazards"
figdir = Path(config["figs"])
resdir = Path(config["results"])

figdir.mkdir(exist_ok=True, parents=True)

hazards = []
epochs = []
scenarios = []
rps = []
totals = []
means = []
mins = []
maxs = []
stds = []
counts = []

# calculate key stats for all rasters
for hazard in config["hazards"]:
    rasters = glob(str(hazdir / f"{hazard}*.tif"))

    for tif in tqdm(rasters):
        name = Path(tif).stem.split('_')
        epoch, scen, rp = name[1], name[2], name[3]

        with rasterio.open(tif) as src:
            stats = src.statistics(1)
            count = src.width * src.height
        
        totals.append(stats.mean * count)
        means.append(stats.mean)
        mins.append(stats.min)
        maxs.append(stats.max)
        stds.append(stats.std)
        counts.append(count)
        hazards.append(hazard)
        epochs.append(int(epoch))
        scenarios.append(scen)
        rps.append(int(rp.replace('rp', '')))


df = pd.DataFrame({
    "hazard": hazards,
    "epoch": epochs,
    "scenario": scenarios,
    "rp": rps,
    "total": totals,
    "mean": means,
    "min": mins,
    "max": maxs,
    "std": stds,
    "n": counts
})

df["scenario"] = pd.Categorical(df["scenario"], categories=["historical", "rcp26", "ssp126", "rcp45", "ssp245", "rcp85", "ssp585"], ordered=True)
df = df.sort_values(by=["epoch", "scenario", "rp"])

# plotting settings
rps = [100, 250]
stat = "mean"
lw = 3
alpha = 0.3
palette = {
    "low": "#f4a582",
    "mid": "#d6604d",
    "high": "#8b0000"
}

fig, axs = plt.subplots(4, 2, figsize=(12, 8))

for ax in axs.flat[len(config["hazards"]):]:
    ax.set_visible(False)

for ax, hazard in zip(axs.flat, config["hazards"]):
    df_haz = df[df["hazard"] == hazard].copy()
    df_haz = df_haz[df_haz["rp"].apply(lambda x: x in rps)].copy()

    scens = list(df_haz["scenario"].unique())
    scens = [s for s in scens if s != "historical"]

    y = []
    labels = []
    for i, scen in enumerate(scens):
        df_haz_scen = df_haz[df_haz["scenario"].apply(lambda x: x in ["historical", scen])]
        x = df_haz_scen.epoch
        y.append(df_haz_scen[stat])
        labels.append(scen)

    ax.plot(x, y[0], lw=lw, color=palette["low"], label=labels[0])
    ax.plot(x, y[1], lw=lw, color=palette["mid"], label=labels[1])
    ax.plot(x, y[2], lw=lw, color=palette["high"], label=labels[2])

    ax.fill_between(x, y[0], y[1], alpha=alpha, color=palette["mid"])
    ax.fill_between(x, y[1], y[2], alpha=alpha, color=palette["high"])
        
    ax.legend(loc="best", frameon=False)
    ax.set_title(hazard, fontsize=11)
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.subplots_adjust(bottom=0.08)

fig.savefig(figdir / "hazard_stats.png", transparent=True, dpi=300)
df.to_csv(resdir / "hazard_stats.csv")
# %%
