"""Plot the deltas for each hazard."""

"""Plot the hazard statistics as maps."""
#%%
from glob import glob
from pathlib import Path
import matplotlib.pyplot as plt
import geopandas as gpd
import xarray as xr
from direct_damages.utils import read_config


config = read_config("../workflow/config.yaml")
hazdir = Path(config["input"]) / "hazards"
admdir = Path(config["input"]) / "admin"
figdir = Path(config["figs"])

figdir.mkdir(exist_ok=True, parents=True)

# choose
hazard = "coastal"
scenario = "ssp245"
returnp  = "rp00100"
baseyr = 2020
years = [2030, 2050, 2080]
units = "m"
coarsen = 100 # set to 1 for identity

# glob list of rasters
rasters = glob(str(hazdir / f"{hazard}*.tif"))
baseline = f"{hazard}_{baseyr}_historical_{returnp}.tif"

admin = gpd.read_parquet(admdir / "level01.geoparquet").to_crs(4326)

# %%
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature

base = xr.open_dataarray(hazdir / baseline)
base = base.coarsen(x=100, y=100, boundary="trim").max()

# %%
xmin = base.x.min().item()
xmax = base.x.max().item()
ymin = base.y.min().item()
ymax = base.y.max().item()

# %%
fig, axs = plt.subplots(1, 4, figsize=(8, 2), subplot_kw=dict(projection=ccrs.PlateCarree()))

# axs[0].imshow(base, cmap="OrRd")
base_plot = base.plot(ax=axs[0], cmap="summer_r", add_colorbar=False)
axs[0].set_title('')

deltas = []
vmin = float("inf")
vmax = -float("inf")
for i, year in enumerate(years):
    raster = f"{hazard}_{year}_{scenario}_{returnp}.tif"
    raster_da = xr.open_dataarray(hazdir / raster)
    raster_da = raster_da.coarsen(x=100, y=100, boundary="trim").max()
    raster_aligned, base_aligned = xr.align(raster_da, base, join="override")
    delta = raster_aligned - base_aligned

    dmin = delta.min().item()
    dmax = delta.max().item()
    if dmin < vmin:
        vmin = dmin
    if dmax > vmax:
        vmax = dmax
    deltas.append(delta)

for i, delta in enumerate(deltas):
    delta_plot = delta.plot(ax=axs[i+1], cmap="coolwarm", add_colorbar=False, vmin=vmin, vmax=vmax)
    axs[i+1].set_title(years[i], fontsize=9)

for ax in axs.flat:
    ax.add_feature(cfeature.COASTLINE, linewidth=0.2, zorder=10)
    ax.add_feature(cfeature.OCEAN, color="#d4d4dc", zorder=9)
    admin.boundary.plot(ax=ax, linewidth=0.2, color='k')
    ax.add_feature(cfeature.LAND, color="#fafaf8")
    ax.set_extent([xmin, xmax, ymin, ymax], crs=ccrs.PlateCarree())

fig.colorbar(
    base_plot, 
    ax=axs[0], 
    orientation="horizontal", 
    pad=0.05, 
    label=f"baseline ({units})",
    shrink=0.9
)

# Shared colorbar for all deltas (spanning axs[1:])
fig.colorbar(
    delta_plot,  # reuse last delta plot handle (all share same vmin/vmax)
    ax=axs[1:],  # span all delta axes
    orientation="horizontal",
    pad=0.05,
    label=f"shift ({units})",
    shrink=0.3
)

fig.savefig(figdir / f"{hazard}_deltas.png", transparent=True, dpi=300)
# %%
fig.savefig(figdir / f"{hazard}_deltas.pdf", transparent=True, dpi=300)
# %%
