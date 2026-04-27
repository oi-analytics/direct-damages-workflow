#%%
import pandas as pd
import geopandas as gpd
from pathlib import Path
from glob import glob
from tqdm import tqdm
import matplotlib.pyplot as plt

from direct_damages.utils import read_config

config = read_config("../workflow/config.yaml")

indir = Path(config["input"])

    
admin = gpd.read_parquet(indir / "admin" / "level01.geoparquet")
# %%
admin.plot()
# %%
exclude = ['kaskazini_unguja', 'kusini_unguja', 'mjini_magharibi', 'kaskazini_pemba', 'kusini_pemba']
admin = admin.set_index('subregion')
# %%
fig, ax = plt.subplots(figsize=(10, 10))
admin.boundary.plot(ax=ax, color='k', linewidth=0.5)
admin.loc[exclude].plot(ax=ax, color='crimson')

# %%
