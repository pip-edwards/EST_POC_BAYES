"""


Author: Pippa Edwards (UoS, NOCS)

This code:
- Creates maps of supplementary figure

"""
#%%

#%%
#import packages
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import seaborn as sns
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcolors
from cmcrameri import cm
from scipy.stats import ks_2samp
from scipy.stats import pearsonr
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pyproj import Geod
from cmcrameri import cm
#%%
#DATASETS

fp = ""
poc = pd.read_csv(f"{fp}input_data/merged_POC_all_ugchla_1103.csv")
allpoc = pd.read_csv(f"{fp}input_data/Global_POC_Database_2025-07-01.csv")
allpoc = allpoc[["poc", "latitude", "longitude", "depth", "date_formatted"]]
allpoc = allpoc.dropna()
sns.set_theme(context = "paper", style = "ticks", font = "arial")
#%%

#set up figure
fig = plt.figure(figsize=(8, 6), dpi = 300)
ax = plt.axes(projection=ccrs.PlateCarree())
land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                    ,edgecolor="#1B1B1B",facecolor="#D2D2D2") 
ax.add_feature(land, zorder=0)
sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                   edgecolor='face',facecolor="#99C2E8")
ax.add_feature(sea, zorder=0)
ax.coastlines(resolution='110m', alpha=0.3)
ax.gridlines(draw_labels = True, linestyle = "--", color = "#454545", alpha  = 0.5)
#add points
ax.scatter(allpoc["longitude"], allpoc["latitude"],transform=ccrs.PlateCarree(), 
            s=15, color="#00203C", linewidth = 0.2, label = "Available Data")
ax.scatter(poc["Longitude"], poc["Latitude"],transform=ccrs.PlateCarree(), 
             s=10, color="#f8a436", linewidth = 0.2, label = "Matched Data")
plt.legend(ncols =2, bbox_to_anchor=(0.7, -0.05))
plt.savefig(f"{fp}figs/overall_map.png")
plt.show()
# %%
