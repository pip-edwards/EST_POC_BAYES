#%%
#Code to make poc and b maps 
#NOTE: ADDING FEATURES AND COLOURBARS/LABELS WAS TROUBLESHOOTED WITH CHATGPT

#packages
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import os
from cmcrameri import cm
import glob
import xarray as xr
import cartopy.feature as cfeature
from cmcrameri import cm
from scipy.stats import linregress
sns.set_context("paper")
#set up fp
fp = "C:/Users/pe1n24/OneDrive - University of Southampton/EST_POC_BAYES/"
os.chdir(fp)

#only list usable files
files = []
for f in os.listdir(fp):
    if "_run" in f:
        #print(f)
        sumstats = glob.glob(f"{f}/*summary_stats.csv")
        #print(sumstats)
        if len(sumstats) != 1:
            print(f, "SOMETHING IS WRONG")
        else:
            rhats = pd.read_csv(sumstats[0])
            if 0.99 < np.mean(rhats["R_hat"]) < 1.01:
                files.append(f)
                print("ANALYSE:", f, np.mean(rhats["R_hat"]))
            else:
                print("DROPPED:", f, np.mean(rhats["R_hat"]))

print(files)

#mpa of where there is data (for the white line on maps)
allsum = xr.open_dataset(f"{fp}input_data/months_of_data.nc")["Months of Data"]
#%%

for file in files:
    print(file)
    #get naming convention
    parts = f.split("_")
    met = parts[0]
    chl = parts[1]

    #check if there is a folder to put the figures into 
    if not os.path.isdir(f"{fp}/figs/{file}"):
        os.makedirs(f"{fp}/figs/{file}")

    #import datasets
    betas = pd.read_csv(glob.glob(f"{fp}{f}/*_beta_vals.csv")[0])
    gammas = pd.read_csv(glob.glob(f"{fp}{f}/*_gamma_vals.csv")[0])
    
    #set up vars for each gamma (to make code lines shorter)
    gamma1 = np.mean(gammas["gamma.1"])
    gamma2 = np.mean(gammas["gamma.2"])
    gamma3 = np.mean(gammas["gamma.3"])
    gamma4 = np.mean(gammas["gamma.4"])
    gamma5 = np.mean(gammas["gamma.5"])
    gamma6 = np.mean(gammas["gamma.6"])
    gamma7 = np.mean(gammas["gamma.7"])

    beta1 = np.mean(betas["beta.1"])
    beta2 = np.mean(betas["beta.2"])
    beta3 = np.mean(betas["beta.3"])
    beta4 = np.mean(betas["beta.4"])
    beta5 = np.mean(betas["beta.5"])
    beta6 = np.mean(betas["beta.6"])
    beta7 = np.mean(betas["beta.7"])

    #import climatologies
    #import chla
    chla_ntd = xr.open_dataset(f"{fp}input_data/occci/occci_overall_climatology.nc")["chlor_a"]
    #transform chla depending on dataset
    if chl == "ugchla":
        chla_ntd = np.log(chla_ntd * 1000)
    else:
        chla_ntd = np.log(chla_ntd)
    #import sst and depth (Chla will be later)
    #import sst and depth and trasnform
    depth = xr.open_dataset(f"{fp}input_data/bathy/depth100_map.nc")["depth"]
    depth = depth.sortby("lat")
    depth = np.log(depth)
    sst_ntd = xr.open_dataset(f"{fp}input_data/sst/SST_overall_climatology.nc")["SST"]
    sst_ntd = np.log(sst_ntd + 1.79)

    #get lat and lon for plotting/mapping purposes
    lats = sst_ntd["lat"]
    lons = sst_ntd["lon"]

    mumap = beta1*np.ones((180,360)) + sst_ntd*beta2 + chla_ntd*beta3 + depth*beta4 + sst_ntd*chla_ntd*beta5 + sst_ntd*depth*beta6 + chla_ntd*depth*beta7
    sigmap =  gamma1*np.ones((180,360)) + sst_ntd*gamma2 + chla_ntd*gamma3 + depth*gamma4 + sst_ntd*chla_ntd*gamma5 + sst_ntd*depth*gamma6 + chla_ntd*depth*gamma7
    pocmap = mumap + (sigmap**2)/2


    fig = plt.figure(figsize=(10, 6), dpi = 300)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines('110m', alpha=0.1)
    ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")

    #MEAN POC map
    pocrange = np.arange(2.8, 6, 0.4)
    contour = ax.contourf(lons, lats, pocmap, transform=ccrs.PlateCarree(),
                        cmap= cm.lipari,extend = "both", levels = pocrange)

    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', 
                        shrink=0.8, pad=0.05, extend = "both")
    cbar.set_label(f"Mean ln(POC) (mg/m²/day)", size = 14)#, weight = "bold")
    land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                        ,edgecolor="#1B1B1B",facecolor="#555555") 
    contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                        linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
    ax.add_feature(land, zorder=2)
    sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                    edgecolor='face',facecolor="#ADADAD")
    ax.add_feature(sea, zorder=0)
    cbar.ax.tick_params(labelsize=14)
    plt.savefig(f"{fp}/figs/{file}/climatology_poc_map.png", transparent = True)
    plt.show()


    #MU MAP
    murange = np.arange(1.8, 5, 0.4)
    fig = plt.figure(figsize=(10, 6), dpi = 300)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines('110m', alpha=0.1)
    ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")
    contour = ax.contourf(lons, lats, mumap, transform=ccrs.PlateCarree(),
                        cmap= cm.lipari,extend = "both", levels = murange)
    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', 
                        shrink=0.8, pad=0.05, extend = "both")
    cbar.set_label(f"μ", size = 14)#, weight = "bold")
    land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                        ,edgecolor="#1B1B1B",facecolor="#555555") 
    contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                        linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
    ax.add_feature(land, zorder=2)
    sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                    edgecolor='face',facecolor="#ADADAD")
    ax.add_feature(sea, zorder=0)
    cbar.ax.tick_params(labelsize=14)
    plt.savefig(f"{fp}/figs/{file}/climatology_mu_map.png", transparent = True)
    plt.show()

    #SIG map
    fig = plt.figure(figsize=(10, 6), dpi = 300)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines('110m', alpha=0.1)
    ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")
    sigrange = np.arange(1.2,2.2,0.1)
    contour = ax.contourf(lons, lats, sigmap, transform=ccrs.PlateCarree(),
                        cmap= cm.lipari,extend = "both", levels = sigrange)
    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', 
                        shrink=0.8, pad=0.05, extend = "both")
    cbar.set_label(f"σ", size = 14)#, weight = "bold")
    land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                        ,edgecolor="#1B1B1B",facecolor="#555555") 
    contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                        linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
    ax.add_feature(land, zorder=2)
    sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                    edgecolor='face',facecolor="#ADADAD")
    ax.add_feature(sea, zorder=0)
    cbar.ax.tick_params(labelsize=14)
    plt.savefig(f"{fp}/figs/{file}/climatology_sig_map.png", transparent = True)
    plt.show()

    #BetaZ and GammaZ term generation
    betaz = beta4 + sst_ntd*beta6 + chla_ntd*beta7
    gammaz =  gamma4 + sst_ntd*gamma6 + chla_ntd*gamma7

    #Bz Map
    fig = plt.figure(figsize=(10, 6), dpi = 300)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines('110m', alpha=0.1)
    ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")
    crange = np.arange(-1, 0.1, 0.1)
    contour = ax.contourf(lons, lats, betaz, transform=ccrs.PlateCarree(),
                        cmap= cm.lipari_r,extend = "min", levels = crange)
    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', 
                        shrink=0.8, pad=0.05, extend = "min")
    land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                        ,edgecolor="#1B1B1B",facecolor="#555555") 
    contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                        linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
    ax.add_feature(land, zorder=2)
    sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                    edgecolor='face',facecolor="#ADADAD")
    ax.add_feature(sea, zorder=0)
    cbar.set_label(r"$β_z$", size = 14)
    cbar.ax.tick_params(labelsize=14)
    plt.savefig(f"{fp}/figs/{file}/betaz_map.png", transparent = True)
    plt.show()

    #Yz Map
    fig = plt.figure(figsize=(10, 6), dpi = 300)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines('110m', alpha=0.1)
    ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")
    crange = np.arange(-0.4, 0.05, 0.05)
    contour = ax.contourf(lons, lats, gammaz, transform=ccrs.PlateCarree(),
                        cmap= cm.lipari_r,extend = "both", levels = crange)
    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', 
                        shrink=0.8, pad=0.05, extend = "both")
    land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                        ,edgecolor="#1B1B1B",facecolor="#555555") 
    contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                        linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
    ax.add_feature(land, zorder=2)
    sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                    edgecolor='face',facecolor="#ADADAD")
    ax.add_feature(sea, zorder=0)
    cbar.set_label(r'$γ_z$', size = 14)
    cbar.ax.tick_params(labelsize=14)
    plt.savefig(f"{fp}/figs/{file}/gammaz_map.png", transparent = True)
    plt.show()

    #Map b
    #for each depth:
    for x in np.arange(100,1600,100):
        #print(x)

        #open depth map 
        depth = xr.open_dataset(f"{fp}input_data/bathy/depth{x}_map.nc")["depth"]
        depth = depth.sortby("lat")
        depth = np.log(depth)

        #make new calc for mean poc
        mumap = beta1*np.ones((180,360)) + sst_ntd*beta2 + chla_ntd*beta3 + depth*beta4 + sst_ntd*chla_ntd*beta5 + sst_ntd*depth*beta6 + chla_ntd*depth*beta7
        sigmap =  gamma1*np.ones((180,360)) + sst_ntd*gamma2 + chla_ntd*gamma3 + depth*gamma4 + sst_ntd*chla_ntd*gamma5 + sst_ntd*depth*gamma6 + chla_ntd*depth*gamma7
        pocmap = mumap + (sigmap**2)/2

        #add to a dataframe
        if x == 100:
            depthpoc = pocmap.expand_dims(depth = [x])
        else:
            depthpoc1 = pocmap.expand_dims(depth = [x]) #make a new one for every point that is not the first
            depthpoc = xr.concat([depthpoc, depthpoc1], dim = "depth") #add to the first one

    #set up empty data array
    bmap = np.ones((180,360))
    #set up list of log depths for the linear regression
    depths = np.log(np.arange(100,1600,100))

    #for each lat and lon
    for lat in range(180):
        for lon in range(360):
            #get the pocs value at each depth
            pocvals = depthpoc[:,lat,lon]
            pocvals = np.array(pocvals.values)
            #make sure it has data within it
            i = len(pocvals) - np.isnan(pocvals).sum()
            #if it doesnt, put a nan in
            if i == 0:
                bmap[lat,lon] = np.nan
            #otherwise linear regress for b
            else:
                b,_,_,_,_ = linregress(depths[:i],pocvals[:i])
                bmap[lat,lon] = b

    #plot b map
    fig = plt.figure(figsize=(10, 6), dpi = 300)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines('110m', alpha=0.1)
    ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")
    #bmap = bmap.where(bmap<0, np.nan)
    crange = np.arange(-1, -0.1, 0.1)
    contour = ax.contourf(lons, lats, bmap, transform=ccrs.PlateCarree(),
                        cmap= cm.lipari_r,extend = "both", levels = crange)
    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', 
                        shrink=0.8, pad=0.05, extend = "both")
    land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                        ,edgecolor="#1B1B1B",facecolor="#555555") 
    sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                    edgecolor='face',facecolor="#ADADAD")
    ax.add_feature(sea, zorder=0)
    contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                        linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
    ax.add_feature(land, zorder=2)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label(r'b', size = 14)
    plt.savefig(f"{fp}/figs/{file}/bmap_pos.png", transparent = True)
    plt.show()
# %%
