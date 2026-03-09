#%%

#Make the matches between the data frames

#import packages
import pandas as pd
import numpy as np
import xarray as xr
import datetime as dt
import os
import warnings
warnings.filterwarnings("ignore")

#%%

fp = "/iridisfs/scratch/pe1n24/EST_POC_BAYES/"
fp = "C:/Users/pe1n24/OneDrive - University of Southampton/EST_POC_BAYES/"

def temp_match(l):
    l = int(l)
    if l >= 0:
        l = l + 0.5
    else:
        l = l - 0.5
    return l

#load POC base database depending on which one is being used
pocdf = pd.read_csv(f"{fp}input_data/Global_POC_Database_2025-07-01.csv")

print(pocdf.columns)

#remove if it is on land
pocdf = pocdf[pocdf["on_land"] == False]

#format date for matching
pocdf["date"] = pd.to_datetime(pocdf["date_formatted"], format = "%Y-%m")

#the sst finishes in 2022 and the chla starts in 1997.
#bound the dataset to this.
pocdf = pocdf[pocdf["date"] >= dt.datetime(1997, 9, 1)]
pocdf = pocdf[pocdf["date"] < dt.datetime(2022, 7, 1)]

#set up lists to append data to
ssts = []
chlas = []

#inport SST data for matching
sst = xr.open_dataset(f"{fp}input_data/correct_lon_SST.nc")["SST"]

#slow loop for peace of mind it is working correctly
#for each row in the pocdf daatframe
for i, r in  pocdf.iterrows():
    print(i)
    #get lat and lon position
    lat = r["latitude"]
    lon = r["longitude"]

    #set the day as the middate of the month to match to SST
    sday = dt.datetime(int(r["year"]), int(r["month"]), 15)

    #set the latitude and longitude to match the middle of the grid cell
    if lon == 180: #does not read in properly if not
        lon = 179
    
    #make new lat and lon
    tlat = temp_match(lat)
    tlon = temp_match(lon)

    #extract the sst value for this datapoint
    lattemp = sst.sel(lat = tlat) #select the values at the latitude
    lontemp = lattemp.sel(lon = tlon) #from this select the values at the longitude
    daytemp = lontemp.sel(time = sday, method = "nearest") #from this select the values at the date
                                        #nearest is used because the date can vary from th 14-16th
    
    #add the mean of this to the sst list
    print(daytemp.shape())
    ssts.append(np.nanmean(daytemp.values))

    #open the correct chla file
    chla = xr.open_dataset(f"{fp}occci/occci_chla_{int(r['year'])}.nc")["chlor_a"]

    #set the date to match the chla days
    #chla is always the first of the month apart from the first month
    cday = dt.datetime(int(r["year"]), int(r["month"]), 1)
    
    #find mean oc-cci for the degree cell it is in
    latchla = chla.sel(lat = slice(lat+ 0.5, lat- 0.5))
    lonchla = latchla.sel(lon = slice(lon- 0.5, lon+ 0.5))
    daychla = lonchla.sel(time = cday, method = "nearest") #one month is on the fourth 

    #add the mean to the dataset
    chlas.append(np.nanmean(daychla.values))

#add to overall daatframe
pocdf["SST"] = ssts
pocdf["Chla"] = chlas

#select wanted columns and rename them
pocdf = pocdf[["date", "year", "month", "season", "latitude", "longitude", "New_category",
            "depth", "SST", "Chla", "poc_converted"]]

pocdf= pocdf.rename(columns= {"date":"Date", "year" :"Year", "month":"Month",
                            "season":"Season", "latitude":"Latitude", "longitude":"Longitude", 
                            "New_category":"Method", "depth":"Depth",  "poc_converted":"POC"})

#put bounds on the data
pocdf = pocdf[pocdf["Depth"] > 100]
pocdf = pocdf[pocdf["Chla"] >= 0.02]
pocdf = pocdf[pocdf["POC"] >= 0.1]

methodruns = ["all", "traps"]
methods = ["sediment trap","Thorium","Marine snow catcher"]

chlaruns = ["mgchla", "ugchla"]

for met in methodruns:
    for chl in chlaruns:

        #drop nan values
        poc = pocdf.dropna()

        #Make variables log

        #depth
        poc["log_Depth"] = np.log(poc["Depth"])
        #for a dataset where chla in in micrograms instead of miligrams
        if chl == "ugchla":
            poc["log_Chla"] = np.log(poc["Chla"]*1000)
        else:
            poc["log_Chla"] = np.log(poc["Chla"])

        #sst
        poc = poc[poc["SST"] > -1.79] #bound to make it all positive when log transformed (also below ice)
        poc["log_SST"] = np.log(poc["SST"] + 1.79)

        #POC    
        poc["log_POC"] = np.log(poc["POC"])
    
        #filter for only sediment traps and thorium
        if met == "traps":
            poc = poc[poc["Method"].isin(methods)]

        #print shape as a failsafe that it has run properly
        print(poc.shape)

        #save as a csv
        poc.to_csv(f"{fp}input_data/merged_POC_{met}_{chl}.csv", index = False)

