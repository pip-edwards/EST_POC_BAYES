#%%
#This makes the histograms in figure 2 :D
#the scatter plots are made in R

#Packages
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import glob
from scipy.stats import ks_2samp
from scipy.stats import pearsonr
sns.set_context("paper")
from sklearn.metrics import r2_score
#colour scheme :)
cs = ["#e99871", "#416c99"]

#%%
#filepaths
fp = ""

#make list of available files that have an rhat score close to one
os.chdir(fp)
files = []
for f in os.listdir(fp):
    if "1103" in f:
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

#%%
#make an empty dataframe to save the stats into 
comparison_stats = pd.DataFrame({"File":[],
                                 "r^2 sigma":[],
                                 "r2 no sigma":[],
                                 "r sigma":[],
                                 "r pval sigma":[],
                                 "r no sigma":[],
                                 "r pval no sigma":[],
                                 "ks sigma":[],
                                 "ks no sigma":[]})

for f in files:
    #split to get input data
    parts = f.split("_")
    met = parts[0]
    chl = parts[1]

    #### HAsH OUT THE FOLLOWING IF ALREADY GENERATED
    #(this is to save time in rerunning so mu and sigma don't have to be loaded each time)
    #get input data 

    # pocdf = pd.read_csv(f"{fp}input_data/merged_POC_{met}_{chl}_1103.csv")
    # sigmas = pd.read_csv(glob.glob(f"{fp}{f}/*_sigma_vals.csv")[0])
    # mus = pd.read_csv(glob.glob(f"{fp}{f}/*_mu_vals.csv")[0])
    ##SET UP GENERATED LOG POC DATASETS
    # mmus = mus.mean(axis = 0).to_list()
    # msigmas = sigmas.mean(axis = 0).to_list()
    # ygen = [] #with sigma
    # ygen0 = [] #without sigma
    # for x in range(len(mmus)):
    #     y = np.random.normal(mmus[x], msigmas[x])
    #     ygen.append(y)
    #     y0 = np.random.normal(mmus[x], 0)
    #     ygen0.append(y0)
    # pocdf["ygen"] = ygen
    # pocdf["ygen0"] = ygen0
    # print(np.min(ygen), np.max(ygen))
    # pocdf.to_csv(f"{fp}{f}/{f}_input_withygen.csv")
    ####

    #load in the data 
    pocdf = pd.read_csv(f"{fp}{f}/{f}_input_withygen.csv")
    ygen = pocdf["ygen"]
    ygen0 = pocdf["ygen0"]

    #make histogram
    #this is where they don't share axis
    fig, ax = plt.subplots(1,2,figsize = (9,4), sharey= True, sharex = True, dpi = 300)
    plt.tight_layout()

    #without sigma
    sns.histplot(x = pocdf["log_POC"], alpha = 0.5, fill = True, linewidth  = 1.5,kde = True,
                color = "black", edgecolor ="black",label = "Observed POC Flux",  binwidth =0.2, ax = ax[1])
    sns.histplot(x = ygen0, alpha = 0.5, fill = True, edgecolor = cs[1], linewidth  = 1.5, kde = True,
                color = cs[1], label = "Modelled POC, no σ", binwidth =0.2, ax = ax[1])

    #with sigma
    sns.histplot(x = pocdf["log_POC"], alpha = 0.5, fill = True, linewidth  = 1.5,kde = True,
                color = "black", edgecolor ="black",label = "Observed POC Flux",  binwidth =0.2, ax = ax[0])
    sns.histplot(x = ygen, alpha = 0.5, fill = True, edgecolor = cs[0], linewidth  = 1.5,kde = True,
                color = cs[0], label = "Modelled POC with σ", binwidth =0.2, ax = ax[0])

    ks,p1 = ks_2samp(pocdf["log_POC"],ygen)
    ks0,p0 = ks_2samp(pocdf["log_POC"],ygen0)
    #ax[0].title(f"KS with sigma = {round(ks,3)}, KS without sigma = {round(ks0,3)}")

    #set axis ticks/labels
    for i in range(2):
        ax[i].set_xlabel("ln(POC Flux) (mg/m²/day)")
        ax[i].legend(loc = "upper left")
    #make a folder if not already exisiting and save file
    if not os.path.isdir(f"{fp}/figs/{f}"):
        os.makedirs(f"{fp}/figs/{f}")
    plt.savefig(f"{fp}/figs/{f}/POC_histograms")
    plt.show()

    #print r2, ks and pcc
    print(f)
    print("r2 yegn ygen0", r2_score(pocdf["log_POC"], pocdf["ygen"]), r2_score(pocdf["log_POC"], pocdf["ygen0"])) 
    print("pcc: ygen", pearsonr(pocdf["log_POC"], pocdf["ygen"]), "ygen0",pearsonr(pocdf["log_POC"], pocdf["ygen0"]))
    print("ks ygen ygen0", ks, ks0)

    r0, p0 = pearsonr(pocdf["log_POC"], pocdf["ygen0"])
    r1, p1 = pearsonr(pocdf["log_POC"], pocdf["ygen"])

    comparison_stats.loc[len(comparison_stats)] = [f,
                                                   r2_score(pocdf["log_POC"], pocdf["ygen"]), 
                                                   r2_score(pocdf["log_POC"], pocdf["ygen0"]),
                                                   r1, p1, r0, p0, ks, ks0]

comparison_stats.to_csv("histogram_comparison_stats_1103runs.csv", index = False)

#make the outlines for the scatter graphs so thye match style of python
#only needs to be made once
fig, ax = plt.subplots(figsize = (4,4), dpi = 300)
sns.scatterplot(x = [0], y = [0], alpha = 0)
plt.yticks([])
plt.xticks([])
plt.savefig(f"{fp}/figs/scatterbackground.png", transparent = True)



# %%
