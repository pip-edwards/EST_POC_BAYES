#%%
#This is to check that the methodology doesn't have a large effect on the overall data
import pandas as pd
import seaborn as sns
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
fp = "C:/Users/pe1n24/OneDrive - University of Southampton/EST_POC_BAYES/"
os.chdir(fp)

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

#%%
#repeat for all files
for file in files:
    print(file)
    #get naming convention
    parts = f.split("_")
    met = parts[0]
    chl = parts[1]

    #import datasets
    betas = pd.read_csv(glob.glob(f"{fp}{f}/*_beta_vals.csv")[0])
    gammas = pd.read_csv(glob.glob(f"{fp}{f}/*_gamma_vals.csv")[0])
    sigmas = pd.read_csv(glob.glob(f"{fp}{f}/*_sigma_vals.csv")[0])
    mus = pd.read_csv(glob.glob(f"{fp}{f}/*_mu_vals.csv")[0])
    pocdf = pd.read_csv(f"input_data/merged_POC_{met}_{chl}.csv")


    #check with beta and gamma first

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

    newpocs = []
    zscores = []
    for i, r in pocdf.iterrows():
        m = beta1 + beta2*r["log_SST"] + beta3*r["log_Chla"] + beta4*r["log_Depth"] + beta5*r["log_SST"]*r["log_Chla"] + beta6*r["log_SST"]*r["log_Depth"] + beta7*r["log_Depth"]*r["log_Chla"]
        s = gamma1 + gamma2*r["log_SST"] + gamma3*r["log_Chla"] + gamma4*r["log_Depth"] + gamma5*r["log_SST"]*r["log_Chla"] + gamma6*r["log_SST"]*r["log_Depth"] + gamma7*r["log_Depth"]*r["log_Chla"]
        newpoc = m + (s**2)/2
        newpocs.append(newpoc)
        zscore = (newpoc-np.mean(pocdf["log_POC"]))/np.std(pocdf["log_POC"])
        zscores.append(zscore)
    pocdf["Modelled_POC"] = newpocs
    pocdf["Zscore"] = zscores

    sns.boxplot(data = pocdf, x = "Method", y = "Zscore")
    plt.xticks(rotation = 90)
    plt.axline((0,0), (1,0), ls = "--", color = "black")
    plt.title(file)
    plt.show()


    newpocs = []
    zscores = []
    for i, r in pocdf.iterrows():
        m = mus[f"mu.{i+1}"]
        s = sigmas[f"sigma.{i+1}"]
        newpoc = np.mean(m) + (np.mean(s)**2)/2
        newpocs.append(newpoc)
        zscore = (newpoc-np.mean(pocdf["log_POC"]))/np.std(pocdf["log_POC"])
        zscores.append(zscore)
    pocdf["Modelled_POC"] = newpocs
    pocdf["Zscore"] = zscores

    fig, ax = plt.subplots(dpi = 300)
    sns.boxplot(data = pocdf, x = "Method", y = "Zscore")
    plt.xticks(rotation = 90)
    plt.axline((0,0), (1,0), ls = "--", color = "black")
    plt.title(file)
    if not os.path.isdir(f"{fp}/figs/{file}"):
        os.makedirs(f"{fp}/figs/{file}")
    plt.savefig(f"{fp}/figs/{file}/method_zscore_boxplot.png")
    plt.show()
    # %%
