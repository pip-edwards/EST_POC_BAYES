#%%
#Plotting the caterpillars

#Packages
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
sns.set_context("paper")
cs = ["#00203C","#df4e38","#62a87c","#3AAAEB","#724e91"]

#set fp
fp = ""
os.chdir(fp)


#%%
#get usable files
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
for f in files:
    print(f)
    #open data
    betas = pd.read_csv(glob.glob(f"{fp}{f}/*_beta_vals.csv")[0])
    gammas = pd.read_csv(glob.glob(f"{fp}{f}/*_gamma_vals.csv")[0])

    #rename the columns
    betas["T"] = betas["beta.2"] 
    betas["C"] = betas["beta.3"] 
    betas["Z"] = betas["beta.4"] 
    betas["T*C"] = betas["beta.5"]
    betas["T*Z"] = betas["beta.6"] 
    betas["C*Z"] = betas["beta.7"]

    gammas["T"] = gammas["gamma.2"] 
    gammas["C"] = gammas["gamma.3"] 
    gammas["Z"] = gammas["gamma.4"] 
    gammas["T*C"] = gammas["gamma.5"]
    gammas["T*Z"] = gammas["gamma.6"] 
    gammas["C*Z"] = gammas["gamma.7"]

    
    #set up new data frames with the nicer names
    bbetas = betas[["T","C","Z","T*C","T*Z","C*Z"]]
    ggammas = gammas[["T","C","Z","T*C","T*Z","C*Z"]]
    #make a list of column names
    columns = ["T","C","Z","T*C","T*Z","C*Z"]

    #make plots of the betas and gammes on their own
    #(little test plots)
    fig,ax = plt.subplots(len(columns), 1, sharex = True)
    fig.suptitle(f)
    for i in range(len(columns)):
        c = columns[i]
        ax[i].hist(bbetas[c], label = c)
        ax[i].axline((0,0), (0,2000),  ls = "--", color = "black")
        ax[i].legend()
    ax[5].set_xlabel("Effect Size on μ")
    if not os.path.isdir(f"{fp}/figs/{f}"):
        os.makedirs(f"{fp}/figs/{f}")
    plt.savefig(f"{fp}/figs/{f}/beta_ES_testplot.png")
    plt.show()

    fig,ax = plt.subplots(len(columns), 1, sharex = True)
    fig.suptitle(f)
    for i in range(len(columns)):
        c = columns[i]
        ax[i].hist(ggammas[c], label = c)
        ax[i].axline((0,0), (0,2000),  ls = "--", color = "black")
        ax[i].legend()
    ax[5].set_xlabel("Effect Size on σ")
    plt.savefig(f"{fp}/figs/{f}/gamma_ES_testplot.png")
    plt.show()

#%%

# %%
#Plot summative effect sizes and the complicated caterpillar plot
#Transform the data by multiply by SD

#COLOUR PALETTE
cs = ["#00203C","#df4e38","#62a87c","#3AAAEB","#724e91"]
cs1 = ["#00203C","#df4e38","#62a87c","#3AAAEB","#724e91"]#["#133859","#000000",  "#62a87c","#000000","#000000",]
aas = [0.7, 0.4, 0.4, 0.4, 0.4]

#THIS IS FOR BETA
for f in files:
    print(f)

    #split to get the file naming convention
    parts = f.split("_")
    met = parts[0]
    chl = parts[1]

#REMOVED FROM LOOP AS NOW ONLY WORKING WITH ONE FILE

#load in betas
betas = pd.read_csv(glob.glob(f"{fp}{f}/*_beta_vals.csv")[0])

#rename the columns
betas["T"] = betas["beta.2"] 
betas["C"] = betas["beta.3"] 
betas["Z"] = betas["beta.4"] 
betas["T*C"] = betas["beta.5"]
betas["T*Z"] = betas["beta.6"] 
betas["C*Z"] = betas["beta.7"]

#set up new data frames with the nicer names
bbetas = betas[["T","C","Z","T*C","T*Z","C*Z"]]
ggammas = gammas[["T","C","Z","T*C","T*Z","C*Z"]]
#make a list of column names
columns = ["T","C","Z","T*C","T*Z","C*Z"]

#load in original data
pocdf = pd.read_csv(f"{fp}input_data/merged_POC_{met}_{chl}_1103.csv")
#set arrays of predictor variables 
sst = np.array(pocdf["log_SST"])
chla = np.array(pocdf["log_Chla"])
depth = np.array(pocdf["log_Depth"])

q = "multiply" #this is how we transform by SD

#make datasets
call =  (bbetas["C"] + np.median(np.multiply.outer(np.array(bbetas["T*C"]),sst), axis = 1) + np.median(np.multiply.outer(np.array(bbetas["C*Z"]),depth), axis = 1))*np.std(pocdf["log_Chla"])
cnos = (bbetas["C"] + np.median(np.multiply.outer(np.array(bbetas["C*Z"]),depth), axis = 1))*np.std(pocdf["log_Chla"])
cnod = (bbetas["C"] + np.median(np.multiply.outer(np.array(bbetas["T*C"]),sst), axis = 1))*np.std(pocdf["log_Chla"])
cnone = np.array(betas["beta.3"]*np.std(pocdf["log_Chla"]))
chlas = pd.DataFrame({"All terms": call,
                    "No SST crossterm": cnos,
                    "No Chl-a crossterm" : np.ones(8000)*np.nan,
                    "No Depth crossterm": cnod})#,
                    #"No crossterms": cnone})

sall = (bbetas["T"] + np.mean(np.multiply.outer(np.array(bbetas["T*C"]),chla), axis = 1) + np.mean(np.multiply.outer(np.array(bbetas["T*Z"]),depth), axis = 1))*np.std(pocdf["log_SST"])
snoc =(bbetas["T"]  + np.mean(np.multiply.outer(np.array(bbetas["T*Z"]),depth), axis = 1))*np.std(pocdf["log_SST"])
snod =(bbetas["T"] + np.mean(np.multiply.outer(np.array(bbetas["T*C"]),chla), axis = 1))*np.std(pocdf["log_SST"])
snone = np.array(betas["beta.2"]*np.std(pocdf["log_SST"]))
ssts = pd.DataFrame({"All terms": sall,
                    "No SST crossterm" : np.ones(8000)*np.nan,
                    "No Chl-a crossterm": snoc,    
                    "No Depth crossterm": snod})#,
                    # "v": snone})

dall = (bbetas["Z"] +  np.mean(np.multiply.outer(np.array(bbetas["C*Z"]),chla), axis = 1) +  np.mean(np.multiply.outer(np.array(bbetas["T*Z"]),sst),axis = 1))*np.std(pocdf["log_Depth"])
dnos = (bbetas["Z"] +  np.mean(np.multiply.outer(np.array(bbetas["C*Z"]),chla), axis = 1) )*np.std(pocdf["log_Depth"])
dnoc = (bbetas["Z"] +  np.mean(np.multiply.outer(np.array(bbetas["T*Z"]),sst),axis = 1))*np.std(pocdf["log_Depth"])
dnone = np.array(betas["beta.4"]*np.std(pocdf["log_Depth"]))
depths = pd.DataFrame({"All terms": dall,
                    "No SST crossterm": dnos,
                    "No Chl-a crossterm": dnoc,
                    "No Depth crossterm" : [100]*4000 + [200]*4000})#, #this is different beause i need the term in the legend
                    #"No crossterms": dnone})

#summarise all betas
balls = pd.DataFrame({"T":sall,"C":call, "Z":dall})

#load in gammas
gammas = pd.read_csv(glob.glob(f"{fp}{f}/*_gamma_vals.csv")[0])

#rename the columns
gammas["T"] = gammas["gamma.2"] 
gammas["C"] = gammas["gamma.3"] 
gammas["Z"] = gammas["gamma.4"] 
gammas["T*C"] = gammas["gamma.5"]
gammas["T*Z"] = gammas["gamma.6"] 
gammas["C*Z"] = gammas["gamma.7"]

#set up new data frames with the nicer names
bgammas = gammas[["T","C","Z","T*C","T*Z","C*Z"]]
ggammas = gammas[["T","C","Z","T*C","T*Z","C*Z"]]
#make a list of column names
columns = ["T","C","Z","T*C","T*Z","C*Z"]

#load in original data
pocdf = pd.read_csv(f"{fp}input_data/merged_POC_{met}_{chl}_1103.csv")
#set arrays of predictor variables 
sst = np.array(pocdf["log_SST"])
chla = np.array(pocdf["log_Chla"])
depth = np.array(pocdf["log_Depth"])

q = "multiply"
call =  (ggammas["C"] + np.median(np.multiply.outer(np.array(ggammas["T*C"]),sst), axis = 1) + np.median(np.multiply.outer(np.array(ggammas["C*Z"]),depth), axis = 1))*np.std(pocdf["log_Chla"])
cnos = (ggammas["C"] + np.median(np.multiply.outer(np.array(ggammas["C*Z"]),depth), axis = 1))*np.std(pocdf["log_Chla"])
cnod = (ggammas["C"] + np.median(np.multiply.outer(np.array(ggammas["T*C"]),sst), axis = 1))*np.std(pocdf["log_Chla"])
cnone = np.array(gammas["gamma.3"]*np.std(pocdf["log_Chla"]))
chlas = pd.DataFrame({"All terms": call,
                    "No SST crossterm": cnos,
                    "No Chl-a crossterm" : np.ones(8000)*np.nan,
                    "No Depth crossterm": cnod})#,
                    #"No crossterms": cnone})

sall = (ggammas["T"] + np.mean(np.multiply.outer(np.array(ggammas["T*C"]),chla), axis = 1) + np.mean(np.multiply.outer(np.array(ggammas["T*Z"]),depth), axis = 1))*np.std(pocdf["log_SST"])
snoc =(ggammas["T"]  + np.mean(np.multiply.outer(np.array(ggammas["T*Z"]),depth), axis = 1))*np.std(pocdf["log_SST"])
snod =(ggammas["T"] + np.mean(np.multiply.outer(np.array(ggammas["T*C"]),chla), axis = 1))*np.std(pocdf["log_SST"])
snone = np.array(gammas["gamma.2"]*np.std(pocdf["log_SST"]))
ssts = pd.DataFrame({"All terms": sall,
                    "No SST crossterm" : np.ones(8000)*np.nan,
                    "No Chl-a crossterm": snoc,    
                    "No Depth crossterm": snod})#,
                    #"No crossterms": snone})

dall = (ggammas["Z"] +  np.mean(np.multiply.outer(np.array(ggammas["C*Z"]),chla), axis = 1) +  np.mean(np.multiply.outer(np.array(ggammas["T*Z"]),sst),axis = 1))*np.std(pocdf["log_Depth"])
dnos = (ggammas["Z"] +  np.mean(np.multiply.outer(np.array(ggammas["C*Z"]),chla), axis = 1) )*np.std(pocdf["log_Depth"])
dnoc = (ggammas["Z"] +  np.mean(np.multiply.outer(np.array(ggammas["T*Z"]),sst),axis = 1))*np.std(pocdf["log_Depth"])
dnone = np.array(gammas["gamma.4"]*np.std(pocdf["log_Depth"]))
depths = pd.DataFrame({"All terms": dall,
                    "No SST crossterm": dnos,
                    "No Chl-a crossterm": dnoc,
                    "No Depth crossterm" : [100]*4000 + [200]*4000})#,
                    #"No crossterms": dnone})

galls = pd.DataFrame({"T":sall,"C":call, "Z":dall})

print(np.mean(sall),np.mean(call),np.mean(dall),)
print(np.std(sall),np.std(call),np.std(dall),)
print(np.mean(dall), np.mean(dnos), np.mean(dnoc))

#%%
betas = pd.read_csv(glob.glob(f"{fp}{f}/*_beta_vals.csv")[0])
gammas = pd.read_csv(glob.glob(f"{fp}{f}/*_gamma_vals.csv")[0])

#rename the columns
betas["T"] = betas["beta.2"]*np.std(pocdf["log_SST"]) 
betas["C"] = betas["beta.3"] *np.std(pocdf["log_Chla"]) 
betas["Z"] = betas["beta.4"] *np.std(pocdf["log_Depth"]) 
betas["T*C"] = betas["beta.5"]*np.std(pocdf["log_SST"]*pocdf["log_Chla"]) 
betas["T*Z"] = betas["beta.6"] *np.std(pocdf["log_SST"]*pocdf["log_Depth"]) 
betas["C*Z"] = betas["beta.7"]*np.std(pocdf["log_Depth"]*pocdf["log_Chla"]) 

gammas["T"] = gammas["gamma.2"]*np.std(pocdf["log_SST"]) 
gammas["C"] = gammas["gamma.3"] *np.std(pocdf["log_Chla"]) 
gammas["Z"] = gammas["gamma.4"] *np.std(pocdf["log_Depth"]) 
gammas["T*C"] = gammas["gamma.5"]*np.std(pocdf["log_SST"]*pocdf["log_Chla"]) 
gammas["T*Z"] = gammas["gamma.6"] *np.std(pocdf["log_SST"]*pocdf["log_Depth"]) 
gammas["C*Z"] = gammas["gamma.7"]*np.std(pocdf["log_Depth"]*pocdf["log_Chla"]) 
    
#set up new data frames with the nicer names
bbetas = betas[["T","C","Z","T*C","T*Z","C*Z"]]
ggammas = gammas[["T","C","Z","T*C","T*Z","C*Z"]]
#%%

#make plots
#do in two halves for ease
cs1 = ["#ff7433","#62a87c","#3969d1"]
cs2 = ["#ffbe33","#9ad4b0","#8CC1E0","#cabf59","#beabf3","#CC89A8"]
yticks1 = [r"$\beta_{T}$",r"$\beta_{C}$",r"$\beta_{Z}$"]
yticks2 = [r"$\beta_{T0}$",r"$\beta_{C0}$",r"$\beta_{Z0}$",r"$\beta_{TC}$",r"$\beta_{TZ}$",r"$\beta_{CZ}$",]

fig,ax = plt.subplots(2,1, figsize = (6,7), dpi = 300)
sns.violinplot(data = balls, palette = cs1, inner = None, orient = "h", gap = -0.3,
               linewidth = 0.8, linecolor = "#333333", alpha = 0.8, ax = ax[0])
sns.violinplot(data = bbetas, palette = cs2, inner = None, orient = "h", gap = -0.3,
               linewidth = 0.8, linecolor = "#333333", alpha = 0.8)

ax[0].scatter(np.median(balls, axis = 0),np.arange(0,3,1), c= "#333333", s = 10, marker  = "|", alpha = 0.8)
ax[1].scatter(np.median(bbetas, axis = 0),np.arange(0,6,1), c= "#333333", s = 10, marker  = "|", alpha = 0.8)

for i in range(2):
    ax[i].axline((0,0), (0,1),  ls = "--", color = "black")
    ax[i].set_xlabel("Effect Size on μ")
ax[0].set_yticks(np.arange(0,3,1), yticks1)
ax[1].set_yticks(np.arange(0,6,1), yticks2)

ax[0].set_xticks(np.arange(-0.8, 0.9,0.2))
ax[1].set_xlim([-1.45,1.25])
ax[1].set_xticks(np.arange(-1.4, 1.3,0.2))
plt.savefig(f"{fp}/figs/{f}/effect_sizes_beta.png")#, transparent = True)

cs1 = ["#ff7433","#62a87c","#3969d1"]
cs2 = ["#ffbe33","#9ad4b0","#8CC1E0","#cabf59","#beabf3","#CC89A8"]
yticks3 = [r"$\gamma_{T}$",r"$\gamma_{C}$",r"$\gamma_{Z}$"]
yticks4 = [r"$\gamma_{T0}$",r"$\gamma_{C0}$",r"$\gamma_{Z0}$",r"$\gamma_{TC}$",r"$\gamma_{TZ}$",r"$\gamma_{CZ}$",]

fig,ax = plt.subplots(2,1, figsize = (6,7), dpi = 300)
sns.violinplot(data = galls, palette = cs1, inner = None, orient = "h", gap = -0.3,
               linewidth = 0.8, linecolor = "#333333", alpha = 0.8, ax = ax[0])

sns.violinplot(data = ggammas, palette = cs2, inner = None, orient = "h",gap = -0.3,
               linewidth = 0.8, linecolor = "#333333", alpha = 0.8)#,  

ax[0].scatter(np.median(galls, axis = 0),np.arange(0,3,1), c= "#333333", s = 10, marker  = "|", alpha = 0.8)
ax[1].scatter(np.median(ggammas, axis = 0),np.arange(0,6,1), c= "#333333", s = 10, marker  = "|", alpha = 0.8)

for i in range(2):
    ax[i].axline((0,0), (0,1),  ls = "--", color = "black")
    ax[i].set_xlabel("Effect Size on σ")
ax[0].set_yticks(np.arange(0,3,1), yticks3)
ax[1].set_yticks(np.arange(0,6,1), yticks4)

ax[0].set_xticks(np.arange(-0.2, 0.06,0.05))
ax[1].set_xlim([-1.2,0.8])
ax[1].set_xticks(np.arange(-1.2, 0.85,0.2))

plt.savefig(f"{fp}/figs/{f}/effect_sizes_gamma.png")
