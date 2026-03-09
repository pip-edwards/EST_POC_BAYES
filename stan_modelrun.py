#%%

#run the stan datasets
import pandas as pd
import numpy as np
import cmdstanpy
import glob
import os

fp = "/iridisfs/scratch/pe1n24/EST_POC_BAYES/"
#fp = "C:/Users/pe1n24/OneDrive - University of Southampton/EST_POC_BAYES/"

#import stan file
mod = cmdstanpy.CmdStanModel(stan_file= f"{fp}stan_noy.stan")

methodruns = ["all", "traps"]
chlaruns = ["mgchla", "ugchla"]

#main_datas = ["run11","run12", "run13", "run14", "run15","run16","run17","run18","run19", "run10"]
runs = ["run1","run2", "run3", "run4", "run5","run6","run7","run8","run9", "run0"]

for r in runs:
    for m in methodruns:
        for c in chlaruns:
            #set filename
            filename = m + "_" + c + "_" + r

            poc = pd.read_csv(f"{fp}input_data/merged_POC_{m}_{c}.csv")

            poc["TxC"] = poc["log_SST"] * poc["log_Chla"]
            poc["CxD"] = poc["log_Depth"] * poc["log_Chla"]
            poc["TxD"] = poc["log_SST"] * poc["log_Depth"]

            #set up dataframe
            xvars = ["log_SST", "log_Chla", "log_Depth", "TxC","TxD","CxD"]
            x = np.column_stack([np.ones(poc.shape[0]), poc["log_SST"],
                                 poc["log_Chla"], poc["log_Depth"],
                                 poc["TxC"], poc["TxD"], poc["CxD"]])
            data = {"N": poc.shape[0],  #number of observations
                "p": x.shape[1],    #number of x variables
                "y": poc["log_POC"],#dependent variable
                "x": x}             #independent variables
            
            #run the mcmc
            mcmc = None
            try:
                mcmc = mod.sample(data = data, chains = 4, iter_sampling=2000, parallel_chains= 4)#, outputdir = "/iridisfs/scratch/pe1n24/run1310/run1612/stan_output")
            except Exception as e:
                with open(f"{fp}failed_terms.log", "a") as logf:
                    logf.write(f"{filename}: {e}\n")

            if mcmc is not None:        
                sumstats =  mcmc.summary()
                sumstats.to_csv(f"{fp}{filename}_summary_stats.csv", index = False)
                
                #print("saving separate draws...")
                os.mkdir(f"{fp}stan_output/{filename}")
                os.mkdir(f"{fp}{filename}/")
                mcmc.save_csvfiles(f"{fp}stan_output/{filename}")
                
                #for the draw files:
                #add them all into one list
                df_files = []
                for file in glob.glob(f"{fp}stan_output/{filename}/*.csv"):
                        df_files.append(file)
                        #print(file)
                
                
                #for these files, compile them into dataframes of mu, sigma, beta and gamma
                for file in df_files:
                    #print(file)
                
                    #this is so that in concats easy 
                    #(I have not looked into a diff quicker way to do this yet as this works fine for now)
                    if file == df_files[0]:
                        headers = pd.read_csv(file, skiprows=47, nrows = 1).columns
                        df = pd.read_csv(file, skiprows=52, skipfooter=5, names=headers) #this is where the data starts and ends in the files
                        betas = df.filter(like = "beta")
                        gammas = df.filter(like = "gamma")
                        mus = df.filter(like = "mu")
                        sigmas = df.filter(like = "sigma")
                        
                    else:
                        df = pd.read_csv(file, skiprows=52, skipfooter=5, names=headers)
                        mus1 = df.filter(like = "mu")
                        sigmas1 = df.filter(like = "sigma")
                        betas1 = df.filter(like = "beta")
                        gammas1 = df.filter(like = "gamma")
                       
                
                        mus = pd.concat([mus, mus1])
                        sigmas = pd.concat([sigmas, sigmas1])
                        betas = pd.concat([betas, betas1])
                        gammas = pd.concat([gammas, gammas1])
                        
                #save as csv files
                mus.to_csv(f"{fp}{filename}/{filename}_mu_vals.csv", index = False)
                sigmas.to_csv(f"{fp}{filename}/{filename}_sigma_vals.csv", index = False)
                gammas.to_csv(f"{fp}{filename}/{filename}_gamma_vals.csv", index = False)
                betas.to_csv(f"{fp}{filename}/{filename}_beta_vals.csv", index = False)     


            




        


