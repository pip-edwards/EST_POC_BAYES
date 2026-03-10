#%%
#Plotting the caterpillars

#Packages
sns.set_context("paper")

#Check your sticky notes gal!

#%%
#get usable files
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