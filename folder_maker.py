#%%
import glob
import os
import shutil

for file in glob.glob("*summary_stats.csv"):
    print(file[:-18])
    f = file[:-18] 
    fp = f"{f}_0903"
    if not os.path.isdir(f) or os.path.isdir(fp):
        os.mkdir(fp)
        f_files = glob.glob(f"{f}*.csv")
        for ff in f_files:
            shutil.move(ff, f"{fp}/{ff}")
