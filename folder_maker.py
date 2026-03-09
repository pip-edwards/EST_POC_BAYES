#%%
import glob
import os
import shutil

for file in glob.glob("*summary_stats.csv"):
    print(file[:-18])
    f = file[:-18] + "_"
    if not os.path.isdir(f):
        os.mkdir(f)
    f_files = glob.glob(f"{f}*.csv")
    for ff in f_files:
        shutil.move(ff, f"{f}/{ff}")
