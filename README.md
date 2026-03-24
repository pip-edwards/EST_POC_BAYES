Code to run analysis and methodology for paper submitted 2026

Running order
make_climatology.py -> this makes climatologies for SST and Chl-a aswell as depth arrays
match_poc.py -> matches POC flux dataset to SST and Chl-a
stan_model_all.py -> runs the stan model on the match data (has loops as may potentially fail out/ have high Rhat scores)
site_map.py -> make supplementary figure S1
histogramplots.py -> make histograms of F2
scatter.r -> make scatter plots of F2
effect_sizes.py -> plot the histograms of effect sizes
poc_figures.py -> make poc maps
