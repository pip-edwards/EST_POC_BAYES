Code to run analysis and methodology for paper submitted 2026

Running order
1. make_climatology.py -> this makes climatologies for SST and Chl-a aswell as depth arrays
2. match_poc.py -> matches POC flux dataset to SST and Chl-a
3. stan_model_all.py -> runs the stan model on the match data (has loops as may potentially fail out/ have high Rhat scores)
4. site_map.py -> make supplementary figure S1
5. histogramplots.py -> make histograms of F2
6. scatter.r -> make scatter plots of F2
7. effect_sizes.py -> plot the histograms of effect sizes F3
8. poc_figures.py -> make poc maps F4
