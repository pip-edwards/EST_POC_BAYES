#####

#setwd("")

library(ggplot2)   
library(fields)
library(MASS)   
library(dplyr) 
library(patchwork)
#####
#input data
allpoc <- read.csv("all_ugchla_run1_1103/all_ugchla_run1_1103_input_withygen.csv")

#make the three variables to plot against eachother
poc <- as.list(allpoc$log_POC)
ygen <-as.list(allpoc$ygen)
ygen0 <- as.list(allpoc$ygen0)

poc <- unlist(poc)
ygen <- unlist(ygen)
ygen0 <- unlist(ygen0)

#####
#for with sigma
#make a density plot for the map
dens <- kde2d(poc, ygen, n = 200)

#map these back to their original points and add to dataframe
dens_vals <- fields::interp.surface(obj = dens, loc = cbind(poc, ygen))
df <- data.frame(x = poc, y = ygen, density = dens_vals)

#plot on its own with colour changing for density
#ab line = 1:! line
#match the axis to that of the histograms
ggplot(df, aes(x = x, y = y, color = density)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black",
              linewidth  = 1) +
  scico::scale_color_scico(palette = "lipari") +
  theme_classic() +
  labs(color = "Density") +
  xlim(-2.6, 8)+ ylim(-4.15, 8.15) + 
  scale_x_continuous(name ="Actual ln(POC Flux)(mg/m²/day) ", 
                     breaks=c(-2,0,2,4,6,8)) +
  scale_y_continuous(name ="Modelled ln(POC Flux)(mg/m²/day) ", 
                     breaks=c(-4,-2,0,2,4,6,8))


#####

#repeat for without sigma
dens <- kde2d(poc, ygen0, n = 200)
dens_vals <- fields::interp.surface(obj = dens, loc = cbind(poc, ygen0))
df <- data.frame(x = poc, y = ygen0, density = dens_vals)

# Plot: points colored by density
ggplot(df, aes(x = x, y = y, color = density)) +
  geom_point(size = 1) +
  scico::scale_color_scico(palette = "lipari") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black",
              linewidth  = 1) +
  theme_classic() +
  labs(color = "Density")+
  xlim(-2.6, 8)+ ylim(-4.15, 8.15) + 
  scale_x_continuous(name ="Actual ln(POC Flux)(mg/m²/day) ", 
                     breaks=c(-2,0,2,4,6,8)) +
  scale_y_continuous(name ="Modelled ln(POC Flux)(mg/m²/day) ", 
                     breaks=c(-4,-2,0,2,4,6,8))

###

#now do this so they have a matching density function
dens1 <- kde2d(poc, ygen, n = 200)
dens_vals1 <- fields::interp.surface(obj = dens1, loc = cbind(poc, ygen))
df1 <- data.frame(x = poc, y = ygen, density = dens_vals1)

dens2 <- kde2d(poc, ygen0, n = 200)
dens_vals2 <- fields::interp.surface(obj = dens2, loc = cbind(poc, ygen0))
df2 <- data.frame(x = poc, y = ygen0, density = dens_vals2)

dens_min <- min(c(df1$density, df2$density), na.rm = TRUE)
dens_max <- max(c(df1$density, df2$density), na.rm = TRUE)

#repeat plot as above but this time set limits onto colour map
p1 <- ggplot(df1, aes(x = x, y = y, color = density)) +
    geom_point(size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black",
                linewidth  = 1) +
    scico::scale_color_scico(palette = "lipari", limits = c(dens_min, dens_max)) +
    theme_classic() +
    labs(color = "Density") +
    xlim(-2.6, 8)+ ylim(-4.15, 8.15) + 
    scale_x_continuous(name = "Observed ln(POC Flux)(mg/m²/day)", breaks = c(-2,0,2,4,6,8)) +
    scale_y_continuous(name = "Modelled ln(POC Flux)(mg/m²/day)", breaks = c(-4,-2,0,2,4,6,8))
p1
p2 <- ggplot(df2, aes(x = x, y = y, color = density)) +
    geom_point(size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black",
                linewidth  = 1) +
    scico::scale_color_scico(palette = "lipari", limits = c(dens_min, dens_max)) +
    theme_classic() +
    labs(color = "Density") +
    xlim(-2.6, 8)+ ylim(-4.15, 8.15) + 
    scale_x_continuous(name = "Observed ln(POC Flux)(mg/m²/day)", breaks = c(-2,0,2,4,6,8)) +
    scale_y_continuous(name = "Modelled ln(POC Flux)(mg/m²/day)", breaks = c(-4,-2,0,2,4,6,8))
p2
#plot together
p1 + p2

#save plot with same dimensions as histograms
ggsave(filename = "figs/scatters.png", width = 10, height = 4, units = "in", dpi =300)


