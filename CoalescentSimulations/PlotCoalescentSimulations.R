rm( list = ls())

library(reshape2)
library(ggplot2)
library(scales)
library(ggpubr)

fst <- read.csv('~/UBC/LinkageCoalescent/CoalescentSimulations/island_1000sims.csv') ## You could swap this out for any of the demographic models you care to examine
#fst <- read.csv('~/UBC/LinkageCoalescent/CoalescentSimulations/RING_1000sims.csv')
#fst <- read.csv('~/UBC/LinkageCoalescent/CoalescentSimulations/CLINE_1000sims.csv')
#fst <- read.csv('~/UBC/LinkageCoalescent/CoalescentSimulations/2PATCH_1000sims.csv')


mean_fst_hist <- ggplot(data = fst, aes( x = weighted_Fst, col = as.factor(rec), fill = as.factor(rec))) + 
  geom_histogram( alpha = 0.65, lwd = 0, bins = 75, position = 'identity')+
  scale_x_continuous( expression(italic(F["ST - Weir and Cockerham"])) ) +
  scale_fill_brewer("Recombination\nRate\n(M/bp)", palette = "Set1")+
  scale_color_brewer("Recombination\nRate\n(M/bp)",palette = "Set1")+
  scale_y_continuous("Frequency")+
  theme_bw()

mean_fst_dens <- ggplot(data = fst, aes( x = weighted_Fst, col = as.factor(rec), fill = as.factor(rec))) + 
  geom_density(adjust = 2, alpha = 0.65, lwd = 0, position = 'identity')+
  scale_x_continuous( expression(italic(F["ST - Weir and Cockerham"])) ) +
  scale_fill_brewer("Recombination\nRate\n(M/bp)", palette = "Set1")+
  scale_color_brewer("Recombination\nRate\n(M/bp)",palette = "Set1")+
  scale_y_continuous("Frequency")+
  theme_bw()

se_fst_hist <- ggplot(data = fst, aes( x = mean_SE, col = as.factor(rec), fill = as.factor(rec))) + 
  scale_x_continuous( expression( "Standard Error ["* italic(F["ST - Weir and Cockerham"])*"] - Jackknife") ) +
  geom_histogram( alpha = 0.65, lwd = 0, bins = 75, position = 'identity')+
  scale_fill_brewer("Recombination\nRate\n(M/bp)", palette = "Set1")+
  scale_color_brewer("Recombination\nRate\n(M/bp)",palette = "Set1")+
  scale_y_continuous("Frequency")+
  theme_bw()

se_fst_dens <- ggplot(data = fst, aes( x = mean_SE, col = as.factor(rec), fill = as.factor(rec))) + 
  scale_x_continuous( expression( "Standard Error ["* italic(F["ST - Weir and Cockerham"])*"] - Jackknife") , limits = c(0,0.03)) +
  geom_density(adjust = 2, alpha = 0.65, lwd = 0, position = 'identity')+
  scale_fill_brewer("Recombination\nRate\n(M/bp)", palette = "Set1")+
  scale_color_brewer("Recombination\nRate\n(M/bp)",palette = "Set1")+
  scale_y_continuous("Frequency")+
  theme_bw()

true_se_fst_dens <- ggplot(data = fst, aes( x = sqrt(variance), col = as.factor(rec), fill = as.factor(rec))) + 
  scale_x_continuous( expression( "Standard Error ["* italic(F["ST - Weir and Cockerham"])*"] - True") , limits = c(0,0.03)) +
  geom_density(adjust = 2, alpha = 0.65, lwd = 0, position = 'identity')+
  scale_fill_brewer("Recombination\nRate\n(M/bp)", palette = "Set1")+
  scale_color_brewer("Recombination\nRate\n(M/bp)",palette = "Set1")+
  scale_y_continuous("Frequency")+
  theme_bw()

trees_hist <- ggplot(data = fst, aes( x = num_trees, col = as.factor(rec), fill = as.factor(rec))) + 
  geom_histogram( alpha = 0.65, lwd = 0, bins = 75, position = 'identity')+
  scale_fill_brewer("Recombination\nRate\n(M/bp)", palette = "Set1")+
  scale_color_brewer("Recombination\nRate\n(M/bp)",palette = "Set1")+
  scale_x_log10( expression( "Number of Trees") ) +
  scale_y_continuous("Frequency")+
  theme_bw()

trees_dens <- ggplot(data = fst, aes( x = num_trees, col = as.factor(rec), fill = as.factor(rec))) + 
  geom_density(adjust = 2, alpha = 0.65, lwd = 0, position = 'identity')+
  scale_fill_brewer("Recombination\nRate\n(M/bp)", palette = "Set1")+
  scale_color_brewer("Recombination\nRate\n(M/bp)",palette = "Set1")+
  scale_x_log10( expression( "Number of Trees") ) +
  scale_y_continuous("Density")+
  theme_bw()



ggarrange(mean_fst_hist, se_fst_hist, trees_hist, ncol = 1, nrow = 3, common.legend = T, legend = 'right', labels = 'AUTO')

ggarrange(mean_fst_dens, se_fst_dens, trees_dens, ncol = 1, nrow = 3, common.legend = T, legend = 'right', labels = 'AUTO')

ggarrange(mean_fst_dens, se_fst_dens, trees_hist, ncol = 1, nrow = 3, common.legend = T, legend = 'right', labels = 'AUTO', align = 'v')

ggarrange(se_fst_dens, true_se_fst_dens, ncol = 1, nrow = 2, common.legend = T, legend = 'right', labels = 'AUTO', align = 'v')


unshuffled_comp<-ggplot(data = fst, aes(y = mean_SE, x = variance, col = as.factor(rec), fill = as.factor(rec)))+
  geom_abline(aes(intercept = 0, slope = 1), lty = 2, col = 'black')+
  geom_point(alpha = 0.5)+
  geom_smooth(method="lm") +
  scale_x_continuous("True Error", limits = c(0,0.02))+
  scale_y_continuous("Estimated Error (Block Jackknife - Unshuffled)", limits = c(0,0.03))+
  scale_color_brewer("Recombination\nRate\n(M/bp)",palette = "Set1")+
  scale_fill_brewer("Recombination\nRate\n(M/bp)",palette = "Set1")+
  theme_bw()

shuffled_comp<-ggplot(data = fst, aes(y = mean_SE_shuf, x = variance, col = as.factor(rec), fill = as.factor(rec)))+
  geom_abline(aes(intercept = 0, slope = 1), lty = 2, col = 'black')+
  geom_point(alpha = 0.5)+
  geom_smooth(method="lm") +
  scale_x_continuous("True Error", limits = c(0,0.02))+
  scale_y_continuous("Estimated Error (Block Jackknife - Shuffled)", limits = c(0,0.03))+
  scale_color_brewer("Recombination\nRate\n(M/bp)",palette = "Set1")+
  scale_fill_brewer("Recombination\nRate\n(M/bp)",palette = "Set1")+
  theme_bw()

ggarrange(unshuffled_comp, shuffled_comp, ncol = 1, nrow = 2, common.legend = T, legend = 'right', labels = 'AUTO', align = 'v')


fst_means<-tapply(fst$weighted_Fst, fst$rec, mean)

FstMeanDF <- data.frame( weightedFst = fst_means, rec = names(fst_means))
## Figure for the paper
cbPalette <- c( "#009E73", "#0072B2", "#D55E00", "#CC79A7")

pal<- c('#a6611a','#dfc27d','#80cdc1','#018571')


mean_fst_dens <- ggplot(data = fst, aes( x = weighted_Fst, fill = as.factor(rec), col = as.factor(rec))) + 
  geom_vline(data = FstMeanDF, aes(xintercept = weightedFst, col = as.factor(rec)), lwd = 1.5 )+
  geom_vline(aes(xintercept = 1/9), lwd = 1, lty = 2, alpha = 0.5)+
  geom_line(adjust = 1, stat = 'density',lwd = 1, position = 'identity')+
  
#    geom_line(adjust = 1, stat = 'density',lwd = 1, position = 'identity')+
   # geom_density(adjust = 2, lwd = 1, alpha = 0.15,  position = 'identity')+
  scale_x_continuous( expression(italic(hat(F)["ST - Weir and Cockerham"])) ) +
  scale_fill_manual("Recombination\nRate (M/bp)", values = pal)+
  scale_color_manual("Recombination\nRate (M/bp)",values = pal)+
  scale_y_continuous("Density", expand = c(0, 0), limits = c(0,10))+
  theme_bw()+
  theme(
    axis.text.x  = element_text(colour = 'black'),
    axis.text.y  = element_text(colour = 'black'),
    #    axis.text.y  = element_blank(),
#    axis.title.y  = element_blank(),
#    axis.ticks.y  = element_blank(),
    legend.position = c(0.8, 0.50),
    legend.text.align = 0,
#    legend.background = element_rect(size=0.2, linetype="solid", 
#                                     colour ="grey20"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

mean_fst_dens

pdf("~/UBC/LinkageCoalescent/writeUp/Fst_Plot.pdf", height = 4.5, width = 5)
print(mean_fst_dens)
dev.off()

plot(density(fst[fst$rec == 0,]$weighted_Fst, bw = 0.01), xlim = c(0, 1 ), xaxs="i")

dxy_means<-tapply(fst$d_xy, fst$rec, mean)

DxyMeanDF <- data.frame( Dxy_mean = dxy_means, rec = names(dxy_means))

mean_dxy_dens <- ggplot(data = fst, aes( x = d_xy, fill = as.factor(rec), col = as.factor(rec))) + 
  geom_vline(data = DxyMeanDF, aes(xintercept = dxy_means, col = as.factor(rec)), lwd = 1.5 )+
  geom_line( adjust = 1, lwd = 1, stat = 'density',  position = 'identity')+
  scale_x_continuous( expression(italic(D[XY])), limits = c(0,0.15)) +
  scale_fill_manual("Recombination\nRate (M/bp)", values = pal)+
  scale_color_manual("Recombination\nRate (M/bp)",values = pal)+
  scale_y_continuous("Density", expand = c(0, 0), limits = c(0,130))+
  guides(colour = F, fill = F)+
         #  scale_y_continuous("Frequency", expand = c(0, 0), limits = c(0,14))+
  theme_bw()+
  theme(
    axis.text.x  = element_text(colour = 'black'),
    axis.text.y  = element_text(colour = 'black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

pdf("../writeUp/Fst_Dxy_Plot.pdf", height = 10, width = 5)
ggarrange(mean_fst_dens, mean_dxy_dens, ncol = 1, nrow = 2)
dev.off()


png("../writeUp/Fst_Dxy_Plot.png", height = 10, width = 5, units = 'in', res = 300)
ggarrange(mean_fst_dens, mean_dxy_dens, ncol = 1, nrow = 2)
dev.off()


  
  
