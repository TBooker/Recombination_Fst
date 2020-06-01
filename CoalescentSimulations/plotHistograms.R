
rm( list = ls())
options(scipen=6)


library(reshape2)
library(ggplot2)
library(scales)
library(ggpubr)

fst <- read.csv('finishedRuns/Island.csv')

fst_means<-tapply(fst$weighted_Fst, fst$rec, mean)

FstMeanDF <- data.frame( weightedFst = fst_means, rec = as.numeric(names(fst_means))*1e8)
## Figure for the paper

cbPalette <- c( "#009E73", "#0072B2", "#D55E00", "#CC79A7")

pal<- c('#a6611a','#dfc27d','#80cdc1','#018571')


mean_fst_dens_main_text <- ggplot(data = fst, aes( x = weighted_Fst, col = as.factor(rec*1e8))) + 
  geom_vline(data = FstMeanDF, aes(xintercept = weightedFst, col = as.factor(rec)), lwd = 0.7, alpha = 0.7 , show.legend = F, lty = 2)+
  geom_vline(aes(xintercept = 1/9), lwd = 0.7, lty = 2)+
  geom_line(adjust = 1, stat = 'density',lwd = 1, position = 'identity')+
  scale_x_continuous( expression(italic(hat(F)["ST - Weir and Cockerham"])) ) +
  scale_fill_manual("Recombination\nRate (cM/Mbp)", values = pal)+
  scale_color_manual("Recombination\nRate (cM/Mbp)",values = pal)+
  scale_y_continuous("Density", expand = c(0, 0), limits = c(0,30))+
  theme_bw()+
  theme(
    axis.text.x  = element_text(colour = 'black'),
    axis.text.y  = element_text(colour = 'black'),
    legend.position = c(0.8, 0.50),
    legend.text.align = 0,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

### For figure 1

pdf("Plots/Island_Fst_Plot.pdf", height = 4.5, width = 5)
print(mean_fst_dens_main_text)
dev.off()


mean_fst_dens <- ggplot(data = fst, aes( x = weighted_Fst, col = as.factor(rec*1e8))) + 
  geom_vline(data = FstMeanDF, aes(xintercept = weightedFst, col = as.factor(rec)), lwd = 0.7, alpha = 0.7 , show.legend = F, lty = 2)+
  geom_line(adjust = 1, stat = 'density',lwd = 1, position = 'identity')+
  scale_x_continuous( expression(italic(hat(F)["ST - Weir and Cockerham"])) ) +
  scale_fill_manual("Recombination\nRate (cM/Mbp)", values = pal)+
  scale_color_manual("Recombination\nRate (cM/Mbp)",values = pal)+
  scale_y_continuous("Density", expand = c(0, 0), limits = c(0,30))+
 # guides(vline= FALSE)+
    theme_bw()+
  theme(
    axis.text.x  = element_text(colour = 'black'),
    axis.text.y  = element_text(colour = 'black'),
    legend.text.align = 0,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )


dxy_means<-tapply(fst$d_xy, fst$rec, mean)
DxyMeanDF <- data.frame( Dxy_mean = dxy_means, rec= as.numeric(names(dxy_means))*1e8)

mean_dxy_dens <- ggplot(data = fst, aes( x = d_xy, fill = as.factor(rec), col = as.factor(rec*1e8))) + 
  geom_vline(data = DxyMeanDF, aes(xintercept = dxy_means, col = as.factor(rec)), lwd = 0.7, alpha = 0.7 , show.legend = F, lty = 2)+
  geom_line( adjust = 1, lwd = 1, stat = 'density',  position = 'identity')+
  scale_x_continuous( expression(italic(D[XY]))) +
  scale_fill_manual("Recombination\nRate (cM/Mbp)", values = pal)+
  scale_color_manual("Recombination\nRate (cM/Mbp)",values = pal)+
  scale_y_continuous("Density", expand = c(0, 0), limits = c(0,800))+
  guides(colour = F, fill = F)+
         #  scale_y_continuous("Frequency", expand = c(0, 0), limits = c(0,14))+
  theme_bw()+
  theme(
    axis.text.x  = element_text(colour = 'black'),
    axis.text.y  = element_text(colour = 'black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )




pdf("Plots/Fst_Dxy_Plot.pdf", height = 10, width = 5)
ggarrange(mean_fst_dens, mean_dxy_dens, ncol = 1, nrow = 2)
dev.off()


png("Plots/Fst_Dxy_Plot.png", height = 10, width = 5, units = 'in', res = 300)
ggarrange(mean_fst_dens, mean_dxy_dens, ncol = 1, nrow = 2)
dev.off()

#### Density plots for each summary stat 

diversity_means<-tapply(fst$diversity, fst$rec, mean)
diversityMeanDF <- data.frame( diversity_mean = diversity_means, rec = as.numeric(names(diversity_means))*1e8)

mean_diversity_dens <- ggplot(data = fst, aes( x = diversity, fill = as.factor(rec), col = as.factor(rec*1e8))) + 
  geom_vline(data = diversityMeanDF, aes(xintercept = diversity_means, col = as.factor(rec)), lwd = 0.7, alpha = 0.7 , show.legend = F, lty = 2)+
  geom_line( adjust = 1, lwd = 1, stat = 'density',  position = 'identity')+
  scale_x_continuous( expression(italic(pi[W,1]))) +
  scale_fill_manual("Recombination\nRate (cM/Mbp)", values = pal)+
  scale_color_manual("Recombination\nRate (cM/Mbp)",values = pal)+
  scale_y_continuous("Density", expand = c(0, 0), limits = c(0,800))+
#  guides(colour = F, fill = F)+
  #  scale_y_continuous("Frequency", expand = c(0, 0), limits = c(0,14))+
  theme_bw()+
  theme(
    axis.text.x  = element_text(colour = 'black'),
    axis.text.y  = element_text(colour = 'black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

tajima_means<-tapply(fst$Tajima_D, fst$rec, mean)
TajimaMeanDF <- data.frame( tajima_mean = tajima_means, rec = as.numeric(names(tajima_means))*1e8)

mean_tajima_dens <- ggplot(data = fst, aes( x = Tajima_D, fill = as.factor(rec), col = as.factor(rec*1e8))) + 
  geom_vline(data = TajimaMeanDF, aes(xintercept = tajima_mean, col = as.factor(rec)), lwd = 0.7, alpha = 0.7 , show.legend = F, lty = 2)+
  geom_line( adjust = 1, lwd = 1, stat = 'density',  position = 'identity')+
  scale_x_continuous( expression(italic("Tajima's D"))) +
  scale_fill_manual("Recombination\nRate (cM/Mbp)", values = pal)+
  scale_color_manual("Recombination\nRate (cM/Mbp)",values = pal)+
#  scale_y_continuous("Density", expand = c(0, 0), limits = c(0,800))+
#  guides(colour = F, fill = F)+
  #  scale_y_continuous("Frequency", expand = c(0, 0), limits = c(0,14))+
  theme_bw()+
  theme(
    axis.text.x  = element_text(colour = 'black'),
    axis.text.y  = element_text(colour = 'black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

H12_means<-tapply(fst$H12, fst$rec, mean)
H12MeanDF <- data.frame( H12_mean = H12_means, rec = as.numeric(names(H12_means))*1e8)

mean_H12_dens <- ggplot(data = fst, aes( x = H12, fill = as.factor(rec), col = as.factor(rec*1e8))) + 
  geom_vline(data = H12MeanDF, aes(xintercept = H12_mean, col = as.factor(rec)), lwd = 0.7, alpha = 0.7 , show.legend = F, lty = 2)+
  geom_line( adjust = 1, lwd = 1, stat = 'density',  position = 'identity')+
  scale_x_continuous( expression(italic("Garud's H12")), limits = c(0,0.15)) +
  scale_fill_manual("Recombination\nRate (cM/Mbp)", values = pal)+
  scale_color_manual("Recombination\nRate (cM/Mbp)",values = pal)+
#  scale_y_continuous("Density", expand = c(0, 0), limits = c(0,800))+
#  guides(colour = F, fill = F)+
  #  scale_y_continuous("Frequency", expand = c(0, 0), limits = c(0,14))+
  theme_bw()+
  theme(
    axis.text.x  = element_text(colour = 'black'),
    axis.text.y  = element_text(colour = 'black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

# Combine into one 5 panel plot with a common legend
all_summary_stats <- ggarrange(mean_fst_dens, mean_dxy_dens, mean_diversity_dens, mean_tajima_dens, mean_H12_dens, ncol = 5, nrow = 1, common.legend = T)

pdf("Plots/summary_stat_densities.pdf", width = 18, height = 4)
print(all_summary_stats)
dev.off()


##############################################
##############################################
##### THERE IS NO SPOON#######################
##############################################
##############################################
### MAKE THE MATRIX PLOT #####################
##############################################

fst_v_dxy <- ggplot(data = fst, aes(x = weighted_Fst, y = d_xy, col = as.factor(rec*1e8)))+
    geom_point(alpha = 0.5)+
    facet_wrap(~rec)+
    scale_y_continuous(expression(italic(D[XY])))+
    scale_x_continuous(expression(italic(hat(F)["ST - Weir and Cockerham"])))+
    theme_bw()+
    scale_colour_manual("Recombination \nRate (cM/Mbp)", values = pal)+
    theme(
      strip.text = element_blank(),
      axis.text.x = element_text( size = 9),
      axis.title.x = element_blank( ),
      axis.text.y = element_text( size = 9),
      axis.title.y = element_text( size = 9)
      #    legend.position = c(0.85, 0.85)
    )
  
  
fst_v_diversity <- ggplot(data = fst, aes(x = weighted_Fst, y = diversity, col = as.factor(rec*1e8)))+
  geom_point(alpha = 0.5)+
  facet_wrap(~rec)+
  scale_y_continuous(expression(italic(pi)))+
  scale_x_continuous(expression(italic(hat(F)["ST - Weir and Cockerham"])))+
  theme_bw()+
  scale_colour_manual("Recombination \nRate (cM/Mbp)", values = pal)+
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text( size = 9),
    axis.title.x = element_blank( ),
    axis.text.y = element_text( size = 9),
    axis.title.y = element_text( size = 9)
    #    legend.position = c(0.85, 0.85)
  )

fst_v_tajima <- ggplot(data = fst, aes(x = weighted_Fst, y = Tajima_D, col = as.factor(rec*1e8)))+
  geom_point(alpha = 0.5)+
  facet_wrap(~rec)+
  scale_y_continuous(expression(italic("Tajima's D")))+
  scale_x_continuous(expression(italic(hat(F)["ST - Weir and Cockerham"])))+
  theme_bw()+
  scale_colour_manual("Recombination \nRate (cM/Mbp)", values = pal)+
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text( size = 9),
    axis.title.x = element_blank( ),
    axis.text.y = element_text( size = 9),
    axis.title.y = element_text( size = 9)
    #    legend.position = c(0.85, 0.85)
  )


fst_v_H12 <- ggplot(data = fst, aes(x = weighted_Fst, y = H12, col = as.factor(rec*1e8)))+
  geom_point(alpha = 0.5)+
  facet_wrap(~rec)+
  scale_y_continuous(expression(italic("Garud's H12")))+
  scale_x_continuous(expression(italic(hat(F)["ST - Weir and Cockerham"])))+
  theme_bw()+
  scale_colour_manual("Recombination \nRate (cM/Mbp)", values = pal)+
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text( size = 9),
    axis.title.x = element_text( size = 9),
    axis.text.y = element_text( size = 9),
    axis.title.y = element_text( size = 9)
    #    legend.position = c(0.85, 0.85)
  )

dxy_v_diversity <- ggplot(data = fst, aes(x = d_xy, y = diversity, col = as.factor(rec*1e8)))+
  geom_point(alpha = 0.5)+
  facet_wrap(~rec)+
  scale_y_continuous(expression(italic(pi["W,1"])))+
  scale_x_continuous(expression(italic(D[XY])))+
  theme_bw()+
  scale_colour_manual("Recombination \nRate (cM/Mbp)", values = pal)+
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text( size = 9),
    axis.title.x = element_blank( ),
    axis.text.y = element_text( size = 9),
    axis.title.y = element_blank( )
    #    legend.position = c(0.85, 0.85)
  )

dxy_v_tajima <- ggplot(data = fst, aes(x = d_xy, y = Tajima_D, col = as.factor(rec*1e8)))+
  geom_point(alpha = 0.5)+
  facet_wrap(~rec)+
  scale_y_continuous(expression(italic("Tajima's D")))+
  scale_x_continuous(expression(italic(D[XY])))+
  theme_bw()+
  scale_colour_manual("Recombination \nRate (cM/Mbp)", values = pal)+
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text( size = 9),
    axis.title.x = element_blank( ),
    axis.text.y = element_text( size = 9),
    axis.title.y = element_blank( )
    #    legend.position = c(0.85, 0.85)
  )


dxy_v_H12 <- ggplot(data = fst, aes(x = d_xy, y = H12, col = as.factor(rec*1e8)))+
  geom_point(alpha = 0.5)+
  facet_wrap(~rec)+
  scale_y_continuous(expression(italic("Garud's H12")))+
  scale_x_continuous(expression(italic(D[XY])))+
  theme_bw()+
  scale_colour_manual("Recombination \nRate (cM/Mbp)", values = pal)+
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text( size = 9),
    axis.title.x = element_text( size = 9),
    axis.text.y = element_text( size = 9),
    axis.title.y = element_blank( )
    #    legend.position = c(0.85, 0.85)
  )

diversity_v_tajima <- ggplot(data = fst, aes(x = diversity, y = Tajima_D, col = as.factor(rec*1e8)))+
  geom_point(alpha = 0.5)+
  facet_wrap(~rec)+
  scale_y_continuous(expression(italic("Tajima's D")))+
  scale_x_continuous(expression(italic(pi)))+
  theme_bw()+
  scale_colour_manual("Recombination \nRate (cM/Mbp)", values = pal)+
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text( size = 9),
#    axis.title.x = element_text( size = 9),
    axis.text.y = element_text( size = 9),
#    axis.title.y = element_text( size = 9),
    axis.title.x = element_blank( ),
    axis.title.y = element_blank( )
    #    legend.position = c(0.85, 0.85)
  )


diversity_v_H12 <- ggplot(data = fst, aes(x = diversity, y = H12, col = as.factor(rec*1e8)))+
  geom_point(alpha = 0.5)+
  facet_wrap(~rec)+
  scale_y_continuous(expression(italic("Garud's H12")))+
  scale_x_continuous(expression(italic(pi)))+
  theme_bw()+
  scale_colour_manual("Recombination \nRate (cM/Mbp)", values = pal)+
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text( size = 9),
    axis.title.x = element_text( size = 9),
    axis.text.y = element_text( size = 9),
 #   axis.title.y = element_text( size = 9),
    #    axis.title.y = element_text( size = 9),
    axis.title.y = element_blank( )
  )


tajima_v_H12 <- ggplot(data = fst, aes(x = Tajima_D, y = H12, col = as.factor(rec*1e8)))+
  geom_point(alpha = 0.5)+
  facet_wrap(~rec)+
  scale_y_continuous(expression(italic("Garud's H12")))+
  scale_x_continuous(expression(italic("Tajima's D")))+
  theme_bw()+
  scale_colour_manual("Recombination \nRate (cM/Mbp)", values = pal)+
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text( size = 9),
    axis.title.x = element_text( size = 9),
    axis.text.y = element_text( size = 9),
#    axis.title.y = element_text( size = 12),
    axis.title.y = element_blank( )
    #    legend.position = c(0.85, 0.85)
  )


library(ggpubr)
library(gridExtra)
library(grid)
ng = nullGrob()
matrixPlot <- ggarrange( fst_v_dxy, ng, ng, ng,
           fst_v_diversity, dxy_v_diversity, ng, ng,
           fst_v_tajima, dxy_v_tajima, diversity_v_tajima, ng,
           fst_v_H12, dxy_v_H12, diversity_v_H12, tajima_v_H12, common.legend = T, nrow = 4, ncol = 4, legend = 'bottom')

pdf('Plots/summary_stat_matrix.pdf', height  = 20, width = 20)
print(matrixPlot)
dev.off()

png('Plots/summary_stat_matrix.png', height  = 20, width = 20, res = 300, units = 'in')
print(matrixPlot)
dev.off()

#### Plot the Fst density for each demographic model

two_patch <- read.csv('finishedRuns/2Patch.csv')
two_patch$model <- 'Two Deme Model'
island <- read.csv('finishedRuns/Island.csv')
island$model <- '100 Deme Island Model'
ring <- read.csv('finishedRuns/Ring.csv')
ring$model <- '20 Deme Ring Model'
cline <- read.csv('finishedRuns/cline.csv')
cline$model <- '20 Deme 1-Dimensional Cline'

models <- rbind ( two_patch, island, ring, cline )

#fst_means<-tapply(fst$weighted_Fst, fst$rec, mean)

FstMeanDF <- data.frame( weightedFst = fst_means, rec = names(fst_means))
## Figure for the paper
cbPalette <- c( "#009E73", "#0072B2", "#D55E00", "#CC79A7")

pal<- c('#a6611a','#dfc27d','#80cdc1','#018571')


model_plots <- ggplot(data = models, aes( x = weighted_Fst, col = as.factor(rec*1e8))) + 
  geom_line(stat = 'density',lwd = 1, position = 'identity')+
  facet_wrap(~model)+
  scale_x_continuous( expression(italic(hat(F)["ST - Weir and Cockerham"])) ) +
  scale_fill_manual("Recombination\nRate (cM/Mbp)", values = pal)+
  scale_color_manual("Recombination\nRate (cM/Mbp)",values = pal)+
  scale_y_continuous("Density", expand = c(0, 0), limits = c(0,25))+
  theme_bw()+
  theme(
    axis.text.x  = element_text(colour = 'black'),
    axis.text.y  = element_text(colour = 'black'),
    legend.text.align = 0,
    strip.text.x = element_text(face = 'bold'), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

pdf("Plots/Fst_by_model.pdf", height = 8, width = 10)
print(model_plots)
dev.off()

