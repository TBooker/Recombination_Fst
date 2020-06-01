rm (list = ls())
library(moments)
library(ggplot2)
library(dplyr)
library(ggpubr)
### Read in the data and include only the "normally" recombining Drosophila regions as defined in the Reinhardt paper..

d2L <- read.csv('../Reinhardt_et_al_data/namerica.alleles.chr2L.MAF0.05.fst', sep = ' ', head = 0)
# Exclude the last row as it is often screwy...
d2L_n <- dim(d2L)[1]
d2L <- d2L[1:(d2L_n-1),]

## Now read in the recombination rate estimates, add headers to each CSV and add the recombination rates 
d2L_rates <- read.csv('../Reinhardt_et_al_data/RRC_input/RRC_Input_empirical.chr2L.WeirCockerham.MAF0.05.txt.rrc', sep = '\t', head = 0)
names( d2L_rates ) <- c('window','nothing','RRC_start','RRC_mid','RRC_end','Comeron_start','Comeron_mid','Comeron_end')
names(d2L) <- c('chrom','start','end','variants','T1_sum','T2_sum','Fst')
d2L$rec <- d2L_rates$Comeron_mid

## Repeat for the other chromosomes...
d2R <- read.csv('../Reinhardt_et_al_data/namerica.alleles.chr2R.MAF0.05.fst', sep = ' ', head = 0)
d2R_n <- dim(d2R)[1]
d2R <- d2R[1:(d2R_n-1),]
d2R_rates <- read.csv('../Reinhardt_et_al_data/RRC_input/RRC_Input_empirical.chr2R.WeirCockerham.MAF0.05.txt.rrc', sep = '\t', head = 0)
names( d2R_rates ) <- c('window','nothing','RRC_start','RRC_mid','RRC_end','Comeron_start','Comeron_mid','Comeron_end')
names(d2R) <- c('chrom','start','end','variants','T1_sum','T2_sum','Fst')
d2R$rec <- d2R_rates$Comeron_mid


d3L <- read.csv('../Reinhardt_et_al_data/namerica.alleles.chr3L.MAF0.05.fst', sep = ' ', head = 0)
d3L_n <- dim(d3L)[1]
d3L <- d3L[1:(d3L_n-1),]
d3L_rates <- read.csv('../Reinhardt_et_al_data/RRC_input/RRC_Input_empirical.chr3L.WeirCockerham.MAF0.05.txt.rrc', sep = '\t', head = 0)
names( d3L_rates ) <- c('window','nothing','RRC_start','RRC_mid','RRC_end','Comeron_start','Comeron_mid','Comeron_end')
names(d3L) <- c('chrom','start','end','variants','T1_sum','T2_sum','Fst')
d3L$rec <- d3L_rates$Comeron_mid


d3R <- read.csv('../Reinhardt_et_al_data/namerica.alleles.chr3R.MAF0.05.fst', sep = ' ', head = 0)
d3R_n <- dim(d3R)[1]
d3R <- d3R[1:(d3R_n-1),]
d3R_rates <- read.csv('../Reinhardt_et_al_data/RRC_input/RRC_Input_empirical.chr3R.WeirCockerham.MAF0.05.txt.rrc', sep = '\t', head = 0)
names( d3R_rates ) <- c('window','nothing','RRC_start','RRC_mid','RRC_end','Comeron_start','Comeron_mid','Comeron_end')
names(d3R) <- c('chrom','start','end','variants','T1_sum','T2_sum','Fst')
d3R$rec <- d3R_rates$Comeron_mid

## Now combine all chromosomes into one file
dros <- na.omit(rbind(d2L, d2R, d3L , d3R))
# What is the mean Fst across all windows?
mean(dros$Fst)

# Mkae the plot of Fst against recombination for the empirical data...
rec_v_fst_dros <- ggplot( data = dros[dros$rec > 0,], aes( x = rec, y = Fst) )+
  geom_hline(aes(yintercept = quantile(dros$Fst, 0.95)), lty = 2)+
  scale_x_continuous(expression('Recombination Rate ('*italic("cM/Mbp")*")"))+
  geom_point(alpha = 0.3, col = '#fc8d59', shape = 16)+
#  annotate("text", label = expression(italic("Drosophila melanogaster")), x = 7.5, y = 0.35, size = 4)+
    scale_y_continuous(expression(italic(hat(F)["ST"])), limits = c(-0.05,0.35))+
  theme_bw()+
  theme(
    strip.text = element_text( size = 12),
    axis.text.x = element_text( size = 12),
    axis.title.x = element_text( size = 12),
    axis.text.y = element_text( size = 12),
    axis.title.y = element_text( size = 12)
  )

# Make the plot of Fst against recombination for the empirical data...
rec_v_fst_dros <- ggplot( data = dros[dros$rec > 0,], aes( x = rec, y = Fst) )+
  geom_hline(aes(yintercept = quantile(dros$Fst, 0.95)), lty = 2)+
  scale_x_continuous(expression('Recombination Rate ('*italic("cM/Mbp")*")"))+
  geom_point(alpha = 0.3, col = '#fc8d59', shape = 16)+
  #  annotate("text", label = expression(italic("Drosophila melanogaster")), x = 7.5, y = 0.35, size = 4)+
  scale_y_continuous(expression(italic(hat(F)["ST"])), limits = c(-0.05,0.35))+
  theme_bw()+
  theme(
    strip.text = element_text( size = 12),
    axis.text.x = element_text( size = 12),
    axis.title.x = element_text( size = 12),
    axis.text.y = element_text( size = 12),
    axis.title.y = element_text( size = 12)
  )

TommyTheme <-   theme_bw() +
  theme( panel.spacing = unit(0, 'lines'),
         strip.background = element_blank(),
         strip.text.x =  element_text( face = 'italic', size = 20),
         axis.title.y = element_text( face = 'italic', size = 20),
         axis.title.x = element_text( face = 'italic', size = 20),
         axis.text.y = element_text(  size = 15),
         axis.text.x = element_text(  size = 15)
  )


Comparison <-ggplot( data = dros[dros$rec > 0,],  aes( x =rec, y = Fst) )+
#  geom_hline(aes(yintercept = quantile(sims$WEIGHTED_FST, 0.99)), col = 'red', lwd = 2, lty = 2)+
  geom_point(alpha = 0.4, col = 'blue')+
  scale_y_continuous(expression(italic(F[ST])), limits = c(-0.05,0.3))+
  scale_x_continuous("Recombination Rate (cM/Mbp)", limits = c(0,17))+
  TommyTheme

png("../Plots/Fst_Recombination_Dros.png",  height = 8, width = 8, units = 'in', res = 300)
print(Comparison)
dev.off()

## remove all windows with 0 recombination...
dros <- dros[ dros$rec != 0 ,]

# Assign each observation into a approximately equally sized recombination bins 
dros$rec_rank <- ntile( dros$rec, 100)

# Calculate the skewness for each bin
dros_skewness <-  tapply(dros$Fst, dros$rec_rank, skewness) 
# Calculate the variance for each bin
dros_variance <-  tapply(dros$Fst, dros$rec_rank, var) 
# Calculate the mean recombination rate for each bin
dros_recs <- tapply(dros$rec, dros$rec_rank, mean) 

## Combine all of those stats into a single data frame...
dros_summaryDF = data.frame(Skewness=dros_skewness, RecombinationRate=dros_recs, Variance=dros_variance  )

## Look at the Kendall correlation for recombination and variance in Fst
cor.test( dros_summaryDF$Variance, dros_summaryDF$RecombinationRate, method = 'kendall'  )
cor.test( dros_summaryDF$Skewness, dros_summaryDF$RecombinationRate, method = 'kendall'  )


## Test for enrichment of outliers in low v high recombination rate regions
## Order the observed data by recombination rate
dros_orders <- dros[order(dros$rec), ]
## Grab the Fst estimates from the top and bottom 1000 recombination rate  
dros_orders_lower <- dros_orders[1:1000, ]
dros_orders_upper <- dros_orders[(nrow(dros_orders)-999):nrow(dros_orders), ]
## Determine how many outliers there are in these sets - using the 95% percentile as a cut-off
dros_lower_hits = dros_orders_lower[dros_orders_lower$Fst > quantile(dros$Fst, 0.95), ]
dros_upper_hits = dros_orders_upper[dros_orders_upper$Fst > quantile(dros$Fst, 0.95), ]
## Format as a matrix and run Fisher's exact test
dros_enrichTest_0.95 <- matrix(c(nrow(dros_orders_lower) - nrow(dros_lower_hits), nrow(dros_lower_hits), nrow(dros_orders_upper)- nrow(dros_upper_hits), nrow(dros_upper_hits)), ncol = 2)
fisher.test(dros_enrichTest_0.95)

dros_lower_hits = dros_orders_lower[dros_orders_lower$Fst > quantile(dros$Fst, 0.975), ]
dros_upper_hits = dros_orders_upper[dros_orders_upper$Fst > quantile(dros$Fst, 0.975), ]

dros_enrichTest_0.975 <- matrix(c(nrow(dros_orders_lower) - nrow(dros_lower_hits), nrow(dros_lower_hits), nrow(dros_orders_upper)- nrow(dros_upper_hits), nrow(dros_upper_hits)), ncol = 2)
fisher.test(dros_enrichTest_0.975)

dros_lower_hits = dros_orders_lower[dros_orders_lower$Fst > quantile(dros$Fst, 0.999), ]
dros_upper_hits = dros_orders_upper[dros_orders_upper$Fst > quantile(dros$Fst, 0.999), ]

dros_enrichTest_0.999 <- matrix(c(nrow(dros_orders_lower) - nrow(dros_lower_hits), nrow(dros_lower_hits), nrow(dros_orders_upper)- nrow(dros_upper_hits), nrow(dros_upper_hits)), ncol = 2)
fisher.test(dros_enrichTest_0.999)


dros$source <- 'Drosophila'

################################################
################################################
### Simulated data
###
###

s2L <- read.csv('../Simulated/drosophilaSimulated.chr2L.fst.csv')
s2L<- s2L[(s2L$BIN_START> 844225)&(s2L$BIN_END< 19946732),]
s2R <- read.csv('../Simulated/drosophilaSimulated.chr2R.fst.csv')
s2R<- s2R[(s2R$BIN_START> 6063980)&(s2R$BIN_END< 20322335),]
s3L <- read.csv('../Simulated/drosophilaSimulated.chr3L.fst.csv')
s3L<- s3L[(s3L$BIN_START> 447386)&(s3L$BIN_END< 18392988),]
s3R <- read.csv('../Simulated/drosophilaSimulated.chr3R.fst.csv')
s3R<- s3R[(s3R$BIN_START> 7940899)&(s3R$BIN_END< 27237549),]

sims <- na.omit(rbind(s2L, s2R, s3L , s3R))

sims$source <- 'simulation'

library(ggplot2)

rec_v_fst_sims <- ggplot( data = sims, aes( x = rec*1e8, y = WEIGHTED_FST) )+
  geom_hline(aes(yintercept = quantile(sims$WEIGHTED_FST, 0.95)), lty = 2)+
  scale_x_continuous(expression('Recombination Rate ('*italic("cM/Mbp")*")"))+
  geom_point(alpha = 0.3, col = '#91bfdb', shape = 17)+
#  annotate("text", label = expression("Simulated Data"), x = 7.5, y = 0.35, size = 4)+
  scale_y_continuous(expression(italic(hat(F)["ST"])), limits = c(-0.05,0.35))+
  theme_bw()+
  theme(
    strip.text = element_text( size = 12),
    axis.text.x = element_text( size = 12),
    axis.title.x = element_text( size = 12),
    axis.text.y = element_text( size = 12),
    axis.title.y = element_text( size = 12)
  )



sims <- sims[ sims$rec != 0 ,]

cor.test(sims$WEIGHTED_FST , sims$rec)

sims$rec_rank <- ntile( sims$rec, 100)
sims_skewness <-  tapply(sims$WEIGHTED_FST, sims$rec_rank, skewness) 
sims_variance <-  tapply(sims$WEIGHTED_FST, sims$rec_rank, var) 
sims_recs <- tapply(sims$rec, sims$rec_rank, mean)*1e8



sims_summaryDF = data.frame(sims_skewness, sims_recs, sims_variance  )
names(sims_summaryDF) = c('Skewness','RecombinationRate','Variance')

cor.test( sims_summaryDF$Variance, sims_summaryDF$RecombinationRate, method = 'kendall'  )
cor.test( sims_summaryDF$Skewness, sims_summaryDF$RecombinationRate, method = 'kendall'  )
plot(sims_summaryDF)

sims_summaryDF$source <- 'Simulated'
dros_summaryDF$source <- 'Empirical'

## Test for enrichment of outliers in low v high recombination rate regions
sims_orders <- sims[order(sims$rec), ]
sims_orders_lower <- sims_orders[1:1000, ]
sims_orders_upper <- sims_orders[(nrow(sims_orders)-999):nrow(sims_orders), ]

sims_lower_hits = sims_orders_lower[sims_orders_lower$WEIGHTED_FST > quantile(sims$WEIGHTED_FST, 0.95), ]
sims_upper_hits = sims_orders_upper[sims_orders_upper$WEIGHTED_FST > quantile(sims$WEIGHTED_FST, 0.95), ]

sims_enrichTest_0.95 <- matrix(c(nrow(sims_orders_lower) - nrow(sims_lower_hits), nrow(sims_lower_hits), nrow(sims_orders_upper)- nrow(sims_upper_hits), nrow(sims_upper_hits)), ncol = 2)
fisher.test(sims_enrichTest_0.95)

sims_lower_hits = sims_orders_lower[sims_orders_lower$WEIGHTED_FST > quantile(sims$WEIGHTED_FST, 0.975), ]
sims_upper_hits = sims_orders_upper[sims_orders_upper$WEIGHTED_FST > quantile(sims$WEIGHTED_FST, 0.975), ]

sims_enrichTest_0.975 <- matrix(c(nrow(sims_orders_lower) - nrow(sims_lower_hits), nrow(sims_lower_hits), nrow(sims_orders_upper)- nrow(sims_upper_hits), nrow(sims_upper_hits)), ncol = 2)
fisher.test(sims_enrichTest_0.975)

sims_lower_hits = sims_orders_lower[sims_orders_lower$WEIGHTED_FST > quantile(sims$WEIGHTED_FST, 0.999), ]
sims_upper_hits = sims_orders_upper[sims_orders_upper$WEIGHTED_FST > quantile(sims$WEIGHTED_FST, 0.999), ]

sims_enrichTest_0.999 <- matrix(c(nrow(sims_orders_lower) - nrow(sims_lower_hits), nrow(sims_lower_hits), nrow(sims_orders_upper)- nrow(sims_upper_hits), nrow(sims_upper_hits)), ncol = 2)
fisher.test(sims_enrichTest_0.999)

### Now let's plot the variance in both simulated and empirical datasets

summaryDF<- rbind(sims_summaryDF, dros_summaryDF)

variance <- ggplot( )+
  geom_point(data = sims_summaryDF, aes( x = RecombinationRate*1e8, y =Variance, col  = 'Simulated'), size = 2,alpha = 1, shape = 17)+
  geom_point(data = dros_summaryDF, aes( x = RecombinationRate, y =Variance, col = 'Empirical'), size = 2 , alpha = 0.8, shape = 16)+
  scale_colour_manual( values = c('#fc8d59','#91bfdb'), labels = c('Empirical', 'Simulated'))+
  scale_x_continuous(expression('Recombination Rate ('*italic("cM/Mbp")*")"))+
  scale_y_continuous( expression("Variance ["*hat(italic(F[ST]))*"]" ))+
#  guides(shape = F)+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    strip.text = element_text( size = 12),
    axis.text.x = element_text( size = 12),
    axis.title.x = element_text( size = 12),
    axis.text.y = element_text( size = 12),
    axis.title.y = element_text( size = 12),
    legend.position = c(0.8, 0.80)
  )

variance <- ggplot( )+
#  geom_point(data = sims_summaryDF, aes( x = RecombinationRate*1e8, y =Variance, col  = 'Simulated'), size = 2,alpha = 1, shape = 17)+
  geom_point(data = summaryDF, aes( x = RecombinationRate, y =Variance*10e2, col = source, shape = source), size = 2 , alpha = 0.8)+
  scale_colour_manual( values = c('#fc8d59','#91bfdb'), labels = c('Empirical', 'Simulated'))+
  scale_x_continuous(expression('Recombination Rate ('*italic("cM/Mbp")*")"))+
  scale_y_continuous( expression("Variance ["*hat(italic(F[ST]))*"] " %*% 100  ))+
  #  guides(shape = F)+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    strip.text = element_text( size = 12),
    axis.text.x = element_text( size = 12),
    axis.title.x = element_text( size = 12),
    axis.text.y = element_text( size = 12),
    axis.title.y = element_text( size = 12),
    legend.position = c(0.8, 0.80)
  )


threePanelFigure <- ggarrange(rec_v_fst_dros, rec_v_fst_sims, variance, nrow = 1, ncol = 3, labels = 'AUTO')

png("../Plots/threePanelFigure.png",width = 12, height = 4, res = 300, units = 'in')
print(threePanelFigure)
dev.off()




