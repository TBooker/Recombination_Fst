rm (list = ls())
library(moments)
library(ggplot2)
library(dplyr)
library(ggpubr)


################################################
################################################
### Simulated data
### 10000bp windows
###

s2L <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/VCF/drosophilaSimulated.chr2L.MAF0.05.fst.csv')
s2L<- s2L[(s2L$BIN_START> 844225)&(s2L$BIN_END< 19946732),]
s2R <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/VCF/drosophilaSimulated.chr2R.MAF0.05.fst.csv')
s2R<- s2R[(s2R$BIN_START> 6063980)&(s2R$BIN_END< 20322335),]
s3L <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/VCF/drosophilaSimulated.chr3L.MAF0.05.fst.csv')
s3L<- s3L[(s3L$BIN_START> 447386)&(s3L$BIN_END< 18392988),]
s3R <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/VCF/drosophilaSimulated.chr3R.MAF0.05.fst.csv')
s3R<- s3R[(s3R$BIN_START> 7940899)&(s3R$BIN_END< 27237549),]

phys_sims <- na.omit(rbind(s2L, s2R, s3L , s3R))
phys_sims$source <- '10,000bp Windows'


################################################
################################################
### Simulated data
### 50 SNP windows
###

snp50_s2L <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/SNP_windows/drosophilaSimulated.chr2L.MAF0.05.50SNPs.fst.csv')
snp50_s2L<- snp50_s2L[(snp50_s2L$BIN_START> 844225)&(snp50_s2L$BIN_END< 19946732),]
snp50_s2R <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/SNP_windows/drosophilaSimulated.chr2R.MAF0.05.50SNPs.fst.csv')
snp50_s2R<- snp50_s2R[(snp50_s2R$BIN_START> 6063980)&(snp50_s2R$BIN_END< 20322335),]
snp50_s3L <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/SNP_windows/drosophilaSimulated.chr3L.MAF0.05.50SNPs.fst.csv')
snp50_s3L<- snp50_s3L[(snp50_s3L$BIN_START> 447386)&(snp50_s3L$BIN_END< 18392988),]
snp50_s3R <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/SNP_windows/drosophilaSimulated.chr3R.MAF0.05.50SNPs.fst.csv')
snp50_s3R<- snp50_s3R[(snp50_s3R$BIN_START> 7940899)&(snp50_s3R$BIN_END< 27237549),]

snp50_sims <- na.omit(rbind(snp50_s2L, snp50_s2R, snp50_s3L , snp50_s3R))

snp50_sims$source <- '50 SNP Windows'

################################################
################################################
### Simulated data
### 250 SNP windows
###

head(snp250_s2R)
snp250_s2L <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/SNP_windows/drosophilaSimulated.chr2L.MAF0.05.250SNPs.fst.csv')
snp250_s2L<- snp250_s2L[(snp250_s2L$BIN_START> 844225)&(snp250_s2L$BIN_END< 19946732),]
snp250_s2R <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/SNP_windows/drosophilaSimulated.chr2R.MAF0.05.250SNPs.fst.csv')
snp250_s2R<- snp250_s2R[(snp250_s2R$BIN_START> 6063980)&(snp250_s2R$BIN_END< 20322335),]
snp250_s3L <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/SNP_windows/drosophilaSimulated.chr3L.MAF0.05.250SNPs.fst.csv')
snp250_s3L<- snp250_s3L[(snp250_s3L$BIN_START> 447386)&(snp250_s3L$BIN_END< 18392988),]
snp250_s3R <- read.csv('~/work/LinkageCoalescent/Drosophila/Simulations/SNP_windows/drosophilaSimulated.chr3R.MAF0.05.250SNPs.fst.csv')
snp250_s3R<- snp250_s3R[(snp250_s3R$BIN_START> 7940899)&(snp250_s3R$BIN_END< 27237549),]

snp250_sims <- na.omit(rbind(snp250_s2L, snp250_s2R, snp250_s3L , snp250_s3R))

snp250_sims$source <- '250 SNP Windows'

mean(snp250_sims$BIN_END - snp250_sims$BIN_START)
mean(snp50_sims$BIN_END - snp50_sims$BIN_START)

sims <- rbind( phys_sims, snp50_sims, snp250_sims)

library(ggplot2)
mean(sims$BIN_END - sims$BIN_START)

rec_v_fst_sims <- ggplot( data = sims, aes( x = rec*1e8, y = WEIGHTED_FST) )+
#  geom_hline(aes(yintercept = quantile(sims$WEIGHTED_FST, 0.95)), lty = 2)+
  scale_x_continuous(expression('Recombination Rate ('*italic("cM/Mbp")*")"))+
  geom_point(alpha = 0.2, shape = 17)+
  #  annotate("text", label = expression("Simulated Data"), x = 7.5, y = 0.35, size = 4)+
  scale_y_continuous(expression(italic(hat(F)["ST"])), limits = c(-0.05,0.35))+
  facet_grid(~source)+
    theme_bw()+
  theme(
    strip.text = element_text( size = 12),
    axis.text.x = element_text( size = 12),
    axis.title.x = element_text( size = 12),
    axis.text.y = element_text( size = 12),
    axis.title.y = element_text( size = 12)
  )



phys_sims <- phys_sims[ phys_sims$rec != 0 ,]

phys_sims$rec_rank <- ntile( phys_sims$rec, 100)

phys_sims_skewness <-  tapply(phys_sims$WEIGHTED_FST, phys_sims$rec_rank, skewness) 
phys_sims_variance <-  tapply(phys_sims$WEIGHTED_FST, phys_sims$rec_rank, var) 
phys_sims_recs <- tapply(phys_sims$rec, phys_sims$rec_rank, mean)*1e8

phys_sims_summaryDF = data.frame(phys_sims_skewness, phys_sims_recs, phys_sims_variance  )
names(phys_sims_summaryDF) = c('Skewness','RecombinationRate','Variance')
phys_sims_summaryDF$source <- '10,000bp Windows'

cor.test( phys_sims_summaryDF$Variance, phys_sims_summaryDF$RecombinationRate, method = 'kendall'  )
cor.test( phys_sims_summaryDF$Skewness, phys_sims_summaryDF$RecombinationRate, method = 'kendall'  )


snp50_sims <- snp50_sims[ snp50_sims$rec != 0 ,]

snp50_sims$rec_rank <- ntile( snp50_sims$rec, 100)

snp50_sims_skewness <-  tapply(snp50_sims$WEIGHTED_FST, snp50_sims$rec_rank, skewness) 
snp50_sims_variance <-  tapply(snp50_sims$WEIGHTED_FST, snp50_sims$rec_rank, var) 
snp50_sims_recs <- tapply(snp50_sims$rec, snp50_sims$rec_rank, mean)*1e8

snp50_sims_summaryDF = data.frame(snp50_sims_skewness, snp50_sims_recs, snp50_sims_variance  )
names(snp50_sims_summaryDF) = c('Skewness','RecombinationRate','Variance')

cor.test( snp50_sims_summaryDF$Variance, snp50_sims_summaryDF$RecombinationRate, method = 'kendall'  )
cor.test( snp50_sims_summaryDF$Skewness, snp50_sims_summaryDF$RecombinationRate, method = 'kendall'  )


snp250_sims <- snp250_sims[ snp250_sims$rec != 0 ,]

snp250_sims$rec_rank <- ntile( snp250_sims$rec, 100)

snp250_sims_skewness <-  tapply(snp250_sims$WEIGHTED_FST, snp250_sims$rec_rank, skewness) 
snp250_sims_variance <-  tapply(snp250_sims$WEIGHTED_FST, snp250_sims$rec_rank, var) 
snp250_sims_recs <- tapply(snp250_sims$rec, snp250_sims$rec_rank, mean)*1e8

snp250_sims_summaryDF = data.frame(snp250_sims_skewness, snp250_sims_recs, snp250_sims_variance  )
names(snp250_sims_summaryDF) = c('Skewness','RecombinationRate','Variance')

cor.test( snp250_sims_summaryDF$Variance, snp250_sims_summaryDF$RecombinationRate, method = 'kendall'  )
cor.test( snp250_sims_summaryDF$Skewness, snp250_sims_summaryDF$RecombinationRate, method = 'kendall'  )


snp250_sims_summaryDF$source <- '250 SNP Windows'
snp50_sims_summaryDF$source <- '50 SNP Windows'

summary_df <- rbind(phys_sims_summaryDF, snp250_sims_summaryDF, snp50_sims_summaryDF)

variance <- ggplot( data = summary_df, aes( x = RecombinationRate, y = Variance) )+
  geom_point()+
#  geom_point(data = summaryDF, aes( x = RecombinationRate, y =Variance*10e2, col = source, shape = source), size = 2 , alpha = 0.8)+
  scale_colour_manual( values = c('#fc8d59','#91bfdb'), labels = c('Empirical', 'Simulated'))+
  scale_x_continuous(expression('Recombination Rate ('*italic("cM/Mbp")*")"))+
  scale_y_continuous( expression("Variance ["*hat(italic(F[ST]))*"] " %*% 100  ))+
  facet_grid(~source)+
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




sixPanelFigure <- ggarrange(rec_v_fst_sims, variance , nrow = 2, ncol = 1, align = 'v')

png("~/work/LinkageCoalescent/writeUp/sixPanelFigure.png",width = 12, height = 8, res = 300, units = 'in')
print(sixPanelFigure)
dev.off()
 

