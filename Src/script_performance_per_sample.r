library(ggplot2)
library(grid)
library(gridExtra)

file_perf_perSample = '../Results/ReconstructionAndPredictionMicrobialComposition/performance_per_sample.tsv'

# Read table with one row per sample, and different columns with different performance metrics (pearson_reconstructed,	pearson_predicted,	braycurtis_reconstructed,	braycurtis_predicted)
perf_perSample <- read.table(file_perf_perSample,sep='\t',row.names=1,header=TRUE,check.names=FALSE)
  
#plot(perf_perSample[2:1], main='Pearson Correlation versus original', xlab='Predicted from environment', ylab='Reconstructed')

#plot(perf_perSample[4:3], main='Bray-Curtis dissimilarity versus original', xlab='Predicted from environment', ylab='Reconstructed')

# Plot
p_pearson <- ggplot(perf_perSample, aes(x = pearson_predicted , y = pearson_reconstructed)) +
  geom_point(size = 0.5) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  xlab('') +
  ylab('Reconstructed') +
  ggtitle('Pearson correlation')

p_braycurtis <- ggplot(perf_perSample, aes(x = braycurtis_predicted , y = braycurtis_reconstructed)) +
  geom_point(size = 0.5) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  xlab('') +
  ylab('') +
  xlim(0, 1) +
  ylim(0, 1) +
  ggtitle('Bray-Curtis dissimilarity')

grid.newpage()
my_layout<-rbind(c(1,2))
combined<-grid.arrange(arrangeGrob(p_pearson),arrangeGrob(p_braycurtis), layout_matrix=my_layout,bottom = textGrob('Predicted from environment'))


pdf('../Results/ReconstructionAndPredictionMicrobialComposition/performance_per_sample_scatterplot.pdf', height=3, width=5)
grid.draw(combined)
dev.off()

