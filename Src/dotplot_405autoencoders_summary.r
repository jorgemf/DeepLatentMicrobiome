library(ggplot2)
library(reshape2)
library(plyr)
library(grid)
library(gridExtra)


df=read.csv('../Results/results_experiments_hyperparameters.csv', header=TRUE)
names(df) <- gsub("\\.", "_", names(df))
names(df) <- gsub("domain_BrayCurtis", "BrayCurtis", names(df))
names(df) <- gsub("domain_pearson_corr", "Pearson", names(df))
names(df) <- gsub("Bioma_Autoencoder", "OTU_AE_architecture", names(df))
names(df) <- gsub("Domain_Autoencoder", "EnvFeatures_AE_architecture", names(df))
names(df) <- gsub("Latent_Space", "Latent_Size", names(df))
gsub("\\.", "_", names(df))
df$OTU_AE_architecture <- gsub("b", "OTU", df$OTU_AE_architecture)
df$EnvFeatures_AE_architecture <- gsub("b", "OTU", df$EnvFeatures_AE_architecture)
df$EnvFeatures_AE_architecture <- gsub("d", "Env", df$EnvFeatures_AE_architecture)
df$Learning_Rate <- gsub("constant = ", "k=", df$Learning_Rate)

# Select subset of columns of interest:
hyper_columns = c('Input_transform', 'Output_transform', 'Reconstruction_Loss', 
                  'Latent_Size', 'OTU_AE_architecture', 'EnvFeatures_AE_architecture',
                  'Activation_Encoder', 'Activation_Decoder', 'Activation_Latent',
                  'Batch_Size', 'Learning_Rate', 'Optimizer','Transformations')
metrics_columns = c('Pearson','BrayCurtis')
df=subset(df,select=c('experiment_n_',c(hyper_columns,metrics_columns)))

# Mean validation results averaging by hyperparameter
# Error bars code adapted from: https://gist.github.com/tomhopper/9674890
dotplot.per.hyperparameter <- function(df,id_hyper,hyper_columns,metrics_columns){
  df_hyper = subset(df,select=c(hyper_columns[id_hyper],metrics_columns))
  # Change from wide to tall format for ggplot2
  melted <- melt(df_hyper, id.vars=hyper_columns[id_hyper])
  # Calculate means
  # #Without melting, it would be:
  # means <- ddply(df, hyper_columns[id_hyper], summarise,
  #               mean_pear=mean(domain_pearson_corr),
  #               mean_bray=mean(domain_BrayCurtis))
  means <- ddply(melted, c(hyper_columns[id_hyper],'variable'), summarise,
                 mean=mean(value))
  # Calculate standard deviation (standard error mean)
  means.sem <- ddply(melted, c(hyper_columns[id_hyper],'variable'),
                     summarise, mean=mean(value), 
                     sem=sd(value)/sqrt(length(value)))
  # Dataframe for the error bars
  means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)
  
  # Curate rows names 
  if(length(id_hyper)>1){
    #means
    means$Transformations = paste(as.character(means[,1]),as.character(means[,2]),as.character(means[,3]),sep="_") 
    means[,hyper_columns[id_hyper]]<-NULL
    means$Transformations = gsub("none_|Wrapper|Loss|Categorical", "", means$Transformations)
    means$Transformations = gsub("Percentage", "TotalSumNormalization", means$Transformations)
    #means.sem
    means.sem$Transformations = paste(as.character(means.sem[,1]),as.character(means.sem[,2]),as.character(means.sem[,3]),sep="_") 
    means.sem[,hyper_columns[id_hyper]]<-NULL
    means.sem$Transformations = gsub("none_|Wrapper|Loss|Categorical", "", means.sem$Transformations)
    means.sem$Transformations = gsub("Percentage", "TotalSumNormalization", means.sem$Transformations)
    id_hyper=13
  }
  means[,c(hyper_columns[id_hyper])] <- factor(means[,c(hyper_columns[id_hyper])])      
  means.sem[,c(hyper_columns[id_hyper])] <- factor(means.sem[,c(hyper_columns[id_hyper])])
  
  # Plot
  p <- ggplot(data = means, aes_string(x = 'mean', y = c(hyper_columns[id_hyper]), group = 'variable'), colour = colsVector) + 
    geom_segment(aes_string(yend=c(hyper_columns[id_hyper])), xend=0, colour="grey75", linetype= 2) +
    geom_errorbarh(data=means.sem, aes(xmax=upper, xmin=lower), height =0.5,  colour = "grey50") + 
    geom_point(size = 2) + 
    theme_bw() +
    facet_grid(facets= . ~ variable) + #, scales="free_x") +
    theme(panel.grid.major.y = element_blank()) +
    theme(text = element_text(size=10)) +
    xlab('Mean validation') + #  'by hyperparameter value'
    ylab(gsub("\\_", " ", hyper_columns[id_hyper])) +
    theme(axis.text.x = element_text(angle = 270, hjust = 1)) +
    theme(legend.position="none")
  #pdf(paste('../Results/performance_hyperparameters_',hyper_columns[id_hyper],'.pdf',sep=''),height=3) 
  #print(p)
  #dev.off()
  return(p)
} # end-function

plots<- list()
p <- dotplot.per.hyperparameter(df,c(1:3),hyper_columns,metrics_columns)
plots=c(plots,list(p))
for (id_hyper in c(4,5,6,7,8,9,10,11)){
  p <- dotplot.per.hyperparameter(df,id_hyper,hyper_columns,metrics_columns)
  plots=c(plots,list(p))
} # end-for

# align plots: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
# p1 <- ggarrange(plots[[1]], plots[[3]], plots[[4]], ncol=1, nrow=3)
# p2 <- ggarrange(plots[[5]], plots[[6]], plots[[7]],
#                 plots[[2]], plots[[8]], plots[[9]], ncol=3, nrow=2)
#                #widths = c(2, 1))
# p <- ggarrange(p1, p2, ncol = 1, nrow = 2)

combined <- grid.arrange(plots[[1]], 
              plots[[3]], 
              plots[[4]],
              arrangeGrob(plots[[5]], plots[[6]], plots[[7]], ncol = 3),
              arrangeGrob(plots[[2]], plots[[8]], plots[[9]], ncol = 3),
              nrow = 5, heights= c(0.27,0.4,0.4,0.23,0.23))
pdf('../Results/performance_hyperparameters_all.pdf',height=9,width=7)
grid.draw(combined)
dev.off()
