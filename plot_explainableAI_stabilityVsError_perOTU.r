library(ggplot2)
library(RColorBrewer)
myPalette <- brewer.pal(12, "Set3")

tax_rank="Phylum"
file_error_perOTU = 'errors_perOTU.tsv'
file_tax = 'data/tax_table_all_80.csv'

taxmat <- read.table(file_tax,sep='\t',row.names=1,header=TRUE)
colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

### Function name_from_OTU_id ###
name_from_OTU_id <- function(dic, id, tax_level){
  return(unlist(strsplit(as.character(dic[id,tax_level]),'__'))[2])
}
###
  
# Read table with one row per OTU, and different columns with different errors (RSD, RSE, RRSE)
error_perOTU <- read.table(file_error_perOTU,sep='\t',row.names=1,header=TRUE,check.names=FALSE)
  
# Add taxonomic info to each OTU
for(level in colnames(taxmat)){
  for (id in row.names(error_perOTU)){
    error_perOTU[id,level] = name_from_OTU_id(taxmat,id,level)
  }
}

# Sort df to obtain best predicted OTUs
error_perOTU_sorted <- error_perOTU[order(error_perOTU$RRSE),]
write.table(error_perOTU_sorted, 'errors_perOTU_withTaxaNames.tsv', sep='\t', row.names = TRUE, col.names = NA)


# Compute % of each Phylum in all 717 OTUs
numOTUs=nrow(error_perOTU)
df_perc_OTU = data.frame(perc=double())
for(phy in sort(unique(error_perOTU$Phylum))){
  perc=nrow(subset(error_perOTU[error_perOTU$Phylum==phy,]))/numOTUs
  df_perc_OTU[phy,'perc']=perc
  print(paste(phy,perc))
}
df_perc_OTU
#                  perc
# Acidobacteria    0.092050209
# Actinobacteria   0.188284519
# Armatimonadetes  0.001394700
# Bacteroidetes    0.133891213
# Chlorobi         0.001394700
# Chloroflexi      0.033472803
# Crenarchaeota    0.005578801
# Elusimicrobia    0.001394700
# Fibrobacteres    0.001394700
# Firmicutes       0.013947001
# Gemmatimonadetes 0.026499303
# Nitrospirae      0.006973501
# Planctomycetes   0.022315202
# Proteobacteria   0.391910739
# Verrucomicrobia  0.075313808
# WS3              0.004184100
pdf('piechart_perc_all_OTUs.pdf')
pie(df_perc_OTU$perc,labels=row.names(df_perc_OTU),col=myPalette,main='OTUs distribution (by Phylum)')
dev.off()


# Compute 5% best predicted OTUs
error_perOTU_sorted_by_RRSE <- error_perOTU[order(error_perOTU$RRSE),]
#717*0.05=35.8
best_predicted_OTUs=head(error_perOTU_sorted_by_RRSE[,c('RRSE','Phylum')],35)
best_predicted_OTUs[order(best_predicted_OTUs$Phylum),]
# Percentage of phylums with OTUs with the best 5% stability
numOTUs=35
df_best_predicted = data.frame(perc=double())
for(phy in sort(unique(best_predicted_OTUs$Phylum))){
  perc=nrow(subset(best_predicted_OTUs[best_predicted_OTUs$Phylum==phy,]))/numOTUs
  df_best_predicted[phy,'perc']=perc
  #print(paste(phy,perc))
}
df_best_predicted
# ==> Acidobacteria, Actinobacteria and Proteobacteria are the best predicted, but in a different order

#                 perc
# Acidobacteria   0.22857143
# Actinobacteria  0.31428571
# Bacteroidetes   0.05714286
# Chloroflexi     0.02857143
# Planctomycetes  0.05714286
# Proteobacteria  0.22857143
# Verrucomicrobia 0.08571429    

# OLD values: 
#                 perc
# Acidobacteria   0.08571429
# Actinobacteria  0.20000000
# Bacteroidetes   0.05714286
# Nitrospirae     0.02857143
# Planctomycetes  0.08571429
# Proteobacteria  0.45714286
# Verrucomicrobia 0.08571429


pdf('piechart_5perc_best_predicted_OTUs.pdf')
pie(df_best_predicted$perc,labels=row.names(df_best_predicted),col=myPalette,main='5% best predicted OTUs (by Phylum)')
dev.off()



##################################
# If Relative Standard Deviation is computed (RSD), with 5CV model
# Plot scatterplot
ggplot(error_perOTU, aes(x=RRSE, y=RSD, color=Phylum)) +
  geom_point(size=2) #+
  #geom_smooth(method=lm, aes(fill=Phylum))


# Check Relative Standard Deviation values
error_perOTU_sorted_by_RSD <- error_perOTU[order(error_perOTU$RSD),]
head(error_perOTU_sorted_by_RSD[,c('RSD','RRSE')])
min(error_perOTU[,'RSD'])
#[1] 0.2711851
max(error_perOTU[,'RSD'])
#[1] 1.873393
mean(error_perOTU[,'RSD'])
#[1} 0.7035187
pdf('hist_RSD_5CV_perOTU.pdf')
hist(error_perOTU[,'RSD'])
dev.off()

# Compute 5% most stable OTUs
#717*0.05=35.8
best_stable_OTUs=head(error_perOTU_sorted_by_RSD[,c('RSD','Phylum')],35)
best_stable_OTUs[order(best_stable_OTUs$Phylum),]
# Percentage of phylums with OTUs with the best 5% stability
numOTUs=35
df_best_stable = data.frame(perc=double())
for(phy in sort(unique(best_stable_OTUs$Phylum))){
  perc=nrow(subset(best_stable_OTUs[best_stable_OTUs$Phylum==phy,]))/numOTUs
  df_best_stable[phy,'perc']=perc
  #print(paste(phy,perc))
}
df_best_stable
        # [1] "Acidobacteria 0.0285714285714286"
        # [1] "Actinobacteria 0.2"
        # [1] "Bacteroidetes 0.0571428571428571"
        # [1] "Chloroflexi 0.0571428571428571"
        # [1] "Crenarchaeota 0.0571428571428571"
        # [1] "Planctomycetes 0.0571428571428571"
        # [1] "Proteobacteria 0.4"
        # [1] "Verrucomicrobia 0.142857142857143"
pdf('piechart_5perc_most_stable_OTUs.pdf')
pie(df_best_stable$perc,labels=row.names(df_best_stable),col=myPalette,main='5% most stable OTUs (by Phylum)')
dev.off()


