library(phyloseq)

otumat <- read.table('../Datasets/otu_table_all_80.csv',sep='\t',row.names=1,header=TRUE,check.names=FALSE)
OTU = otu_table(otumat, taxa_are_rows = TRUE)

mapmat <- read.table('../Datasets/metadata_table_all_80.csv',sep='\t',row.names=1,header=TRUE)
MAP = sample_data(mapmat)


taxmat <- read.table('../Datasets/tax_table_all_80.csv',sep='\t',row.names=1,header=TRUE)
colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat=as.matrix(taxmat)
TAX = tax_table(taxmat)

tSpecies = phyloseq(OTU, MAP, TAX) 
save(tSpecies,file='table_all_80_phyloseqObject.RData')

for (taxLevel in colnames(tax_table(tSpecies))[2:6]){
  print(taxLevel)
  tAgg=tax_glom(tSpecies,taxrank=taxLevel)
  # Save phyloseq object
  save(tAgg,file=paste('../Datasets/Aggregated/table',taxLevel,'phyloseqObject.RData',sep='_'))
  # Save .tsv
  df.otu=as.data.frame(otu_table(tAgg))
  otuids=rownames(df.otu)
  data=cbind(otuids, df.otu)
  write.table(data, paste('../Datasets/Aggregated/otu_table_',taxLevel,'.csv',sep=''), sep="\t", row.names = FALSE, col.names = TRUE)
  df.tax=as.matrix(tax_table(tAgg))
  print('after')
  otuids=row.names(tax_table(tAgg))
  df.tax=cbind(otuids, df.tax)
  write.table(df.tax, paste('../Datasets/Aggregated/tax_table_',taxLevel,'.csv',sep=''), sep="\t", row.names = FALSE, col.names = TRUE)
  #df.sample=as.data.frame(sample_data(tAgg))
  #write.table(df.sample, '../Datasets/Aggregated/metadata_table_all_80.csv', sep="\t", row.names = FALSE, col.names = TRUE)
}
