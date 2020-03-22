library(ggplot2)

tax_rank="Phylum"
file_error_perOTU = 'errors_perOTU.tsv'
file_tax = 'data/tax_table_all_80.csv'

taxmat <- read.table(file_tax,sep='\t',row.names=1,header=TRUE)
colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


# Read table with one row per OTU, and different columns with different errors (RSD, RSE, RRSE)
error_perOTU <- read.table(file_error_perOTU,sep='\t',row.names=1,header=TRUE,check.names=FALSE)

# Add taxonomic info to each OTU
name_from_OTU_id <- function(dic, id, tax_level){
  return(unlist(strsplit(as.character(dic[id,tax_level]),'__'))[2])
}

for (id in row.names(error_perOTU)){
  error_perOTU[id,tax_rank] = name_from_OTU_id(taxmat,id,tax_rank)
}

# Plot scatterplot
ggplot(error_perOTU, aes(x=RRSE, y=RSD, color=Phylum)) +
  geom_point(size=2) #+
  #geom_smooth(method=lm, aes(fill=Phylum))