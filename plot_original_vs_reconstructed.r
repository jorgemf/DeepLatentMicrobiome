library(phyloseq)
library(biomformat)
library(ggplot2)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
`%not_in%` <- purrr::negate(`%in%`)

# pre-built function (public online) to use only one color legend in a composed plot
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

# Author: Sara Cabello Pinedo
plot_taxa_subset = function (data_or, data_tr, trait, trait_level, tax_rank){
  #function to plot comparative barplots between original (upper row) and transformed data by autoencoder (bottom row)
  #color according to tax_rank
  #samples are shown individually and only the ones that belong to a specific class ('trait_level') of a variable ('trait) are plotted
  
  #tax rank name definition
  ranks=c("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES", "none")
  rank_name=ifelse(grepl("2", toString(tax_rank)), "PHYLUM", ifelse(grepl("3", toString(tax_rank)),"CLASS", ifelse(grepl("4", toString(tax_rank)), "ORDER", 
                                                                                                                   ifelse(grepl("5", toString(tax_rank)), "FAMILY", ifelse(grepl("6", toString(tax_rank)), "GENUS", ifelse(grepl("7", toString(tax_rank)),"SPECIES", "OTUs"))))))
  
  if (rank_name != "OTUs") {
    #group data according to tax rank (as a parameter)
    data_or.rank <- data_or %>% 
      tax_glom (taxrank = tax_rank) 
    data_tr.rank <- data_tr %>% 
      tax_glom (taxrank = tax_rank) 
  }
  
  else{
    data_or.rank=data_or 
    data_tr.rank=data_tr 
  }
  
  #show all the possible values that the indicated variable can take
  print(unique(sample_data(data_or.rank)[,trait]))
  
  data_or.rank_df=data_or.rank %>%
    psmelt()
  
  data_tr.rank_df=data_tr.rank %>%
    psmelt()
  
  samples_id=sort(unique(data_or.rank_df$Sample))
  
  #check that every sample is in both original and transformed data
  samples2delete=NULL
  for (i in (1: length(samples_id))){
    sample=samples_id[i]
    selected_rows=data_or.rank_df[data_or.rank_df$Sample == sample,]
    suma = sum(selected_rows$Abundance)
    if (suma==0){
      samples2delete=c(samples2delete, sample)
    }
  }
  
  #select intereseting samples
  samples_subset_or=subset(data_or.rank_df, data_or.rank_df[[trait]] == trait_level)
  samples_subset_or=samples_subset_or[samples_subset_or$Sample %not_in% samples2delete,]
  samples_subset_tr=subset(data_tr.rank_df, data_tr.rank_df[[trait]] == trait_level)
  samples_subset_tr=samples_subset_tr[samples_subset_tr$Sample %in% samples_subset_or$Sample, ]
  
  #plot
  x_text=element_text(size = 0)
  #original plot
  p1 <- ggplot(samples_subset_or, aes(x = samples_subset_or$Sample, y = Abundance, fill = samples_subset_or[,tax_rank]))+
    geom_bar(aes(), stat="identity", position="fill") + 
    guides(fill = guide_legend(title=rank_name)) +
    scale_fill_hue(c=80, l=70) +
    ylab("Relative abundance") +
    #xlab("Samples")+
    ggtitle("Original data") +
    theme(axis.title.x = x_text ,axis.text.x = element_text(angle=60,hjust=1,size='8'))
  
  #transformed plot
  p2 <- ggplot(samples_subset_tr, aes(x = samples_subset_tr$Sample, y = Abundance, fill = samples_subset_tr[,tax_rank]))+
    geom_bar(aes(), stat="identity", position="fill") + 
    guides(fill = guide_legend(title=rank_name)) +
    scale_fill_hue(c=80, l=70) +
    ylab("Relative abundance") +
    xlab("Samples")+
    ggtitle("Transformed data") +
    theme(axis.text.x = element_text(angle=60,hjust=1,size='8'))
  
  
  text_title=paste(rank_name, " COMPOSITION: ", trait, " = ", trait_level, sep = "")
  title=textGrob(text_title, gp=gpar(fontface="bold", fontsize=20))
  plot=grid_arrange_shared_legend(p1, p2)
  # BGJ: Update to save in .pdf file
  pdf(paste('plot_original_vs_reconstructed',trait, trait_level,tax_rank,'.pdf',sep='_'))  
  grid.arrange(plot, top = title)
  dev.off()
}


# Author: Beatriz Garcia-Jimenez
build_physeq_object <- function(fotu,fmap,ftax){
  otumat <- read.table(fotu,sep='\t',row.names=1,header=TRUE,check.names=FALSE)
  OTU = otu_table(otumat, taxa_are_rows = TRUE)

  mapmat <- read.table(fmap,sep='\t',row.names=1,header=TRUE)
  MAP = sample_data(mapmat)
  
  taxmat <- read.table(ftax,sep='\t',row.names=1,header=TRUE)
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxmat=as.matrix(taxmat)
  TAX = tax_table(taxmat)

  phyObj = phyloseq(OTU, MAP, TAX) 
 
  return(phyObj) 
}

file_tax = 'data/tax_table_all_80.csv'
physeq_orig = build_physeq_object('data/otu_table_all_80.csv','data/metadata_table_all_80.csv',file_tax)
physeq_rec = build_physeq_object('otus_predicted_biome.tsv','data/metadata_table_all_80.csv',file_tax)

# Testing with a subset of 4500 samples
physeq_orig_subset=subset_samples(physeq_orig,Maize_Line=='Popcorn')
physeq_rec_subset=subset_samples(physeq_rec,Maize_Line=='Popcorn')


#plot example

plot_taxa_subset(physeq_orig_subset, physeq_rec_subset,  'Maize_Line',  'Popcorn', 'Phylum')
