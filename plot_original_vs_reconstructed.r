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

# Author: Sara Cabello Pinedo. Adapted by Beatriz Garcia Jimenez
plot_taxa_subset = function (data_or, data_tr, data_pred, variable, value, tax_rank){
  #function to plot comparative barplots between original (middle row), reconstructed by autoencoder (top row) and predicted microbial composition by environmental features (bottom row)
  #color according to tax_rank
  #samples are shown individually and only the ones that belong to a specific value of a variable are plotted
  
  if (tax_rank != "Species") {
    #group data according to tax rank (as a parameter)
    data_or.rank <- data_or %>% 
      tax_glom(taxrank = tax_rank) 
    data_tr.rank <- data_tr %>% 
      tax_glom(taxrank = tax_rank) 
    data_pred.rank <- data_pred %>% 
      tax_glom(taxrank = tax_rank)
  }else{
    data_or.rank=data_or 
    data_tr.rank=data_tr 
    data_pred.rank=data_pred
  }
  
  data_or.rank_df=data_or.rank %>%
    psmelt()
  data_tr.rank_df=data_tr.rank %>%
    psmelt()
  data_pred.rank_df=data_pred.rank %>%
    psmelt()
  
  samples_id=sort(unique(data_or.rank_df$Sample))
  
  #select subset of samples
  samples_subset_or=subset(data_or.rank_df, data_or.rank_df[[variable]] == value)
  samples_subset_tr=subset(data_tr.rank_df, data_tr.rank_df[[variable]] == value)
  samples_subset_pred=subset(data_pred.rank_df, data_pred.rank_df[[variable]] == value)
  
  #plot
  # reconstructed 
  p1 <- ggplot(samples_subset_tr, aes(x = samples_subset_tr$Sample, y = Abundance, fill = samples_subset_tr[,tax_rank]))+
    geom_bar(aes(), stat="identity", position="fill") + 
    guides(fill = guide_legend(title=tax_rank)) +
    scale_fill_hue(c=80, l=70) +
    ylab("") +
    #xlab("Samples")+
    ggtitle("Reconstructed") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  # original 
  p2 <- ggplot(samples_subset_or, aes(x = samples_subset_or$Sample, y = Abundance, fill = samples_subset_or[,tax_rank]))+
    geom_bar(aes(), stat="identity", position="fill") + 
    guides(fill = guide_legend(title=tax_rank)) +
    scale_fill_hue(c=80, l=70) +
    ylab("Relative abundance") +
    #xlab("Samples")+
    ggtitle("Original") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  # predicted
  p3 <- ggplot(samples_subset_pred, aes(x = samples_subset_pred$Sample, y = Abundance, fill = samples_subset_pred[,tax_rank]))+
    geom_bar(aes(), stat="identity", position="fill") + 
    guides(fill = guide_legend(title=tax_rank)) +
    scale_fill_hue(c=80, l=70) +
    ylab("") +
    xlab("Samples")+
    ggtitle("Predicted from environment") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
      #axis.text.x=element_text(angle=60,hjust=1,size='8'))
  
  text_title=paste(tax_rank, " composition: ", variable, " = ", value, sep = "")
  title=textGrob(text_title, gp=gpar(fontface="bold", fontsize=20))
  plot=grid_arrange_shared_legend(p1, p2, p3)
  # BGJ: Update to save in .pdf file
  pdf(paste('barplot_original_vs_reconstructed_vs_predicted',variable, value,tax_rank,'.pdf',sep='_'))  
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
file_meta = 'data/metadata_table_all_80.csv'
# Train samples
#physeq_orig = build_physeq_object('data/otu_table_all_80.csv',file_meta,file_tax)
#physeq_rec = build_physeq_object('otus_predicted_biome.tsv',file_meta,file_tax)
# Test samples
physeq_orig = build_physeq_object('otus_original_test.tsv',file_meta,file_tax)
physeq_rec = build_physeq_object('otus_reconstAEfromBiome.tsv',file_meta,file_tax)
physeq_pred = build_physeq_object('otus_predFromDomain.tsv',file_meta,file_tax)

#plot example
#plot_taxa_subset(physeq_orig, physeq_rec,  'Maize_Line',  'Popcorn', 'Phylum')
for (var in c('INBREDS','Maize_Line')){
  for (value in levels(get_variable(physeq_orig,var))){
    plot_taxa_subset(physeq_orig, physeq_rec, physeq_pred, var, value, 'Phylum')
  }
}


#'age','Temperature','Precipitation3Days'

# Testing with a subset of samples
#physeq_orig_subset=subset_samples(physeq_orig,Maize_Line=='Popcorn')
#physeq_rec_subset=subset_samples(physeq_rec,Maize_Line=='Popcorn')

#plot example
#plot_taxa_subset(physeq_orig_subset, physeq_rec_subset,  'Maize_Line',  'Popcorn', 'Phylum')


#plot_bar(physeq_orig_subset, fill="Phylum")
#plot_bar(physeq_rec_subset, fill="Phylum")

