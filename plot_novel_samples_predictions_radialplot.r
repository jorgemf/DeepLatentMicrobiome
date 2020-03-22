library(phyloseq)
file_tax = 'data/tax_table_all_80.csv'


physeq_pred = build_physeq_object('otus_predFromDomain_novelSamples.tsv','data/metadata_novel_samples_only5envFeatures.csv',file_tax)

# To select samples to represent in radar graph
physeq_pred_sub = prune_samples(sample_names(physeq_pred) %in% c('new00','new01'), physeq_pred)


tax_rank="Phylum"
plot_bar(physeq_pred_sub, fill=tax_rank)

physeq_rank <- tax_glom(physeq_pred_sub, taxrank = tax_rank)

# Translate taxa ID to taxa understandable name
dic_taxa_names=tax_table(physeq_rank)[,tax_rank]
for(i in seq(1,length(taxa_names(physeq_rank)))){
  id=taxa_names(physeq_rank)[i]
  taxa=unlist(strsplit(dic_taxa_names[id],'__'))[2]
  taxa_names(physeq_rank)[i]=taxa
}

df=as.data.frame(t(otu_table(physeq_rank)))
suffix=paste('_new00-01_',tax_rank)
write.table(df , paste('otus_predFromDomain_novelSamples_withTaxaNames',suffix,'.tsv',sep=''), sep="\t", row.names = TRUE, col.names = NA)

library(plotrix)
cols=c(rgb(253,174,97,maxColorValue=255),rgb(215,25,28,maxColorValue=255),rgb(171,221,164,maxColorValue=255),rgb(43,131,186,maxColorValue=255))
pdf(paste('radialPlot_otus_predFromDomain_novelSamples',suffix,'.pdf',sep=''))
radial.plot(as.matrix(df),labels=colnames(df),rp.type="p",start=1,lwd=4,line.col=cols,main='')
legend(x=-0.6,y=-0.55,legend=rownames(df),bty='n',col=cols,horiz=TRUE,lty=c(1,1),lwd=4)
dev.off()
