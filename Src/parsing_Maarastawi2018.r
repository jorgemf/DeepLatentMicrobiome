library(phyloseq)
library(biomformat)
library(hash)

##############
# FUNCTIONS
##############
phyloSubset <- function(phyloOb, Var, Value){
  remove_idx <- as.character(get_variable(phyloOb, Var)) == as.character(Value)
  filteredOb <- prune_samples(remove_idx, phyloOb)
  return(filteredOb)
}
##############


otufile="../Datasets/Maarastawi2018/OTUtable_ERP105451_taxonomy_abundances_SSU_v4.1_manualParse_ConsensusLineage_noSuperKingdom.tsv"
#otufile_prefiltered="otu_table_Maarastawi2018_allBacteria_v1.csv"

# Loading OTU table and tax table in a phyloseq object
OTU_TAX <- import_qiime(otufilename=otufile)

mapmat <- read.table('../Datasets/Maarastawi2018/metadata_SRARunTable.csv',sep=',',row.names=1,header=TRUE)
MAP = sample_data(mapmat)

data = merge_phyloseq(OTU_TAX,MAP)
# otu_table()   OTU Table:         [ 1966 taxa and 322 samples ]
# sample_data() Sample Data:       [ 322 samples by 47 sample variables ]
# tax_table()   Taxonomy Table:    [ 1966 taxa by 7 taxonomic ranks ]

#save(data,file='table_1946taxa_Maarastawi2018_phyloseqObject.RData')

######
# 0.- To complete sample_data and filter samples
# Add new variables by blocks to the phyloseq object
subject.data=phyloSubset(data,'geographic_location_.country_and.or_sea.','Italy')
subject.data=phyloSubset(subject.data,'crop_rotation','RR')
sample_data(data)[sample_names(subject.data),'pH'] <- 4.9
sample_data(data)[sample_names(subject.data),'Nmin'] <- 24.97
sample_data(data)[sample_names(subject.data),'N'] <- 0.07
sample_data(data)[sample_names(subject.data),'C'] <- 0.98
sample_data(data)[sample_names(subject.data),'C:N'] <- 12.71
sample_data(data)[sample_names(subject.data),'Corg'] <- 0.98
sample_data(data)[sample_names(subject.data),'soil_type'] <- 'loam'
sample_data(data)[sample_names(subject.data),'clay_fration'] <- 9.51
sample_data(data)[sample_names(subject.data),'water_holding_capacity'] <- 41.8

subject.data=phyloSubset(data,'geographic_location_.country_and.or_sea.','Italy')
subject.data=phyloSubset(subject.data,'crop_rotation','MM')
sample_data(data)[sample_names(subject.data),'pH'] <- 4.2
sample_data(data)[sample_names(subject.data),'Nmin'] <- 27.63
sample_data(data)[sample_names(subject.data),'N'] <- 0.06
sample_data(data)[sample_names(subject.data),'C'] <- 0.74
sample_data(data)[sample_names(subject.data),'C:N'] <- 11.11
sample_data(data)[sample_names(subject.data),'Corg'] <- 0.75
sample_data(data)[sample_names(subject.data),'soil_type'] <- 'sandy_loam'
sample_data(data)[sample_names(subject.data),'clay_fration'] <- 13.22
sample_data(data)[sample_names(subject.data),'water_holding_capacity'] <- 43.9

subject.data=phyloSubset(data,'geographic_location_.country_and.or_sea.','Philippines')
subject.data=phyloSubset(subject.data,'soil_source','IRRI')
subject.data=phyloSubset(subject.data,'crop_rotation','RR')
sample_data(data)[sample_names(subject.data),'pH'] <- 5.7
sample_data(data)[sample_names(subject.data),'Nmin'] <- 6.81
sample_data(data)[sample_names(subject.data),'N'] <- 0.14
sample_data(data)[sample_names(subject.data),'C'] <- 1.73
sample_data(data)[sample_names(subject.data),'C:N'] <- 11.89
sample_data(data)[sample_names(subject.data),'Corg'] <- 1.74
sample_data(data)[sample_names(subject.data),'soil_type'] <- 'silty_clay'
sample_data(data)[sample_names(subject.data),'clay_fration'] <- 59.56
sample_data(data)[sample_names(subject.data),'water_holding_capacity'] <- 79.9 

subject.data=phyloSubset(data,'geographic_location_.country_and.or_sea.','Philippines')
subject.data=phyloSubset(subject.data,'soil_source','IRRI')
subject.data=phyloSubset(subject.data,'crop_rotation','MR')
sample_data(data)[sample_names(subject.data),'pH'] <- 5.7
sample_data(data)[sample_names(subject.data),'Nmin'] <- 4.71
sample_data(data)[sample_names(subject.data),'N'] <- 0.15
sample_data(data)[sample_names(subject.data),'C'] <- 1.81
sample_data(data)[sample_names(subject.data),'C:N'] <- 11.94
sample_data(data)[sample_names(subject.data),'Corg'] <- 1.82
sample_data(data)[sample_names(subject.data),'soil_type'] <- 'silty_clay'
sample_data(data)[sample_names(subject.data),'clay_fration'] <- 60.17
sample_data(data)[sample_names(subject.data),'water_holding_capacity'] <- 72.7 

subject.data=phyloSubset(data,'geographic_location_.country_and.or_sea.','Philippines')
subject.data=phyloSubset(subject.data,'soil_source','Tarlac')
subject.data=phyloSubset(subject.data,'crop_rotation','RR')
sample_data(data)[sample_names(subject.data),'pH'] <- 5.8
sample_data(data)[sample_names(subject.data),'Nmin'] <- 3.89
sample_data(data)[sample_names(subject.data),'N'] <- 0.06
sample_data(data)[sample_names(subject.data),'C'] <- 0.77
sample_data(data)[sample_names(subject.data),'C:N'] <- 11.67
sample_data(data)[sample_names(subject.data),'Corg'] <- 0.77
sample_data(data)[sample_names(subject.data),'soil_type'] <- 'loam'
sample_data(data)[sample_names(subject.data),'clay_fration'] <- 10.19
sample_data(data)[sample_names(subject.data),'water_holding_capacity'] <- 60.3 

subject.data=phyloSubset(data,'geographic_location_.country_and.or_sea.','Philippines')
subject.data=phyloSubset(subject.data,'soil_source','Tarlac')
subject.data=phyloSubset(subject.data,'crop_rotation','MR')
sample_data(data)[sample_names(subject.data),'pH'] <- 5.2
sample_data(data)[sample_names(subject.data),'Nmin'] <- 4.4
sample_data(data)[sample_names(subject.data),'N'] <- 0.06
sample_data(data)[sample_names(subject.data),'C'] <- 0.95
sample_data(data)[sample_names(subject.data),'C:N'] <- 13.96
sample_data(data)[sample_names(subject.data),'Corg'] <- 0.96
sample_data(data)[sample_names(subject.data),'soil_type'] <- 'silty_loam'
sample_data(data)[sample_names(subject.data),'clay_fration'] <- 12.87
sample_data(data)[sample_names(subject.data),'water_holding_capacity'] <- 53.9 

# Remove samples from Marburg (with geoographic_location=Philippines) --> It doesn't make sense!!
remove_idx <- as.character(get_variable(data, 'soil_source')) != as.character('Marburg')
data.final <- prune_samples(remove_idx, data)
data <- data
#save(data,file='table_1946taxa_Maarastawi2018_phyloseqObject_withMetadata.RData')

data.rhiz <- phyloSubset(data,'compartment','rhizosphere')
#save(data.rhiz,file='table_1946taxa_Maarastawi2018_phyloseqObject_withMetadata_onlyRhizosphere.RData')

data <- data.rhiz

####################
### Create sequential otu ID, according to order in Walters et al. 2018 whole dataset, as the model is trained.
# To use these hash's every time I save a otu_table in .csv
# Second and futher times
load('../Datasets/dictsSortedTaxa_byWalters2018.RData')
# # Just run the first time!!!
# load('../Datasets/table_all_80_phyloseqObject.RData')
# dictWalters2Seq <- hash()
# dictSeq2Walters <- hash()
# count=1
# for(taxa in taxa_names(tSpecies)){
#   id=paste("w", formatC(count, width=3, flag="0"), sep="")
#   dictWalters2Seq[[taxa]]=id
#   dictSeq2Walters[[id]]=taxa
#   count=count+1
# }
# save(dictWalters2Seq,dictSeq2Walters,file='../Datasets/dictsSortedTaxa_byWalters2018.RData')
#df.dict=as.data.frame(values(dictWalters2Seq))
#colnames(df.dict)='sequentialTaxaID'
#write.table(df.dict,file='../Datasets/dictSortedTaxa_byWalters2018.txt', sep="\t", row.names = TRUE, col.names = TRUE)
####################

data.backup=data
for (taxLevel in colnames(tax_table(data.backup))[2:7]){
#for (taxLevel in c('Class')){
  print(taxLevel)
  # 1.- To aggregate
  tAggTL=tax_glom(data.backup,taxrank=taxLevel)
  # Save phyloseq object
  save(tAggTL,file=paste('../Datasets/Maarastawi2018/table',taxLevel,'Maarastawi2018_phyloseqObject.RData',sep='_'))
  
  # 2.- To Compare tax_table(tAggTL)
  load(paste('../Datasets/Aggregated/table_',taxLevel,'_phyloseqObject.RData',sep=''))
  tax=tax_table(tSpecies)
  taxTL=tax_table(tAggTL)
  taxaIDsToPreserveTL=c()
  taxaIDsToPreserve=c()
  taxaIDsToAdd=c()
  dictTaxaID <- hash() 
  for(taxaID in taxa_names(tAgg)){
    taxaName=unlist(strsplit(tax[taxaID,taxLevel],'__'))[2]
    taxaName=gsub('\\[|]','',taxaName)
    if(!is.na(taxaName)){
      #print(taxaName)
      if(taxaName %in% taxTL[,taxLevel]){
        # Possible improvement in the comparison, not limited to the current taxLevel, or fill-in the previous levels in Maarostaski when it goes from Species to Kingdom whitout intermediate values!!
        taxaID.TL=taxa_names(taxTL[taxTL[,taxLevel]==taxaName])
        dictTaxaID[[taxaID.TL]] = taxaID
        taxaIDsToPreserve=c(taxaIDsToPreserve,taxaID)
        taxaIDsToPreserveTL=c(taxaIDsToPreserveTL,taxaID.TL)
      }else{
        taxaIDsToAdd=c(taxaIDsToAdd,taxaID)
      } # end-if
    }else{
      taxaIDsToAdd=c(taxaIDsToAdd,taxaID)
    } # end-if-isna
  } # end-for taxa names
  print('Taxa in common:')
  print(length(taxaIDsToPreserveTL))
  print('Taxa with zero abundance:')
  print(length(taxaIDsToAdd))

  # 3.- Replace id's TL by Walters et al. 
  tAggTL.new=prune_taxa(taxaIDsToPreserveTL,tAggTL)
  for(i in 1:ntaxa(tAggTL.new)){
    idTL = taxa_names(tAggTL.new)[i]
    taxa_names(tAggTL.new)[i] = dictTaxaID[[idTL]] 
  } # end-for replace otu IDs
  
  # 4.- Add rows with 0 abundances, by creating new phyloseq object, with the size I want, to replace with zeros, and to merge_phyloseq.
  tAggTL.add=prune_taxa(taxa_names(tAggTL)[1:length(taxaIDsToAdd)],tAggTL)
  len=nsamples(tAggTL.add)
  for(i in 1:ntaxa(tAggTL.add)){
    otu_table(tAggTL.add)[i,]=rep(0,len)
  }
  taxa_names(tAggTL.add)=taxaIDsToAdd
  tAggTL.final=merge_phyloseq(tAggTL.new,tAggTL.add)
  
  # 5.- IMP! Sort by otuID.seq
  df.otu=as.data.frame(otu_table(tAggTL.final))
  otuids=rownames(df.otu)
  otuids.seq=sapply(otuids, function(x){dictWalters2Seq[[x]]})
  data=cbind(otuids.seq, otuids, df.otu)
  data=data[order(data$otuids.seq),]
  data$otuids.seq=NULL
  
  # 6.- Save .tsv  
  write.table(data, paste('../Datasets/Maarastawi2018/otu_table_',taxLevel,'_Maarastawi2018.csv',sep=''), sep="\t", row.names = FALSE, col.names = TRUE)
  df.tax=as.matrix(tax_table(tAggTL.final))
  otuids=row.names(tax_table(tAggTL.final))
  df.tax=cbind(otuids, df.tax)
  write.table(df.tax, paste('../Datasets/Maarastawi2018/tax_table_',taxLevel,'_Maarastawi2018.csv',sep=''), sep="\t", row.names = FALSE, col.names = TRUE)
  
  df.sample=as.data.frame(sample_data(tAggTL.final))
  n=names(df.sample)
  n[1]='X.SampleID'
  names(df.sample)=n
  write.table(df.sample, paste('../Datasets/Maarastawi2018/metadata_table_Maarastawi2018.csv',sep=''), sep="\t", row.names = FALSE, col.names = TRUE)
} # end-for taxLevel  
  
    
# [1] "Phylum"
# [1] "Taxa in common:" 15
# [1] "Taxa with zero abundance:" 1
# [1] "Class" 
# [1] "Taxa in common:" 26
# [1] "Taxa with zero abundance:" 19
# [1] "Order"
# [1] "Taxa in common:" 36
# [1] "Taxa with zero abundance:" 47
# [1] "Family"
# [1] "Taxa in common:" 70
# [1] "Taxa with zero abundance:" 74
# [1] "Genus"
# [1] "Taxa in common:" 84
# [1] "Taxa with zero abundance:" 138
# [1] "Species"
# [1] "Taxa in common:" 0
# [1] "Taxa with zero abundance:" 222


