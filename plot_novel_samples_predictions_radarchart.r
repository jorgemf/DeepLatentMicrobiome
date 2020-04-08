library(fmsb)
library(phyloseq)
library(extrafont)

####
# Taken from: https://stackoverflow.com/questions/54185029/change-labels-colors-in-r-radarchart
# Adapted to avoid overlap in top and bottom of peripherical labels
radarchart2 <- function (df, axistype = 0, seg = 4, pty = 16, pcol = 1:8, plty = 1:6, 
                         plwd = 1, pdensity = NULL, pangle = 45, pfcol = NA, cglty = 3, 
                         cglwd = 1, cglcol = "navy", axislabcol = "blue", vlabcol = "black", title = "", 
                         maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlabels = NULL, 
                         vlcex = NULL, caxislabels = NULL, calcex = NULL, paxislabels = NULL, 
                         palcex = NULL, ...) 
{
  if (!is.data.frame(df)) {
    cat("The data must be given as dataframe.\n")
    return()
  }
  if ((n <- length(df)) < 3) {
    cat("The number of variables must be 3 or more.\n")
    return()
  }
  if (maxmin == FALSE) {
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type = "n", frame.plot = FALSE, 
       axes = FALSE, xlab = "", ylab = "", main = title, asp = 1, 
       ...)
  theta <- seq(90, 450, length = n + 1) * pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) {
    polygon(xx * (i + CGap)/(seg + CGap), yy * (i + CGap)/(seg + 
                                                             CGap), lty = cglty, lwd = cglwd, border = cglcol)
    if (axistype == 1 | axistype == 3) 
      CAXISLABELS <- paste(i/seg * 100, "(%)")
    if (axistype == 4 | axistype == 5) 
      CAXISLABELS <- sprintf("%3.2f", i/seg)
    if (!is.null(caxislabels) & (i < length(caxislabels))) 
      CAXISLABELS <- caxislabels[i + 1]
    if (axistype == 1 | axistype == 3 | axistype == 4 | 
        axistype == 5) {
      if (is.null(calcex)) 
        text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS, 
             col = axislabcol)
      else text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS, 
                col = axislabcol, cex = calcex)
    }
  }
  if (centerzero) {
    arrows(0, 0, xx * 1, yy * 1, lwd = cglwd, lty = cglty, 
           length = 0, col = cglcol)
  }
  else {
    arrows(xx/(seg + CGap), yy/(seg + CGap), xx * 1, yy * 
             1, lwd = cglwd, lty = cglty, length = 0, col = cglcol)
  }
  PAXISLABELS <- df[1, 1:n]
  if (!is.null(paxislabels)) 
    PAXISLABELS <- paxislabels
  if (axistype == 2 | axistype == 3 | axistype == 5) {
    if (is.null(palcex)) 
      text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol)
    else text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol, 
              cex = palcex)
  }
  #
  #print("theta:")
  #print(theta)
  #print("xx:")
  #print(xx)
  #print("yy:")
  #print(yy)
  #
  VLABELS <- colnames(df)
  polars1=(abs(yy)==1)
  polars2=((abs(yy)>=9.8e-01) & (abs(yy)!=1))
  notPolars=!as.logical(polars1 + polars2)
  if (!is.null(vlabels)) 
    VLABELS <- vlabels
  if (is.null(vlcex)){ # BGJ: 2020.04.07 (inspired by https://stackoverflow.com/questions/43403700/controlling-srt-in-text-to-work-selectively-in-r)
    text(1.2 * xx[polars1], 1.2 * (yy[polars1]+(yy[polars1]*0.05)), VLABELS[polars1], col = vlabcol, font=4, family='Times')
    if(sum(polars2)>0)
      text(1.2 * (xx[polars2]+xx[polars2]*0.25), 1.2 * (yy[polars2]+(yy[polars2]*0.02)), VLABELS[polars2], col = vlabcol, font=4, family='Times') # if rotate: ,srt=15
    text(1.2 * xx[notPolars], 1.2 * yy[notPolars], VLABELS[notPolars], col = vlabcol, font=4, family='Times')
    # text(xx * 1.2, yy * 1.2, VLABELS, col = vlabcol) # original
  }else{ #BGJ: 2020.04.07
    text(1.2 * xx[polars1], 1.2 * (yy[polars1]+(yy[polars1]*0.05)), VLABELS[polars1], cex = vlcex, col = vlabcol, font=4, family='Times')
    if(sum(polars2)>0)
      text(1.2 * (xx[polars2]+xx[polars2]*0.25), 1.2 * (yy[polars2]+(yy[polars2]*0.02)), VLABELS[polars2], cex = vlcex, col = vlabcol, font=4, family='Times') # if rotate: ,srt=15
    text(1.2 * xx[notPolars], 1.2 * yy[notPolars], VLABELS[notPolars], cex = vlcex, col = vlabcol, font=4, family='Times')
    #text(xx * 1.2, yy * 1.2, VLABELS, cex = vlcex, col = vlabcol)  # original
  }
  series <- length(df[[1]])
  SX <- series - 2
  if (length(pty) < SX) {
    ptys <- rep(pty, SX)
  }
  else {
    ptys <- pty
  }
  if (length(pcol) < SX) {
    pcols <- rep(pcol, SX)
  }
  else {
    pcols <- pcol
  }
  if (length(plty) < SX) {
    pltys <- rep(plty, SX)
  }
  else {
    pltys <- plty
  }
  if (length(plwd) < SX) {
    plwds <- rep(plwd, SX)
  }
  else {
    plwds <- plwd
  }
  if (length(pdensity) < SX) {
    pdensities <- rep(pdensity, SX)
  }
  else {
    pdensities <- pdensity
  }
  if (length(pangle) < SX) {
    pangles <- rep(pangle, SX)
  }
  else {
    pangles <- pangle
  }
  if (length(pfcol) < SX) {
    pfcols <- rep(pfcol, SX)
  }
  else {
    pfcols <- pfcol
  }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap/(seg + CGap) + (df[i, ] - df[2, ])/(df[1, 
                                                         ] - df[2, ]) * seg/(seg + CGap)
    if (sum(!is.na(df[i, ])) < 3) {
      cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n", i, 
                  df[i, ]))
    }
    else {
      for (j in 1:n) {
        if (is.na(df[i, j])) {
          if (na.itp) {
            left <- ifelse(j > 1, j - 1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left > 1, left - 1, n)
            }
            right <- ifelse(j < n, j + 1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right < n, right + 1, 
                              1)
            }
            xxleft <- xx[left] * CGap/(seg + CGap) + 
              xx[left] * (df[i, left] - df[2, left])/(df[1, 
                                                         left] - df[2, left]) * seg/(seg + CGap)
            yyleft <- yy[left] * CGap/(seg + CGap) + 
              yy[left] * (df[i, left] - df[2, left])/(df[1, 
                                                         left] - df[2, left]) * seg/(seg + CGap)
            xxright <- xx[right] * CGap/(seg + CGap) + 
              xx[right] * (df[i, right] - df[2, right])/(df[1, 
                                                            right] - df[2, right]) * seg/(seg + 
                                                                                            CGap)
            yyright <- yy[right] * CGap/(seg + CGap) + 
              yy[right] * (df[i, right] - df[2, right])/(df[1, 
                                                            right] - df[2, right]) * seg/(seg + 
                                                                                            CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft
              yytmp <- yyleft
              xxleft <- xxright
              yyleft <- yyright
              xxright <- xxtmp
              yyright <- yytmp
            }
            xxs[j] <- xx[j] * (yyleft * xxright - yyright * 
                                 xxleft)/(yy[j] * (xxright - xxleft) - 
                                            xx[j] * (yyright - yyleft))
            yys[j] <- (yy[j]/xx[j]) * xxs[j]
          }
          else {
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j] * CGap/(seg + CGap) + xx[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, 
                                                 j]) * seg/(seg + CGap)
          yys[j] <- yy[j] * CGap/(seg + CGap) + yy[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, 
                                                 j]) * seg/(seg + CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i - 
                                                            2], border = pcols[i - 2], col = pfcols[i - 
                                                                                                      2])
      }
      else {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i - 
                                                            2], border = pcols[i - 2], density = pdensities[i - 
                                                                                                              2], angle = pangles[i - 2], col = pfcols[i - 
                                                                                                                                                         2])
      }
      points(xx * scale, yy * scale, pch = ptys[i - 2], 
             col = pcols[i - 2])
    }
  }
}
####

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
###

plot_radarchart_several_taxaLevels <- function(physeq_pred_sub,suffix0){
  for(tax_rank in c("Phylum","Class")){
    suffix=paste(suffix0,tax_rank,sep='_')
    #plot_bar(physeq_pred_sub, fill=tax_rank)
    physeq_rank <- tax_glom(physeq_pred_sub, taxrank = tax_rank)
    
    # Translate taxa ID to taxa understandable name
    dic_taxa_names=tax_table(physeq_rank)[,tax_rank]
    for(i in seq(1,length(taxa_names(physeq_rank)))){
      id=taxa_names(physeq_rank)[i]
      #taxa=unlist(strsplit(dic_taxa_names[id],'__'))[2]
      taxa=gsub('bacteria-','bact.',dic_taxa_names[id])
      #taxa=dic_taxa_names[id]
      taxa_names(physeq_rank)[i]=taxa
    }
    
    # # Define color scale at Phylum level
    # colorsPhylum = hue_pal(c=80,l=70)(16)
    # taxaPhylum = c(sort(unique(tax_table(physeq_rank)[,"Phylum"])))
    # # HEREEEEEEEEEEEEE
    # #https://stackoverflow.com/questions/37862979/create-color-palette-function-from-named-list-or-vector
    # # HEREE: Averiguar como generar esta lista a partir de pares <taxaPhylum,colorsPhylum>
    # col_universe <- list(Acidobacteria = colorsPhylum[1], etc.)
    #      #"Actinobacteria"   "Armatimonadetes"  "Bacteroidetes"    "Chloroflexi" "Crenarchaeota"    "Elusimicrobia"    "Fibrobacteres"    "Firmicutes"       "Gemmatimonadetes" "Nitrospirae"      "Planctomycetes"   "Proteobacteria"   "Verrucomicrobia"  "WS3"   
    # #col_universe <- list(dark_blue = "#034772", med_blue = "#2888BC", light_blue = "#73B7CE", green = "#699D46", orange = "#EA8936", gold = "#F9C347", dark_grey = "#58595B", medium_grey = "#7D7E81",light_grey = "#C1C2C4")
    # pal1 <- c(tax_table(physeq_rank)[,"Phylum"])
    # vlabcol_values=c()
    # vlabcol_values=c(vlabcol_values,unname(col_universe[pal1]))

    #########
    
    df=as.data.frame(t(otu_table(physeq_rank)))
    df=df[, order(names(df))]
    write.table(df , paste('otus_predFromDomain_novelSamples_withTaxaNames',suffix,'.tsv',sep=''), sep="\t", row.names = TRUE, col.names = NA)
    
    # Plot Radar Chart
    # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
    min=0
    max=round(max(df),2)
    max_perc=max*100
    max_axis=(max_perc+(5-(max_perc%%5)))/100 # To avoid rare numbers in graph
    labels_central_axis=seq(max_axis,0,-(round(max_axis/4,2)))
    labels_central_axis[5]=round(0,0)
    data <- rbind(rep(min,length(df)) , rep(max_axis,length(df)) , df)

    colors_v=c(rgb(0.8,0.2,0.5,0.9), rgb(0.2,0.5,0.5,0.9), rgb(0.7,0.5,0.1,0.9))
    pdf(paste('radarChart_otus_predFromDomain_novelSamples',suffix,'.pdf',sep=''))
    par(xpd=NA)         # Allow plotting outside the plot region
    #par(family='Times',font=3,font.main=2,font.axis=1)
    radarchart2(data,
               axistype=1,
               #custom polygon
               pcol=colors_v,
               plwd=3, plty=3, pty=19,
               #custom the grid
               cglcol="grey88", cglty=1, axislabcol="black",
               # Reverse axis labeling
               caxislabels=labels_central_axis,   
               cglwd=0.8,
               vlcex=0.9,
               # custom color labels: radarchart2
               # TODO: according to higher taxonomic group
               #vlabcol=c(rep('red',15),rep('blue',15),rep('green',15),rep('purple',15)),
               title='Relative abundances in different environmental conditions'
    )
    legend(x=-1, y=1.48, legend = c('actual','hot and dry','cold and wet'), bty = "n", horiz=TRUE, pch=20 , col=colors_v, cex=1.2, pt.cex=3)
    dev.off()
    embed_fonts(paste('radarChart_otus_predFromDomain_novelSamples',suffix,'.pdf',sep=''))
  } # end-for
} # end-function
###





file_tax = 'data/tax_table_all_80_cleanNames.csv'
physeq_pred = build_physeq_object('otus_predFromDomain_novelSamples.tsv','data/metadata_novel_samples_only3envFeatures.csv',file_tax)

# To select samples to represent in radar graph
physeq_pred_sub = prune_samples(sample_names(physeq_pred) %in% c('new00','new01','new02'), physeq_pred)
suffix0='_new00-01-02_age01'

physeq_pred_sub = prune_samples(sample_names(physeq_pred) %in% c('new03','new04','new05'), physeq_pred)
suffix0='_new03-04-05_age10'

plot_radarchart_several_taxaLevels(physeq_pred_sub,suffix0)


# Alternative: To rotate labels: ggradar: ggplot extension for radar graph: https://rpubs.com/updragon/ggradar