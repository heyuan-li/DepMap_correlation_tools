## functions for DepMap analysis and visualizations ##

# required libs
library(ggplot2)
library(ggpubr)
library(pheatmap)

## calculating correlations ##
single_cor = function(gene, data1, data2, method = "pearson", split_tissue = F){
  # This function calculate correlations between gene of interest in data1 and all the genes in data2
  
  # test existence of genes in data1
  if (!(gene %in% colnames(data1))) stop('Gene not existing in data1')
  
  # overlapping cell lines in two datasets
  cl = intersect(rownames(data1), rownames(data2))
  if(length(cl) < 10) warning('number of cell lines less than 10')
  
  
  # correlation tests
  if(split_tissue){
    f = table(ccle_anno[cl,"primary_disease"])
    tt = names(f)[f>=10]
    
    ls = lapply(tt, function(tissue) {
      subcl = cl[ccle_anno[cl,"primary_disease"] == tissue]
      d = data1[subcl, gene]
      c = t(sapply(data2[subcl,], function(x){
        cor = cor.test(d, x, method = method)
        c(cor$estimate, cor$p.value,cor$parameter + 2)
      }))
      colnames(c) = c("cor","pval","n")
      c = data.frame(c, padj = p.adjust(c[,"pval"], method = "BH"),
                     "tissue_type" = tissue, 
                     "gene" = rownames(c), 
                     row.names = seq(1:nrow(c)))[,c(6,1:2,4,5,3)]
    })
    
    c = Reduce(rbind, ls)
  }
  
  else {
    d = data1[cl, gene]
    c = t(sapply(data2[cl,], function(x){
      cor = cor.test(d, x, method = method)
      c(cor$estimate, cor$p.value, cor$parameter + 2)
    }))
    colnames(c) = c("cor","pval","n")
    c = data.frame(c, padj = p.adjust(c[,"pval"], method = "BH"),
                   "tissue_type" = "Pan-cancer", 
                   "gene" = rownames(c), 
                   row.names = seq(1:nrow(c)))[,c(6,1:2,4,5,3)]
  }
  
  return(c)
}
group_cor = function(genelist, data1, data2,method = "pearson"){
  # test existence of genes in data1
  flag = genelist %in% colnames(data1)
  if (sum(flag) < length(genelist)) warning('Not all genes existing in data1, ignoring')
  
  # overlapping cell lines in two datasets
  cl = intersect(rownames(data1), rownames(data2))
  if(length(cl) < 10) warning('number of cell lines less than 10')
  
  # correlation test
  c = lapply(genelist[flag], function(gene){
    d = data1[cl, gene]
    sapply(data2[cl,], function(x){cor(d,x,use = "na.or.complete", method = method)})
  })
  c = Reduce(cbind,c)
  colnames(c) = genelist[flag]
  return(c)
}
group_cor_limited = function(genelist1, genelist2, data1, data2, method = "pearson"){
  # test existence of genes in data1
  flag = genelist1 %in% colnames(data1)
  if (sum(flag) < length(genelist1)) warning('Not all genes existing in data1, ignoring')
  
  # overlapping cell lines in two datasets
  cl = intersect(rownames(data1), rownames(data2))
  if(length(cl) < 10) warning('number of cell lines less than 10')
  
  # correlation test
  c = lapply(genelist1[flag], function(gene){
    d = data1[cl, gene]
    sapply(data2[cl,genelist2], function(x){cor(d,x,use = "na.or.complete", method = method)})
  })
  c = Reduce(cbind,c)
  colnames(c) = genelist1[flag]
  return(c)
}

## plotting ##


tissue_dep = function(genelist, data, dep_cutoff = -0.5){
  # This function plot a heatmap of fraction of cell lines that shows dependency for the genelist of interest for each tissue type, the default dependency cutoff is -0.5 and can be customized
  
  # test existence of genes in data
  flag = genelist %in% colnames(data)
  if (sum(flag) < length(genelist)) warning('Not all genes existing in data, ignoring')
  genelist = genelist[flag]
  
  cl = intersect(rownames(data), rownames(ccle_anno))
  f = table(ccle_anno[cl,"primary_disease"])
  tt = names(f)[f>=10]
  
  # framing #
  ds = matrix(ncol = length(tt), nrow = length(genelist))
  colnames(ds) = tt
  rownames(ds) = genelist
  
  # absolute threshold for dependency score #
  cutoff = dep_cutoff
  
  for (i in 1:length(tt)) {
    d = data[cl[ccle_anno[cl,"primary_disease"] == tt[i]],]
    ds[,i] = sapply(genelist, function(x){sum (na.omit(d[[x]]) < cutoff) / length(na.omit(d[[x]]))})
  }
  
  pheatmap(ds, border_color = "white")
  
}
tissue_dep_binary = function(genelist, data, dep_cutoff = -0.5, perc_cutoff = 0.6){
  # This function is a binary heatmap version of tissue_dep(). the default dependency cutoff is -0.5, the default perc_cutoff is 0.6 and can be customized

  # test existence of genes in data
  flag = genelist %in% colnames(data)
  if (sum(flag) < length(genelist)) warning('Not all genes existing in data, ignoring')
  genelist = genelist[flag]
  
  cl = intersect(rownames(data), rownames(ccle_anno))
  f = table(ccle_anno[cl,"primary_disease"])
  tt = names(f)[f>=10]
    
  # framing #
  ds = matrix(ncol = length(tt), nrow = length(genelist))
  colnames(ds) = tt
  rownames(ds) = genelist
  
  # absolute threshold for dependency score #
  cutoff = dep_cutoff
  perc = perc_cutoff
  
  for (i in 1:length(tt)) {
    d = data[cl[ccle_anno[cl,"primary_disease"] == tt[i]],]
    ds[,i] = sapply(genelist, function(x){
      if (sum (na.omit(d[[x]]) < cutoff) / length(na.omit(d[[x]])) >= perc) {1}
      else {0}
    })
  }
  
  #ds_f = ds[which(sapply(1:nrow(ds), function(x){sum(ds[x,])}) > 0),]
  
  pheatmap(ds, 
               border_color = "white",
               #cellwidth = 5,
               #cellheight = 3,
               col = c('white', 'red'), 
               #breaks = c(min(ds), cutoff_rho, max(ds)),
               #cluster_rows=FALSE, 
               legend = F,
               #show_rownames = F,
               #width = 900, height = 600,
               #fontsize_col = 6,
               #fontsize_row = 6
               
               #filename = paste(paste('heatmap',cutoff_rho, sep = '_'), 'png', sep = '.')
  )
  
}
single_tissue = function(gene,data) {
  # check gene existence 
  if (!(gene %in% colnames(data))) stop('Gene not existing in data')
  
}
single_cor_plot = function(gene1, gene2, data1, data2, split_tissue = F){
  
  # check gene existence 
  if (!(gene1 %in% colnames(data1))) stop('Gene1 not existing in data1')
  if (!(gene2 %in% colnames(data2))) stop('Gene2 not existing in data2')
  
  # overlapping cell lines in two datasets and removing NAs
  cl = intersect(rownames(data1), rownames(data2))
  cl = cl[which(!is.na(data1[cl,gene1]) & !is.na(data2[cl,gene2]))]
  if(length(cl) < 10) warning('number of cell lines less than 10')
  
  # construct dataframe for plotting
  df = data.frame(gene1 = data1[cl,gene1], gene2 = data2[cl,gene2], "tissue_type" = ccle_anno[cl,"primary_disease"])
  
  # scatter plots
  if(split_tissue) {
    f = table(df$tissue_type)
    tt = names(f)[f>=10]
    df = df[df$tissue_type %in% tt,]
    p = ggplot(df, aes(gene1, gene2)) +
      geom_point() +
      geom_smooth(method = "lm") +
      facet_wrap(~tissue_type, scale = "free") +
      stat_cor() +
      xlab(gene1) +
      ylab(gene2) +
      theme_bw() +
      theme(#axis.title.x=element_blank(),
        #axis.text.x=element_text(angle = 45, hjust = 1),
        #axis.ticks.x=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
      )
  }
  else {
    p = ggplot(as.data.frame(df), aes(gene1, gene2)) +
      geom_point() +
      geom_smooth(method = "lm") +
      stat_cor() +
      xlab(gene1) +
      ylab(gene2) +
      theme_bw() +
      theme(#axis.title.x=element_blank(),
        #axis.text.x=element_text(angle = 45, hjust = 1),
        #axis.ticks.x=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
      )
  }
  return(p)
}
single_box = function(gene1, gene2, data1, data2, cutoff = 0.33, abs = F){
  
  # check gene existence 
  if (!(gene1 %in% colnames(data1))) stop('Gene1 not existing in data1')
  if (!(gene2 %in% colnames(data2))) stop('Gene2 not existing in data2')
  
  # overlapping cell lines in two datasets and removing NAs
  cl = intersect(rownames(data1), rownames(data2))
  cl = cl[which(!is.na(data1[cl,gene1]) & !is.na(data2[cl,gene2]))]
  if(length(cl) < 10) warning('number of cell lines less than 10')
  
  # construct dataframe for plotting
  df = data.frame("dep" = data2[cl,gene2],"group" = paste(gene1, "independent"))
  
  if(abs) {df$group[data1[cl,gene1] < cutoff] = paste(gene1, "dependent")}
  else {
    cutoff_abs = quantile(data1[cl,gene1], prob = c(cutoff, 1-cutoff))
    df$group[data1[cl,gene1] < cutoff_abs[1]] = paste(gene1, "dependent")
    df = df[data1[cl,gene1] < cutoff_abs[1] | data1[cl,gene1] > cutoff_abs[2],]
  }
  
  # boxplot
  p = ggplot(df, aes(group, dep, fill=group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5, color = "grey") +
    ylab(paste(gene2, "dependency")) +
    scale_fill_manual(values=c("red", "white")) +
    stat_compare_means(method = "wilcox.test") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle = 45, hjust = 1),
          #axis.ticks.x=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
    )
  
  return(p)
}
cor_heatmap = function(genelist1, genelist2, data1, data2){
  c = group_cor_limited(genelist1, genelist2, data1, data2, method = method)
  pheatmap(t(c),
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           breaks = seq(-max(abs(c)), max(abs(c)), length.out = 100)
  )
}
gsva_cor = function(gene, data, subset = NULL, n = 10){
  # for subset, the options includes "WP","REACTOME","KEGG"
  
  # check gene existence 
  if (!(gene %in% colnames(data))) stop('Gene not existing in data')
  
  # correlation tests
  cl = intersect(rownames(gsva), rownames(data))
  data_goi = data[cl,gene]
  gsva_goicor = apply(gsva[cl,], 2, function(x){cor(x,data_goi,use = "na.or.complete")})
  
  pdata = -sort(gsva_goicor, decreasing = T)
  pdata = pdata[grep(paste(subset, "_", sep = ""), names(pdata))]
  pdata = as.data.frame(c(head(pdata,n), tail(pdata,n)))
  colnames(pdata) = "Correlation"
  pdata = transform(pdata, "gene_set" = rownames(pdata), "color" = c(rep("neg",n), rep("pos",n)))
  
  ggplot(pdata, aes(reorder(gene_set, Correlation), Correlation, fill = color)) + 
    geom_bar(stat = "identity") + 
    coord_flip() +
    xlab("") +
    scale_fill_manual(values=c("#9eadc8", "#af3e4d")) +
    theme(legend.position = "none",
          #axis.title.x=element_blank(),
          #axis.text.x=element_text(angle = 45, hjust = 1),
          #axis.ticks.x=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
    )
  
}
plot_density_cor = function(genelist, data){
  rn = colnames(data)
  data_genelist = crispr_ceres[rn,genelist]
  gsva_glcor = sapply(genelist, function(goi){
    apply(data, 1, function(x){cor(x,data_genelist[,goi],use = "na.or.complete")})
  })
  
  gsva_glcor = transform(gsva_glcor, "rn" = rownames(gsva_glcor))
  pdata = reshape(gsva_glcor, idvar = 'gene', varying = list(1:ncol(gsva_glcor)-1),
                  v.names = 'Correlation', direction = 'long',
                  timevar = 'candidates', times = genelist,
                  new.row.names = seq(1,nrow(gsva_glcor)*ncol(gsva_glcor)))
  
  
  ggplot(pdata, aes(x = Correlation, y = candidates)) +
    geom_density_ridges(aes(fill = candidates, scale = 1.5)) +
    #scale_x_continuous(limits = c(-0.3,0.7)) +
    theme_ridges()
}

