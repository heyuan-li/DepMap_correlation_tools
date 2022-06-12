## data table clean up and formatting ##

setwd("/Users/heyuan/Documents/GitHub/DepMap_correlation_tools/")

## CCLE files

# cell line metadata
ccle_anno = read.csv("data/sample_info.csv", header = T, row.names = 1)

# copy number
ccle_cn = read.csv("data/CCLE_gene_cn.csv", header = T, row.names = 1)
temp = colnames(ccle_cn)
colnames(ccle_cn) = sapply(temp,function(x){strsplit(x,"[.][.]")[[1]][1]})

# rnaseq
ccle_rnaseq = read.csv("data/CCLE_expression.csv", header = T, row.names = 1)
temp = colnames(ccle_rnaseq)
colnames(ccle_rnaseq) = sapply(temp,function(x){strsplit(x,"[.][.]")[[1]][1]})

# protein
ccle_protein = read.csv("data/protein_quant_current_normalized.csv", header = T, stringsAsFactors = F)
ccle_protein = ccle_protein[c(2, 49:426)]
symbols = ccle_protein$Gene_Symbol
dup = which(duplicated(symbols))
ind = sapply(dup, function(dup){
  first = match(symbols[dup], symbols)
  if (sum(is.na(as.numeric(ccle_protein[first,]))) < sum(is.na(as.numeric(ccle_protein[dup,]))))
  {dup}
  else {first}
})
ccle_protein = ccle_protein[-ind,]
symbols = ccle_protein$Gene_Symbol
ccle_protein = t(ccle_protein[,-1])

lines = sapply(rownames(ccle_protein),function(x){strsplit(x,"_")[[1]][1]})
lines_correctx = lines[grep("^X", lines)]
lines[grep("^X", lines)] = unlist(strsplit(lines_correctx, "X"))[seq(2, 2*length(lines_correctx), 2)]
DepMap_ID = sapply(lines, function(x){
  rownames(ccle_anno)[match(x, ccle_anno$stripped_cell_line_name)]
})

dup = which(duplicated(DepMap_ID))
means = t(sapply(dup, function(x){colMeans(ccle_protein[grep(DepMap_ID[x], DepMap_ID),], na.rm = T)}))
rownames(means) = DepMap_ID[dup]
ccle_protein_r = ccle_protein[-as.numeric(sapply(dup, function(x){grep(DepMap_ID[x], DepMap_ID)})),]
ccle_protein = data.frame(rbind(ccle_protein_r,means))
DepMap_ID = c(DepMap_ID[-as.numeric(sapply(dup, function(x){grep(DepMap_ID[x], DepMap_ID)}))], DepMap_ID[dup])
rownames(ccle_protein) = DepMap_ID
colnames(ccle_protein) = symbols

# GSVA output
gsva <- read.table("data/gsva_ccle_output.txt", row.names = 1)
colnames(gsva) = sapply(colnames(gsva), function(str){str_replace(str, "[.]","-")})
gsva = data.frame(t(gsva))

## Dependency data

# CRISPR chronos score
crispr_chronos = read.csv("data/CRISPR_gene_effect.csv", header = T, row.names = 1)
temp = colnames(crispr_chronos)
colnames(crispr_chronos) = sapply(temp,function(x){strsplit(x,"[.][.]")[[1]][1]})

# shRNA DEMETER2 score
shrna_d2 = read.csv("data/D2_combined_gene_dep_scores.csv", header = T, row.names = 1)
symbols = sapply(rownames(shrna_d2),function(x){strsplit(x," [(]")[[1]][1]})
lines = colnames(shrna_d2)
lines_correctx = lines[grep("^X", lines)]
lines[grep("^X", lines)] = unlist(strsplit(lines_correctx, "X"))[seq(2, 2*length(lines_correctx), 2)]
DepMap_ID = sapply(lines, function(x){rownames(ccle_anno)[match(x, ccle_anno$CCLE_Name)]})
DepMap_ID[which(is.na(DepMap_ID))] = rownames(ccle_anno)[c(1004, 1744, 1071, 618, 1419)]
shrna_d2 = data.frame(t(shrna_d2)[c(-301,-356),], row.names = DepMap_ID[c(-301,-356)])
colnames(shrna_d2) = symbols

# shRNA RSA score
drive_rsa = read.table("data/RSA_output_compiled_r.txt", header = T, stringsAsFactors = F)
genes = drive_rsa$gene
names_rsa = c(unlist(strsplit(colnames(drive_rsa)[2:9], "X"))[c(2,4,6,8,10,12,14,16)], 
              colnames(drive_rsa)[10:ncol(drive_rsa)])
drive_rsa = as.data.frame(t(drive_rsa[,-1]))
colnames(drive_rsa) = genes
DepMap_ID = sapply(toupper(names_rsa), function(x){
  rownames(ccle_anno)[match(x, ccle_anno$stripped_cell_line_name)]
})
drive_rsa = cbind(DepMap_ID, drive_rsa)
drive_rsa = transform(drive_rsa, row.names = 1)


## save objects

save(ccle_anno, ccle_cn, ccle_protein, ccle_rnaseq, crispr_chronos, drive_rsa, shrna_d2, gsva, file = "data/global.RData")
