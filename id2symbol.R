# transform all entrez id in regulon object to gene symbol
library(biomaRt)
# object is downloaded from https://figshare.com/articles/Human_breast_carcinoma_signalome_regulatory_network/695962
load('~/Dropbox/_ChrisProject/workspace/viper/regulons/brca_tcga_rnaseq851_signalomeregulon.rda')
all_genes <- unique(c(names(regul),unlist(lapply(regul,function(x)names(x$tfmode))),unlist(lapply(regul,function(x)names(x$likelihood)))))

#mapping, only genes with 1:1 mapping are kept
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("hgnc_symbol", "entrezgene"),filters = "entrezgene", values = all_genes,mart = mart)
write.csv(results,'id2symbol.csv')
genes <- results$hgnc_symbol
names(genes) <- results$entrezgene

#change gene names in regulon object
change_symbol <- function(x){
  ind_keep <- names(x$tfmode) %in% names(genes)
  x$tfmode <- x$tfmode[ind_keep]
  x$likelihood <- x$likelihood[ind_keep]
  names(x$tfmode) <- genes[names(x$tfmode)]
  x
}

regul_symbol <- regul[names(regul) %in% names(genes)]
names(regul_symbol) <- genes[names(regul_symbol)]
regul_symbol <- lapply(regul_symbol,change_symbol)
save(regul_symbol,file = 'regul_symbol.rda')
