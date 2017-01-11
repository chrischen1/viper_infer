library(viper)

#inport log2 fold change dataset
logFC <- read.csv("resources/logFC.csv", row.names=1, stringsAsFactors=FALSE)

# load regulon object
load("resources/regulon_symbol.rdata")

# VIPER for BT20 abe 3uM 6h
x <- logFC[,1]
names(x) <- rownames(logFC)
mrs <- msviper(x, regulon1, verbose = FALSE,minsize = 5)
result <- summary(mrs,mrs = length(mrs$es$nes))
write.csv(result,'viper_prediction.csv') # inferred activity and p values
