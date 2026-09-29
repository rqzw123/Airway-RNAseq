# script to get data from airway package from Bioconductor
# code is slightly modified from bioinformagician

library(airway)

data(airway)
airway

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]

# renaming the categories
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')

write.table(sample_info, file = "data/sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "data/counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

row_data <- as.data.frame(rowData(airway))
write.table(row_data, file = "data/row_data.csv", sep = ',', col.names = T, row.names = T, quote = F)
