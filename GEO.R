options(stringsAsFactors = F)##避免将character转换为因子
library(GEOquery)
# 下载表达矩阵
gse <- getGEO("GSE7765", GSEMatrix = TRUE) 
show(gse)
class(gse)
a<-gse[[1]]
head(a)
b<-gse[[2]]
class(a)
