# 构建假数据
data.matrix <- matrix(nrow=100, ncol=10)
colnames(data.matrix) <- c(
  paste("wt", 1:5, sep=""),
  paste("ko", 1:5, sep="")
)
rownames(data.matrix) <- c(
  paste("gene", 1:100, sep="")
)
for (i in 1:100) {
  wt.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  ko.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  data.matrix[i,] <- c(wt.values, ko.values)
}

# procomp默认样本为行，基因为列，因此在这里需要转置
pca <- prcomp(t(data.matrix), scale=TRUE)
str(pca) # return x, sdev and rotation
# 取第一主成分和第二主成分
plot(pca$x[,1], pca$x[,2]) 
# 第一主成分解释了原始数据（所有10个样本的基因表达）中最大的方差
# 第二主成分解释了...第二大方差，以此类推

# 让我们看看PC1能解释原始数据多大的方差，因此我们使用sdev(标准差)的平方
# 标准差的平方是方差
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, 
        main="Scree Plot", xlab="Principal Component",
        ylab="Percent Variation")

library(ggplot2)
pca.data <- data.frame(
  Sample=rownames(pca$x),
  X=pca$x[,1],
  Y=pca$x[,2]
)

ggplot(data=pca.data, mapping = aes(x=X, y=Y, label=Sample)) +
  geom_text() + 
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) + 
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) + 
  theme_bw() +
  ggtitle("My PCA Graph")

# 下面我们可以通过pca$rotation来看什么基因对结果影响最大
loading_scores <- pca$rotation[,1] # 这里使用第一列因为PC1的对原始数据的影响达到91%
gene_scores <- abs(loading_scores)
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])
top_10_genes
# 查看造成wt与ko不一样的关键基因
pca$rotation[top_10_genes,1]
# 得到的结果说明，PC1下得分为正的是将ko推向PCA Graph右边的关键基因
# 为负的则是将wt推向左边的关键基因
# 说明这些基因是造成wt和ko不同的差异性表达基因

# 实战-单细胞RNA分析
a <- read.table("./GSE111229_Mammary_Tumor_fibroblasts_768samples_rawCounts.txt.gz",
                sep="\t", header=T)
head(a)
colnames(a)
rownames(a)
a[1:6, 1:4]
floor(ncol(a)/50)
dat <- a[apply(a,1,function(x){sum(x>1) > floor(ncol(a)/50)}),]
dim(dat)

dat <- log2(edgeR::cpm(dat)+1)
dat[1:4,1:4]

x <- 1:10
y <- 2*x
z <- rnorm(10)
tmp <- data.frame(x,y,z)
dist(t(tmp))
hc <- hclust(dist(t(dat)))

clus <- cutree(hc, 4)
table(clus)
group_list <- as.factor(clus)

colnames(dat)
library(stringr)
plate <- str_split(colnames(dat), "_", simplify = T)[,3]
head(plate)
table(plate)

n_g <- apply(a,2,function(x) sum(x>1))
table(n_g)
head(n_g)
dim(n_g)
class(n_g)
df <- data.frame(g=group_list, plate=plate, n_g=n_g)

hist(n_g, breaks=30)
save(a,dat,df,file='./input.Rdata')
