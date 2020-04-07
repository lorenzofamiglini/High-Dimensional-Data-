---
title: "Hdda2"
output: html_document
---

```{r setup, include=FALSE}
library("rbiom")
library("biomformat")
```

```{r}
biom <- read_biom("D:/Datasets/Hdda_project/221_otu_table.biom")
summary(biom)
show(biom)
```


```{r}
fecal_samples = as(biom_data(biom), "matrix")
## taxonomic rank: 
# species -genus -family-class- phylum 
bacteria_taxa = observation_metadata(biom)

bacteria_taxa$taxonomy2 <- gsub("\\[|\\]", "", bacteria_taxa$taxonomy2)
bacteria_taxa$taxonomy3 <- gsub("\\[|\\]", "", bacteria_taxa$taxonomy3)
bacteria_taxa$taxonomy4 <- gsub("\\[|\\]", "", bacteria_taxa$taxonomy4)
bacteria_taxa$taxonomy5 <- gsub("\\[|\\]", "", bacteria_taxa$taxonomy5)
bacteria_taxa$taxonomy6 <- gsub("\\[|\\]", "", bacteria_taxa$taxonomy6)
```


```{r}
tax2 <- bacteria_taxa[,c("taxonomy1","taxonomy2")]
tax2["id"] <- rownames(tax2)
tax2["merge"] <- tidyr::unite(tax2[,-ncol(tax2)],"merge_type",sep = ",")

tax3 <- bacteria_taxa[,c("taxonomy1","taxonomy2", "taxonomy3")]
tax3["id"] <- rownames(tax3)
tax3["merge"] <- tidyr::unite(tax3[,-ncol(tax3)],"merge_type",sep = ",")

tax4 <- subset(bacteria_taxa[,c(1:4)], taxonomy4 != "o__")
tax4["id"] <- rownames(tax4)
tax4["merge"] <- tidyr::unite(tax4[,-ncol(tax4)],"merge_type",sep = ",")

tax5 <- subset(bacteria_taxa[,c(1:5)], taxonomy5 != "f__")
tax5["id"] <- rownames(tax5)
tax5["merge"] <- tidyr::unite(tax5[,-ncol(tax5)],"merge_type",sep = ",")

tax6 <- subset(bacteria_taxa[,c(1:6)], taxonomy6 != "g__")
tax6["id"] <- rownames(tax6)
tax6["merge"] <- tidyr::unite(tax6[,-ncol(tax6)],"merge_type",sep = ",")
```

```{r}
tax2[,c("taxonomy1","taxonomy2")] <- list(NULL)
tax3[,c("taxonomy1","taxonomy2","taxonomy3")] <- list(NULL)
tax4[,c("taxonomy1","taxonomy2","taxonomy3","taxonomy4")] <- list(NULL)
tax5[,c("taxonomy1","taxonomy2","taxonomy3","taxonomy4","taxonomy5")] <- list(NULL)
tax6[,c("taxonomy1","taxonomy2","taxonomy3","taxonomy4","taxonomy5","taxonomy6")] <- list(NULL)
```


```{r}
df_fecal <- as.data.frame(fecal_samples)
df_fecal$id <- rownames(df_fecal)
```


```{r}
library(dplyr)
prepare_data <- function(tax2){
  fec_tax2 <- as.data.frame(merge(df_fecal,tax2,by="id", sort = F))
  fec_tax2$id <- NULL
  fec_tax2 <- fec_tax2 %>% group_by(merge) %>% summarize_all(sum) 
  fec_tax2 <- as.data.frame(t(fec_tax2))
  colnames(fec_tax2) <- as.character(unlist(fec_tax2[1,]))
  fec_tax2 <- fec_tax2[-1,]
  return(fec_tax2)
}

```

```{r}
final_tax2 <- prepare_data(tax2)
final_tax3 <- prepare_data(tax3)
final_tax4 <- prepare_data(tax4)
final_tax5 <- prepare_data(tax5)
final_tax6 <- prepare_data(tax6)
final <- cbind(final_tax2, final_tax3, final_tax4, final_tax5, final_tax6)
as.data.frame(t(final))
```


```{r}
write.csv(final_tax6,"D:/Datasets/Hdda_project/otus_tab.csv")
```


```{r}
final_tax6 <- as.data.frame(lapply(final_tax6, as.numeric))
clu <- cluster::pam(final_tax6,2,diss = F, metric = "manhattan")
cluster::plot.partition(clu)
clu$silinfo$avg.width
```

```{r}
fecal_samples <- apply(fecal_samples,2,function(x){x/sum(x)})


dist_jsd <- philentropy::JSD(fecal_samples)


for (i in 2:30){
  clu <- cluster::pam(dist_jsd,i,diss = T, metric = "manhattan")
 #cluster::plot.partition(clu)
  print(clu$silinfo$avg.width)
}




heatmap(dist_jsd)
```

```{r}
data <- read.csv("C:/Users/loren/OneDrive/Desktop/MetaHIT_SangerSamples.genus.txt", header=T, row.names=1, dec=".", sep="\t")
data=data[-1,]
```


```{r}
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
	KLD <- function(x,y) sum(x *log(x/y))
	JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
	matrixColSize <- length(colnames(inMatrix))
	matrixRowSize <- length(rownames(inMatrix))
	colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
        
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

	for(i in 1:matrixColSize) {
		for(j in 1:matrixColSize) { 
			resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
			as.vector(inMatrix[,j]))
		}
	}
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
 }
```


```{r}
data.dist=dist.JSD(fecal_samples)
data.dist
```


```{r}
 pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
                         require(cluster)
                         cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
                         return(cluster)
                        }
```


```{r}
data.cluster=pam.clustering(data.dist, k=2)
data.cluster
```


```{r}
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
```


```{r}
obs.silhouette
```

```{r}
library(ade4)
 obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
 obs.bet=bca(obs.pca, fac=as.factor(data.cluster)) 
 s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F)
```

