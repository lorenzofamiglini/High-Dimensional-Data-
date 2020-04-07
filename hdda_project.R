---
title: "R Notebook"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 

```{r}
#install.packages("devtools")
#devtools::install_github("joey711/biomformat")
library("biomformat")

#install.packages("remotes")
remotes::install_github("cmmr/rbiom")

library("rbiom")
```
```{r}
map_nutriens <- read.csv("D:/Datasets/Hdda_project/46764_mapping_file.txt", sep = "\t")
numeric_cols <- unlist(lapply(map_nutriens, is.numeric))  

nut = map_nutriens[ , numeric_cols]

nut[,c("lib_reads_seqd","run_date","qiita_prep_id","age","anonymized_name","data_gen_ncc_database_ave",
                     "data_ncc_database_version","data_software_version_ave","day_of_intake","dri_life_stage_rda_cat",
                     "elevation","energy_kcal_ave","energy_kj_ave","host_subject_id","host_taxid","intake_amount",
                     "intake_reliability","latitude","longitude","perc_cal_from_protein_ave","perc_cal_from_pufa_ave",
                     "perc_cal_from_sfa_ave","polyunsatsat_fat_ratio_ave","qiita_study_id","taxon_id","visit_number")] <- list(NULL)

nutriens = cbind(map_nutriens$X.SampleID,nut)
names(nutriens)[1] = "id"
colnames(nutriens)
```


```{r}
library(biomformat)
biom <- read_biom("D:/Datasets/Hdda_project/221_otu_table.biom")

summary(biom)
show(biom)

# type: OTU table 
# matrix_type: dense 
# 3119 rows and 100 columns 

## taxonomic rank: 
# species - genus -family-class- phylum 

dati <- biom_data(biom)
dati
str(dati)
#sort(colnames(x))
dati@Dimnames[[2]] #da qui possiamo prendere l'id primario 
dati@x
mat <- dati@x

print(length(unique(observation_metadata(biom)[,6])))

data.frame(table(observation_metadata(biom)[,5])) 


fecal_samples = as(biom_data(biom), "matrix")

#phy_tree(biom)

## taxonomic rank: 
# species -genus -family-class- phylum 
bacteria_taxa = as(observation_metadata(biom), "matrix")
#bacteria_taxa$taxonomy7 <- NULL
#table(bacteria_taxa$taxonomy6)

```

```{r}
library("phyloseq")
fecal_samples <- apply(fecal_samples,2,function(x){x/sum(x)})
OTU <- otu_table(fecal_samples, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(bacteria_taxa))
physeq <- phyloseq(OTU, TAX)

plot_bar(physeq, fill = "taxonomy6")
```
```{r}
library("ape")
random_tree = rtree(3393, rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree,show.tip.label = F)
```


```{r}
physeq1 = merge_phyloseq(physeq, random_tree)

```

```{r}
plot_tree(physeq1,label.tips=F, ladderize="left", plot.margin=0.3)
```

```{r}
plot_heatmap(physeq1)
plot_heatmap(physeq1, taxa.label="taxonomy2")
```

```{r}
distanza <- UniFrac(physeq1, weighted=FALSE, normalized=TRUE, parallel=TRUE, fast=TRUE)
#distance(physeq1, "uunifrac") 

distanza_df <- as.data.frame(as.matrix(distanza))
distanza_df$id <- rownames(distanza_df)


permanova_df <- as.data.frame(merge(distanza_df,nutriens,by='id', sort = F))

perm_bacteria = permanova_df[,2:101]
perm_nutriens = permanova_df[,102:ncol(permanova_df)]


perm_nutriens = perm_nutriens[vapply(perm_nutriens, function(x) length(unique(x)) > 1, logical(1L))]
# all columns to grams
for (i in 1:ncol(perm_nutriens)){
  if (grepl("mg",(colnames(perm_nutriens)[i]))){
    perm_nutriens[,i] = perm_nutriens[,i]/1000
  }
  else if (grepl("mcg",(colnames(perm_nutriens)[i]))){
    perm_nutriens[,i] = perm_nutriens[,i]/1000000
  }
}

standardize = function(col){
  col <- (col - mean(col)) / sd(col)
  return (col)
}

standardized_nutriens = apply(perm_nutriens,2,standardize)
#perm_nutriens = as.data.frame(scale(perm_nutriens))
plot(distanza)
```

```{r}
lista <- c()
for (i in 1:ncol(perm_nutriens)){
  permanova <- vegan::adonis(distanza ~ perm_nutriens[,i],data = perm_nutriens, permutations=1000, method = "bray")
  if (unlist(permanova)$`aov.tab.Pr(>F)1` >= 0){
    lista = append(lista, unlist(permanova)$`aov.tab.Pr(>F)1`)
  }
}
```

```{r}
set.seed(20)
p_lista <- c()
for (i in 1:ncol(perm_nutriens)){
  tests <- lapply(1:3, function(i) vegan::adonis(distanza ~ perm_nutriens[,i], permutations = 100, method = "bray"))
  get_list <- sapply(1:3, function(i) tests[[i]]$aov.tab$`Pr(>F)`[1])
  get_list_adj <- p.adjust(get_list, method = "fdr")
  p_val_adj <- mean(get_list_adj)
  p_lista <- append(p_lista, p_val_adj)
}

p_lista <-p_lista[p_lista <= 0.25]
```


```{r}

conta <- 0
list_idx <- c()
for(i in p_lista){
  conta = conta + 1
  if (i <= 0.25){
    list_idx = append(list_idx, conta)
  }
}


sele_col <- names(perm_nutriens)[list_idx]
perm_nutriens[,sele_col]

permanova$f.perms

```

```{r}
for (i in 2:30){
  clu <- cluster::pam(dist_jsd,i,diss = T, metric = "manhattan")
 #cluster::plot.partition(clu)
  print(clu$silinfo$avg.width)
}

otu <- as.matrix(otu_table(physeq1))

dist_jsd <- philentropy::JSD(t(fecal_samples))
heatmap(dist_jsd)
```


```{r}
UniFrac()
```

