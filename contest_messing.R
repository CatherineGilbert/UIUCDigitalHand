library(tidyverse)
library(FactoMineR)
library(pheatmap)
library(irlba)
library(factoextra)

setwd("C:/Users/cmg3/Documents/GitHub/UIUCDigitalHand")

ng <- read_csv("C:/Users/cmg3/Downloads/numeric_genotype.csv")
ng_dt <- ng %>% column_to_rownames("Hybrid")
ng_mat <- ng %>% column_to_rownames("Hybrid") %>% as.matrix() %>% t()

#ng_cov <- cov(ng_mat, method = c("pearson"), use = "pairwise.complete.obs") #kinship matrix
#ng_cov_csv <- as.data.frame(ng_cov)
#write_csv(ng_cov_csv, "ng_cov.csv")
ng_cov_csv <- read_csv("ng_cov.csv")
ng_cov <- as.matrix(ng_cov_csv)

#ng_cor <- cor(ng_mat, method = c("pearson"), use = "pairwise.complete.obs") #kinship matrix
#ng_cor_csv <- as.data.frame(ng_cor)
#write_csv(ng_cor_csv, "ng_cor.csv")
ng_cor_csv <- read_csv("ng_cor.csv")
ng_cor <- as.matrix(ng_cor_csv)

ng_pca <- irlba::prcomp_irlba(ng_cor, n = 20)

summary(ng_pca)

get_eig(ng_pca)

fviz_screeplot(ng_pca, ncp = 20)

fviz_pca_ind(ng_pca, label = "none")

ng_pca_kmeans <- hkmeans(ng_pca$rotation, k = 10)
fviz_cluster(ng_pca_kmeans)

library(viridis)

ng_cor2 <- round(ng_cor, 2)
image(ng_cor2, useRaster = TRUE, col = magma(100))
