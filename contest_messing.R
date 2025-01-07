library(tidyverse)
library(FactoMineR)
library(pheatmap)
library(irlba)
library(factoextra)
library(viridis)

# add inbred genetics to the PCA, then mark then on the ind. PCA biplot
# the hybrids should show up inbetween the marked inbreds at the corners of the PCA biplot


setwd("C:/Users/cmg3/Documents/GitHub/UIUCDigitalHand")

#GENERATE KINSHIP ---------
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

  #mostly useless but nice looking image of the kinship matrix
ng_cor2 <- round(ng_cor, 2)
image(ng_cor2, useRaster = TRUE, col = magma(100))

#RUN PCA
ng_pca <- irlba::prcomp_irlba(ng_cor, n = 20)
summary(ng_pca)
get_eig(ng_pca)

  #plots
fviz_screeplot(ng_pca, ncp = 20)
fviz_pca_ind(ng_pca, label = "none")
ng_pca_kmeans <- hkmeans(ng_pca$rotation, k = 10)
fviz_cluster(ng_pca_kmeans, geom = "point")

library(NbClust)
nb <- NbClust(ng_cor2, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "kmeans", index ="ccc")


# Create Hybrid parentage design matrix ---------
ng_cor_dt <- mutate(ng_cor_csv, Hybrid = names(ng_cor_csv)) %>% select(Hybrid)
#ng_cor_dt <- separate_wider_delim(ng_cor_dt, Hybrid, delim = "/", names = c("P1","P2"), cols_remove = FALSE) 

ng_cor_dt <- ng_cor_dt %>% mutate(P1 = str_split(Hybrid, "/"))
ng_cor_dt <- ng_cor_dt %>% mutate(P0 = str_split(Hybrid, "[/]|_(?=[A-Z])"))
ng_cor_dt$extra <- mapply(intersect, ng_cor_dt$P1, ng_cor_dt$P0, SIMPLIFY = FALSE)
ng_cor_dt$P0 <- mapply(c, ng_cor_dt$P0, ng_cor_dt$extra, SIMPLIFY = FALSE)
ng_cor_dt$P0 <- mapply(sort, ng_cor_dt$P0, SIMPLIFY = FALSE)
ng_cor_dt <- select(ng_cor_dt, -extra)

pedigree_long <- unnest(ng_cor_dt, cols = P0) %>% group_by(Hybrid) %>% count(P0) %>% mutate(n = n/4)
pedigree <- arrange(pedigree_long, P0) %>% pivot_wider(id_cols = Hybrid, names_from = P0, values_from = n)

write.csv(pedigree, "pedigree.csv", na = "0")
pedigree <- read_csv("pedigree.csv", col_types = cols(...1 = col_skip()))

contribution <- pedigree_long %>% group_by(P0) %>% summarise(sum(n))


# Add maturities to env inv ------------ 
merged_dt <- read_csv("C:/Users/cmg3/Downloads/merged_data5.csv")

merge_mat <- select(merged_dt, Env, Year, Field_Location, Experiment, Hybrid, Date_Planted, Date_Harvested,
       Pollen_DAP_days, Silk_DAP_days, Latitude, Longitude, TT_pGerEme, TT_pEmeEnJ, TT_pEnJFlo, TT_pFloFla,       
       TT_pFlaFlw, TT_pFlwStG, TT_pStGEnG, TT_pEnGMat, TT_pMatHar) %>% 
  mutate(TT_total = TT_pEmeEnJ + TT_pEnJFlo + TT_pFloFla + 
           TT_pFlaFlw + TT_pFlwStG + TT_pStGEnG + TT_pEnGMat + TT_pMatHar) %>%
  mutate(Date_Planted = as.Date(Date_Planted, tryFormats = c("%m/%d/%y")),
         Date_Harvested = as.Date(Date_Harvested, tryFormats = c("%m/%d/%y")))
merge_mat <- unique(merge_mat)

# TT is by env, not hybrid in the given SC data
est_hybrid_mats <- select(merge_mat, Hybrid, Year, Env, Latitude:TT_total) %>% unique() 
# we'll just take the median of those TT to estimate a maturity for each hybrid. 
est_hybrid_mats <- est_hybrid_mats %>% group_by(Hybrid) %>% summarize(Genetics = median(TT_total))

# Assigning a Genetics value for the simulated cultivar to be used by APSIM
est_hybrid_mats <- mutate(est_hybrid_mats, lett = str_to_upper(str_extract(Genetics,"^[A-Za-z]")), 
                                  num = as.numeric(str_extract(Genetics,"\\d+")))
est_hybrid_mats <- mutate(est_hybrid_mats, lett = ifelse(is.na(lett), "B", lett))
corn_mats <- c(80,90,95,100,103,105,108,110,112,115,120,130)
est_hybrid_mats <- est_hybrid_mats %>% rowwise() %>%
  mutate(num = corn_mats[which.min(abs(corn_mats - num))[1]]) %>%
  mutate(Mat = paste0(lett,"_",as.character(num)))
est_hybrid_mats <- select(est_hybrid_mats, -lett, -num)
write_csv(est_hybrid_mats, "est_hybrid_mats.csv")

envirotype_input <- read_csv("C:/Users/cmg3/Downloads/envirotype_input.csv", 
                             +     col_types = cols(...1 = col_skip()))


# Creating a new SC input file ------
env_input <- select(merged_dt, Env, Year, Date_Planted, Latitude, Longitude, Hybrid)
env_input <- mutate(env_input, Date_Planted = as.Date(Date_Planted, tryFormats = c("%m/%d/%y")))

# Joining the SC input with the new mat estimations
env_input <- left_join(env_input, est_hybrid_mats)
env_input <- select(env_input, -Hybrid, -Genetics) %>% rename(Genetics = Mat, Location = Env, Planting = Date_Planted)
env_input <- unique(env_input)

write_csv(env_input, "input_with_mats.csv")

#After running through the SC Tool, compare to what's known to look at tool accuracy ----
daily_charact_x <- read_csv("C:/Users/cmg3/OneDrive/Documents/gtf_output/output/daily_charact_x.csv")
charact_x <- read_csv("C:/Users/cmg3/OneDrive/Documents/gtf_output/output/charact_x.csv")
trials_x <- read_csv("C:/Users/cmg3/OneDrive/Documents/gtf_output/output/trials_x.csv")


# get silk dates (period 6 end / 7 start)

sim_silk_dates <- filter(charact_x, Period == 7) %>% ungroup() %>% select(id_trial, SilkDate_Sim = Period_Start_Date)
silk_comp <- left_join(sim_silk_dates, trials_x) %>% select(-Genetics, -PlantingDate_Sim)
merge_mat2 <- select(merge_mat, Location = Env, Year, Hybrid : TT_total) %>% left_join(est_hybrid_mats) %>% 
  rename(PlantingDate = Date_Planted) %>% mutate(Year = as.character(Year))
silk_comp2 <- left_join(silk_comp, merge_mat2)
library(lubridate)
silk_comp2 <- mutate(silk_comp2, Silk_Sim_days = as.numeric(as_date(SilkDate_Sim) - PlantingDate))


filter(silk_comp2, !is.na(Silk_DAP_days)) %>% 
  select(Silk_Sim_days, Silk_DAP_days) %>%
  ggplot() +
  aes(x = Silk_Sim_days, y = Silk_DAP_days) +
  geom_jitter(alpha = 0.05) +
  geom_line(
    aes(x = Silk_Sim_days,
        y = Silk_Sim_days),
    colour = "#CE4260"
  ) + 
  labs(
    title = "Days from Planting to Silk", 
    x = "Simulated # of Days",
    y = "Observed # of Days"
  ) +
  theme_minimal()
