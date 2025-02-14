library(tidyverse)
select <- dplyr::select
library(FactoMineR)
library(pheatmap)
library(irlba)
library(factoextra)
library(viridis)
library(caret)
library(esquisse)
library(janitor)
library(misty)
library(CCA)
library(pheatmap)


# add inbred genetics to the PCA, then mark then on the ind. PCA biplot
# the hybrids should show up inbetween the marked inbreds at the corners of the PCA biplot


#get TT from planting to harvest 

setwd("C:/Users/cmg3/Documents/GitHub/UIUCDigitalHand")

# Generate a kinship matrix ---------
ng <- read_csv("numeric_genotype.csv")
ng_mat <- ng %>% column_to_rownames("Hybrid") %>% as.matrix() %>% t()

#ng_cor <- cor(ng_mat, method = c("pearson"), use = "pairwise.complete.obs") #kinship matrix
#ng_cor_csv <- as.data.frame(ng_cor)
#write_csv(ng_cor_csv, "ng_cor.csv")
ng_cor_csv <- read_csv("ng_cor.csv")
ng_cor <- as.matrix(ng_cor_csv)

  #mostly useless but nice looking image of the kinship matrix
ng_cor2 <- round(ng_cor, 2)
image(ng_cor2, useRaster = TRUE, col = magma(100))

# Create a hybrid parentage design matrix ---------
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


# Run genetics PCA (reducing genetics dimensionality)----
#with genetics as p (in n x p)
big_matrix <- ng %>% column_to_rownames("Hybrid") %>% as.matrix()

library(tidyverse)
big_matrix <- read_csv("ng_cor.csv")

library(pcaMethods)
object <- pca(big_matrix, nPcs = 20, method = "nipals")
write_rds(object, "genetics_PCA.rds")

hybrid_genetics_PCs <- as_tibble(genetics_PCA@scores) %>% mutate(Hybrid = ng$Hybrid) %>% relocate(Hybrid)
write_csv(hybrid_genetics_PCs, "hybrid_genetics_PCs.csv")


# Assign maturities to hybrids (using mean of total thermal time between environments) ------------ 
merged_dt <- read_csv("C:/Users/cmg3/Downloads/merged_data5.csv")

merge_mat <- select(merged_dt, Env, Year, Field_Location, Experiment, Hybrid, Date_Planted, Date_Harvested,
       Pollen_DAP_days, Silk_DAP_days, Latitude, Longitude, TT_pGerEme, TT_pEmeEnJ, TT_pEnJFlo, TT_pFloFla,       
       TT_pFlaFlw, TT_pFlwStG, TT_pStGEnG, TT_pEnGMat, TT_pMatHar) %>% 
  mutate(TT_total = TT_pEmeEnJ + TT_pEnJFlo + TT_pFloFla + 
           TT_pFlaFlw + TT_pFlwStG + TT_pStGEnG + TT_pEnGMat + TT_pMatHar) %>%
  mutate(Date_Planted = as.Date(Date_Planted, tryFormats = c("%m/%d/%y")),
         Date_Harvested = as.Date(Date_Harvested, tryFormats = c("%m/%d/%y")))
merge_mat <- unique(merge_mat)


# Assigning a Genetics value for the simulated cultivar to be used by APSIM 

# TT is by env, not hybrid in the given SC data
est_hybrid_mats <- select(merge_mat, Hybrid, Year, Env, Latitude:TT_total) %>% unique() 
# we'll just take the median of those TT to estimate a maturity for each hybrid. 
est_hybrid_mats <- est_hybrid_mats %>% group_by(Hybrid) %>% summarize(Genetics = median(TT_total))

est_hybrid_mats <- mutate(est_hybrid_mats, lett = str_to_upper(str_extract(Genetics,"^[A-Za-z]")), 
                                  num = as.numeric(str_extract(Genetics,"\\d+")))
est_hybrid_mats <- mutate(est_hybrid_mats, lett = ifelse(is.na(lett), "B", lett))
corn_mats <- c(80,90,95,100,103,105,108,110,112,115,120,130)
est_hybrid_mats <- est_hybrid_mats %>% rowwise() %>%
  mutate(num = corn_mats[which.min(abs(corn_mats - num))[1]]) %>%
  mutate(Mat = paste0(lett,"_",as.character(num)))
est_hybrid_mats <- select(est_hybrid_mats, -lett, -num)
write_csv(est_hybrid_mats, "est_hybrid_mats.csv")


# Creating a new SC input file
env_input <- select(merged_dt, Env, Year, Date_Planted, Latitude, Longitude, Hybrid)
env_input <- mutate(env_input, Date_Planted = as.Date(Date_Planted, tryFormats = c("%m/%d/%y")))

# Joining the SC input with the new mat estimations
env_input <- left_join(env_input, est_hybrid_mats)
env_input <- select(env_input, -Hybrid, -Genetics) %>% rename(Genetics = Mat, Location = Env, Planting = Date_Planted)
env_input <- unique(env_input)

write_csv(env_input, "input_with_mats.csv")

#After running through the SC Tool, compare to what's known to look at tool accuracy
daily_charact_x <- read_csv("C:/Users/cmg3/OneDrive/Documents/gtf_output/output/daily_charact_x.csv")
charact_x <- read_csv("C:/Users/cmg3/OneDrive/Documents/gtf_output/output/charact_x.csv")
trials_x <- read_csv("C:/Users/cmg3/OneDrive/Documents/gtf_output/output/trials_x.csv")


# get silk dates (period 6 end / 7 start)
sim_silk_dates <- filter(charact_x, Period == 7) %>% ungroup() %>% select(id_trial, SilkDate_Sim = Period_Start_Date)
silk_comp <- left_join(sim_silk_dates, trials_x) %>% select(-Genetics, -PlantingDate_Sim)
merge_mat2 <- select(merge_mat, Location = Env, Year, Hybrid : TT_total) %>% left_join(est_hybrid_mats) %>% 
  rename(PlantingDate = Date_Planted)
silk_comp2 <- left_join(silk_comp, merge_mat2)
silk_comp2 <- mutate(silk_comp2, Silk_Sim_days = as.numeric(as_date(SilkDate_Sim) - PlantingDate))

hmmm <- select(silk_comp2, Silk_Sim_days, Silk_DAP_days, Mat, Genetics) %>% 
  filter(!is.na(Silk_DAP_days)) %>% mutate(Silk_Diff = Silk_Sim_days - Silk_DAP_days)

silk_plt <- filter(silk_comp2, !is.na(Silk_DAP_days)) %>% mutate(Silk_Diff = Silk_Sim_days - Silk_DAP_days)

  ggplot(silk_plt) +
  aes(x = Silk_Sim_days, y = Silk_DAP_days) +
  geom_jitter(alpha = 0.10, size = 0.05) +
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

#HMMMM
  ggplot(silk_plt) +
    aes(x = Silk_DAP_days, y = Silk_Diff) +
    geom_jitter() +
    theme_minimal()
  
# Run environmental PCA -------

setwd("C:/Users/cmg3/Documents/GitHub/UIUCDigitalHand")

#get the seasonal characterization outputs, tie to the SC inputs
charact_x <- read_csv("C:/Users/cmg3/OneDrive/Documents/gtf_output/output/charact_x.csv")
charact_x <- filter(charact_x, !Period %in% c(6,9)) # removing 6 because it's almost always 1 day long, 9 because missing values
charact_x <- pivot_wider(charact_x, id_cols = id_trial, names_from = Period, values_from = Rain:Period_End_DOY)
trials_x <- read_csv("C:/Users/cmg3/OneDrive/Documents/gtf_output/output/trials_x.csv")
tc_x <- left_join(trials_x, charact_x)
tc_x <- select(tc_x, -Genetics) #remove dummy genetics used for shorter input file

#get the merged dataset, tie to genetics used for the SC inputs
merged_dt <- read_csv("merged_data5.csv")
est_hybrid_mats <- read_csv("est_hybrid_mats.csv")
env_merge <- left_join(merged_dt, est_hybrid_mats)

#tie the merged dataset to the seasonal characterization 
env_merge <- mutate(env_merge, Planting = as.Date(Date_Planted, tryFormats = c("%m/%d/%y")))
env_merge <- left_join(env_merge, tc_x, 
                       by = c("Year", "Latitude", "Longitude", 
                              "Planting", "Mat", "Env" = "Location"))

env <- select(env_merge, Hybrid, Env, Planting, HI30_pGerEme:LL__10, Latitude:`% Clay`, 
              Rain_1:WaterStress_11, Length_1:Length_11)

# var_check <- sapply(names(env), FUN = function(x){var(env[x], na.rm = TRUE)}) %>%
#   as.data.frame() %>% rownames_to_column("var") 
# names(var_check) <- c("var","variance")
# arrange(var_check, variance)

env <- select(env, !c(AccEmTT_1,AccEmTT_2,NutrientStress_1,NutrientStress_2,
                 WaterStress_1,WaterStress_2,WaterStress_3,WaterStress_5,
                 WaterStress_10,WaterStress_11))

env <- remove_constant(env)
env <- remove_empty(env, which = c("rows","cols"))
env <- distinct(env)

env_PCA <- pca(env, nPcs = 20, method = "nipals", scale = "pareto")
write_rds(env_PCA, "env_PCA.rds")
env_PCs <- as_tibble(env_PCA@scores) %>% cbind(env[1:3],.)
write_csv(env_PCs, "env_PCs.csv")

# na_prop <- sapply(names(env), FUN = function(x){mean(is.na(env[x]))}) %>%
#   as.data.frame() %>% rownames_to_column("var") 
# names(na_prop) <- c("var","na_prop")
# 
# arrange(na_prop, desc(na_prop)) %>% head()

# #run PCA for hybrids
# env_hybrid <- env %>% group_by(Hybrid) %>% select(-Env) %>% summarize(across(everything(), mean))
# env_hybrid_pca <- pca(env_hybrid, nPcs = 20, method = "nipals", scale = "pareto")
# write_rds(env_hybrid_pca, "env_hybrid_PCA.rds")
# env_hybrid_PCs <- as_tibble(env_pca_hybrid@scores) %>% mutate(Hybrid = env_hybrid$Hybrid) %>% relocate(Hybrid)
# write_csv(env_hybrid_PCs, "hybrid_env_PCs.csv")
# 
# #run PCA for locations
# env_loc <- env %>% group_by(Env) %>% select(-Hybrid) %>% summarize(across(everything(), mean))
# env_loc_pca <- pca(env_loc, nPcs = 20, method = "nipals", scale = "pareto")
# write_rds(env_loc_pca, "env_loc_PCA.rds")



# Plot the PCAs -------

env_PCA <- read_rds("env_PCA.rds")
genetics_PCA <- read_rds("genetics_PCA.rds")

slplot(env_PCA)
pheatmap(loadings(env_PCA), cluster_cols = FALSE,
         breaks = seq(from = -0.8, to = 0.8, length.out = 101))
env_screeplt <- data.frame("Eigenvalues" = env_PCA@R2, "PC" = paste0("Env_PC",c(1:10)))
ggplot(env_screeplt) +
  aes(x = fct_inorder(PC), y = Eigenvalues) +
  geom_bar(stat = "summary", fun = "sum", fill = "#112446") +
  theme_minimal()

slplot(genetics_PCA, lcex = 0.1, scex = 0.1)
pheatmap(loadings(genetics_PCA), cluster_cols = FALSE)
genetics_screeplt <- data.frame("Eigenvalues" = genetics_PCA@R2, "PC" = paste0("Genetics_PC",c(1:10)))
ggplot(genetics_screeplt) +
  aes(x = fct_inorder(PC), y = Eigenvalues) +
  geom_bar(stat = "summary", fun = "sum", fill = "#112446") +
  theme_minimal()



# Create training data with yield and PCs ---------

merged_dt <- read_csv("merged_data5.csv")
mergedPC <- mutate(merged_dt, Planting = as.Date(Date_Planted, tryFormats = c("%m/%d/%y")))
#mergedPC <- select(mergedPC, Env:Yield_Mg_ha, Planting)

env_PCs <- read_csv("env_PCs.csv")
names(env_PCs) <- c("Hybrid", "Env", "Planting", paste0("Env_PC",c(1:20)))

hybrid_genetics_PCs <- read_csv("hybrid_genetics_PCs.csv")
names(hybrid_genetics_PCs) <- c("Hybrid", paste0("Genetics_PC",c(1:20)))

mergedPC <- left_join(mergedPC, env_PCs)
mergedPC <- left_join(mergedPC, hybrid_genetics_PCs)

write_csv(mergedPC, "merged_w_PCs.csv")


# Create a model


# Predict PCs for testing data ------

testing_data <- read_csv("testing_data.csv")
merged_w_PCs <- read_csv("merged_w_PCs.csv")

#ENV PCs
env_PCA <- readRDS("env_PCA.rds")
required_var <- loadings(env_PCA) %>% rownames()
missing_var <- required_var[!required_var %in% names(testing_data)]

#Soil data
  #for the testing sites missing soil values, which we have from the training data, join the soil values we have
  
  merged_dt <- read_csv("merged_data5.csv")
  
  train_soils <- select(merged_dt, Field_Location, Latitude, Longitude, Env, `%K Sat`:`% Clay`) %>% distinct()
  train_soils <- filter(train_soils, !is.na(`%K Sat`))
  
  #one canadian location's soil data is unavailable
  unique(testing_data$Field_Location)[!unique(testing_data$Field_Location) %in% train_soils$Field_Location]
  
  test_soils <- select(testing_data, Field_Location, Latitude, Longitude) %>% distinct()
  
  library(terra)
  test_soils_v <- vect(test_soils, geom = c("Longitude","Latitude"))
  train_soils_v <- vect(train_soils, geom = c("Longitude","Latitude"))
  
  near_soil <- nearest(test_soils_v, train_soils_v) %>% values()
  near_soil
  
  test_soils <- train_soils[near_soil$to_id,] %>% select(-c(Latitude, Longitude, Field_Location)) %>% 
    rename(Ref_Env = Env) %>% cbind(distance = near_soil$distance) %>% cbind(test_soils, .)
  
  #join the soil data
  testing_data2 <- left_join(testing_data, test_soils)

#Seasonal parameters from SCT

  testing_data2 <- mutate(testing_data2, TT_total = TT_pEmeEnJ + TT_pEnJFlo + TT_pFloFla + 
           TT_pFlaFlw + TT_pFlwStG + TT_pStGEnG + TT_pEnGMat + TT_pMatHar,
         Planting = as.Date(paste0(Date_Planted,"/",Year), tryFormats = c("%j/%Y"))) 

  # assign maturities
  est_hybrid_mats <- select(testing_data2, Hybrid, Year, Env, Latitude, Longitude, TT_total) %>% unique() 
  # we'll just take the median of those TT to estimate a maturity for each hybrid. 
  est_hybrid_mats <- est_hybrid_mats %>% group_by(Hybrid) %>% summarize(Specific_Genetics = median(TT_total, na.rm=TRUE))
  corn_mats <- c(80,90,95,100,103,105,108,110,112,115,120,130)
  est_hybrid_mats <- est_hybrid_mats %>% rowwise() %>%
    mutate(Genetics = corn_mats[which.min(abs(corn_mats - Specific_Genetics))[1]])
  
  #make the testing input, run through
  testing_data2 <- left_join(testing_data2, est_hybrid_mats)
  testing_input <- select(testing_data2, Location = Env, Latitude, Longitude, Year, Planting, Genetics) %>% distinct()
  
  charact_x <- read_csv("C:/Users/cmg3/OneDrive/Documents/testing_gtf/output/charact_x.csv")
  trials_x <- read_csv("C:/Users/cmg3/OneDrive/Documents/testing_gtf/output/trials_x.csv")
  
  tc_x <- pivot_wider(charact_x, id_cols = "id_trial", names_from = "Period", values_from = Rain:Period_End_DOY) %>% 
    left_join(trials_x,.) %>% relocate(c(Year, Latitude, Longitude, Planting, Genetics))
  
  testing_data2 <- relocate(testing_data2, c(Year, Latitude, Longitude, Planting, Genetics))
  
  testing_data3 <- left_join(testing_data2, tc_x) %>% select(all_of(required_var))

  env_PCs <- predict(env_PCA, testing_data3) %>% .$scores %>% as_tibble()
  names(env_PCs) <- c(paste0("Env_PC",c(1:20)))
  
  testing_data_full <- select(testing_data, Env:Longitude) %>% cbind(env_PCs)

# add Genetics PCs to testing
  genetics_PCA <- readRDS("genetics_PCA.rds")
  
  big_matrix <- ng %>% column_to_rownames("Hybrid") %>% as.matrix()
  
  gen_PCs <- as_tibble(genetics_PCA@scores) %>% mutate(Hybrid = rownames(big_matrix)) %>% 
    relocate(Hybrid)
  names(gen_PCs) <- c("Hybrid",paste0("Genetics_PC",c(1:20)))
  
  testing_data_full <- left_join(testing_data_full, gen_PCs)

  write_csv(testing_data_full, "testing_data_w_PCs.csv")  
  
  