# PM25 analysis from the full spatio-temporal random forest models
# -- see cv_rf_st.R, cv_rf_st_m5.R
# Create tables and plots of results
# -- This file has been simplified to only report the results we care about.

# Helper functions of obtaining results
getRMSE <- function(true, pred) {
  return(sqrt(mean( (true - pred)^2 )))
}
getR2 <- function(true, pred) {
  return(summary( lm( true ~ pred ) )$r.squared)
}
getIntSlope <- function(true, pred) {
  return(as.numeric(coef( lm( true ~ pred) )))
}
getInt <- function(true, pred) {
  return(as.numeric(coef( lm( true ~ pred) ))[1])
}
getSlope <- function(true, pred) {
  return(as.numeric(coef( lm( true ~ pred) ))[2])
}

library(testthat) # Helps me run various checks
library(dplyr) # Manipulates data more easily
packageVersion("dplyr") # MUST BE AT LEAST 1.0.0
library(xtable) # Produces basic latex tables -- I modify later by hand
sessionInfo()

# Read in the results from cv, b1, b2
all_results <- readRDS("./pm25/cvresults/rftemp_comb.rds")
df_full_cv <- all_results$df_save_cv 
df_full_b1 <- all_results$df_save_b1
df_full_b2 <- all_results$df_save_b2
param <- all_results$param
rm(all_results)
gc()

# Read in additional results
extra_results <- readRDS("./pm25/cvresults/rf_st_results2.rds")
df_full_cv_2 <- extra_results$df_save_cv
df_full_b1_2 <- extra_results$df_save_b1
df_full_b2_2 <- extra_results$df_save_b2
expect_identical(df_full_cv_2$fold_cv, df_full_cv$fold_cv)
expect_equivalent(df_full_cv_2$cmaq_x, df_full_cv$cmaq_x)
expect_equivalent(df_full_cv_2$pm25_value, df_full_cv$pm25_value)
expect_identical(df_full_b1_2$fold_b1, df_full_b1$fold_b1)
expect_equivalent(df_full_b1_2$cmaq_x, df_full_b1$cmaq_x)
expect_equivalent(df_full_b1_2$pm25_value, df_full_b1$pm25_value)
expect_identical(df_full_b2_2$fold_b2, df_full_b2$fold_b2)
expect_equivalent(df_full_b2_2$cmaq_x, df_full_b2$cmaq_x)
expect_equivalent(df_full_b2_2$pm25_value, df_full_b2$pm25_value)

# Combine datasets
df_cv <- df_full_cv %>% 
  left_join(df_full_cv_2 %>% select(day, cmaq_id, p5a_1:p5c_4),
            by = c("cmaq_id" = "cmaq_id", "day" = "day")) %>%
  rename(fold = fold_cv) %>%
  mutate(fold_type = "Random")

df_b1 <- df_full_b1 %>%
  left_join(df_full_b1_2 %>% select(day, cmaq_id, p5a_1:p5c_4),
            by = c("cmaq_id" = "cmaq_id", "day" = "day")) %>%
  rename(fold = fold_b1) %>%
  mutate(fold_type = "Spatial Cluster 1")

df_b2 <- df_full_b2 %>%
  left_join(df_full_b2_2 %>% select(day, cmaq_id, p5a_1:p5c_4),
            by = c("cmaq_id" = "cmaq_id", "day" = "day")) %>%
  rename(fold = fold_b2) %>%
  mutate(fold_type = "Spatial Cluster 2")

df_stack <- rbind(df_cv, df_b1, df_b2) %>%
  # New functionality from dplyr 1.0
  relocate(fold_type) %>%
  mutate(aod_missing = ifelse(is.na(aod_value), "Missing", "Observed")) %>%
  select(-contains(c("b_", "c_", "d_1", "d_2", "d_3", "d_4")))
nrow(df_stack) # 63,136

# Summarize by fold type -- overall
all_r2 <- df_stack %>%
  group_by(fold_type) %>%
  summarize(
    across(contains(c("_1", "_2", "_3", "_4")), ~ getR2(pm25_value, .x))    
  ) %>%
  ungroup(.)

all_rmse <- df_stack %>%
  group_by(fold_type) %>%
  summarize(
    across(contains(c("_1", "_2", "_3", "_4")), ~ getRMSE(pm25_value, .x))    
  ) %>%
  ungroup(.)

all_int <- df_stack %>%
  group_by(fold_type) %>%
  summarize(
    across(contains(c("_1", "_2", "_3", "_4")), ~ getInt(pm25_value, .x))    
  ) %>%
  ungroup(.)

all_slope <- df_stack %>%
  group_by(fold_type) %>%
  summarize(
    across(contains(c("_1", "_2", "_3", "_4")), ~ getSlope(pm25_value, .x))    
  ) %>%
  ungroup(.)

# Summarize by fold type -- by missingness
miss_r2 <- df_stack %>%
  group_by(fold_type, aod_missing) %>%
  summarize(
    across(contains(c("_1", "_2", "_3", "_4")), ~ getR2(pm25_value, .x))    
  ) %>%
  ungroup(.)

miss_rmse <- df_stack %>%
  group_by(fold_type, aod_missing) %>%
  summarize(
    across(contains(c("_1", "_2", "_3", "_4")), ~ getRMSE(pm25_value, .x))    
  ) %>%
  ungroup(.)

miss_int <- df_stack %>%
  group_by(fold_type, aod_missing) %>%
  summarize(
    across(contains(c("_1", "_2", "_3", "_4")), ~ getInt(pm25_value, .x))    
  ) %>%
  ungroup(.)

miss_slope <- df_stack %>%
  group_by(fold_type, aod_missing) %>%
  summarize(
    across(contains(c("_1", "_2", "_3", "_4")), ~ getSlope(pm25_value, .x))    
  ) %>%
  ungroup(.)


# For model (a), (1, 2, 3, 4, 5)
#   pick the best mtry
# For the time being, ignoring the other models (b, c, d)
mnum <- seq(5)
lnum <- c("a")
mgrid <- expand.grid(mnum = mnum, lnum = lnum, 
                     stringsAsFactors = F) %>%
  mutate(mname = paste0(mnum, lnum))
attributes(mgrid)$out.attrs <- NULL
mgrid$cv = NA; mgrid$b1 = NA; mgrid$b2 <- NA
mgrid$cv_mtry = NA; mgrid$b1_mtry = NA; mgrid$b2_mtry <- NA

for (mi in seq(nrow(mgrid))) {
  use_names <- paste0("p", mgrid$mnum[mi], mgrid$lnum[mi], "_", seq(4))
  mgrid$cv[mi] <- which.min(as.vector(all_rmse[1, use_names]))
  mgrid$b1[mi] <- which.min(as.vector(all_rmse[2, use_names]))
  mgrid$b2[mi] <- which.min(as.vector(all_rmse[3, use_names]))
  mgrid$cv_mtry[mi] <- param$mtry[mgrid$cv[mi]]
  mgrid$b1_mtry[mi] <- param$mtry[mgrid$b1[mi]]
  mgrid$b2_mtry[mi] <- param$mtry[mgrid$b2[mi]]
} 

saveRDS(
  list(
    all_r2 = all_r2,
    all_rmse = all_rmse,
    all_int = all_int,
    all_slope = all_slope,
    miss_r2 = miss_r2,
    miss_rmse = miss_rmse,
    miss_int = miss_int,
    miss_slope = miss_slope,
    mgrid = mgrid,
    param = param,
    descriptions = c("a = random forest with convolution",
                     "b = extra trees", "c = convolution and force (x,y)",
                     "d = convolution and force day, day of week"),
    file_source = c("st_results.R")
  ),
  "pm25/cvresults/cvbest_st.rds"
)

print(mgrid)
# CV: 4 for m1, m2; 8 for m3-m5
# B1: 8 for all m1-m5
# B2: 4 for m1-m4; 8 for m5

#### Scatter plots ####
# -- Construct for each of 5 models for both mtry = 4 and mtry = 8
png(filename = "./pm25/tables/scatter_cvr_m3_mtry4.png")
plot(df_cv$p3a_1, df_cv$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()

png(filename = "./pm25/tables/scatter_cvr_m3_mtry8.png")
plot(df_cv$p3a_2, df_cv$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()

png(filename = "./pm25/tables/scatter_cvb1_m3_mtry4.png")
plot(df_b1$p3a_1, df_b1$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()

png(filename = "./pm25/tables/scatter_cvb1_m3_mtry8.png")
plot(df_b1$p3a_2, df_b1$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()

png(filename = "./pm25/tables/scatter_cvb2_m3_mtry4.png")
plot(df_b2$p3a_1, df_b2$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()

png(filename = "./pm25/tables/scatter_cvb2_m3_mtry8.png")
plot(df_b2$p3a_2, df_b2$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()

#### Create tables ####
r2_best <- 
  rbind(
    all_r2 %>% filter(fold_type == "Random") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    miss_r2 %>% filter(fold_type == "Random") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    all_r2 %>% filter(fold_type == "Spatial Cluster 1") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    miss_r2 %>% filter(fold_type == "Spatial Cluster 1") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    all_r2 %>% filter(fold_type == "Spatial Cluster 2") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2),
    miss_r2 %>% filter(fold_type == "Spatial Cluster 2") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2)
  ) %>%
  rename(Method = fold_type) %>%
  arrange(Method, AOD)

rmse_best <- 
  rbind(
    all_rmse %>% filter(fold_type == "Random") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    miss_rmse %>% filter(fold_type == "Random") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    all_rmse %>% filter(fold_type == "Spatial Cluster 1") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    miss_rmse %>% filter(fold_type == "Spatial Cluster 1") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    all_rmse %>% filter(fold_type == "Spatial Cluster 2") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2),
    miss_rmse %>% filter(fold_type == "Spatial Cluster 2") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2)
  ) %>%
  rename(Method = fold_type) %>%
  arrange(Method, AOD)

int_best <- 
  rbind(
    all_int %>% filter(fold_type == "Random") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    miss_int %>% filter(fold_type == "Random") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    all_int %>% filter(fold_type == "Spatial Cluster 1") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    miss_int %>% filter(fold_type == "Spatial Cluster 1") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    all_int %>% filter(fold_type == "Spatial Cluster 2") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2),
    miss_int %>% filter(fold_type == "Spatial Cluster 2") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2)
  ) %>%
  rename(Method = fold_type) %>%
  arrange(Method, AOD)

slope_best <- 
  rbind(
    all_slope %>% filter(fold_type == "Random") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    miss_slope %>% filter(fold_type == "Random") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    all_slope %>% filter(fold_type == "Spatial Cluster 1") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    miss_slope %>% filter(fold_type == "Spatial Cluster 1") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    all_slope %>% filter(fold_type == "Spatial Cluster 2") %>% mutate(AOD = "All") %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2),
    miss_slope %>% filter(fold_type == "Spatial Cluster 2") %>% rename(AOD = aod_missing) %>%
      select(fold_type, AOD, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2)
  ) %>%
  rename(Method = fold_type) %>%
  arrange(Method, AOD)

print(xtable(r2_best, digits = 3), file = "./pm25/tables/rf_st_r2.txt",
      include.rownames = F)
print(xtable(rmse_best, digits = 2), file = "./pm25/tables/rf_st_rmse.txt",
      include.rownames = F)
print(xtable(int_best, digits = 2), file = "./pm25/tables/rf_st_int.txt",
      include.rownames = F)
print(xtable(slope_best, digits = 2), file = "./pm25/tables/rf_st_slope.txt",
      include.rownames = F)

saveRDS(
  list(r2_table = r2_best,
       rmse_table = rmse_best,
       int_table = int_best,
       slope_table = slope_best),
  "./pm25/tables/st_tables.rds")


#### Regional Summaries ####
library(sf)
cmaq_region <- readRDS("./cmaq/region_assignment.rds") %>%
  st_drop_geometry(.) %>%
  select(cmaq_id, state_name, region)

region_stack <- left_join(df_stack, cmaq_region,
                          by = c("cmaq_id" = "cmaq_id")) %>%
  select(fold_type:fold, state_name, region,
         contains(c("a_1", "a_2")))

region_all_r2 <- region_stack %>%
  group_by(fold_type, region) %>%
  summarize(
    across(p1a_1:p5a_2, ~ getR2(pm25_value, .x))    
  ) %>%
  ungroup(.) %>%
  rename(Method = fold_type) %>%
  arrange(Method, region)

region_all_rmse <- region_stack %>%
  group_by(fold_type, region) %>%
  summarize(
    across(p1a_1:p5a_2, ~ getRMSE(pm25_value, .x))    
  ) %>%
  ungroup(.) %>%
  rename(Method = fold_type) %>%
  arrange(Method, region)

region_all_int <- region_stack %>%
  group_by(fold_type, region) %>%
  summarize(
    across(p1a_1:p5a_2, ~ getInt(pm25_value, .x))    
  ) %>%
  ungroup(.) %>%
  rename(Method = fold_type) %>%
  arrange(Method, region)

region_all_slope <- region_stack %>%
  group_by(fold_type, region) %>%
  summarize(
    across(p1a_1:p5a_2, ~ getSlope(pm25_value, .x))    
  ) %>%
  ungroup(.) %>%
  rename(Method = fold_type) %>%
  arrange(Method, region)

# Summarize -- we focus on the "best" fits
print(mgrid)

region_r2 <- 
  rbind(
    region_all_r2 %>% filter(Method == "Random") %>%
      select(Method, region, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    region_all_r2 %>% filter(Method == "Spatial Cluster 1") %>%
      select(Method, region, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    region_all_r2 %>% filter(Method == "Spatial Cluster 2") %>%
      select(Method, region, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2)
  ) %>%
  arrange(Method, region)

region_rmse <- 
  rbind(
    region_all_rmse %>% filter(Method == "Random") %>%
      select(Method, region, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    region_all_rmse %>% filter(Method == "Spatial Cluster 1") %>%
      select(Method, region, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    region_all_rmse %>% filter(Method == "Spatial Cluster 2") %>%
      select(Method, region, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2)
  ) %>%
  arrange(Method, region)


region_int <- 
  rbind(
    region_all_int %>% filter(Method == "Random") %>%
      select(Method, region, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    region_all_int %>% filter(Method == "Spatial Cluster 1") %>%
      select(Method, region, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    region_all_int %>% filter(Method == "Spatial Cluster 2") %>%
      select(Method, region, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2)
  ) %>%
  arrange(Method, region)

region_slope <- 
  rbind(
    region_all_slope %>% filter(Method == "Random") %>%
      select(Method, region, p1a_1, p2a_1, p3a_2, p4a_2, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    region_all_slope %>% filter(Method == "Spatial Cluster 1") %>%
      select(Method, region, contains("a_2")) %>%
      rename(p1a = p1a_2, p2a = p2a_2, p3a = p3a_2, p4a = p4a_2, p5a = p5a_2),
    region_all_slope %>% filter(Method == "Spatial Cluster 2") %>%
      select(Method, region, p1a_1, p2a_1, p3a_1, p4a_1, p5a_2) %>%
      rename(p1a = p1a_1, p2a = p2a_1, p3a = p3a_1, p4a = p4a_1, p5a = p5a_2)
  ) %>%
  arrange(Method, region)

print(xtable(region_r2, digits = 3), file = "./pm25/tables/region_st_r2.txt",
      include.rownames = F)

print(xtable(region_rmse, digits = 2), file = "./pm25/tables/region_st_rmse.txt",
      include.rownames = F)

print(xtable(region_int, digits = 2), file = "./pm25/tables/region_st_int.txt",
      include.rownames = F)

print(xtable(region_slope, digits = 2), file = "./pm25/tables/region_st_slope.txt",
      include.rownames = F)

saveRDS(list(reg_r2 = region_r2, 
             reg_rmse = region_rmse,
             reg_int = region_int,
             reg_slope = region_slope),
        "./pm25/tables/regional_st_tables.rds")
