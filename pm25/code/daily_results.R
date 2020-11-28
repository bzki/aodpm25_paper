# PM25 analysis
# -- see cv_rf_daily.R, cv_rf_daily_m5.R
# Create tables and plots of results

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

library(testthat)
library(dplyr)
library(xtable)
sessionInfo()

comb_list <- readRDS("./pm25/cvresults/cv_rf_results.rds")
param <- comb_list$param
cvdf_list <- comb_list$cvdf_list
b1df_list <- comb_list$b1df_list
b2df_list <- comb_list$b2df_list

cv_stack1 <- do.call("rbind", cvdf_list)
b1_stack1 <- do.call("rbind", b1df_list)
b2_stack1 <- do.call("rbind", b2df_list)

# Second set of additional model results
comb_extra3 <- readRDS("./pm25/cvresults/cv_rf_results_m5.rds")
cv.stack3 <- do.call("rbind", comb_extra3$cvdf_list) %>%
  select(cmaq_id, day, p5a_1:p5c_4)
b1.stack3 <- do.call("rbind", comb_extra3$b1df_list) %>%
  select(cmaq_id,  day, p5a_1:p5c_4)
b2.stack3 <- do.call("rbind", comb_extra3$b2df_list) %>%
  select(cmaq_id,  day, p5a_1:p5c_4)

# Combine all predictions together
cv_stack <- cv_stack1 %>% 
  left_join(cv.stack3, by = c("cmaq_id" = "cmaq_id", "day" = "day"))
b1_stack <- b1_stack1 %>% 
  left_join(b1.stack3, by = c("cmaq_id" = "cmaq_id", "day" = "day"))
b2_stack <- b2_stack1 %>% 
  left_join(b2.stack3, by = c("cmaq_id" = "cmaq_id", "day" = "day"))
rm(comb_list, comb_extra3, cv.stack3, b1.stack3, b2.stack3,
   cvdf_list, b1df_list, b2df_list, cv_stack1, b1_stack1, b2_stack1)

#### Updated with additional results, simplified presentation ####
# - Focus on convolution results only 
# -- none of the other models for the time being
# This is model (b) and there are specifications 1 to 5.

combined_data <-
  rbind(cv_stack %>% rename(fold = fold_cv) %>% 
          mutate(fold_type = "Random"), 
        b1_stack %>% rename(fold = fold_b1) %>% 
          mutate(fold_type = "Spatial Cluster 1"), 
        b2_stack %>% rename(fold = fold_b2) %>%
          mutate(fold_type = "Spatial Cluster 2")) %>%
  mutate(AOD = ifelse(is.na(aod_value), "Missing", "Observed")) %>%
  select(cmaq_id:day, fold_type, AOD, 
         contains("b_"))

orig_df <- combined_data

combined_data <- combined_data %>%
  select(cmaq_id:fold_type, contains(c("_1", "_2", "_3", "_4")))

x_rmse <- combined_data %>%
  group_by(fold_type) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getRMSE(pm25_value, .x))    
  ) %>%
  ungroup(.)

x_r2 <- combined_data %>%
  group_by(fold_type) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getR2(pm25_value, .x))    
  ) %>%
  ungroup(.)

x_int <- combined_data %>%
  group_by(fold_type) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getInt(pm25_value, .x))    
  ) %>%
  ungroup(.)

x_slope <- combined_data %>%
  group_by(fold_type) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getSlope(pm25_value, .x))    
  ) %>%
  ungroup(.)

xm_rmse <- combined_data %>% 
  mutate(AOD = ifelse(is.na(aod_value), "Missing", "Observed")) %>%
  group_by(fold_type, AOD) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getRMSE(pm25_value, .x))    
  ) %>%
  ungroup(.)

xm_r2 <- combined_data %>% 
  mutate(AOD = ifelse(is.na(aod_value), "Missing", "Observed")) %>%
  group_by(fold_type, AOD) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getR2(pm25_value, .x))    
  ) %>%
  ungroup(.)

xm_int <- combined_data %>% 
  mutate(AOD = ifelse(is.na(aod_value), "Missing", "Observed")) %>%
  group_by(fold_type, AOD) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getInt(pm25_value, .x))    
  ) %>%
  ungroup(.)

xm_slope <- combined_data %>% 
  mutate(AOD = ifelse(is.na(aod_value), "Missing", "Observed")) %>%
  group_by(fold_type, AOD) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getSlope(pm25_value, .x))    
  ) %>%
  ungroup(.)

# Summarize the best mtry for each feature list and cross-validation type
mnum <- seq(5)
lnum <- c("b")
mgrid <- expand.grid(mnum = mnum, lnum = lnum, stringsAsFactors = F)
mgrid <- mgrid %>%
  mutate(mname = paste0(mnum, lnum))
attributes(mgrid)$out.attrs <- NULL
mgrid$cv = NA; mgrid$b1 = NA; mgrid$b2 <- NA
mgrid$cv_mtry = NA; mgrid$b1_mtry = NA; mgrid$b2_mtry = NA

for (mi in seq(nrow(mgrid))) {
  use_names <- paste0(mgrid$mnum[mi], mgrid$lnum[mi])
  print(use_names)
  mgrid$cv[mi] <- which.min(as.vector(x_rmse[1, ] %>% 
                                        select(contains(use_names))))
  mgrid$b1[mi] <- which.min(as.vector(x_rmse[2, ] %>% 
                                        select(contains(use_names))))
  mgrid$b2[mi] <- which.min(as.vector(x_rmse[3, ] %>% 
                                        select(contains(use_names))))
  
  mgrid$cv_mtry[mi] <- ifelse(mgrid$cv[mi] == 1, 4, 
                              ifelse(mgrid$cv[mi] == 2, 8,
                                     ifelse(mgrid$cv[mi] == 3, 12, 16)))
  mgrid$b1_mtry[mi] <- ifelse(mgrid$b1[mi] == 1, 4, 
                              ifelse(mgrid$b1[mi] == 2, 8,
                                     ifelse(mgrid$b1[mi] == 3, 12, 16)))
  mgrid$b2_mtry[mi] <- ifelse(mgrid$b2[mi] == 1, 4, 
                              ifelse(mgrid$b2[mi] == 2, 8,
                                     ifelse(mgrid$b2[mi] == 3, 12, 16)))
} 

print(mgrid_sub <- mgrid %>%
        mutate(mname = paste0(mnum, lnum)) %>%
        select(mname, cv_mtry, b1_mtry, b2_mtry, cv, b1, b2))

saveRDS(list(mgrid = mgrid, 
             x_rmse = x_rmse, x_r2 = x_r2, x_int = x_int, x_slope = x_slope), 
        "./pm25/cvresults/cvbest_daily.rds")

#### Scatter plots ####
# -- Construct for each of 5 models for both mtry = 4 and mtry = 8
df_cv <- combined_data %>% filter(fold_type == "Random")
df_b1 <- combined_data %>% filter(fold_type == "Spatial Cluster 1")
df_b2 <- combined_data %>% filter(fold_type == "Spatial Cluster 2")

png(filename = "./pm25/tables/scatter_daily_cvr_m3_mtry4.png")
plot(df_cv$p3b_1, df_cv$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()

png(filename = "./pm25/tables/scatter_daily_cvr_m3_mtry8.png")
plot(df_cv$p3b_2, df_cv$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()

png(filename = "./pm25/tables/scatter_daily_cvb1_m3_mtry4.png")
plot(df_b1$p3b_1, df_b1$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()

png(filename = "./pm25/tables/scatter_daily_cvb1_m3_mtry8.png")
plot(df_b1$p3b_2, df_b1$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()

png(filename = "./pm25/tables/scatter_daily_cvb2_m3_mtry4.png")
plot(df_b2$p3b_1, df_b2$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()

png(filename = "./pm25/tables/scatter_daily_cvb2_m3_mtry8.png")
plot(df_b2$p3b_2, df_b2$pm25_value, xlim = c(0, 70), ylim = c(0, 70),
     xlab = "CV PM2.5", ylab = "Observed PM2.5")
abline(a = 0, b = 1, col = "red")
dev.off()


# We pick mtry = 8 for all model b (convolution PM2.5) 
#   for all folds/models except model 5 -- set this to mtry = 16 for all folds
r2x_main <- 
  rbind(
    x_r2 %>% mutate(AOD = "All") %>% 
      select(fold_type, AOD, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4),
    xm_r2 %>% 
      select(fold_type, AOD, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4)
  ) %>%
  rename(Method = fold_type) %>%
  arrange(Method, AOD)

rmsex_main <- 
  rbind(
    x_rmse %>% mutate(AOD = "All") %>% 
      select(fold_type, AOD, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4),
    xm_rmse %>% 
      select(fold_type, AOD, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4)
  ) %>%
  rename(Method = fold_type) %>%
  arrange(Method, AOD)

intx_main <- 
  rbind(
    x_int %>% mutate(AOD = "All") %>% 
      select(fold_type, AOD, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4),
    xm_int %>% 
      select(fold_type, AOD, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4)
  ) %>%
  rename(Method = fold_type) %>%
  arrange(Method, AOD)

slopex_main <- 
  rbind(
    x_slope %>% mutate(AOD = "All") %>% 
      select(fold_type, AOD, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4),
    xm_slope %>% 
      select(fold_type, AOD, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4)
  ) %>%
  rename(Method = fold_type) %>%
  arrange(Method, AOD)

print(xtable(r2x_main, digits = 3), file = "./pm25/tables/rfdailyx_r2.txt",
      include.rownames = F)
print(xtable(rmsex_main, digits = 2), file = "./pm25/tables/rfdailyx_rmse.txt",
      include.rownames = F)
print(xtable(intx_main, digits = 2), file = "./pm25/tables/rfdailyx_int.txt",
      include.rownames = F)
print(xtable(slopex_main, digits = 2), file = "./pm25/tables/rfdailyx_slope.txt",
      include.rownames = F)

saveRDS(
  list(r2_table = r2x_main,
       rmse_table = rmsex_main,
       int_table = intx_main,
       slope_table = slopex_main),
  "./pm25/tables/daily_tables.rds")

#### Results by region ####
# Read in the file I prepared that has region based on point
library(sf)
cmaq_region <- readRDS("./cmaq/region_assignment.rds") %>%
  st_drop_geometry(.) %>%
  select(cmaq_id, state_name, region)

region_stack <- left_join(combined_data, cmaq_region,
                          by = c("cmaq_id" = "cmaq_id"))

reg_r2 <- region_stack %>%
  group_by(fold_type, region) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getR2(pm25_value, .x))    
  ) %>%
  ungroup(.) %>%
  select(fold_type, region, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4) %>%
  rename(Method = fold_type) %>%
  arrange(Method, region)

reg_rmse <- region_stack %>%
  group_by(fold_type, region) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getRMSE(pm25_value, .x))    
  ) %>%
  ungroup(.) %>%
  select(fold_type, region, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4) %>%
  rename(Method = fold_type) %>%
  arrange(Method, region)

reg_int <- region_stack %>%
  group_by(fold_type, region) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getInt(pm25_value, .x))    
  ) %>%
  ungroup(.) %>%
  select(fold_type, region, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4) %>%
  rename(Method = fold_type) %>%
  arrange(Method, region)

reg_slope <- region_stack %>%
  group_by(fold_type, region) %>%
  summarize(
    across(p1b_1:p5b_4, ~ getSlope(pm25_value, .x))    
  ) %>%
  ungroup(.) %>%
  select(fold_type, region, p1b_2, p2b_2, p3b_2, p4b_2, p5b_4) %>%
  rename(Method = fold_type) %>%
  arrange(Method, region)

print(xtable(reg_r2, digits = 3), file = "./pm25/tables/region_daily_r2.txt",
      include.rownames = F)

print(xtable(reg_rmse, digits = 2), file = "./pm25/tables/region_daily_rmse.txt",
      include.rownames = F)

print(xtable(reg_int, digits = 2), file = "./pm25/tables/region_daily_int.txt",
      include.rownames = F)

print(xtable(reg_slope, digits = 2), file = "./pm25/tables/region_daily_slope.txt",
      include.rownames = F)

saveRDS(list(reg_r2 = reg_r2, 
             reg_rmse = reg_rmse,
             reg_int = reg_int,
             reg_slope = reg_slope),
        "./pm25/tables/regional_daily_tables.rds")
