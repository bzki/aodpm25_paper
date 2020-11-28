library(testthat)
library(dplyr)
library(xtable)

#### Overall and by AOD missingness ####
daily_overall_list <- readRDS("./pm25/tables/daily_tables.rds")
st_overall_list <- readRDS("./pm25/tables/st_tables.rds")

do_r2 <- daily_overall_list$r2_table
do_rmse <- daily_overall_list$rmse_table
do_int <- daily_overall_list$int_table
do_slope <- daily_overall_list$slope_table

sto_r2 <- st_overall_list$r2_table
sto_rmse <- st_overall_list$rmse_table
sto_int <- st_overall_list$int_table
sto_slope <- st_overall_list$slope_table

expect_identical(do_r2$Method, sto_r2$Method)
expect_identical(do_r2$AOD, sto_r2$AOD)
expect_identical(do_rmse$Method, sto_rmse$Method)
expect_identical(do_rmse$AOD, sto_rmse$AOD)
expect_identical(do_int$Method, sto_int$Method)
expect_identical(do_int$AOD, sto_int$AOD)
expect_identical(do_slope$Method, sto_slope$Method)
expect_identical(do_slope$AOD, sto_slope$AOD)

co_r2 <- left_join(do_r2, sto_r2)
co_rmse <- left_join(do_rmse, sto_rmse)
co_int <- left_join(do_int, sto_int)
co_slope <- left_join(do_slope, sto_slope)

colnames(co_r2)[3:ncol(co_r2)] <- c(paste0("M", seq(5), "a"),
                                    paste0("M", seq(5), "b"))
colnames(co_rmse)[3:ncol(co_rmse)] <- c(paste0("M", seq(5), "a"),
                                    paste0("M", seq(5), "b"))
colnames(co_int)[3:ncol(co_int)] <- c(paste0("M", seq(5), "a"),
                                    paste0("M", seq(5), "b"))
colnames(co_slope)[3:ncol(co_slope)] <- c(paste0("M", seq(5), "a"),
                                    paste0("M", seq(5), "b"))

# Multiply co_r2 by 100
mx100 <- function(x) { return(x*100) }
co_r2_x100 <- co_r2 %>%
  mutate(across(M1a:M5b, mx100))

print(xtable(co_r2_x100, digits = 1), file = "./pm25/tables/overall_combine_r2.txt",
      include.rownames = F)
print(xtable(co_rmse, digits = 2), file = "./pm25/tables/overall_combine_rmse.txt",
      include.rownames = F)
print(xtable(co_int, digits = 2), file = "./pm25/tables/overall_combine_int.txt",
      include.rownames = F)
print(xtable(co_slope, digits = 2), file = "./pm25/tables/overall_combine_slope.txt",
      include.rownames = F)

#### Regional ####
rm(list=ls())
daily_reg_list <- readRDS("./pm25/tables/regional_daily_tables.rds")
st_reg_list <- readRDS("./pm25/tables/regional_st_tables.rds")

dr_r2 <- daily_reg_list$reg_r2
dr_rmse <- daily_reg_list$reg_rmse
dr_int <- daily_reg_list$reg_int
dr_slope <- daily_reg_list$reg_slope

str_r2 <- st_reg_list$reg_r2
str_rmse <- st_reg_list$reg_rmse
str_int <- st_reg_list$reg_int
str_slope <- st_reg_list$reg_slope

testthat::expect_identical(dr_r2$Method, str_r2$Method)
testthat::expect_identical(dr_r2$region, str_r2$region)
testthat::expect_identical(dr_rmse$Method, str_rmse$Method)
testthat::expect_identical(dr_rmse$region, str_rmse$region)
testthat::expect_identical(dr_int$Method, str_int$Method)
testthat::expect_identical(dr_int$region, str_int$region)
testthat::expect_identical(dr_slope$Method, str_slope$Method)
testthat::expect_identical(dr_slope$region, str_slope$region)

colnames(dr_r2)
colnames(str_r2)
cr_r2 <- left_join(dr_r2, str_r2)
cr_rmse <- left_join(dr_rmse, str_rmse)
cr_int <- left_join(dr_int, str_int)
cr_slope <- left_join(dr_slope, str_slope)

colnames(cr_r2)[3:ncol(cr_r2)] <- c(paste0("M", seq(5), "a"),
                                    paste0("M", seq(5), "b"))
colnames(cr_rmse)[3:ncol(cr_rmse)] <- c(paste0("M", seq(5), "a"),
                                        paste0("M", seq(5), "b"))
colnames(cr_int)[3:ncol(cr_int)] <- c(paste0("M", seq(5), "a"),
                                      paste0("M", seq(5), "b"))
colnames(cr_slope)[3:ncol(cr_slope)] <- c(paste0("M", seq(5), "a"),
                                          paste0("M", seq(5), "b"))

# Multiply co_r2 by 100
mx100 <- function(x) { return(x*100) }
cr_r2_x100 <- cr_r2 %>%
  mutate(across(M1a:M5b, mx100))

print(xtable(cr_r2_x100, digits = 1), file = "./pm25/tables/region_combine_r2.txt",
      include.rownames = F)
print(xtable(cr_rmse, digits = 2), file = "./pm25/tables/region_combine_rmse.txt",
      include.rownames = F)
print(xtable(cr_int, digits = 2), file = "./pm25/tables/region_combine_int.txt",
      include.rownames = F)
print(xtable(cr_slope, digits = 2), file = "./pm25/tables/region_combine_slope.txt",
      include.rownames = F)
