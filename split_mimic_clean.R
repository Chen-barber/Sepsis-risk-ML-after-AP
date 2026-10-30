rm(list = ls())
gc()

library(readr)
library(caret)

data_clean <- "data_clean"
data_sets <- "data_sets"
dir.create(data_sets, showWarnings = FALSE, recursive = TRUE)

path_in <- file.path(data_clean, "mimiciv_imputed.csv")
path_wmu <- file.path(data_clean, "wmu_imputed.csv")
path_train <- file.path(data_sets, "mimic_train.csv")
path_val <- file.path(data_sets, "mimic_internal_val.csv")
path_all <- file.path(data_sets, "mimic_all.csv")
path_ext <- file.path(data_sets, "wmu_external_val.csv")

outcome <- "Sepsis3"
train_frac <- 0.7
split_seed <- 341L

stopifnot(file.exists(path_in), file.exists(path_wmu))
data <- read_csv(path_in, show_col_types = FALSE)
y <- as.integer(as.numeric(data[[outcome]]))
stopifnot(!anyNA(y), all(y %in% c(0L, 1L)))
data[[outcome]] <- y

set.seed(split_seed)
idx_train <- createDataPartition(
  factor(y, levels = c(0, 1)),
  p = train_frac,
  list = FALSE
)[, 1]
train <- data[idx_train, , drop = FALSE]
val <- data[-idx_train, , drop = FALSE]

write.csv(train, path_train, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(val, path_val, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data, path_all, row.names = FALSE, fileEncoding = "UTF-8")
wmu <- read_csv(path_wmu, show_col_types = FALSE)
write.csv(wmu, path_ext, row.names = FALSE, fileEncoding = "UTF-8")
