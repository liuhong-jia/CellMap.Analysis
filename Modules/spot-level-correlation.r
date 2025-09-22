library(dplyr)
library(tidyr)

str_remove_spe <- function(string) gsub("[^[:alnum:]]", "_", string)


spot_correlation <- function(pred, truth, norm = TRUE) {
  pred  <- as.data.frame(pred,  check.names = FALSE)
  truth <- as.data.frame(truth, check.names = FALSE)
  
  if (!identical(rownames(pred), rownames(truth))) {
    common_rows <- intersect(rownames(pred), rownames(truth))
    if (length(common_rows) == 0) stop("No overlapping row names.")
    message("⚠ Row names differ — aligning by common rows.")
    pred  <- pred[common_rows, , drop = FALSE]
    truth <- truth[common_rows, , drop = FALSE]
  }
  
  if (!identical(colnames(pred), colnames(truth))) {
    common_cols <- intersect(colnames(pred), colnames(truth))
    if (length(common_cols) == 0) stop("No overlapping column names.")
    message("⚠ Column names differ — aligning by common columns.")
    pred  <- pred[, common_cols, drop = FALSE]
    truth <- truth[, common_cols, drop = FALSE]
  }
  
  cor_vec <- vapply(seq_len(nrow(pred)), function(i)
    cor(as.numeric(pred[i, ]), as.numeric(truth[i, ]), use = "pairwise.complete.obs"),
    numeric(1))
  
  if (norm) {
    rng <- range(cor_vec, na.rm = TRUE)
    if (diff(rng) != 0) {
      cor_vec <- (cor_vec - rng[1]) / diff(rng)
    }
  }
  
  list(mean = mean(cor_vec, na.rm = TRUE),
       median = median(cor_vec, na.rm = TRUE))
}

input_root  <- "../data"
result_root <- "../results"
out_root    <- "../score"

datasets <- list.dirs(input_root, recursive = FALSE, full.names = FALSE)

method_patterns <- list(
  CellMap       = "CellMap_result.txt$",
  RCTD          = "RCTD_result.txt$",
  DestVI        = "DestVI_result.txt$",
  Stereoscope   = "Stereoscope_result.txt$",
  CARD          = "CARD_result.txt$",
  Redeconve     = "Redeconve.result.txt$",
  Seurat        = "Seurat_result.txt$",
  Tangram       = "Tangram_result.txt$",
  Cell2location = "Cell2location_result.txt$",
  novoSpaRc     = "novoSpaRc_result.txt$",
  SpatialDWLS   = "SpatialDWLS_result.txt$",
  SpaOTsc       = "SpaOTsc.result.csv$",
  CytoSPACE     = "fractional_abundances_by_spot\\.csv$",
  SpatialDecon = "spatialdecon_result.csv"
)

result_list <- lapply(names(method_patterns), function(x) list())
names(result_list) <- names(method_patterns)

for (dataset in datasets) {
  cat("\n=== Dataset:", dataset, "===\n")
  
  input_dir  <- file.path(input_root,  dataset)
  result_dir <- file.path(result_root, dataset)
  
  truth_file <- file.path(input_dir, "true.matrix.rds")
  if (!file.exists(truth_file)) {
    warning("Missing true.matrix.rds in ", dataset); next
  }
  truth <- readRDS(truth_file)
  colnames(truth) <- str_remove_spe(colnames(truth))

  for (method in names(method_patterns)) {
    
    pattern <- method_patterns[[method]]
    files   <- list.files(result_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) next
    
    for (f in files) {
      cat("  ", method, "->", basename(f), "\n")
      
      df <- tryCatch({
		if (grepl("\\.csv$", f)) {
			read.csv(f, row.names = 1)
		} 
		} else {
			sep_used <- if (grepl("\t$", pattern)) "\t" else ","
		}		
		read.table(f, sep = sep_used, header = TRUE, row.names = 1)
		}
		}, error = function(e) { 
		warning("    read error"); 
		return(NULL) 
		})
      if (is.null(df)) next
      
      if (method == "Redeconve") df <- t(df)
      if (method == "Seurat") {
		colnames(df) <- gsub("^prediction\\.score\\.", "", colnames(df))
		}
      if (method == "SpatialDWLS") { rownames(df) <- df$cell_ID; df <- df[,-1] }
      if (method == "Cell2location") {
        colnames(df) <- gsub("^q05cell_abundance_w_sf_", "", colnames(df))
        df <- df / rowSums(df)
      }
	  
      colnames(df) <- str_remove_spe(colnames(df))
      common <- intersect(colnames(df), colnames(truth))
	  truth_sub <- truth[rownames(df),common, drop = FALSE]
      df_sub <- df[, common, drop = FALSE]
      
      if (ncol(df_sub) == 0) { warning("    no overlap"); next }
      
      corr <- tryCatch(spot_correlation(df_sub, truth_sub), error = function(e) NULL)
      if (is.null(corr)) next
      
      result_list[[method]][[length(result_list[[method]]) + 1]] <-
        data.frame(Dataset = dataset,
                   File    = basename(f),
                   MeanCorrelation   = corr$mean,
                   MedianCorrelation = corr$median)
    }
  }
}



for (method in names(result_list)) {
  df <- bind_rows(result_list[[method]])
  if (nrow(df) == 0) next
  
  out_file <- file.path(out_root, paste0(method, "_correlation.csv"))
  write.csv(df, out_file, row.names = FALSE)
  cat("✔ Saved:", out_file, "\n")
}


combined_df <- do.call(rbind, lapply(names(result_list), function(method) {
  df <- bind_rows(result_list[[method]])
  if (nrow(df) == 0) return(NULL)
  
  df %>%
    group_by(Dataset) %>%
    slice(1) %>%  
    ungroup() %>%
    mutate(Method = method) %>%
	select(Dataset, Method, MeanCorrelation)
    ##select(Dataset, Method, MedianCorrelation)
	
}))

wide_df <- combined_df %>%
  tidyr::pivot_wider(names_from = Method, values_from = MeanCorrelation)

final_outfile <- file.path(out_root, "MeanCorrelation.csv")
write.csv(wide_df, final_outfile, row.names = FALSE)
cat("✔ Final combined file saved:", final_outfile, "\n")



