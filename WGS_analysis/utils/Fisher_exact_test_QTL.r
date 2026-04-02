
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(rlang)
  library(stringr)
  library(tibble)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -------- Paths --------
base_path <- "~/Dropbox/Jiawan/JS_paper_summary/iPSC_paper/iPSC data summary/11_SUBMISSION/Nature communications submission/Reviews/analysis/common_effects/isogenic"
cell_types <- c("ipsc","npc","imn","mn")

# -------- File resolvers (flexible names) --------

# DEGs (padj0.05) per cell type
deg_file_for <- function(ct) {
  # canonical name first
  cand <- file.path(base_path, sprintf("%s_isogenic_deseq2_padj0.05.csv", ct))
  if (file.exists(cand)) return(cand)
  # fallback: any file with ct + 'padj0.05'
  hits <- list.files(
    base_path, full.names = TRUE, ignore.case = TRUE,
    pattern = sprintf("^%s_.*deseq2.*padj0.*\\.(csv|tsv|csn)$", ct)
  )
  if (length(hits)) hits[1] else NA_character_
}

# Background/population gene list per ct (supports .csv / .tsv / .csn and ...deseq2.*)
pop_file_for <- function(ct) {
  pats <- c(
    sprintf("^%s_isogenic_deseq2\\.(csv|tsv|csn)$", ct),
    sprintf("^%s_isogenic.*deseq2.*\\.(csv|tsv|csn)$", ct),
    sprintf("^%s_.*deseq2.*\\.(csv|tsv|csn)$", ct)
  )
  hits <- unlist(lapply(
    pats,
    function(p) list.files(base_path, pattern = p, full.names = TRUE, ignore.case = TRUE)
  ))
  if (length(hits)) hits[1] else NA_character_
}

# Slopes/variant-joined file per ct (contains log2FoldChange + GTEX/Psych slopes)
slopes_file_for <- function(ct) {
  cand <- c(
    file.path(base_path, sprintf("%s_isogenic_DEGs_with_variant.csv", ct)),
    file.path(base_path, sprintf("%s_isogenic_DEGs_with_variant.csv", ct))
  )
  if (!any(file.exists(cand))) {
    extra <- list.files(base_path,
                        pattern = sprintf("^%s.*DEGs_with_variant\\.csv$", ct),
                        full.names = TRUE)
    cand <- c(cand, extra)
  }
  cand[file.exists(cand)][1] %||% NA_character_
}

eqtl_file <- file.path(base_path, "common_variants_in_isogenic.csv")

# -------- Safe readers --------
safe_read_table <- function(fp) {
  # Most of your files are CSV even if extension is .csn, so try CSV first.
  tryCatch(readr::read_csv(fp, show_col_types = FALSE),
           error = function(e1) {
             tryCatch(readr::read_tsv(fp, show_col_types = FALSE),
                      error = function(e2) {
                        stop("Failed to read file as CSV or TSV: ", fp)
                      })
           })
}


strip_ver <- function(x) sub("\\.\\d+$","", trimws(as.character(x)))  # ENSG....13 → ENSG....
numfix <- function(x) as.numeric(suppressWarnings(as.character(x)))
sgn <- function(x) ifelse(is.na(x) | x == 0, NA_real_, ifelse(x > 0, 1, -1))

fix_gene_col <- function(df) {
  df <- df %>% rename_with(tolower)
  cand <- grep("ensembl", names(df), value = TRUE)
  if (length(cand) == 0) cand <- grep("^gene(_id)?$|geneid|gene.id", names(df), value = TRUE)
  if (length(cand) == 0) stop("No gene column found. Columns: ", paste(names(df), collapse=", "))
  cand <- cand[1]
  df %>% rename(ensembl_gene_id = !!sym(cand))
}

add_log2fc <- function(df) {
  nm_raw <- names(df)
  nm_key <- gsub("[^a-z0-9]", "", tolower(nm_raw))
  hits <- which(nm_key %in% c("log2foldchange","log2fc","lfc"))
  if (length(hits) == 0) return(df)
  lfc_col <- nm_raw[hits[1]]
  df %>% mutate(log2FoldChange = suppressWarnings(as.numeric(.data[[lfc_col]])))
}

# -------- Load eQTL set once (robust to column names) --------
if (!file.exists(eqtl_file)) stop("eQTL file not found: ", eqtl_file)
eqtl_raw0 <- safe_read_table(eqtl_file) %>% rename_with(tolower)

if ("gene_id" %in% names(eqtl_raw0)) {
  eqtl_genes <- strip_ver(eqtl_raw0$gene_id)
} else {
  eqtl_genes <- fix_gene_col(eqtl_raw0) %>% pull(ensembl_gene_id) %>% strip_ver()
}
eqtl_genes <- unique(eqtl_genes)
message("Loaded eQTL genes: ", length(eqtl_genes))

# -------- Diagnostics (what exists?) --------
diag <- tibble(
  cell_type     = cell_types,
  deg_file      = sapply(cell_types, deg_file_for),
  deg_exists    = file.exists(sapply(cell_types, deg_file_for)),
  pop_file      = sapply(cell_types, pop_file_for),
  pop_exists    = file.exists(sapply(cell_types, pop_file_for)),
  slopes_file   = sapply(cell_types, slopes_file_for),
  slopes_exists = !is.na(sapply(cell_types, slopes_file_for)) &
                  file.exists(sapply(cell_types, slopes_file_for))
)
print(diag)


compute_one <- function(ct, mode = c("overlap","at_least_one","both")) {
  mode <- match.arg(mode)

  deg_fp <- deg_file_for(ct)
  pop_fp <- pop_file_for(ct)

  if (!file.exists(deg_fp)) { message(sprintf("[%s|%s] Missing DEG file: %s", toupper(ct), mode, deg_fp)); return(NULL) }
  if (!file.exists(pop_fp)) { message(sprintf("[%s|%s] Missing population file: %s", toupper(ct), mode, pop_fp)); return(NULL) }

  # read population + deg
  df_population <- try(safe_read_table(pop_fp), silent = TRUE)
  if (inherits(df_population, "try-error")) { message(sprintf("[%s|%s] Read error (population).", toupper(ct), mode)); return(NULL) }
  df_deg <- try(safe_read_table(deg_fp), silent = TRUE)
  if (inherits(df_deg, "try-error")) { message(sprintf("[%s|%s] Read error (DEG).", toupper(ct), mode)); return(NULL) }

  df_population <- try(fix_gene_col(df_population), silent = TRUE)
  if (inherits(df_population, "try-error")) { message(sprintf("[%s|%s] No gene column (population).", toupper(ct), mode)); return(NULL) }
  df_deg <- try(fix_gene_col(df_deg), silent = TRUE)
  if (inherits(df_deg, "try-error")) { message(sprintf("[%s|%s] No gene column (DEG).", toupper(ct), mode)); return(NULL) }

  all_genes <- unique(strip_ver(df_population$ensembl_gene_id))
  deg_genes <- unique(strip_ver(df_deg$ensembl_gene_id))

  if (length(all_genes) == 0) { message(sprintf("[%s|%s] Population gene list empty.", toupper(ct), mode)); return(NULL) }
  if (length(deg_genes) == 0) { message(sprintf("[%s|%s] DEG gene list empty.", toupper(ct), mode)); return(NULL) }


  if (mode == "overlap") {
    a <- sum(deg_genes %in% eqtl_genes)

  } else {
    s_fp <- slopes_file_for(ct)
    if (is.na(s_fp) || !file.exists(s_fp)) {
      message(sprintf("[%s|%s] Slopes file not found (needed for direction modes).", toupper(ct), mode))
      return(NULL)
    }

    slopes_df <- try(readr::read_csv(s_fp, show_col_types = FALSE), silent = TRUE)
    if (inherits(slopes_df, "try-error")) {
      message(sprintf("[%s|%s] Read error (slopes): %s", toupper(ct), mode, s_fp))
      return(NULL)
    }


    names(slopes_df) <- sub("^hange$", "log2FoldChange", names(slopes_df))

    slopes_df <- slopes_df %>%
      rename_with(tolower) %>%
      fix_gene_col() %>%
      add_log2fc()

    if (!("log2foldchange" %in% names(slopes_df))) {
      message(sprintf("[%s|%s] No log2FoldChange column in slopes file: %s", toupper(ct), mode, s_fp))
      return(NULL)
    }

    # collect slopes with flexible names
    lfc   <- numfix(slopes_df$log2foldchange)
    gtex  <- numfix(slopes_df$gtex_slope %||% slopes_df$gtex_slop)
    psych <- numfix(slopes_df$psychencode_slope %||% slopes_df$psychencode_slop %||%
                    slopes_df$regression_slope %||% slopes_df$regression_slop)

    if (all(is.na(gtex)) & all(is.na(psych))) {
      message(sprintf("[%s|%s] Both slope columns missing/NA in: %s", toupper(ct), mode, s_fp))
      return(NULL)
    }

    s_lfc   <- sgn(lfc)
    s_gtex  <- sgn(gtex)
    s_psych <- sgn(psych)

    if (mode == "at_least_one") {
      agree <- (!is.na(s_lfc)) & ( (!is.na(s_gtex)  & s_gtex  == s_lfc) |
                                    (!is.na(s_psych) & s_psych == s_lfc) )
    } else { # both
      agree <- (!is.na(s_lfc)) & (!is.na(s_gtex)) & (!is.na(s_psych)) &
               (s_gtex == s_lfc) & (s_psych == s_lfc)
    }

    agree_genes <- unique(strip_ver(slopes_df$ensembl_gene_id[which(agree)]))
    agree_genes <- intersect(agree_genes, deg_genes)  # among DEGs only
    a <- length(agree_genes)
  }

  b <- sum(deg_genes %in% setdiff(all_genes, eqtl_genes))
  non_deg <- setdiff(all_genes, deg_genes)
  c <- sum(non_deg %in% eqtl_genes)
  d <- sum(non_deg %in% setdiff(all_genes, eqtl_genes))

  contingency <- matrix(c(a,b,c,d), nrow = 2, byrow = TRUE,
                        dimnames = list(DEG = c("Yes","No"),
                                        eQTL = c("Yes","No")))

  ft_greater  <- fisher.test(contingency, alternative = "greater")
  ft_twosided <- fisher.test(contingency, alternative = "two.sided")

  out <- tibble(
    cell_type          = ct,
    mode               = mode,
    test               = c("greater","two.sided"),
    a = a, b = b, c = c, d = d,
    DEG_total          = length(deg_genes),
    eQTL_total         = length(eqtl_genes),
    background_total   = length(all_genes),
    odds_ratio         = c(unname(ft_greater$estimate),  unname(ft_twosided$estimate)),
    conf_low           = c(ft_greater$conf.int[1],       ft_twosided$conf.int[1]),
    conf_high          = c(ft_greater$conf.int[2],       ft_twosided$conf.int[2]),
    p_value            = c(ft_greater$p.value,           ft_twosided$p.value)
  )

  # save per-ct per-mode
  out_name <- file.path(base_path, sprintf("%s_isogenic_eQTL_DEG_%s.csv", ct, mode))
  readr::write_csv(out, out_name)
  message(sprintf("Saved -> %s", out_name))

  out
}

# -------- Run all --------
all_results <- list()
for (ct in cell_types) {
  for (mode in c("overlap","at_least_one","both")) {
    message(sprintf("=== %s | %s ===", toupper(ct), mode))
    res <- try(compute_one(ct, mode), silent = TRUE)
    if (!inherits(res, "try-error") && !is.null(res)) {
      all_results[[length(all_results)+1]] <- res
    } else {
      message(sprintf("[SKIP] %s | %s -> see messages above.", ct, mode))
    }
  }
}

if (length(all_results) == 0) {
  message(
    "\nNo results produced. Check:\n",
    " - DEG files exist & have a gene column\n",
    " - Population files exist (now supports .csv/.tsv/.csn with 'deseq2' in name)\n",
    " - Direction modes need *_DEGs_with_variant.csv containing log2FoldChange + GTEX/Psych slopes\n",
    " - eQTL IDs vs Ensembl versions (versions are stripped here)\n"
  )
  stop("No results produced.")
}

combined <- dplyr::bind_rows(all_results)

# FDR within (cell_type, mode, test)
combined <- combined %>%
  group_by(cell_type, mode, test) %>%
  mutate(FDR_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()

combined_out <- file.path(base_path, "isogenic_eqtl_deg_enrichment_allmodes.csv")
readr::write_csv(combined, combined_out)
message(sprintf("Saved combined -> %s", combined_out))

# quick console snapshot (two-sided)
print(
  combined %>%
    filter(test == "two.sided") %>%
    arrange(cell_type, mode) %>%
    select(cell_type, mode, odds_ratio, conf_low, conf_high, p_value, FDR_BH)
)



#plot

library(dplyr)
library(ggplot2)
library(dplyr)







path <- "~/Dropbox/Jiawan/JS_paper_summary/iPSC_paper/iPSC data summary/11_SUBMISSION/Nature communications submission/Reviews/analysis/common_effects/fisher_exact/all"

files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)

summary <- files %>%
  lapply(function(f) {
    df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    df <- df %>% dplyr::filter(test == "two.sided")
    df$source_file <- basename(f)
    return(df)
  }) %>%
  bind_rows()

# Add cell type from filename (before first underscore)
summary$cell_type <- sub("_.*", "", summary$source_file)
print(summary)

summary <- summary %>%
  mutate(FDR = p.adjust(p_value, method = "BH"),
         star = case_when(
           FDR < 0.05 ~ "**",
           p_value < 0.05 ~ "*",
           TRUE ~ ""
         ))

summary$metric <- "eQTL enrichment"

p <- ggplot(summary, aes(x = cell_type, y = metric, fill = odds_ratio)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = star), fontface = "bold", vjust = 0.5, size = 5) +
  # center color map at OR = 1
  scale_fill_gradient2(
    low = "#2c7bb6", mid = "white", high = "#d7191c",
    midpoint = 1, name = "Odds ratio"
  ) +
  labs(x = NULL, y = NULL,
       title = "Enrichment of DEGs among eQTL-regulated genes",
       subtitle = "* p < 0.05   ** FDR < 0.05 (BH)") +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p)

ggsave("isogenic_eQTL_DEG_enrichemnt_all.pdf", p, width = 6.18, height = 5.18, dpi = 300)
