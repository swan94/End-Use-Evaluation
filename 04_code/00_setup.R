# ==========================================================================
# 00_setup.R
#
# Project: qFL-Chr.25 End-Use Evaluation
#
# Purpose:
#   1. Define the main project folder
#   2. Load required packages
#   3. Define all project paths
#   4. Create required output folders
#   5. Set global R options
#   6. Define reusable cleaning functions
#   7. Create the genotype classification key
#
# Project location:
#   C:/Users/sjw23099/Desktop/End-Use Evaluation
# ==========================================================================


# --------------------------------------------------------------------------
# 1. Define the project root
# --------------------------------------------------------------------------

project_root <- normalizePath(
  "C:/Users/sjw23099/Desktop/End-Use Evaluation",
  winslash = "/",
  mustWork = TRUE
)

message(
  "\nProject root detected:\n",
  project_root,
  "\n"
)


# --------------------------------------------------------------------------
# 2. Confirm expected project folders exist
# --------------------------------------------------------------------------

required_project_folders <- c(
  "01_manuscript",
  "02_data_raw",
  "03_data_processed",
  "04_code",
  "05_results",
  "06_submission"
)

missing_project_folders <- required_project_folders[
  !dir.exists(
    file.path(
      project_root,
      required_project_folders
    )
  )
]

if (length(missing_project_folders) > 0) {
  stop(
    paste0(
      "\nThe following project folders were not found:\n\n",
      paste(
        missing_project_folders,
        collapse = "\n"
      ),
      "\n\nExpected project root:\n",
      project_root
    )
  )
}


# --------------------------------------------------------------------------
# 3. Required packages
# --------------------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "readxl",
  "janitor",
  "lme4",
  "lmerTest",
  "emmeans",
  "broom.mixed",
  "performance",
  "patchwork",
  "knitr",
  "rmarkdown"
)


# Identify packages that are not installed
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]


# Stop and provide an installation command when packages are missing
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "\nThe following packages are not installed:\n\n",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n\nInstall them using:\n\n",
      'install.packages(c("',
      paste(
        missing_packages,
        collapse = '", "'
      ),
      '"))\n'
    )
  )
}


# --------------------------------------------------------------------------
# 4. Load packages
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  
  library(tidyverse)
  library(readxl)
  library(janitor)
  
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(broom.mixed)
  library(performance)
  
  library(patchwork)
  
  library(knitr)
  library(rmarkdown)
  
})


# --------------------------------------------------------------------------
# 5. Define project paths
# --------------------------------------------------------------------------

paths <- list(
  
  root = project_root,
  
  manuscript = file.path(
    project_root,
    "01_manuscript"
  ),
  
  data_raw = file.path(
    project_root,
    "02_data_raw"
  ),
  
  data_processed = file.path(
    project_root,
    "03_data_processed"
  ),
  
  code = file.path(
    project_root,
    "04_code"
  ),
  
  results = file.path(
    project_root,
    "05_results"
  ),
  
  submission = file.path(
    project_root,
    "06_submission"
  ),
  
  audit = file.path(
    project_root,
    "05_results",
    "audit"
  ),
  
  afis_results = file.path(
    project_root,
    "05_results",
    "AFIS"
  ),
  
  yarn_results = file.path(
    project_root,
    "05_results",
    "yarn"
  ),
  
  figures = file.path(
    project_root,
    "05_results",
    "figures"
  ),
  
  tables = file.path(
    project_root,
    "05_results",
    "tables"
  ),
  
  notebook_outputs = file.path(
    project_root,
    "05_results",
    "notebook_outputs"
  ),
  
  notebook_figures = file.path(
    project_root,
    "05_results",
    "notebook_outputs",
    "figures"
  )
  
)


# --------------------------------------------------------------------------
# 6. Create required output folders
# --------------------------------------------------------------------------

folders_to_create <- c(
  
  paths$data_processed,
  
  paths$audit,
  paths$afis_results,
  paths$yarn_results,
  paths$figures,
  paths$tables,
  
  paths$notebook_outputs,
  paths$notebook_figures
  
)

purrr::walk(
  folders_to_create,
  ~ dir.create(
    path = .x,
    recursive = TRUE,
    showWarnings = FALSE
  )
)


# --------------------------------------------------------------------------
# 7. Global R options
# --------------------------------------------------------------------------

options(
  
  stringsAsFactors = FALSE,
  
  dplyr.summarise.inform = FALSE,
  
  scipen = 999,
  
  contrasts = c(
    "contr.sum",
    "contr.poly"
  )
  
)


# Use Kenward-Roger degrees of freedom for lmer emmeans results
emmeans::emm_options(
  lmer.df = "kenward-roger"
)


# --------------------------------------------------------------------------
# 8. Helper function: clean identifiers
# --------------------------------------------------------------------------

clean_identifier <- function(x) {
  
  x |>
    as.character() |>
    stringr::str_trim() |>
    stringr::str_replace_all(
      "\\s+",
      ""
    ) |>
    dplyr::na_if("")
  
}


# --------------------------------------------------------------------------
# 9. Helper function: safely convert values to numeric
# --------------------------------------------------------------------------

safe_numeric <- function(x) {
  
  readr::parse_number(
    as.character(x),
    na = c(
      "",
      "NA",
      "N/A",
      ".",
      "-"
    )
  )
  
}


# --------------------------------------------------------------------------
# 10. Helper function: standardize genotype names
# --------------------------------------------------------------------------

standardize_genotype <- function(x) {
  
  cleaned <- x |>
    as.character() |>
    stringr::str_trim() |>
    stringr::str_replace_all(
      "[[:space:]_-]+",
      ""
    ) |>
    stringr::str_to_upper()
  
  dplyr::recode(
    
    cleaned,
    
    "AC107NEG" = "AC107Neg",
    "AC107POS" = "AC107Pos",
    
    "AC164NEG" = "AC164Neg",
    "AC164POS" = "AC164Pos",
    
    "ACSJ4" = "ACSJ4",
    "ACALASJ4" = "ACSJ4",
    
    "DP213NEG" = "DP213Neg",
    "DP213POS" = "DP213Pos",
    
    "DP263NEG" = "DP263Neg",
    "DP263POS" = "DP263Pos",
    
    "DP50" = "DP50",
    "DELTAPINE50" = "DP50",
    
    "GA428NEG" = "GA428Neg",
    "GA428POS" = "GA428Pos",
    
    "GA437NEG" = "GA437Neg",
    "GA437POS" = "GA437Pos",
    
    "GA089" = "GA089",
    "GA2004089" = "GA089",
    
    "PM304NEG" = "PM304Neg",
    "PM304POS" = "PM304Pos",
    
    "PM321NEG" = "PM321Neg",
    "PM321POS" = "PM321Pos",
    
    "PMHS26" = "PMHS26",
    "PAYMASTERHS26" = "PMHS26",
    
    "SEALAND883" = "SL883",
    "SL883" = "SL883",
    
    "GA201800" = "GA201800",
    "GA201808" = "GA201808",
    
    .default = NA_character_
    
  )
  
}


# --------------------------------------------------------------------------
# 11. Helper function: standardize AFIS processing stages
# --------------------------------------------------------------------------

standardize_stage <- function(x) {
  
  cleaned <- x |>
    as.character() |>
    stringr::str_trim() |>
    stringr::str_to_lower() |>
    stringr::str_replace_all(
      "[[:punct:]]",
      " "
    ) |>
    stringr::str_squish()
  
  dplyr::case_when(
    
    stringr::str_detect(
      cleaned,
      "bale"
    ) ~ "Bale",
    
    stringr::str_detect(
      cleaned,
      "open|opener"
    ) ~ "Opened",
    
    stringr::str_detect(
      cleaned,
      "card"
    ) ~ "Carded",
    
    stringr::str_detect(
      cleaned,
      "sliver"
    ) ~ "Sliver",
    
    TRUE ~ NA_character_
    
  )
  
}


# --------------------------------------------------------------------------
# 12. Helper function: confirm that a file exists
# --------------------------------------------------------------------------

check_file_exists <- function(file_path) {
  
  if (!file.exists(file_path)) {
    
    stop(
      paste0(
        "\nRequired file was not found:\n\n",
        file_path,
        "\n\nCheck the filename and confirm that the file is stored in ",
        "02_data_raw."
      )
    )
    
  }
  
  invisible(TRUE)
  
}


# --------------------------------------------------------------------------
# 13. Helper function: show Excel worksheet names
# --------------------------------------------------------------------------

show_excel_sheets <- function(file_path) {
  
  check_file_exists(file_path)
  
  workbook_sheets <- readxl::excel_sheets(
    file_path
  )
  
  message(
    "\nWorkbook: ",
    basename(file_path)
  )
  
  print(
    workbook_sheets
  )
  
  invisible(
    workbook_sheets
  )
  
}


# --------------------------------------------------------------------------
# 14. Genotype classification key
# --------------------------------------------------------------------------

genotype_key <- tibble::tribble(
  
  ~genotype,   ~background, ~nil_set, ~qtl_status,    ~material_type,
  
  "AC107Neg",  "Acala",     "AC107",   "Negative",     "NIIL",
  "AC107Pos",  "Acala",     "AC107",   "Positive",     "NIIL",
  
  "AC164Neg",  "Acala",     "AC164",   "Negative",     "NIIL",
  "AC164Pos",  "Acala",     "AC164",   "Positive",     "NIIL",
  
  "ACSJ4",     "Acala",     NA,        "Parent_Check", "Recurrent parent",
  
  
  "DP213Neg",  "Deltapine", "DP213",   "Negative",     "NIIL",
  "DP213Pos",  "Deltapine", "DP213",   "Positive",     "NIIL",
  
  "DP263Neg",  "Deltapine", "DP263",   "Negative",     "NIIL",
  "DP263Pos",  "Deltapine", "DP263",   "Positive",     "NIIL",
  
  "DP50",      "Deltapine", NA,        "Parent_Check", "Recurrent parent",
  
  
  "GA428Neg",  "Georgia",   "GA428",   "Negative",     "NIIL",
  "GA428Pos",  "Georgia",   "GA428",   "Positive",     "NIIL",
  
  "GA437Neg",  "Georgia",   "GA437",   "Negative",     "NIIL",
  "GA437Pos",  "Georgia",   "GA437",   "Positive",     "NIIL",
  
  "GA089",     "Georgia",   NA,        "Parent_Check", "Recurrent parent",
  
  
  "PM304Neg",  "Paymaster", "PM304",   "Negative",     "NIIL",
  "PM304Pos",  "Paymaster", "PM304",   "Positive",     "NIIL",
  
  "PM321Neg",  "Paymaster", "PM321",   "Negative",     "NIIL",
  "PM321Pos",  "Paymaster", "PM321",   "Positive",     "NIIL",
  
  "PMHS26",    "Paymaster", NA,        "Parent_Check", "Recurrent parent",
  
  
  "SL883",     "Donor",     NA,        "Parent_Check", "QTL donor",
  
  "GA201800",  "Check",     NA,        "Parent_Check", "Check",
  "GA201808",  "Check",     NA,        "Parent_Check", "Check"
  
)


# --------------------------------------------------------------------------
# 15. Set useful factor orders in the genotype key
# --------------------------------------------------------------------------

genotype_key <- genotype_key |>
  mutate(
    
    background = factor(
      background,
      levels = c(
        "Acala",
        "Deltapine",
        "Georgia",
        "Paymaster",
        "Donor",
        "Check"
      )
    ),
    
    qtl_status = factor(
      qtl_status,
      levels = c(
        "Negative",
        "Positive",
        "Parent_Check"
      )
    ),
    
    material_type = factor(
      material_type,
      levels = c(
        "NIIL",
        "Recurrent parent",
        "QTL donor",
        "Check"
      )
    )
    
  )


# --------------------------------------------------------------------------
# 16. Save genotype classification key
# --------------------------------------------------------------------------

readr::write_csv(
  genotype_key,
  file.path(
    paths$data_processed,
    "genotype_key.csv"
  ),
  na = ""
)


# --------------------------------------------------------------------------
# 17. Setup summary
# --------------------------------------------------------------------------

setup_summary <- tibble::tibble(
  
  item = c(
    "Project root",
    "Raw data folder",
    "Processed data folder",
    "Code folder",
    "Results folder",
    "Audit folder",
    "Notebook output folder",
    "Notebook figure folder"
  ),
  
  path = c(
    paths$root,
    paths$data_raw,
    paths$data_processed,
    paths$code,
    paths$results,
    paths$audit,
    paths$notebook_outputs,
    paths$notebook_figures
  ),
  
  exists = file.exists(
    c(
      paths$root,
      paths$data_raw,
      paths$data_processed,
      paths$code,
      paths$results,
      paths$audit,
      paths$notebook_outputs,
      paths$notebook_figures
    )
  )
  
)


print(
  setup_summary,
  n = Inf
)


message(
  "\nSetup completed successfully.",
  "\nNotebook outputs will be saved in:",
  "\n",
  paths$notebook_outputs,
  "\n"
)