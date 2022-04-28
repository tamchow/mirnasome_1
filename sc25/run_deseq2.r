options(repos = list(CRAN="http://cran.rstudio.com/"))
if(!suppressMessages(require("DESeq2", quietly = TRUE))) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2", lib = .libPaths()[1])
  # deps <- unlist(tools::package_dependencies('DESeq2', which=c('Depends','Imports')), use.names=T, recursive=TRUE)
  # install.packages(deps, lib = .libPaths()[1], repos = BiocManager::repositories())
}
if (!suppressMessages(require("stringr", quietly = TRUE))) {
  BiocManager::install("stringr", lib = .libPaths()[1])
}
require("DESeq2")
require("stringr")
suppressMessages(library(tools))
my.args <- commandArgs(trailingOnly = TRUE)
if (length(my.args) < 3) {
  stop(
    "No arguments provided.\nUsage:\n'RScript run_deseq2.r --args <count file path> <count metadata file path> <experimental treatment> <control treatment>'"
  )
}
my.countFile <- my.args[2]
my.countMetaFile <- my.args[3]
my.countData <- read.csv(my.countFile, header = TRUE, row.names = 1)
my.countData.meta <- read.csv(
  my.countMetaFile,
  row.names = 1,
  colClasses = c("character", "factor", "factor", "factor")
)
my.countData[] <- lapply(my.countData, function(x)
  if (is.numeric(x))
    as.integer(round(x))
  else
    x)
my.countData <- my.countData[, rownames(my.countData.meta)]
if (length(my.args) > 5) {
  my.normalizer <- my.countData[my.args[6],]
  my.countData[] <-
    apply(my.countData, 1, function(x)
      as.integer(round(x / my.normalizer * 1e4)))
}
my.sample_groups <- unique(my.countData.meta$group)
for(my.sample_group in my.sample_groups) {
  my.sample_group.meta <- my.countData.meta[my.countData.meta$group == my.sample_group, ]
  my.sample_group.meta$set <- droplevels(my.sample_group.meta$set)
  my.countData.group <- my.countData[, names(my.countData) %in% rownames(my.sample_group.meta)]
  my.countData.dds <- DESeqDataSetFromMatrix(
    countData = my.countData.group,
    colData = my.sample_group.meta[, names(my.sample_group.meta) != "group"],
    design = ~ set + treatment
  )
  my.countData.analysis <- DESeq(my.countData.dds, quiet = TRUE)
  my.countData.resultsNames <-resultsNames(my.countData.analysis)
  my.resultsFile <- paste0(tools::file_path_sans_ext(my.countFile), "_", my.sample_group, "_deseq2", ".csv")
  my.groups <- NA
  if (is.na(my.args[4]) || is.na(my.args[5])) {
    my.groups <-
      sort(unique(my.countData.meta$treatment), decreasing = TRUE)
    if (length(my.groups) != 2)
      stop(
        "Contrast groups could not be unambiguously estimated from the metadata.\nPlease provide contrast groups as 'RScript run_deseq2.r --args <count file path> <count metadata file path> <experimental treatment> <control treatment>'"
      )
  } else {
    my.groups <- c(my.args[4], my.args[5])
  }
  my.contrast <- c("treatment", as.character(my.groups))
  cat(
    str_interp(
      "Detected experimental group: ${my.groups[1]}\nDetected control group: ${my.groups[2]}\nIf this is incorrect, please specify groups manually after other args.\n"
    )
  )
  if (!all(
    str_interp("treatment_${my.groups[1]}_vs_${my.groups[2]}") %in% my.countData.resultsNames
  )) {
    stop(
      str_interp(
        "contrast groups are not consistent with resultsNames: ${paste(my.countData.resultsNames, collapse=',')}"
      )
    )
  }
  my.countData.results <- results(my.countData.analysis, contrast = my.contrast)
  write.csv(my.countData.results, my.resultsFile)
  cat(str_interp("Results saved to ${my.resultsFile}\n")
)
}