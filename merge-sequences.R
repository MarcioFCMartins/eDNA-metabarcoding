# Perform sample inference and merge reads in demultiplexed, primer-less
# and filtered fastq files 
# Author: Márcio Martins (marciomartinsred@gmail.com)
# Date: 2022-05-24

# Install and load required packages --------------------------------------
packages <- c("BiocManager", "optparse", "dada2", "stringr", "ggplot2", "dplyr")

installed_packages <- packages %in% rownames(installed.packages())
# Try to install packages from CRAN
if (any(installed_packages == FALSE)) {
    invisible(install.packages(packages[!installed_packages]))
}
# Try to install packages from Bioconductor
if (any(installed_packages == FALSE)) {
    invisible(BiocManager::install(packages[!installed_packages]))
}

if (any(installed_packages == FALSE)) {
    missing_packages <- paste(
        packages[!installed_packages],
        collapse = "\n"
    )
    stop("The following packages could not be installed:\n")
}

invisible(lapply(packages, library, character.only = TRUE))
set.seed(123)
# Read arguments ----------------------------------------------------------
option_list <- list(
    make_option(
        c("-d", "--data-directory"),
        type = "character",
        help = "Path for directory containing your demultiplexed fastq files"
    ),
    make_option(
        c("-o", "--output"),
        type = "character",
        help = "Folder to store your clean files and metadata"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)


# Prepare list of samples -------------------------------------------------

files <- data.frame(
    "path" = list.files(
        path = opts$`data-directory`,
        pattern = "\\.fastq\\.gz",
        recursive = TRUE,
        full.names = TRUE
    )
)

files$filename <- basename(files$path)
# Check if a file is reverse or forward read
files$forward <- str_detect(files$filename, "_R1_")
# Get sample name
files$sample <- str_extract(files$filename, "[^_]*")
# Prepare paths for output file
files$output_file <- file.path(opts$output, files$filename)

# Check if there's 2 files for every sample
# This check is rudimentary and doesn't see if every sample name has a file
if (length(files$forward[files$forward]) != length(files$forward[files$forward])) {
    stop("The number of reverse and forward files does not match.")
}

fwd_samples <- files[files$forward, ]
rev_samples <- files[!files$forward, ]

fwd_samples_vec <- fwd_samples$path
names(fwd_samples_vec) <- fwd_samples$sample

rev_samples_vec <- rev_samples$path
names(rev_samples_vec) <- rev_samples$sample


# Perform sample inference ------------------------------------------------
# Learn forward error rates
fwd_error <- learnErrors(
    fwd_samples_vec, 
    nbases=1e8, 
    multithread=TRUE, 
    randomize = TRUE)

# Learn reverse error rates
rev_error <- learnErrors(
    rev_samples_vec, 
    nbases=1e8, 
    multithread=TRUE,
    randomize = TRUE)

fwd_dada <- dada(
    derep = fwd_samples_vec, 
    err   = fwd_error,
    pool  = "pseudo",
    multithread = TRUE)

rev_dada <- dada(
    derep = rev_samples_vec, 
    err   = rev_error,
    pool  = "pseudo",
    multithread = TRUE)

# Merge sequences ---------------------------------------------------------

merged <- mergePairs(
    fwd_dada,
    fwd_samples_vec, 
    rev_dada,
    rev_samples_vec,
    verbose=TRUE
)


# Sequence table ----------------------------------------------------------
seqtab <- makeSequenceTable(merged)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

