# Perform sample inference and merge reads in demultiplexed, primer-less
# and filtered fastq files 
# Author: Márcio Martins (marciomartinsred@gmail.com)
# Date: 2022-05-24

# Install and load required packages --------------------------------------
packages <- c("BiocManager", "optparse", "dada2", "stringr", 
              "ggplot2", "dplyr", "tidyr")

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
files <- files[order(files$sample), ]

# Check if there's 2 files for every sample
# This check is rudimentary and doesn't see if every sample name has a file
if (length(files$forward[files$forward]) != length(files$forward[files$forward])) {
    stop("The number of reverse and forward files does not match.")
}

fwd_samples <- files[files$forward, ]
rev_samples <- files[!files$forward, ]

# Preview forward read quality and determine cutoff  ----------------------
message("Preparing preview of forward read quality. Please hold.")
message("If you are using this program over ssh, make sure X11 fowarding is on!")
message("The plot might take a while to be shown.")

# Create quality plot
qual_plot <- plotQualityProfile(fwd_samples$path, aggregate = TRUE) +
    labs(
        title = "Forward read sequence quality",
        subtitle = "Select cutoff point"
    )

# Calculate average quality for each sequence length
qual_stats <- qual_plot$data |>
    group_by(Cycle) |>
    summarise(avg = weighted.mean(Score, Count))

# Display quality plot
plot_temp <- tempfile("plot", fileext = ".png") # temporary file for plot
ggsave(plot_temp, qual_plot)
browseURL(plot_temp)

# Wait for user to input the desired length
message("Please input the foward read cutoff point.")
message(
    paste0(
        "The average nucleotide quality dips below 30 for the ",
        "first time in position ",
        qual_stats$Cycle[min(which(qual_stats$avg < 30))],
        "."
    )
)

fwd_cut <- scan("stdin", integer(), n = 1)


# Preview reverse read quality and determine cutoff -----------------------
message("Preparing preview of reverse read quality. Please hold.")
message("If you are using this program over ssh, make sure X11 fowarding is on!")
message("The plot might take a while to be shown.")

# Create quality plot
qual_plot <- plotQualityProfile(rev_samples$path, aggregate = TRUE) +
    labs(
        title = "Reverse read sequence quality",
        subtitle = "Select cutoff point"
    )

# Calculate average quality for each sequence length
qual_stats <- qual_plot$data |>
    group_by(Cycle) |>
    summarise(avg = weighted.mean(Score, Count))

# Display quality plot
plot_temp2 <- tempfile("plot2", fileext = ".png") # temporary file for plot
ggsave(plot_temp2, qual_plot)
browseURL(plot_temp2)

message("Please input the reverse read cutoff point.")
message(
    paste0(
        "The average nucleotide quality dips below 30 for the ",
        "first time in position ",
        qual_stats$Cycle[min(which(qual_stats$avg < 30))],
        "."
    )
)

rev_cut <- scan("stdin", integer(), n = 1)


# Perform sequence filtering ----------------------------------------------
# Have to use matchID - https://github.com/benjjneb/dada2/issues/737

message("Filtering sequences. Please hold")
filtered_output <- dada2::filterAndTrim(
    fwd = fwd_samples$path,
    filt = fwd_samples$output_file,
    rev = rev_samples$path,
    filt.rev = rev_samples$output_file,
    truncLen = c(fwd_cut, rev_cut),
    matchIDs = TRUE,
    multithread = TRUE,
    verbose = TRUE
)

write.csv(
    filtered_output,
    file.path(opts$output, "0000-filtering-stats.csv")
)

# Perform sample inference ------------------------------------------------
fwd_samples_vec <- fwd_samples$output_file
names(fwd_samples_vec) <- fwd_samples$sample

rev_samples_vec <- rev_samples$output_file
names(rev_samples_vec) <- rev_samples$sample

message("Performing sequence inference.")
# Learn forward error rates
fwd_error <- learnErrors(
    fwd_samples_vec, 
    nbases=1e8, 
    multithread=TRUE, 
    randomize = TRUE
)

# Learn reverse error rates
rev_error <- learnErrors(
    rev_samples_vec, 
    nbases=1e8, 
    multithread=TRUE,
    randomize = TRUE
)


# Dereplicate amplicon sequences to improve computational speed of inference
fwd_derep <- derepFastq(fwd_samples_vec)
rev_derep <- derepFastq(rev_samples_vec)

# Perform inference
fwd_dada <- dada(
    derep = fwd_derep, 
    err   = fwd_error,
    pool  = "pseudo",
    multithread = TRUE
)

rev_dada <- dada(
    derep = rev_derep, 
    err   = rev_error,
    pool  = "pseudo",
    multithread = TRUE
)

# Merge sequences ---------------------------------------------------------
message("Merging paired-end sequences")
merged <- mergePairs(
    fwd_dada,
    fwd_derep, 
    rev_dada,
    rev_derep,
    verbose=TRUE
)


# Sequence table ----------------------------------------------------------
seqtab <- makeSequenceTable(merged)



# Assign taxonomy script starts here --------------------------------------
###########################################################################

seqtab <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)


# TEST, MOVE LATER --------------------------------------------------------

tax <- assignTaxonomy(
    seqtab,
    "./Onedrive/700-in-progress-projects/706 - Arrabida restored seagrass meadows/2_data/edna-sequencing-processing/silva_132.18s.99_rep_set.dada2.fa.gz",
    multithread = TRUE,
    verbose = TRUE)
