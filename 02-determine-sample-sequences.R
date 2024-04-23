# Perform sample inference and merge reads in demultiplexed, primer-less
# and filtered fastq files 
# Author: Márcio Martins (marciomartinsred@gmail.com)
# Date: 2022-05-24

# Install and load required packages --------------------------------------
packages <- c("BiocManager", "optparse", "dada2", "stringr", 
              "ggplot2", "dplyr", "tidyr", "parallel", "patchwork")

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
        help = "Path for directory containing your demultiplexed fastq files."
    ),
    make_option(
        c("-e", "--error-function"),
        type = "character",
        default = "loessErrfun",
        help = "Function used to estimate errors for sample inference. For binned quality data, such as NovaSeq, `custom-loess` is recommended."
    ),
    make_option(
        c("-o", "--output"),
        type = "character",
        help = "Folder to store your filtered fastq files, sequence table and diagnostics."
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)


# Prepare list of samples and necessary files -----------------------------
# Create directory to export data
dir.create(opts$output, recursive = TRUE)

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
files$output_file <- file.path(opts$output,"filtered-fastq", files$filename)
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
plot_path <- file.path(opts$output, "00-foward-read-quality.png")
ggsave(plot_path, qual_plot)
browseURL(plot_path)

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

# Save and display quality plot
plot_path <- file.path(opts$output, "00-reverse-read-quality.png")
ggsave(plot_path, qual_plot)
browseURL(plot_path)

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
    maxEE = Inf,
    multithread = parallel::detectCores() - 1,
    verbose = TRUE
)

write.csv(
    filtered_output,
    file.path(opts$output, "00-filtering-stats.csv")
)


# Estimate errors in data - modified due to binned qualities --------------
# NovaSeq bins nucleotide quality to save storage space, which interferes
# with the loess function used to estimate error. See discussion in
# https://github.com/benjjneb/dada2/issues/791#issuecomment-502256869
# https://github.com/benjjneb/dada2/issues/1307
# Solution used here is copy-pasted from Option 4 in :
# https://github.com/benjjneb/dada2/issues/1307#issuecomment-957680971

loessErrfun_modified <- function(trans) {
    qq <- as.numeric(colnames(trans))
    est <- matrix(0, nrow=0, ncol=length(qq))
    for(nti in c("A","C","G","T")) {
        for(ntj in c("A","C","G","T")) {
            if(nti != ntj) {
                errs <- trans[paste0(nti,"2",ntj),]
                tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
                rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
                rlogp[is.infinite(rlogp)] <- NA
                df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
                
                # original
                # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
                # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
                # #        mod.lo <- loess(rlogp ~ q, df)
                
                # jonalim's solution
                # https://github.com/benjjneb/dada2/issues/938
                mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),
                                degree = 1, span = 0.95)
                
                pred <- predict(mod.lo, qq)
                maxrli <- max(which(!is.na(pred)))
                minrli <- min(which(!is.na(pred)))
                pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
                pred[seq_along(pred)<minrli] <- pred[[minrli]]
                est <- rbind(est, 10^pred)
            } # if(nti != ntj)
        } # for(ntj in c("A","C","G","T"))
    } # for(nti in c("A","C","G","T"))
    
    # HACKY
    MAX_ERROR_RATE <- 0.25
    MIN_ERROR_RATE <- 1e-7
    est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
    est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
    
    # enforce monotonicity
    # https://github.com/benjjneb/dada2/issues/791
    estorig <- est
    est <- est %>%
        data.frame() %>%
        mutate_all(funs(case_when(. < X40 ~ X40,
                                  . >= X40 ~ .))) %>% as.matrix()
    rownames(est) <- rownames(estorig)
    colnames(est) <- colnames(estorig)
    
    # Expand the err matrix with the self-transition probs
    err <- rbind(1-colSums(est[1:3,]), est[1:3,],
                 est[4,], 1-colSums(est[4:6,]), est[5:6,],
                 est[7:8,], 1-colSums(est[7:9,]), est[9,],
                 est[10:12,], 1-colSums(est[10:12,]))
    rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
    colnames(err) <- colnames(trans)
    # Return
    return(err)
}

if(opts$`error-function` == "custom-loess") {
    error_func <- loessErrfun_modified
} else {
    error_func <- loessErrfun
}

fwd_samples_vec <- fwd_samples$output_file
names(fwd_samples_vec) <- fwd_samples$sample

rev_samples_vec <- rev_samples$output_file
names(rev_samples_vec) <- rev_samples$sample

message("Performing sequence inference.")
# Learn forward error rates
fwd_error <- learnErrors(
    fwd_samples_vec, 
    errorEstimationFunction = error_func,
    multithread=TRUE, 
    randomize = TRUE
)

# Learn reverse error rates
rev_error <- learnErrors(
    rev_samples_vec,
    errorEstimationFunction = error_func,
    multithread=TRUE,
    randomize = TRUE
)

# Save and display quality plot
plot_path <- file.path(opts$output, "01-error-model-fit.png")
error_plot <- (
    plotErrors(fwd_error) + 
        labs(title = "Foward sequence error model")) / (
    plotErrors(rev_error) +
         labs(title = "Reverse sequence error model")
    )
ggsave(plot_path, error_plot, height = 14)


# Perform sample inference ------------------------------------------------
# Dereplicate amplicon sequences to improve computational speed of inference
fwd_derep <- derepFastq(fwd_samples_vec)
rev_derep <- derepFastq(rev_samples_vec)

# Perform inference
fwd_dada <- dada(
    derep = fwd_derep, 
    err   = fwd_error,
    pool  = TRUE,
    multithread = TRUE
)

rev_dada <- dada(
    derep = rev_derep, 
    err   = rev_error,
    pool  = TRUE,
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
write.csv(
    seqtab,
    file.path(opts$output, "02-sequence-table.csv")
)

seqtab_nochimera <- removeBimeraDenovo(
    seqtab, 
    method="consensus", 
    multithread=TRUE
)

write.csv(
    seqtab_nochimera,
    file.path(opts$output, "02-sequence-table-nochimeras.csv")
)

# Save the sequence table as RDS for further analysis in R
saveRDS(
    seqtab_nochimera,
    file.path(opts$output, "02-sequence-table-nochimeras.rds")
)

# Export assortment of diagnostics ----------------------------------------

# Histogram of sequence length
png(file.path(opts$output, "02-sequence-length-histogram.png"))
hist(nchar(getSequences(seqtab)), main = "Sequence length histogram (chimeras removed)")
dev.off()

# How the number of sequences evolved along the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(
    filtered_output, 
    sapply(fwd_dada, getN), 
    sapply(rev_dada, getN), 
    sapply(merged, getN), 
    rowSums(seqtab))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- fwd_samples$sample

write.csv(
    track,
    file.path(opts$output, "02-sequence-number-evolution.csv")
)


track <- as.data.frame(track)
track$sample <- rownames(track)

track_long <- track |>
    mutate(
        depth = as.numeric(str_replace(
            str_extract(sample, "-[\\d]*"),
            "-", 
            ""
        )
        )) |>
    pivot_longer(
        cols = input:nonchim
    ) |>
    mutate(
        name = factor(name, levels = names(track))
    )

track_plot <- ggplot(track_long) +
    geom_line(
        aes(x = name,
            y = log10(value),
            group = sample)
    ) +
    labs(title = "Evolution of sequence number per sample along pipeline") +
    theme_bw()

ggsave(
    file.path(opts$output, "02-sequence-number-evolution.png"),
    track_plot
)



