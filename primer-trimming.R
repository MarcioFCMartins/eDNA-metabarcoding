# Remove primers from demultiplexed fastq files based on mapping file
# provided by MrDNA
# Author: Márcio Martins (marciomartinsred@gmail.com)
# Date: 2022-05-19

# Install and load required packages --------------------------------------
packages <- c("BiocManager", "optparse", "stringr", 
              "doMC", "foreach", "parallel", "dada2")

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

n_cores <- parallel::detectCores()
n_cores <- ifelse(n_cores > 8, 8L, n_cores)
registerDoMC()

# Read arguments ----------------------------------------------------------
option_list <- list(
    make_option(
        c("-d", "--data-directory"),
        type = "character",
        help = "Path for directory containing your demultiplexed fastq files"
    ),
    make_option(
        c("-m", "--mapping"),
        type = "character",
        help = "File containing primer mapping (provided by MrDNA)"
    ),
    make_option(
        c("-o", "--output"),
        type = "character",
        help = "Folder to store your clean files and metadata"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts = parse_args(opt_parser)


# Prepare table with files ------------------------------------------------
files <- data.frame(
    "path" = list.files(
        path    = opts$`data-directory`,
        pattern = "\\.fastq\\.gz",
        recursive = TRUE,
        full.names = TRUE
    )
)

files$filename <- basename(files$path)

# Check if a file is reverse or foward read
files$foward      <- str_detect(files$filename, "_R1_")
# Get sample name
files$sample <- str_extract(files$filename, "[^_]*")


# Prepare table with primer used per sample -------------------------------
mapping <- read.delim(
    file = opts$mapping
)
names(mapping) <- c("sample", "barcode", "foward_primer", "barcode_name",
                    "reverse_primer", "projectname", "description")

# Primers for forward samples
primers_fwd <- mapping[, c("sample", "foward_primer")]
primers_fwd$foward <- TRUE
names(primers_fwd) <- c("sample", "primer", "foward")
# Primers for reverse samples
primers_rev <- mapping[, c("sample", "reverse_primer")]
primers_rev$foward <- FALSE
names(primers_rev) <- c("sample", "primer", "foward")

primers <- rbind(
    primers_fwd,
    primers_rev
)


# Remove primers ----------------------------------------------------------

files <- merge(
    files,
    primers,
    by = c("sample", "foward")
)

foreach(i = 1:nrow(files)) %dopar% {
    x <- files[i, ]
    input  <- x$path
    # Format output file name according to CAVASA :
    # sample identifier _ barcode sequence _ lane number
    # _ direction of read (i.e. R1 or R2) _ set number
    newname <- str_split(x$filename, "_", simplify = TRUE)
    newname <- paste(newname[1], "S1", newname[3], newname[4], newname[5],
                     sep = "_")

    primer <- x$primer

    dada2::removePrimers(
        fn      = input,
        fout    = file.path(opts$output, newname),
        primer.fwd = primer,
        orient  = FALSE,
        verbose = TRUE
    )
}



