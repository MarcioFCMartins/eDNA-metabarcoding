# Remove primers from demultiplexed fastq files based on mapping file
# provided by MrDNA
# Author: Márcio Martins (marciomartinsred@gmail.com)
# Date: 2022-05-20

# Install and load required packages --------------------------------------
packages <- c("BiocManager", "optparse", "stringr", 
              "doMC", "foreach", "parallel", "dada2")
installed_packages <- packages %in% rownames(installed.packages())

message("Setting up libraries")
sink(nullfile()) # suppress output
# Try to install packages from CRAN
if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
}
installed_packages <- packages %in% rownames(installed.packages())
# Try to install packages from Bioconductor
if (any(installed_packages == FALSE)) {
    BiocManager::install(packages[!installed_packages], update = FALSE)
}
sink() # end suppressing output

# Warn user if any package could not be installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
    missing_packages <- paste(
        packages[!installed_packages],
        collapse = "\n"
    )
    stop(
    	paste("The following packages could not be installed:\n",
    	packages[!installed_packages])
    )
} 

invisible(lapply(packages, library, character.only = TRUE))

# Use as many physical cores as available
n_cores <- parallel::detectCores(logical = FALSE) - 1 
n_cores <- ifelse(n_cores > 8, 8, n_cores)
registerDoMC(n_cores)

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

# Warn user if there are missing arguments
if(any(is.null(opts$output), 
       is.null(opts$mapping), 
       is.null(opts$`data-directory`))
) {
    stop(
        paste("There are missing arguments. This program requires:\n,
              --data-directory --mapping and --output arguments")
    )
}

# Prepare table with files ------------------------------------------------
# opts <- list(
#     "data-directory" = "~/Onedrive/700-in-progress-projects/706 - Arrabida restored seagrass meadows/2_data/raw_data/edna-sequencing/raw-data/",
#     "mapping"        = "~/Onedrive/700-in-progress-projects/706 - Arrabida restored seagrass meadows/2_data/raw_data/edna-sequencing/012522RSeuk1391F-mapping.txt"
# )

files <- data.frame(
    "path" = list.files(
        path    = opts$`data-directory`,
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

# Check if there's 2 files for every sample 
# This check is rudimentary and doesn't see if every sample name has a file
if(length(files$forward[files$forward]) != length(files$forward[files$forward])){
    stop("The number of reverse and forward files does not match.")
}

# Prepare table with primer used per sample -------------------------------
mapping <- read.table(
    file = opts$mapping,
    row.names = NULL,
    sep = "\t", 
    dec = "."
)
names(mapping) <- c("sample", "barcode", "forward_primer", "barcode_name",
                    "reverse_primer", "projectname", "description")

# Primers for forward samples
primers_fwd <- mapping[, c("sample", "forward_primer")]
primers_fwd$forward <- TRUE
names(primers_fwd) <- c("sample", "primer", "forward")
# Primers for reverse samples
primers_rev <- mapping[, c("sample", "reverse_primer")]
primers_rev$forward <- FALSE
names(primers_rev) <- c("sample", "primer", "forward")

primers <- rbind(
    primers_fwd,
    primers_rev
)


# Remove primers ----------------------------------------------------------

# Create final table with all information about each file
files <- merge(
    files,
    primers,
    by = c("sample", "forward")
)

# Iterate trough file table (in parallel)
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
    # Clean memory after the end of each loop
    invisible(gc(verbose = FALSE))
}



