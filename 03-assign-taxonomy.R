# Take a sequence table and assign taxonomies
# Author: Márcio Martins (marciomartinsred@gmail.com)
# Date: 2022-05-30

# Install and load required packages --------------------------------------
packages <- c("BiocManager", "optparse", "dada2", "dplyr")

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
        c("-d", "--sequence-table"),
        type = "character",
        help = "Path to the file containing your sequence table."
    ),
    make_option(
        c("-r", "--reference-sequences"),
        type = "character",
        help = "Path to the file containing the reference fastq file used for the classifier."
    ),
    make_option(
        c("-o", "--output"),
        type = "character",
        help = "Directory to export your taxonomy table."
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Assign taxonomy
seqtab <- readRDS(
    opts$`sequence-table`
)

tax <- assignTaxonomy(
    seqtab,
    opts$`reference-sequences`,
    multithread = TRUE,
    verbose = TRUE)

tax <- data.frame(tax)
tax$sequence <- rownames(tax)

seqtab_transposed <- data.frame(t(seqtab))
seqtab_transposed$sequence <- rownames(seqtab_transposed)

final <- left_join(seqtab_transposed, tax, by = "sequence")

reference_db <- basename(opts$`reference-sequences`)
reference_db <- sub("\\.dada2.*", "", reference_db)

write.csv(
    final,
    file.path(opts$output, paste0("taxonomy-", reference_db, ".csv"))
)
