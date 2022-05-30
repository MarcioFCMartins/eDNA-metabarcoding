# Assign taxonomy to sequences determined via dada2 pipeline
# Author: Márcio Martins (marciomartinsred@gmail.com)
# Date: 2022-05-26

# Assign taxonomy script starts here --------------------------------------
###########################################################################



# TEST, MOVE LATER --------------------------------------------------------

tax <- assignTaxonomy(
    seqtab,
    "./Onedrive/700-in-progress-projects/706 - Arrabida restored seagrass meadows/2_data/edna-sequencing-processing/silva_132.18s.99_rep_set.dada2.fa.gz",
    multithread = TRUE,
    verbose = TRUE)

tax <- data.frame(tax)
tax$sequence <- rownames(tax)

seqtab_transposed <- data.frame(t(seqtab))
seqtab_transposed$sequence <- rownames(seqtab_transposed)


final <- left_join(seqtab_transposed, tax, by = "sequence")