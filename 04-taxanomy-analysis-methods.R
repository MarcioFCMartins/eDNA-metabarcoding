# Contains methods to assist in analysis taxonomy results from DNA metabarcoding
# Source this file in your analysis and follow along examples
# Author: Márcio Martins - September 2021




# Sankey plots ------------------------------------------------------------

# Function to order taxonomic factor levels to improve readability in ggsankey
# Levels are ordered based on observed number of sequences, but ensuring that
# the order of the parent taxonomic groups is maintained. I.e. levels are first
# ordered based on the "parent" taxonomic level and only then by sequence count.
# This minimizes the crossing of flow lines in the diagram and improves
# readability.
# @param df Name of data.frame with taxonomic information
# @param cols Character vector with names of columns containing taxonomic levels
# @param value Character name of the column holding the counts
fct_hierarchical_sort <- function(df, cols, value) {
    factor_levels <- c()
    previous_levels <- NULL
    current_levels <- NULL
    
    for (i in 1:length(cols)) {
        if (i == 1) {
            col1 <- cols[i]
            
            order <- df |>
                group_by(.data[[col1]]) |>
                summarise(total = sum(.data[[value]])) |>
                arrange(total) |>
                as.data.frame()
            
            current_levels <- order[, 1]
            previous_levels <- data.frame(
                "groups" = current_levels,
                "p_order" = c(1:length(current_levels))
            )
            
            factor_levels <- c(factor_levels, current_levels)
        } else {
            col1 <- cols[i - 1]
            col2 <- cols[i]
            
            order <- df |>
                group_by(.data[[col1]], .data[[col2]]) |>
                summarise(total = sum(.data[[value]]))
            
            order <- merge(
                order,
                previous_levels,
                by.x = col1,
                by.y = "groups"
            ) |>
                arrange(p_order, total) |>
                as.data.frame()
            
            current_levels <- order[, 2]
            previous_levels <- data.frame(
                "groups" = current_levels,
                "p_order" = c(1:length(current_levels))
            )
            
            factor_levels <- c(factor_levels, current_levels)
        }
    }
    
    return(unique(factor_levels))
    
}

# Tax is the taxonomy csv exported by 03-assign-taxonomy.R

# library(ggsankey)
# library(dplyr)
# library(tidyr)
# library(forcats)

# tax_counts <- tax |>
#     rowwise() |>
#     mutate(total = sum(c_across(starts_with("ar") | starts_with("ss")))) |>
#     mutate(
#         across(c(kingdom, phylum, class, order, family, genus),
#                .fns = \(x) ifelse(is.na(x), "Unknown", x)
#         )
#     ) |>
#     ungroup()
# 
# # Sort ALL taxonomic groups by number of observations and order of
# # previous taxonomic level's groups
# sankey_factor_levels <- fct_hierarchical_sort(
#     tax_counts,
#     c("kingdom", "phylum", "class", "order", "family", "genus"),
#     "total"
# )
# 
# # Format in long mode with in/out nodes for sankey
# nodes <- tax_counts |>
#     select(sequence:total) |>
#     ggsankey::make_long(kingdom, phylum, class, order, family, genus,
#                         value = total
#     )
# 
# 
# nodes <- mutate(
#     nodes,
#     node = factor(node, levels = sankey_factor_levels),
#     next_node = factor(next_node, levels = sankey_factor_levels)
# ) |>
#     mutate(
#         node = fct_relevel(node, c("Unknown"), after = 0)
#     )
# 
# 
# ggplot(nodes, aes(
#     x = x,
#     next_x = next_x,
#     node = node,
#     next_node = next_node,
#     fill = node,
#     value = value,
#     label = node
# )) +
#     geom_sankey(
#         color = "#4f4f4f",
#         flow.alpha = .5
#     ) +
#     geom_sankey_label(
#         fill = "#FFFFFF70",
#         color = "#00000070",
#         label.size = NA,
#         size = 3
#     ) +
#     theme_sankey(base_size = 16) +
#     theme(legend.position = "none")



# Re-level taxonomic groups into "interest groups" -------------------

# Function to assign a "common" name to all observations belonging
# to arbitrary taxonomic levels. You have to give it a table with
# taxonomic data and one with the new names for taxonomic groups
# it will then see if any of the taxonomic levels given matches
# one of your interest group and return a new vector with the interest names
# @param taxa Data.frame with taxonomic groups
# @param group_names Data.frame with 2 columns: taxon and new_name
taxa_to_common_group <- function(taxa, group_names) {
    # Return the highest level taxon for each observation
    # that is found in the interest taxon data.frame
    # Important that the taxa table follows a left to right nesting
    highest_level_matching_taxon <- apply(
        taxa,
        1,
        \(x) ifelse(
            any(x %in% group_names$taxon),
            x[max(which(x %in% group_names$taxon))],
            NA
        )
    )
    
    new_group_indexes <- sapply(
        highest_level_matching_taxon,
        \(x) ifelse(
            !is.na(x),
            which(group_names$taxon %in% x),
            NA
        ),
        USE.NAMES = FALSE
    )
    
    new_groups <- group_names$new_name[new_group_indexes]
    
    return(new_groups)
}

# interest_groups <- data.frame(
#     "taxon" = c(
#         "Zostera",
#         "Posidonia",
#         "Embryophyta",
#         "Phragmoplastophyta",
#         "Diatomea",
#         "Chlorophyta_ph",
#         "Florideophycidae",
#         "Dinoflagellata",
#         "Cryptomonadales",
#         "Phaeophyceae",
#         "Ochrophyta"
#     ),
#     "new_name" = c(
#         "Zostera",
#         "Cymodocea",
#         "Other plants",
#         "Other plants",
#         "Diatoms",
#         "Green algae",
#         "Red algae",
#         "Other microalgae",
#         "Other microalgae",
#         "Brown algae",
#         "Brown algae"
#     )
# )
# 
# tax$interest_group <- taxa_to_common_group(
#     tax_counts_long[, 2:7],
#     interest_groups
# )


# Perform NMDS to view clustering of smaples ------------------------------
# library(vegan)

# # https://jkzorz.github.io/2020/04/04/NMDS-extras.html
# nmds_loadings <- envfit(tax_mds, tax_vegan_rel)
# 
# 
# nmds_loadings_ggplot <- data.frame(
#     nmds_loadings$vectors$arrows,
#     "r" = nmds_loadings$vectors$r,
#     "p" = nmds_loadings$vectors$pvals
# ) |>
#     rownames_to_column("OTU") |>
#     filter(r > 0.2) |>
#     rowwise() |>
#     mutate(
#         theta = atan2(NMDS2, NMDS1),
#         x = r * cos(theta),
#         y = r * sin(theta)
#     ) |>
#     mutate(
#         OTU = str_replace(OTU, "(\\d+)", ""),
#         OTU = str_replace_all(OTU, "\\.+", " "),
#         OTU = str_trim(OTU),
#         OTU = factor(OTU, levels = unique(OTU))
#     )
# 
# nmds_points <- tax_mds$points |>
#     as_tibble(rownames = "samples") |>
#     mutate(core_id = toupper(substr(samples, 1, 4))) |>
#     left_join(cores, by = "core_id")
# 
# 
# group_polygons <- nmds_points |>
#     group_by(sites_notes) |>
#     # chull will say which rows comprise the convex polygon
#     # and slice will extract those rows
#     slice(chull(MDS1, MDS2))
# 
# site_labels <- c(
#     "Ponta do Adoxe" = "Nearby (PA)",
#     "Ria Formosa"    = "Donor (RF)",
#     "Arrábida 10yr"  = "Restored-10yr (AR10)",
#     "Arrábida 3yr"   = "Restored-3yr"
# )
# 
# ggplot() +
#     geom_point(
#         data = nmds_points,
#         aes(x = MDS1, y = MDS2, fill = sites_notes, shape = sites_notes),
#     ) +
#     geom_polygon(
#         data = group_polygons,
#         aes(x = MDS1, y = MDS2, fill = sites_notes),
#         alpha = .15
#     ) +
#     geom_segment(
#         data = nmds_loadings_ggplot,
#         aes(x = 0, y = 0, xend = x, yend = y, color = OTU),
#         key_glyph = "rect"
#     ) +
#     geom_text_repel(
#         data = nmds_loadings_ggplot,
#         aes(x = x, y = y, label = OTU, color = OTU),
#         max.overlaps = 30
#     ) +
#     scale_shape_manual(
#         labels = site_labels,
#         values = c(21:24),
#         name = "Samples meadows"
#     ) +
#     scale_fill_discrete(
#         labels = site_labels
#     ) +
#     scale_color_manual(
#         values = group_colors
#     ) +
#     labs(
#         color = "Biological group",
#         fill = "Samples meadows"
#     ) 

