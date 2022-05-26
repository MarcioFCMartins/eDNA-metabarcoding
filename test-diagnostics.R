dim(seqtab)
table(nchar(getSequences(seqtab)))

hist(nchar(getSequences(seqtab)))

getN <- function(x) sum(getUniques(x))
track <- cbind(
    filtered_output, 
    sapply(fwd_dada, getN), 
    sapply(rev_dada, getN), 
    sapply(merged, getN), 
    rowSums(seqtab))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- fwd_samples$sample
View(track)

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

ggplot(track_long) +
    geom_line(
        aes(x = name,
            y = log10(value),
            group = sample)
    ) +
    theme_minimal()
