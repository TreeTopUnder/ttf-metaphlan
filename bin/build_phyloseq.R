#!/usr/bin/env Rscript
# Build a phyloseq object from MetaPhlAn 4 species/SGB tables.
# Research use only. Not for clinical or diagnostic use.

parse_args <- function(args) {
    parsed <- list(
        species = NULL,
        sgb = NULL,
        metadata = NULL,
        db = NULL,
        output = "metaphlan_phyloseq.rds"
    )
    i <- 1
    while (i <= length(args)) {
        key <- args[[i]]
        value <- if (i < length(args)) args[[i + 1]] else NA
        if (key == "--species") parsed$species <- value
        else if (key == "--sgb") parsed$sgb <- value
        else if (key == "--metadata") parsed$metadata <- value
        else if (key == "--db") parsed$db <- value
        else if (key == "--output") parsed$output <- value
        else stop("Unknown argument: ", key)
        i <- i + 2
    }
    parsed
}

strip_rank <- function(token) {
    sub("^[kpcofgst]__", "", token)
}

parse_taxonomy <- function(clades) {
    ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB")
    prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__", "t__")
    mat <- matrix(NA_character_, nrow = length(clades), ncol = length(ranks), dimnames = list(clades, ranks))
    for (i in seq_along(clades)) {
        tokens <- unlist(strsplit(clades[[i]], "|", fixed = TRUE))
        for (token in tokens) {
            for (j in seq_along(prefixes)) {
                if (startsWith(token, prefixes[[j]])) {
                    mat[i, j] <- strip_rank(token)
                }
            }
        }
    }
    mat
}

read_abundance <- function(path) {
    if (is.null(path) || !file.exists(path) || file.info(path)$size == 0) {
        return(NULL)
    }
    table <- tryCatch(
        utils::read.delim(path, check.names = FALSE, row.names = 1, comment.char = "#"),
        error = function(e) NULL
    )
    if (is.null(table) || nrow(table) == 0) {
        return(NULL)
    }
    table
}

normalize_tip <- function(x) {
    x <- sub(".*\\|", "", x)
    x <- sub("^[st]__", "", x)
    x
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package phyloseq is required")
}

sgb <- read_abundance(args$sgb)
species <- read_abundance(args$species)
otu_src <- if (!is.null(sgb) && nrow(sgb) > 0) sgb else species
if (is.null(otu_src) || nrow(otu_src) == 0) {
    stop("No species or SGB abundance rows were available to build phyloseq")
}

otu_mat <- as.matrix(otu_src)
storage.mode(otu_mat) <- "numeric"
otu_mat[is.na(otu_mat)] <- 0
otu <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)

tax_mat <- parse_taxonomy(rownames(otu_mat))
tax <- phyloseq::tax_table(tax_mat)

sample_data_obj <- NULL
if (!is.null(args$metadata) && nzchar(args$metadata) && file.exists(args$metadata) && file.info(args$metadata)$size > 0) {
    md <- utils::read.delim(args$metadata, check.names = FALSE, comment.char = "#")
    if (ncol(md) >= 1 && nrow(md) >= 1) {
        id_col <- names(md)[[1]]
        rownames(md) <- as.character(md[[id_col]])
        common <- intersect(rownames(md), colnames(otu_mat))
        if (length(common) > 0) {
            md <- md[common, , drop = FALSE]
            otu <- phyloseq::prune_samples(common, otu)
            sample_data_obj <- phyloseq::sample_data(md)
        }
    }
}

tree_obj <- NULL
if (!is.null(args$db) && dir.exists(args$db)) {
    nwk <- list.files(args$db, pattern = "\\.nwk$", recursive = TRUE, full.names = TRUE)
    if (length(nwk) > 0) {
        tree <- tryCatch(phyloseq::read_tree(nwk[[1]]), error = function(e) NULL)
        if (!is.null(tree)) {
            otu_ids <- phyloseq::taxa_names(otu)
            tip_labels <- tree$tip.label
            mapped <- character(length(otu_ids))
            names(mapped) <- otu_ids
            norm_tips <- normalize_tip(tip_labels)
            names(norm_tips) <- tip_labels
            for (otu_id in otu_ids) {
                if (otu_id %in% tip_labels) {
                    mapped[[otu_id]] <- otu_id
                    next
                }
                hits <- tip_labels[norm_tips == normalize_tip(otu_id)]
                if (length(hits) == 1) {
                    mapped[[otu_id]] <- hits[[1]]
                    next
                }
                hits <- tip_labels[grepl(normalize_tip(otu_id), tip_labels, fixed = TRUE)]
                if (length(hits) == 1) {
                    mapped[[otu_id]] <- hits[[1]]
                }
            }
            keep_otus <- names(mapped)[nzchar(mapped)]
            keep_tips <- unique(mapped[keep_otus])
            if (length(keep_tips) >= 2) {
                drop_tips <- setdiff(tip_labels, keep_tips)
                if (length(drop_tips) > 0 && requireNamespace("ape", quietly = TRUE)) {
                    tree <- ape::drop.tip(tree, drop_tips)
                } else if (length(drop_tips) > 0) {
                    tree <- phyloseq::prune_taxa(keep_tips, phyloseq::phy_tree(tree))
                }
                rename <- mapped[keep_otus]
                # Align OTU names to remaining tree tips when the tree uses SGB-only labels.
                if (!is.null(tree) && !is.null(tree$tip.label)) {
                    new_otu_names <- phyloseq::taxa_names(otu)
                    for (otu_id in keep_otus) {
                        tip <- rename[[otu_id]]
                        if (!is.null(tip) && tip %in% tree$tip.label) {
                            new_otu_names[new_otu_names == otu_id] <- tip
                        }
                    }
                    phyloseq::taxa_names(otu) <- new_otu_names
                    phyloseq::taxa_names(tax) <- phyloseq::taxa_names(otu)
                    common_taxa <- intersect(phyloseq::taxa_names(otu), tree$tip.label)
                    if (length(common_taxa) >= 2) {
                        otu <- phyloseq::prune_taxa(common_taxa, otu)
                        tax <- phyloseq::prune_taxa(common_taxa, tax)
                        tree_obj <- tree
                    }
                }
            }
        }
    }
}

components <- list(otu, tax)
if (!is.null(sample_data_obj)) {
    components <- c(components, list(sample_data_obj))
}
if (!is.null(tree_obj)) {
    components <- c(components, list(phyloseq::phy_tree(tree_obj)))
}

ps <- do.call(phyloseq::phyloseq, components)
dir.create(dirname(args$output), showWarnings = FALSE, recursive = TRUE)
saveRDS(ps, args$output)
message("Wrote ", args$output, " with ", phyloseq::ntaxa(ps), " taxa and ", phyloseq::nsamples(ps), " samples")
