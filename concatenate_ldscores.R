suppressPackageStartupMessages({
    library(purrr)
    library(dplyr)
    # library(vroom)
    library(readr)
    library(tibble)
    library(stringr)
    library(logger)
    library(glue)
    library(arrow)
})

args = commandArgs(trailingOnly = TRUE)

log_info(glue("Reading in {length(args)} results files"))

files = tibble(
    file = args
) 

print(files)
# saveRDS(files, "files.Rds")
# saveRDS(args, "args.Rds")
# files = readRDS("files.Rds")

files = files %>%
    mutate(
        splits = map(file, ~str_split(.x, "/")[[1]]),
        split_len = map_int(splits, length),
        BED = map2_chr(splits, split_len, ~.x[[.y - 1]]),
        S2G = map_chr(file, ~str_split(basename(.x), "\\.")[[1]][[1]]),
        PHENO = map_chr(file, ~str_split(basename(.x), "\\.")[[1]][[2]])
    ) 

print(head(files))

enrichments = pmap_dfr(
    list(w = files$file, x = files$BED, y = files$S2G, z = files$PHENO),
    function(w, x, y, z) {
        log_info(glue("reading in {w}"))
        readr::read_tsv(w, col_types = cols()) %>%
            dplyr::mutate(
                gene_set = x,
                s2g = y,
                pheno = z
            )
    }
)

log_info(glue("Done reading in {length(args)} results files"))

readr::write_tsv(enrichments, "concatenated_ldscore_results_modules_2023_05_18.tsv")
arrow::write_feather(enrichments, "concatenated_ldscore_results_modules_2023_05_18.feather")
