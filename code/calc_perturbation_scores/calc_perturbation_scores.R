suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(glue)
})

args = commandArgs(trailingOnly = TRUE)
N_ARGS = length(args)

# print(args)

input = args[1]
output = args[2]
group = paste(args[3:N_ARGS], collapse = " ")


cat(glue::glue("input = {input}\n"))
cat(glue::glue("output = {output}\n"))
cat(glue::glue("group = {group}\n"))

stopifnot(file.exists(input))

# df = readRDS(input) %>%
# df = qs::qread(input) %>%
df = readr::read_tsv(input) %>%
    dplyr::rename(gene_group = cluster) %>%
    dplyr::filter(gene_group == group) %>%
    # dplyr::select(gene_name, n) %>%
    # dplyr::filter(n >= quantile(n, probs = .50)) %>%
    # dplyr::filter(n >= 5) %>%
    dplyr::mutate(score = 1L) %>%
    dplyr::select(gene_name, score)

stopifnot(nrow(df) >= 10)

readr::write_tsv(df, output, col_names = FALSE)
