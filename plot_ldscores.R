suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(glue)
})

ldscore = vroom::vroom("concatenated_ldscore_results.tsv") 

filter_l2 = function(x) x %>% dplyr::filter(Category == "L2_0")
filter_perturb = function(x) x %>% dplyr::filter(stringr::str_detect(gene_set, "perturbation"))

relabel_pheno = function(x) {
    x %>%
        dplyr::mutate(
            pheno_label = case_when(
                pheno == "PASS_Lupus" ~ "Lupus", 
                pheno == "PASS_Rheumatoid_Arthritis" ~ "Rheumatoid Arthritis", 
                pheno == "PASS_Type_1_Diabetes" ~ "Type 1 Diabetes", 
                pheno == "PASS_Ulcerative_Colitis" ~ "Ulcerative Colitis", 
                pheno == "UKB_460K_disease_AID_SURE" ~ "Autoimmune disease", 
                pheno == "UKB_460K_disease_PSORIASIS" ~ "Psoriasis",
                TRUE ~ pheno
            )
        )
}

relabel_gene_set = function(x) {
    x %>%
        dplyr::mutate(
            gene_set_label = case_when(
                gene_set == "perturbation_indegree_control" ~ "Control TFs",
                gene_set == "perturbation_indegree_IEI" ~ "IEI Targets",
                gene_set == "perturbation_indegree_IL2RA" ~ "IL2RA Regulators"
            )
        )
}

relabel_S2G = function(x) {
    x %>%
        dplyr::mutate(
            s2g_label = case_when(
                s2g == "100kb" ~ "100kb around gene body",
                s2g == "5kb" ~ "5kb around gene body",
                s2g == "ABC" ~ "Activity-By-Contact",
                s2g == "PCHiC" ~ "SNPs with promoter-capture\nHiC connections",
                s2g == "Roadmap_Enhancer" ~ "Roadmap enhancer-gene links",
                s2g == "Yoshida" ~ "SNPs in ATAC-peaks\nassociated with expression",
                s2g == "FinemapBloodeQTL" ~ "Finemapped whole-blood\nGTEx eQTLs"
            )
        )
}

plot = ldscore %>%
        filter_l2 %>%
        relabel_gene_set %>%
        group_by(gene_set) %>%
        mutate(
            mean_enrichment = mean(Enrichment),
            gene_set_label = forcats::fct_reorder(factor(gene_set_label), mean_enrichment, .desc = TRUE),
        ) %>%
        relabel_pheno %>%
        relabel_S2G %>%
        filter_perturb %>%
        ggplot(data = ., aes(x = gene_set_label, y = Enrichment)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(aes(x = gene_set_label, y = Enrichment, colour = s2g_label, shape = Enrichment_p < .05), alpha = .9) +
            cowplot::theme_cowplot(font_size = 12) +
            cowplot::panel_border() +
            labs(colour = "SNP to gene strategy", y = "h2 Enrichment", shape = "Enrichment pvalue < 0.05") +
            facet_wrap(~pheno_label) +
            theme(
                axis.title.x = element_blank(), 
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
            ) +
            guides(colour = guide_legend(keyheight = 2))

ggsave("ldscore_h2.pdf", plot, width = 8, height = 6, units = "in")

plot = ldscore %>%
        filter_l2 %>%
        relabel_gene_set %>%
        group_by(gene_set) %>%
        mutate(
            mean_enrichment = mean(Enrichment),
            gene_set_label = forcats::fct_reorder(factor(gene_set_label), mean_enrichment, .desc = TRUE)
        ) %>%
        relabel_pheno %>%
        relabel_S2G %>%
        filter_perturb %>%
        ggplot(data = ., aes(x = gene_set_label, y = Prop._h2)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(aes(x = gene_set_label, y = Prop._h2, colour = s2g_label, shape = Enrichment_p < .05), alpha = .9) +
            cowplot::theme_cowplot(font_size = 12) +
            cowplot::panel_border() +
            labs(colour = "SNP to gene strategy", y = "h2 Proportion", shape = "Enrichment pvalue < 0.05") +
            facet_wrap(~pheno_label) +
            theme(
                axis.title.x = element_blank(), 
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
            ) +
            guides(colour = guide_legend(keyheight = 2))

ggsave("ldscore_h2_prop.pdf", plot, width = 8, height = 6, units = "in")
