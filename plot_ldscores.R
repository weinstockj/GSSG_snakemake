suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(glue)
})

DATE = "2023_05_18"
# ldscore = vroom::vroom("concatenated_ldscore_results.tsv") 
ldscore = vroom::vroom(glue::glue("concatenated_ldscore_results_modules_{DATE}.tsv"))
S2G = "ABC"
# S2G = "100kb"

filter_l2 = function(x) x %>% dplyr::filter(Category == "L2_0")
filter_perturb = function(x) x %>% dplyr::filter(stringr::str_detect(gene_set, "perturbation"))

# relabel_pheno = function(x) {
#     x %>%
#         dplyr::mutate(
#             pheno_label = case_when(
#                 pheno == "PASS_Lupus" ~ "Lupus", 
#                 pheno == "PASS_Rheumatoid_Arthritis" ~ "Rheumatoid Arthritis", 
#                 pheno == "PASS_Type_1_Diabetes" ~ "Type 1 Diabetes", 
#                 pheno == "PASS_Ulcerative_Colitis" ~ "Ulcerative Colitis", 
#                 pheno == "UKB_460K_disease_AID_SURE" ~ "Autoimmune disease", 
#                 pheno == "UKB_460K_disease_PSORIASIS" ~ "Psoriasis",
#                 TRUE ~ pheno
#             )
#         )
# }
relabel_pheno = function(x) {
    x %>%
        dplyr::mutate(
            pheno_label = case_when(
               pheno == "finngen_R8_AUTOIMMUNE" ~ "All Autoimmune",
               pheno == "finngen_R8_CHRONNAS" ~ "Crohn's Disease",
               pheno == "finngen_R8_G6_MS" ~ "Multiple Sclerosis",
               pheno == "finngen_R8_J10_ASTHMA_EXMORE" ~ "Asthma",
               pheno == "finngen_R8_K11_UC_STRICT2" ~ "Ulcerative Colitis",
               pheno == "finngen_R8_L12_DERMATITISECZEMA" ~ "Dermatitis Eczema",
               pheno == "finngen_R8_L12_LUPUS" ~ "Lupus",
               pheno == "finngen_R8_L12_PSORIASIS" ~ "Psoriasis",
               pheno == "finngen_R8_RHEUMA_NOS" ~ "Rheumatoid Arthritis",
               pheno == "finngen_R8_T1D_WIDE" ~ "Type 1 Diabetes",
               pheno == "RA" ~ "Rheumatoid Arthritis",
               pheno == "LUPUS" ~ "Lupus",
               TRUE ~ pheno
            )
        )
}

# relabel_gene_set = function(x) {
#     x %>%
#         dplyr::mutate(
#             gene_set_label = case_when(
#                 gene_set == "perturbation_indegree_control" ~ "Control TFs",
#                 gene_set == "perturbation_indegree_IEI" ~ "IEI Targets",
#                 gene_set == "perturbation_indegree_IL2RA" ~ "IL2RA Regulators"
#             )
#         )
# }

relabel_gene_set = function(x) {
    x %>%
        dplyr::mutate(
            gene_set_label = stringr::str_remove(gene_set, "perturbation_indegree_")
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

# plot = ldscore %>%
#         filter_l2 %>%
#         # relabel_gene_set %>%
#         group_by(gene_set) %>%
#         mutate(
#             mean_enrichment = mean(Enrichment),
#             gene_set_label = forcats::fct_reorder(factor(gene_set_label), mean_enrichment, .desc = TRUE),
#         ) %>%
#         relabel_pheno %>%
#         relabel_S2G %>%
#         filter_perturb %>%
#         ggplot(data = ., aes(x = gene_set_label, y = Enrichment)) +
#             geom_boxplot(outlier.shape = NA) +
#             geom_jitter(aes(x = gene_set_label, y = Enrichment, colour = s2g_label, shape = Enrichment_p < .05), alpha = .9) +
#             cowplot::theme_cowplot(font_size = 12) +
#             cowplot::panel_border() +
#             labs(colour = "SNP to gene strategy", y = "h2 Enrichment", shape = "Enrichment pvalue < 0.05") +
#             facet_wrap(~pheno_label) +
#             theme(
#                 axis.title.x = element_blank(), 
#                 axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
#             ) +
#             guides(colour = guide_legend(keyheight = 2))

# ggsave("ldscore_h2.pdf", plot, width = 8, height = 6, units = "in")

# plot = ldscore %>%
#         filter_l2 %>%
#         relabel_gene_set %>%
#         group_by(gene_set) %>%
#         mutate(
#             mean_enrichment = mean(Enrichment),
#             gene_set_label = forcats::fct_reorder(factor(gene_set_label), mean_enrichment, .desc = TRUE)
#         ) %>%
#         relabel_pheno %>%
#         relabel_S2G %>%
#         filter_perturb %>%
#         ggplot(data = ., aes(x = gene_set_label, y = Prop._h2)) +
#             geom_boxplot(outlier.shape = NA) +
#             geom_jitter(aes(x = gene_set_label, y = Prop._h2, colour = s2g_label, shape = Enrichment_p < .05), alpha = .9) +
#             cowplot::theme_cowplot(font_size = 12) +
#             cowplot::panel_border() +
#             labs(colour = "SNP to gene strategy", y = "h2 Proportion", shape = "Enrichment pvalue < 0.05") +
#             facet_wrap(~pheno_label) +
#             theme(
#                 axis.title.x = element_blank(), 
#                 axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
#             ) +
#             guides(colour = guide_legend(keyheight = 2))

# ggsave("ldscore_h2_prop.pdf", plot, width = 8, height = 6, units = "in")

plot = ldscore %>%
        filter_l2 %>%
        dplyr::filter(s2g == .env[["S2G"]]) %>%
        relabel_gene_set %>%
        relabel_pheno %>%
        relabel_S2G %>%
        filter_perturb %>%
        dplyr::mutate(
            conf.low = Coefficient - 1.96 * Coefficient_std_error,
            conf.high = Coefficient + 1.96 * Coefficient_std_error,
            Coefficient_pvalue = pnorm(abs(`Coefficient_z-score`), lower.tail = FALSE) * 2,
            sig = Coefficient_pvalue < .05
        )
        # dplyr::mutate(
        #     conf.low = Enrichment - 1.96 * Enrichment_std_error,
        #     conf.high = Enrichment + 1.96 * Enrichment_std_error,
        #     sig = Enrichment_p < .05
        # )

plot = plot %>%
        # ggplot(data = ., aes(y = pheno_label, x = Enrichment, xmin = conf.low, xmax = conf.high, colour = sig)) +
        ggplot(data = ., aes(y = pheno_label, x = Coefficient, xmin = conf.low, xmax = conf.high, colour = sig)) +
            # facet_grid(rows = vars(gene_set_label)) +
            facet_wrap(~gene_set_label, ncol = 2) +
            geom_pointrange() +
            # geom_vline(xintercept = 1.0, linetype = "dashed", colour = "gray", alpha = .7) +
            geom_vline(xintercept = 0.0, linetype = "dashed", colour = "gray", alpha = .7) +
            cowplot::theme_cowplot(font_size = 12) +
            cowplot::panel_border() +
            scale_color_manual(values = c("TRUE" = "#c23121", "FALSE" = "grey")) +
            # labs(colour = "Enrichment pvalue < 0.05") +
            labs(colour = "Tau pvalue < 0.05", x = "Tau (95% CI)") +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 3), labels =scales::label_scientific(digits = 2)) +
            # coord_cartesian(xlim = c(-1e-8, 3e-8)) +
            # coord_cartesian(xlim = c(-4e-8, 5e-7)) +
            theme(
                axis.title.y = element_blank(),
                axis.text.y  = element_text(size = 10),
                legend.position = "bottom"
            ) +
            guides(colour = guide_legend(keyheight = 1.3))

ggsave(
       glue::glue("ldscore_h2_modules_tau_{S2G}_{DATE}.pdf"), 
       plot,
       width = 6,
       height = 10,
       units = "in")

plot = ldscore %>%
        filter_l2 %>%
        dplyr::filter(s2g == .env[["S2G"]]) %>%
        relabel_gene_set %>%
        relabel_pheno %>%
        relabel_S2G %>%
        filter_perturb %>%
        dplyr::filter(gene_set_label == "4A") %>%
        dplyr::mutate(
            conf.low = Coefficient - 1.96 * Coefficient_std_error,
            conf.high = Coefficient + 1.96 * Coefficient_std_error,
            Coefficient_pvalue = pnorm(abs(`Coefficient_z-score`), lower.tail = FALSE) * 2,
            sig = Coefficient_pvalue < .05
        )

plot = plot %>%
        ggplot(data = ., aes(y = pheno_label, x = Coefficient, xmin = conf.low, xmax = conf.high, colour = sig)) +
            # facet_grid(rows = vars(gene_set_label)) +
            geom_pointrange() +
            # geom_vline(xintercept = 1.0, linetype = "dashed", colour = "gray", alpha = .7) +
            geom_vline(xintercept = 0.0, linetype = "dashed", colour = "gray", alpha = .7) +
            cowplot::theme_cowplot(font_size = 12) +
            scale_color_manual(values = c("TRUE" = "#c23121", "FALSE" = "grey")) +
            # labs(colour = "Enrichment pvalue < 0.05") +
            labs(colour = "Tau pvalue < 0.05", x = "Tau (95% CI)") +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 3), labels =scales::label_scientific(digits = 2)) +
            # coord_cartesian(xlim = c(-1e-8, 3e-8)) +
            # coord_cartesian(xlim = c(-4e-8, 5e-7)) +
            theme(
                axis.title.y = element_blank(),
                axis.text.y  = element_text(size = 10),
                legend.position = "bottom"
            )

ggsave(
       glue::glue("ldscore_h2_modules_4A_tau_{S2G}_{DATE}.pdf"), 
       plot,
       width = 5,
       height = 7,
       units = "in")

plot = ldscore %>%
        filter_l2 %>%
        dplyr::filter(s2g == .env[["S2G"]]) %>%
        relabel_gene_set %>%
        relabel_pheno %>%
        relabel_S2G %>%
        filter_perturb %>%
        dplyr::mutate(
            conf.low = Enrichment - 1.96 * Enrichment_std_error,
            conf.high = Enrichment + 1.96 * Enrichment_std_error,
            sig = Enrichment_p < .05
        )

plot = plot %>%
        # ggplot(data = ., aes(y = pheno_label, x = Enrichment, xmin = conf.low, xmax = conf.high, colour = sig)) +
        ggplot(data = ., aes(y = pheno_label, x = Enrichment, xmin = conf.low, xmax = conf.high, colour = sig)) +
            # facet_grid(rows = vars(gene_set_label)) +
            facet_wrap(~gene_set_label, ncol = 2) +
            geom_pointrange() +
            # geom_vline(xintercept = 1.0, linetype = "dashed", colour = "gray", alpha = .7) +
            geom_vline(xintercept = 0.0, linetype = "dashed", colour = "gray", alpha = .7) +
            cowplot::theme_cowplot(font_size = 12) +
            cowplot::panel_border() +
            scale_color_manual(values = c("TRUE" = "#c23121", "FALSE" = "grey")) +
            # labs(colour = "Enrichment pvalue < 0.05") +
            labs(colour = "h2 Enrichment pvalue < 0.05", x = "Tau (95% CI)") +
            # coord_cartesian(xlim = c(-1e-8, 3e-8)) +
            # coord_cartesian(xlim = c(-4e-8, 5e-7)) +
            theme(
                axis.title.y = element_blank(),
                axis.text.y  = element_text(size = 10),
                legend.position = "bottom"
            ) +
            guides(colour = guide_legend(keyheight = 1.3))

ggsave(
       glue::glue("ldscore_h2_modules_marginal_enrichment_{S2G}_{DATE}.pdf"), 
       plot,
       width = 6,
       height = 10,
       units = "in")
