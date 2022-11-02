import os 
OAK = os.environ['OAK_DIR']

log_dir = "logs"
OAK_resource_dir = os.path.join(OAK, "resources", "GSSG")
local_resource_dir = "processed_data"

local_geneset_dir = "genesets"

OAK_bed_dir = os.path.join(OAK, "GSSG_bed")
OAK_annot_dir = os.path.join(OAK, "GSSG_annot")
OAK_LDSC_dir = os.path.join(OAK, "GSSG_ldsc")
OAK_sumstats_dir = os.path.join(OAK, "resources", "LDSC", "sumstats")
local_bed_dir = "bed"

#PRICE_LAB_URL="https://alkesgroup.broadinstitute.org/LDSCORE/Dey_Enhancer_MasterReg/processed_data"
PRICE_LAB_URL="https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/Dey_Enhancer_MasterReg/processed_data/*"
PRICE_LAB_URL="gs://broad-alkesgroup-public/LDSCORE/Dey_Enhancer_MasterReg/processed_data"

downstream_indegree_path = "/oak/stanford/groups/pritch/users/jweinstk/perturbation_data/rnaseq_pipeline/scripts/_targets/objects/downstream_indegree_by_group"

gene_group_prefix = {
    "control" : "Control",
    "IEI"     : "IEI Target",
    "IL2RA"   : "IL2RA Regulators"
}

gene_groups = list(gene_group_prefix.keys())

print("gene groups are: ", gene_groups)

R_VERSION = "4.1.2"
GCC_VERSION = "12.1.0"
BEDTOOLS = "~/downloads/bedtools"
BIMFILE_PATH = OAK_resource_dir
CONDA_INIT = os.path.join(OAK, "miniconda3/etc/profile.d/conda.sh")
BIM_DIR = os.path.join(OAK, "resources", "LDSC", "1KG_EUR")
WEIGHTS_DIR="/oak/stanford/groups/pritch/users/jweinstk/resources/LDSC/weights_hm3_no_hla"
FREQ_DIR="/oak/stanford/groups/pritch/users/jweinstk/resources/LDSC/1000G_Phase3_frq"
hapmap = "/oak/stanford/groups/pritch/users/jweinstk/resources/LDSC/hapmap3_snps"
BASELINE="/oak/stanford/groups/pritch/users/jweinstk/resources/LDSC/baselinev2_1"

LDSC_SCRIPT_DIR="~/downloads/ldsc"
LDSC_RESOURCE_DIR=os.path.join(OAK, "resources", "LDSC")

CHRs = range(1, 23)

PHENOs = (
    "PASS_Lupus",
    "PASS_Rheumatoid_Arthritis",
    "PASS_Type_1_Diabetes",
    "PASS_Ulcerative_Colitis",
    "UKB_460K_disease_AID_SURE",
    "UKB_460K_disease_PSORIASIS"
)

BEDs = (
    "ABC9_G_top10_0_015",
    "ABC9_I_top10_0_015",
    "All_genes",
    "EDS_Binary_top10",
    "eQTL_CTS_prob",
    "Expecto_MVP",
    "HOMOD_prob",
    "master_regulator",
    "PCHiC_binary",
    # "perturbation_indegree",
    "perturbation_indegree_control",
    "perturbation_indegree_IEI",
    "perturbation_indegree_IL2RA",
    "pLI_genes2",
    "PPI_All",
    "PPI_control",
    "PPI_Enhancer",
    "PPI_Master",
    "RegNet_Enhancer",
    "SEG_GTEx_top10",
    "TF_genes_curated",
    "Trans_Reg_genes"
)

S2Gs = (
    "5kb",
    "FinemapBloodeQTL",
    "100kb",
    "Yoshida",
    "PCHiC",
    "Roadmap_Enhancer",
    "ABC"
)

if not os.path.isdir(OAK_resource_dir):
    IOError(f"{OAK_resource_dir} does not exist")

if not os.path.isdir(local_resource_dir):
    IOError(f"{local_resource_dir} does not exist")

if not os.path.isdir(local_bed_dir):
    IOError(f"{local_bed_dir} does not exist")

if not os.path.isdir(hapmap):
    IOError(f"{hapmap} does not exist")

if not os.path.isdir(OAK_sumstats_dir):
    IOError(f"{OAK_sumstats_dir} does not exist")

rule all:
    input:
        os.path.join(OAK_resource_dir, "resources.download.DONE"),
        os.path.join(local_resource_dir, "resources.download.DONE"),
        expand(os.path.join(local_geneset_dir, "perturbation_indegree_{group}.txt"), group = gene_groups),
        expand(os.path.join(local_bed_dir, "{BED}.DONE"), BED = BEDs),
        #expand(os.path.join(local_bed_dir, "{BED}.DONE.ln"), BED = BEDs),
        expand(os.path.join(local_bed_dir, "{BED}.{S2G}.DONE.merged.sorted"), BED = BEDs, S2G = S2Gs),
        expand(os.path.join(OAK_bed_dir, "{BED}", "sorted.merged.{S2G}.bed"), BED = BEDs, S2G = S2Gs),
        expand(os.path.join(OAK_annot_dir, "annot.{S2G}.{BED}.DONE"), BED = BEDs, S2G = S2Gs),
        expand(os.path.join(OAK_annot_dir, "{BED}", "sorted.merged.{S2G}.{CHR}.annot.gz"), BED = BEDs, S2G = S2Gs, CHR = CHRs),
        expand(os.path.join(OAK_LDSC_dir, "{BED}", "{S2G}.{CHR}.l2.ldscore.gz"), BED = BEDs, S2G = S2Gs, CHR = CHRs),
        expand(os.path.join(OAK_LDSC_dir, "{BED}", "{S2G}.{PHENO}.results"), BED = BEDs, S2G = S2Gs, PHENO = PHENOs),
        "concatenated_ldscore_results.tsv"


rule download_resources:
    output:
        OAK = os.path.join(OAK_resource_dir, "resources.download.DONE"),
        local = os.path.join(local_resource_dir, "resources.download.DONE")
    log:
        stderr = os.path.join(log_dir, "stderr.log"),
        stdout = os.path.join(log_dir, "stdout.log")
    shell:
        """
        ml load google-cloud-sdk/338.0.0
        #wget -v -P {OAK_resource_dir} {PRICE_LAB_URL}
        gsutil cp {PRICE_LAB_URL}/* {OAK_resource_dir}
        ln -s {OAK_resource_dir}/* {local_resource_dir}
        touch {output.OAK}
        touch {output.local}
        """
rule create_perturbation_geneset:
    input:
        downstream_indegree_path
    output:
        os.path.join(local_geneset_dir, "perturbation_indegree_{group}.txt")
    params:
        full_group_name = lambda wildcards: gene_group_prefix[wildcards.group]
    shell:
        """
        ml load R/{R_VERSION}
        ml load gcc/{GCC_VERSION}
        Rscript code/calc_perturbation_scores/calc_perturbation_scores.R {input} {output} {params.full_group_name}
        """

rule create_bedgraphs:
    input:
        #os.path.join(local_geneset_dir, "perturbation_indegree.txt")
        expand(os.path.join(local_geneset_dir, "perturbation_indegree_{group}.txt"), group = gene_groups)
    output:
        os.path.join(local_bed_dir, "{BED}.DONE")
    log:
        stderr = os.path.abspath(os.path.join(log_dir, "create_bedgraphs.{BED}.stderr.log")),
        stdout = os.path.abspath(os.path.join(log_dir, "create_bedgraphs.{BED}.stdout.log"))
    threads:
        1
    shell:
        """
        ml load R/{R_VERSION}
        ml load gcc/{GCC_VERSION}
        orig_dir=$(pwd)
        cd code/GeneSet_toS2G
        Rscript geneset_to_bed.R ../../{local_geneset_dir} {OAK_bed_dir} {wildcards.BED} 1>{log.stdout} 2>{log.stderr}
        cd $orig_dir
        touch {output}
        """

rule sort_merge:
    input:
        rules.create_bedgraphs.output
    output:
        touch = os.path.join(local_bed_dir, "{BED}.{S2G}.DONE.merged.sorted"),
        bed = os.path.join(OAK_bed_dir, "{BED}", "sorted.merged.{S2G}.bed")
    params:
        dir = os.path.join(OAK_bed_dir, "{BED}")
    shell:
        """
        for name in {params.dir}/{wildcards.S2G}.bed
        do
            echo "name is $name"
            base=$(basename $name)
            echo "now writing to {params.dir}/sorted.merged.$base"
            #{BEDTOOLS} sort -i $name | {BEDTOOLS} merge -i stdin > {params.dir}/sorted.merged.$base
            {BEDTOOLS} sort -i $name | {BEDTOOLS} merge -i stdin > {output.bed}
        done

        touch {output.touch}
        """
        
        
rule soft_link_bedgraphs:
    input:
        rules.create_bedgraphs.output,
        rules.sort_merge.output.bed
    output:
        os.path.join(local_bed_dir, "{BED}.DONE.ln")
    threads:
        1
    shell:
        """
        mkdir -p {local_bed_dir}/{wildcards.BED}
        ln -sf {OAK_bed_dir}/{wildcards.BED}/* {local_bed_dir}/{wildcards.BED}
        touch {output}
        """
    
rule annot:
    input:
        rules.sort_merge.output.bed
    output:
        touch = os.path.join(OAK_annot_dir, "annot.{S2G}.{BED}.DONE"),
        annot = expand(os.path.join(OAK_annot_dir, "{{BED}}", "sorted.merged.{{S2G}}.{CHR}.annot.gz"), CHR = CHRs)
    params:
        bed_dir = os.path.join(OAK_bed_dir, "{BED}"),
        annot_dir = os.path.join(OAK_annot_dir, "{BED}"),
        base = "sorted.merged.{S2G}"
    shell:
        """
        source {CONDA_INIT}
        conda activate GSSG_annot
        python code/GeneSet_toS2G/bedgraph_to_annot.py --bedname {params.base} \
            --bedfile_path {params.bed_dir} \
            --bimfile_path {BIMFILE_PATH} \
            --annot_path {params.annot_dir}

        touch {output.touch}
        """

rule create_ldscore:
    input:
        annot = os.path.join(OAK_annot_dir, "{BED}", "sorted.merged.{S2G}.{CHR}.annot.gz")
    output:
        os.path.join(OAK_LDSC_dir, "{BED}", "{S2G}.{CHR}.l2.ldscore.gz")
    params:
        BFILE_PREFIX = os.path.join(BIM_DIR, "1000G.EUR.QC.{CHR}"),
        out_prefix = os.path.join(OAK_LDSC_dir, "{BED}", "{S2G}.{CHR}"),
        LDSC_SCRIPT_DIR = LDSC_SCRIPT_DIR,
        hapmap = hapmap
    shell:
        """
        source {CONDA_INIT}
        conda activate ldsc
        mkdir -p $(dirname {output})
        python {params.LDSC_SCRIPT_DIR}/ldsc.py --bfile {params.BFILE_PREFIX} \
            --l2 \
            --ld-wind-cm 1 \
            --yes-really \
            --annot {input.annot} \
            --print-snps {params.hapmap}/hm.{wildcards.CHR}.snp \
            --out {params.out_prefix}
        """

rule ldscore_reg:
    input:
        annot = expand(os.path.join(OAK_annot_dir, "{{BED}}", "sorted.merged.{S2G}.{CHR}.annot.gz"), CHR = CHRs, S2G = S2Gs),
        ldscore = expand(os.path.join(OAK_LDSC_dir, "{{BED}}", "{S2G}.{CHR}.l2.ldscore.gz"), CHR = CHRs, S2G = S2Gs),
        PHENO = os.path.join(OAK_sumstats_dir, "{PHENO}.sumstats") 
    output:
        os.path.join(OAK_LDSC_dir, "{BED}", "{S2G}.{PHENO}.results")
    params:
        out_prefix = os.path.join(OAK_LDSC_dir, "{BED}", "{S2G}.{PHENO}"),
        weight_prefix = os.path.join(WEIGHTS_DIR, "weights."),
        frq_prefix = os.path.join(FREQ_DIR, "1000G.EUR.QC."),
        annot_prefix = os.path.join(OAK_LDSC_dir, "{BED}", "{S2G}."),
        LDSC_SCRIPT_DIR = LDSC_SCRIPT_DIR,
        baseline_prefix = os.path.join(BASELINE, "baselineLD.")
    shell:
        """
        source {CONDA_INIT}
        conda activate ldsc
        mkdir -p $(dirname {output})
        annot_dir=$(dirname {input.annot})

        for chr in {{1..22}}
            do
            orig={OAK_annot_dir}/{wildcards.BED}/sorted.merged.{wildcards.S2G}.$chr.annot.gz 
            new={OAK_LDSC_dir}/{wildcards.BED}/{wildcards.S2G}.$chr.annot.gz 
            ln -sf $orig $new
        done
        
        python {params.LDSC_SCRIPT_DIR}/ldsc.py --h2 {input.PHENO} \
            --ref-ld-chr {params.annot_prefix},{params.baseline_prefix} \
            --frqfile-chr {params.frq_prefix} \
            --w-ld-chr {params.weight_prefix} \
            --overlap-annot \
            --print-coefficients \
            --print-delete-vals \
            --out {params.out_prefix}

        """

rule concat_ldscore:
    input:
        expand(os.path.join(OAK_LDSC_dir, "{BED}", "{S2G}.{PHENO}.results"), BED = BEDs, S2G = S2Gs, PHENO = PHENOs)
    output:
        "concatenated_ldscore_results.tsv"
    shell:
        """
        ml load R/{R_VERSION}
        ml load gcc/{GCC_VERSION}
        Rscript concatenate_ldscores.R {input}
        """
