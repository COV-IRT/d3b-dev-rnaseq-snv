# RNAseq SNV Calling Workflow (beta)
This repo contains snv calling methods using Broad GATK best practices, Strelka2, and VarDict
It also contains a very basic VEP annotation tool

<p align="center">
  <img alt="Logo for The Center for Data Driven Discovery" src="https://raw.githubusercontent.com/d3b-center/handbook/master/website/static/img/chop_logo.svg?sanitize=true" width="400px" />
</p>

For all workflows, input bams should be indexed beforehand.  This tool is provided in `tools/samtools_index.cwl`

## GATK4 v4.1.7.0
The overall [workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) picks up from post-STAR alignment, starting at mark duplicates.
For the most part, tool parameters follow defaults from the GATK Best Practices [WDL](https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl), written in cwl with added optimization for use on the Cavatica platform.
A mild warning, sambamba is used in this workflow to mark dulicates for speed and efficiency instead of picard.
Behavior should be the same except for markig optical duplicates.
`workflows/d3b_gatk_rnaseq_snv_wf.cwl` is the wrapper cwl used to run all tools for GATK4.
Run time (n=782) COV-IRT dataset,  ~3 hours, cost on cavatica ~$1.05 per sample

### Inputs
```yaml
inputs:
  output_basename: string
  scatter_ct: {type: int?, doc: "Number of interval lists to split into", default: 50}
  STAR_sorted_genomic_bam: {type: File, doc: "STAR sorted alignment bam"}
  reference_fasta: {type: File, secondaryFiles: ['^.dict', '.fai'], doc: "Reference genome used"}
  reference_dict: File
  call_bed_file: {type: File, doc: "BED or GTF intervals to make calls"}
  exome_flag: {type: string?, default: "Y", doc: "Whether to run in exome mode for callers. Should be Y or leave blank as default is Y. Only make N if you are certain"}
  knownsites: {type: 'File[]', doc: "Population vcfs, based on Broad best practices"}
  dbsnp_vcf: {type: File, secondaryFiles: ['.idx']}
  tool_name: {type: string, doc: "description of tool that generated data, i.e. gatk_haplotypecaller"}
  mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
```

### Outputs
```yaml
outputs:
  filtered_hc_vcf: {type: File, outputSource: gatk_filter_vcf/filtered_vcf, doc: "Haplotype called vcf with Broad-recommended FILTER values added"}
  pass_vcf: {type: File, outputSource: gatk_pass_vcf/pass_vcf, doc: "Filtered vcf selected for PASS variants"}
  anaylsis_ready_bam: {type: File, outputSource: gatk_applybqsr/recalibrated_bam, doc: "Duplicate marked, Split N trimmed CIGAR BAM, BQSR recalibratede, ready for RNAseq calling"}
  bqsr_table: {type: File, outputSource: gatk_baserecalibrator/output, doc: "BQSR table"}
```

### Docker Pulls
 - `kfdrc/sambamba:0.7.1`
 - `kfdrc/gatk:4.1.7.0R`

### Workflow Diagram

![WF diagram](misc/d3b_gatk_rnaseq_snv_wf.cwl.svg)

### GATK4 simulated bash calls

 | Step                                         | Type         | Num scatter            | Command                                                                                                                                                                                                         |
 | -------------------------------------------- | ------------ | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
 | bedtools_gtf_to_bed                                         | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | bedtools_gtf_to_bed                                         | run step         | NA            | cp /sbgenomics/Projects/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/gencode.v33.primary_assembly.annotation.gtf.bed .; exit 0;                                                                                                                                                                                                         |
 | bedtools_gtf_to_bed                                         | run step         | NA            | cat /sbgenomics/Projects/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/gencode.v33.primary_assembly.annotation.gtf.bed | grep -vE "^#" | cut -f 1,4,5 | awk '{OFS = "\t";a=$2-1;print $1,a,$3; }' | bedtools sort | bedtools merge > gencode.v33.primary_assembly.annotation.gtf.bed.bed                                                                                                                                                                                                         |
 | gatk_intervallisttools                                         | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | gatk_intervallisttools                                         | run step         | NA            | /gatk BedToIntervalList -I /sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/8b37cb5b-ef56-4f3b-a7c9-87a89b07171f/bedtools_gtf_to_bed/gencode.v33.primary_assembly.annotation.gtf.bed -O gencode.v33.primary_assembly.annotation.gtf.interval_list -SD /sbgenomics/Projects/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/GRCh38.primary_assembly.genome.dict; LIST=gencode.v33.primary_assembly.annotation.gtf.interval_list;BANDS=0;                                                                                                                                                                                                         |
 | gatk_intervallisttools                                         | run step         | NA            | /gatk IntervalListTools --java-options "-Xmx2000m" --SCATTER_COUNT=50 --SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW --UNIQUE=true --SORT=true --BREAK_BANDS_AT_MULTIPLES_OF=$BANDS --INPUT=$LIST --OUTPUT=.;CT=`find . -name 'temp_0*' | wc -l`;seq -f "%04g" $CT | xargs -I N -P 4 /gatk IntervalListToBed --java-options -Xmx100m -I temp_N_of_$CT/scattered.interval_list -O temp_N_of_$CT/scattered.interval_list.N.bed;mv temp_0*/*.bed .;                                                                                                                                                                                                         |
 | preprocess_rnaseq_bam_sambamba_md_sorted                                         | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | preprocess_rnaseq_bam_sambamba_md_sorted                                         | run step         | NA            | mkdir TMP                                                                                                                                                                                                         |
 | preprocess_rnaseq_bam_sambamba_md_sorted                                         | run step         | NA            | sambamba markdup --tmpdir TMP -t 4 /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/NASA_HEC_transferred_data/COVIRT_all_data/02-AlignedData/COVHA-20200403-P2-H04-N.all-reads/COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.bam COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.md.bam                                                                                                                                                                                                         |
 | preprocess_rnaseq_bam_sambamba_md_sorted                                         | run step         | NA            | mv COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.md.bam.bai COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.md.bai                                                                                                                                                                                                         |
 | preprocess_rnaseq_bam_gatk_splitntrim                                         | run step         | NA            | /gatk SplitNCigarReads --java-options "-Xmx30G -XX:+PrintFlagsFinal -Xloggc:gc_log.log -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" --seconds-between-progress-updates 30 -R /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/NASA_HEC_transferred_data/genome_files/Homo_sapiens_ensembl_release100/Homo_sapiens.GRCh38.dna.primary_assembly_mainchrs_and_SARS-CoV-2-NC_045512.2.fa -I /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/preprocess_rnaseq_bam_sambamba_md_sorted/COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.md.bam  -OBI -O COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.md.splitn.bam                                                                                                                                                                                                         |
 | gatk_baserecalibrator                                         | run step         | NA            | /gatk BaseRecalibrator --java-options "-Xmx7500m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log" -R /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/NASA_HEC_transferred_data/genome_files/Homo_sapiens_ensembl_release100/Homo_sapiens.GRCh38.dna.primary_assembly_mainchrs_and_SARS-CoV-2-NC_045512.2.fa -I /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/preprocess_rnaseq_bam_gatk_splitntrim/COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.md.splitn.bam --use-original-qualities -O COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.md.splitn.recal_data.csv  --known-sites /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/GATK4_VCF_REFS/1000G_phase1.snps.high_confidence.hg38.chr_converted_canonical.vcf.gz --known-sites /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/GATK4_VCF_REFS/1000G_omni2.5.hg38.chr_converted_canonical.vcf.gz --known-sites /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/GATK4_VCF_REFS/Mills_and_1000G_gold_standard.indels.hg38.chr_converted_canonical.vcf.gz --known-sites /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/GATK4_VCF_REFS/Homo_sapiens_assembly38.known_indels.chr_converted_canonical.vcf.gz                                                                                                                                                                                                         |
 | gatk_applybqsr                                         | run step         | NA            | /gatk ApplyBQSR --java-options "-Xms3000m -Xmx7500m -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" --create-output-bam-md5 --add-output-sam-program-record -R /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/NASA_HEC_transferred_data/genome_files/Homo_sapiens_ensembl_release100/Homo_sapiens.GRCh38.dna.primary_assembly_mainchrs_and_SARS-CoV-2-NC_045512.2.fa -I /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/preprocess_rnaseq_bam_gatk_splitntrim/COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.md.splitn.bam --use-original-qualities -O COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.md.splitn.aligned.duplicates_marked.recalibrated.bam -bqsr /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/gatk_baserecalibrator/COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.md.splitn.recal_data.csv                                                                                                                                                                                                         |
 | gatk_haplotype_rnaseq                                         | scatter         | 50            | /gatk HaplotypeCaller -R /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/NASA_HEC_transferred_data/genome_files/Homo_sapiens_ensembl_release100/Homo_sapiens.GRCh38.dna.primary_assembly_mainchrs_and_SARS-CoV-2-NC_045512.2.fa -I /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/gatk_applybqsr/COVHA-20200403-P2-H04-N.all-reads_Aligned.sortedByCoord.out.md.splitn.aligned.duplicates_marked.recalibrated.bam --standard-min-confidence-threshold-for-calling 20 -dont-use-soft-clipped-bases -L /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/gatk_intervallisttools/scattered.interval_list.0039.bed -O 4ee96c3b-dbf3-4a97-b492-6699ec38adb5.gatk.hc.called.vcf.gz --dbsnp /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/GATK4_VCF_REFS/Homo_sapiens_assembly38.dbsnp138.chr_converted_canonical.vcf                                                                                                                                                                                                         |
 | merge_hc_vcf                                         | run step         | NA            | /gatk MergeVcfs --java-options "-Xmx2000m" --TMP_DIR=./TMP --CREATE_INDEX=true --SEQUENCE_DICTIONARY=/sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/NASA_HEC_transferred_data/genome_files/Homo_sapiens_ensembl_release100/Homo_sapiens.GRCh38.dna.primary_assembly_mainchrs_and_SARS-CoV-2-NC_045512.2.dict --OUTPUT=4ee96c3b-dbf3-4a97-b492-6699ec38adb5.STAR_GATK4.merged.vcf.gz  -I /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/gatk_haplotype_rnaseq_1_s/4ee96c3b-dbf3-4a97-b492-6699ec38adb5.gatk.hc.called.vcf.gz -I /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/gatk_haplotype_rnaseq_2_s/4ee96c3b-dbf3-4a97-b492-6699ec38adb5.gatk.hc.called.vcf.gz -I /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/gatk_haplotype_rnaseq_3_s/4ee96c3b-dbf3-4a97-b492-6699ec38adb5.gatk.hc.called.vcf.gz -I /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/gatk_haplotype_rnaseq_4_s/4ee96c3b-dbf3-4a97-b492-6699ec38adb5.gatk.hc.called.vcf.gz                                                                                                                                                                                                         |
 | gatk_filter_vcf                                         | run step         | NA            | /gatk VariantFiltration -R /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/NASA_HEC_transferred_data/genome_files/Homo_sapiens_ensembl_release100/Homo_sapiens.GRCh38.dna.primary_assembly_mainchrs_and_SARS-CoV-2-NC_045512.2.fa -V /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/merge_hc_vcf/4ee96c3b-dbf3-4a97-b492-6699ec38adb5.STAR_GATK4.merged.vcf.gz --window 35 --cluster 3 --filter-name "FS" --filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0"  -O 4ee96c3b-dbf3-4a97-b492-6699ec38adb5.gatk.hc.filtered.vcf.gz                                                                                                                                                                                                         |
 | gatk_pass_vcf                                         | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | gatk_pass_vcf                                         | run step         | NA            | /gatk SelectVariants --java-options "-Xmx7500m" -V /sbgenomics/workspaces/7c22397d-ea11-4069-9687-3049119a2d37/tasks/4ee96c3b-dbf3-4a97-b492-6699ec38adb5/gatk_filter_vcf/4ee96c3b-dbf3-4a97-b492-6699ec38adb5.gatk.hc.filtered.vcf.gz -O 4ee96c3b-dbf3-4a97-b492-6699ec38adb5.STAR_GATK4.PASS.vcf.gz --exclude-filtered TRUE                                                                                                                                                                                                         |


## Strelka2 v2.9.10
This [workflow](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#rna-seq) is pretty straight forward, with a `PASS` filter step added to get `PASS` calls.
`workflows/d3b_strelka2_rnaseq_snv_wf.cwl` is the wrapper cwl that runs this workflow.
Run time (n=1) ~50 minutes, cost on cavatica ~$0.40

### Inputs
```yaml
inputs:
  reference: { type: File, secondaryFiles: [.fai] }
  input_rna_bam: {type: File, secondaryFiles: [^.bai]}
  strelka2_bed: {type: File?, secondaryFiles: [.tbi], label: gzipped bed file}
  cores: {type: ['null', int], default: 16, doc: "Num cores to use"}
  ram: {type: ['null', int], default: 30, doc: "Max mem to use in GB"}
  output_basename: string
```

### Outputs
```yaml
  strelka2_prepass_vcf: {type: File, outputSource: strelka2_rnaseq/output_vcf, doc: "Strelka2 SNV calls"}
  strelka2_pass_vcf: {type: File, outputSource: gatk_pass_vcf/pass_vcf, doc: "Strelka2 calls filtered on PASS"}
```

### Docker Pulls
 - `kfdrc/strelka2:2.9.10`
 - `kfdrc/gatk:4.1.1.0`

### Workflow Diagram

![WF diagram](misc/d3b_strelka2_rnaseq_snv_wf.cwl.svg)

### Strelka2 simulated bash calls

 | Step                | Type         | Num scatter            | Command                                                                                                                                                                                                         |
 | ------------------- | ------------ | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
 | strelka2_rnaseq                | run step         | NA            | /strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam /sbgenomics/Projects/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/da63df67-62a4-487b-aa68-d7f139809160.Aligned.out.sorted.bam --reference /sbgenomics/Projects/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/GRCh38.primary_assembly.genome.fa  --rna --runDir ./ && ./runWorkflow.py -m local -j 16 -g 30                                                                                                                                                                                                         |
 | strelka2_rnaseq                | run step         | NA            | mv results/variants/variants.vcf.gz STRELKA2_TEST.strelka2.rnaseq.vcf.gz                                                                                                                                                                                                         |
 | strelka2_rnaseq                | run step         | NA            | mv results/variants/variants.vcf.gz.tbi STRELKA2_TEST.strelka2.rnaseq.vcf.gz.tbi                                                                                                                                                                                                         |
 | gatk_pass_vcf                | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | gatk_pass_vcf                | run step         | NA            | /gatk SelectVariants --java-options "-Xmx7500m" -V /sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/5f77306d-9650-4b82-83e1-1623eb07e211/strelka2_rnaseq/STRELKA2_TEST.strelka2.rnaseq.vcf.gz -O STRELKA2_TEST.strelka2.PASS.vcf.gz --exclude-filtered TRUE                                                                                                                                                                                                         |

## VardictJava v1.7.0
This [workflow](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/rnaseq/variation.py) is based on the Vardict run style of BC Bio.
`workflows/d3b_vardict_rnaseq_snv_wf.cwl` is the wrapper cwl that runs this workflow.
Tweaking `vardict_bp_target` and `vardict_intvl_target_size` maybe be needed to improve run time in high coverage areas, by reducing their values from defaults.
Run time (n=1) ~9.5 hours, cost on cavatica ~$5.50.

### Inputs
```yaml
inputs:
  output_basename: string
  STAR_sorted_genomic_bam: {type: File, doc: "STAR sorted alignment bam", secondaryFiles: ['^.bai']}
  sample_name: string
  reference_fasta: {type: File, secondaryFiles: ['.fai', '^.dict'], doc: "Reference genome used"}
  reference_dict: File
  vardict_min_vaf: {type: ['null', float], doc: "Min variant allele frequency for vardict to consider.  Recommend 0.2", default: 0.2}
  vardict_cpus: {type: ['null', int], default: 4}
  vardict_ram: {type: ['null', int], default: 8, doc: "In GB"}
  vardict_bp_target: {type: ['null', int], doc: "Intended max number of base pairs per file.  Existing intervals large than this will NOT be split into another file. Make this value smaller to break up the work into smaller chunks", default: 60000000}
  vardict_intvl_target_size: {type: ['null', int], doc: "For each file, split each interval into chuck of this size", default: 20000}
  call_bed_file: {type: File, doc: "BED or GTF intervals to make calls"}
  tool_name: {type: string, doc: "description of tool that generated data, i.e. gatk_haplotypecaller"}
  padding: {type: ['null', int], doc: "Padding to add to input intervals, recommened 0 if intervals already padded, 150 if not", default: 150}
  mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
```

### Outputs
```yaml
outputs:
  vardict_prepass_vcf: {type: File, outputSource: sort_merge_vardict_vcf/merged_vcf, doc: "VarDict SNV calls"}
  vardict_pass_vcf: {type: File, outputSource: gatk_pass_vcf/pass_vcf, doc: "VarDict calls filtered on PASS"}
```

### Docker Pulls
- `kfdrc/vardict:1.7.0`
- `kfdrc/gatk:4.1.1.0`
- `kfdrc/python:2.7.13`

### Workflow Diagram

![WF diagram](misc/d3b_vardict_rnaseq_snv_wf.cwl.svg)

### Vardict simulated bash calls

 | Step                              | Type         | Num scatter            | Command                                                                                                                                                                                                         |
 | --------------------------------- | ------------ | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
 | bedtools_gtf_to_bed                              | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | bedtools_gtf_to_bed                              | run step         | NA            | cat /sbgenomics/Projects/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/gencode.v33.primary_assembly.annotation.gtf | grep -vE "^#" | cut -f 1,4,5 | awk '{OFS = "\t";a=$2-1;print $1,a,$3; }' | bedtools sort | bedtools merge > gencode.v33.primary_assembly.annotation.gtf.bed                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            | python -c 'def main():                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |     import sys                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |     bp_target = 20000000                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |     intvl_target_size = 20000                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |     bed_file = open("/sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/f85f4ba0-7927-435e-b39a-b8c6571baa4c/bedtools_gtf_to_bed/gencode.v33.primary_assembly.annotation.gtf.bed")                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                                                                                                                                                                                                          |
 | python_vardict_interval_split                              | run step         | NA            |     i=0                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |     intvl_set = {}                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |     cur_size = 0                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |     for cur_intvl in bed_file:                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         f = 0                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         if i not in intvl_set:                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |             intvl_set[i] = []                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         data = cur_intvl.rstrip("\n").split("\t")                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         (chrom, start, end) = (data[0], data[1], data[2])                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         intvl_size = int(end) - int(start)                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         if intvl_size >= bp_target:                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |             if len(intvl_set[i]) != 0:                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                 i += 1                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                 intvl_set[i] = []                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                 f = 1                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         elif cur_size + intvl_size > bp_target:                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |             if len(intvl_set[i]) != 0:                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                 i += 1                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                 intvl_set[i] = []                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                 cur_size = intvl_size                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         else:                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |             cur_size += intvl_size                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         intvl_set[i].append([chrom, start, end])                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         if f == 1:                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |             i += 1                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |             cur_size = 0                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |     bed_file.close()                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                                                                                                                                                                                                          |
 | python_vardict_interval_split                              | run step         | NA            |     for set_i, invtl_list in sorted(intvl_set.items()):                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         set_size = 0                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         out = open("set_" + str(set_i) + ".bed", "w")                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         for intervals in invtl_list:                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |             (chrom, start, end) = (intervals[0], intervals[1], intervals[2])                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |             intvl_size = int(end) - int(start)                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |             set_size += intvl_size                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |             for j in range(int(start), int(end), intvl_target_size):                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                 new_end = j + intvl_target_size                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                 if new_end > int(end):                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                     new_end = end                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                 out.write(chrom + "\t" + str(j) + "\t" + str(new_end) + "\n")                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         sys.stderr.write("Set " + str(set_i) + " size:\t" + str(set_size) + "\n")                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |         out.close()                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |                                                                                                                                                                                                          |
 | python_vardict_interval_split                              | run step         | NA            | if __name__ == "__main__":                                                                                                                                                                                                         |
 | python_vardict_interval_split                              | run step         | NA            |     main()'                                                                                                                                                                                                         |
 | vardict                              | scatter         | 89            | /bin/bash -c set -eo pipefail; export VAR_DICT_OPTS='"-Xms768m" "-Xmx6g"'; /VarDict-1.7.0/bin/VarDict -G /sbgenomics/Projects/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/GRCh38.primary_assembly.genome.fa -f 0.2 -th 4 --nosv --deldupvar -N VARDICT_NEW_SPLIT -b '/sbgenomics/Projects/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/da63df67-62a4-487b-aa68-d7f139809160.Aligned.out.sorted.bam' -z -c 1 -S 2 -E 3 -g 4 -F 0x700 -V 0.01 -x 150 /sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/f85f4ba0-7927-435e-b39a-b8c6571baa4c/python_vardict_interval_split/set_82.bed > vardict_results.txt && cat vardict_results.txt | /VarDict-1.7.0/bin/teststrandbias.R > vardict_r_test_results.txt && /VarDict-1.7.0/bin/var2vcf_valid.pl -N 'BS_5P3CZQV8'-E -f 0.2 -v 50 vardict_r_test_results.txt > VARDICT_NEW_SPLIT.result.vcf && cat VARDICT_NEW_SPLIT.result.vcf | perl -e 'while(<>){if ($_ =~ /^#/){print $_;} else{@a = split /\t/,$_; if($a[3] =~ /[KMRYSWBVHDXkmryswbvhdx]/){$a[3] = "N";} if($a[4] =~ /[KMRYSWBVHDXkmryswbvhdx]/){$a[4] = "N";} if($a[3] ne $a[4]){print join("\t", @a);}}}' > VARDICT_NEW_SPLIT.set_82.vcf && bgzip  VARDICT_NEW_SPLIT.set_82.vcf && tabix  VARDICT_NEW_SPLIT.set_82.vcf.gz                                                                                                                                                                                                         |
 | sort_merge_vardict_vcf                              | run step         | NA            | /gatk SortVcf --java-options "-Xmx6g" -O VARDICT_NEW_SPLIT.vardict.merged.vcf --SEQUENCE_DICTIONARY /sbgenomics/Projects/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/GRCh38.primary_assembly.genome.dict --CREATE_INDEX false -I /sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/f85f4ba0-7927-435e-b39a-b8c6571baa4c/vardict_1_s/VARDICT_NEW_SPLIT.set_0.vcf.gz -I /sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/f85f4ba0-7927-435e-b39a-b8c6571baa4c/vardict_2_s/VARDICT_NEW_SPLIT.set_1.vcf.gz -I /sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/f85f4ba0-7927-435e-b39a-b8c6571baa4c/vardict_3_s/VARDICT_NEW_SPLIT.set_10.vcf.gz -I /sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/f85f4ba0-7927-435e-b39a-b8c6571baa4c/vardict_4_s/VARDICT_NEW_SPLIT.set_11.vcf.gz -I /sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/f85f4ba0-7927-435e-b39a-b8c6571baa4c/vardict_5_s/VARDICT_NEW_SPLIT.set_12.vcf.gz -I /sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/f85f4ba0-7927-435e-b39a-b8c6571baa4c/vardict_6_s/VARDICT_NEW_SPLIT.set_13.vcf.gz -I /sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/f85f4ba0-7927-435e-b39a-b8c6571baa4c/vardict_7_s/VARDICT_NEW_SPLIT.set_14.vcf.gz && cat VARDICT_NEW_SPLIT.vardict.merged.vcf | uniq | bgzip > VARDICT_NEW_SPLIT.vardict.merged.vcf.gz && tabix VARDICT_NEW_SPLIT.vardict.merged.vcf.gz                                                                                                                                                                                                         |
 | gatk_pass_vcf                              | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | gatk_pass_vcf                              | run step         | NA            | /gatk SelectVariants --java-options "-Xmx7500m" -V /sbgenomics/workspaces/598f0ba4-d8a8-45e7-8bf2-1fe004e4979a/tasks/f85f4ba0-7927-435e-b39a-b8c6571baa4c/sort_merge_vardict_vcf/VARDICT_NEW_SPLIT.vardict.merged.vcf.gz -O VARDICT_NEW_SPLIT.vardict.PASS.vcf.gz --exclude-filtered TRUE                                                                                                                                                                                                         |


## Variant Effect Predictor
[Variant Effect Predictor](https://useast.ensembl.org/info/docs/tools/vep/index.html) is an ENSEMBL tool for annotating variants.
The tool built for this repo has very basic and rigid functionality, but can be run on any of the vcf outputs from the worfklows.
`tools/variant_effect_predictor.cwl`.
Run time (n=782) COV-IRT dataset,  ~6 minutes, cost on cavatica ~$0.22 per sample


### VEP simulated bash calls

 | Step                 | Type         | Num scatter            | Command                                                                                                                                                                                                         |
 | -------------------- | ------------ | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
 | vep-1oo-annotate                 | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | vep-1oo-annotate                 | run step         | NA            | tar -xzf /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/homo_sapiens_merged_vep_100_GRCh38.tar.gz                                                                                                                                                                                                       |
 | vep-1oo-annotate                 | run step         | NA            | perl /ensembl-vep/vep --cache --dir_cache $PWD --cache_version 100 --vcf --symbol  --merged  --canonical --variant_class --offline --ccds --uniprot --protein --numbers --hgvs --hgvsg --fork 14 --sift b --vcf_info_field ANN -i /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/3c3c570b-cd69-4106-b169-ffa65411ec5a.STAR_GATK4.PASS.vcf.gz -o STDOUT --stats_file 9aeaeaf3-0e59-4c92-a709-c6bd37431294_stats.txt --stats_text --warning_file 9aeaeaf3-0e59-4c92-a709-c6bd37431294_warnings.txt --allele_number --dont_skip --allow_non_variant --fasta /sbgenomics/Projects/7c22397d-ea11-4069-9687-3049119a2d37/NASA_HEC_transferred_data/genome_files/Homo_sapiens_ensembl_release100/Homo_sapiens.GRCh38.dna.primary_assembly_mainchrs_and_SARS-CoV-2-NC_045512.2.fa | /ensembl-vep/htslib/bgzip -@ 14 -c > 9aeaeaf3-0e59-4c92-a709-c6bd37431294.STAR_GATK4.vep.vcf.gz && /ensembl-vep/htslib/tabix 9aeaeaf3-0e59-4c92-a709-c6bd37431294.STAR_GATK4.vep.vcf.gz                                                                                                                                                                                                         |

```yaml
inputs:
  reference: {type: File,  secondaryFiles: [.fai], label: Fasta genome assembly with index}
  input_vcf:
    type: File
    secondaryFiles: [.tbi]
  output_basename: string
  merged_cache: {type: boolean, doc: "If merged cache being used", default: true}
  tool_name: {type: string, doc: "Name of tool used to generate calls"}
  cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  ```

```yaml
outputs:
  output_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]
  output_txt:
    type: File
    outputBinding:
      glob: '*_stats.txt'
  warn_txt:
    type: ["null", File]
    outputBinding:
      glob: '*_warnings.txt'
```


## Kraken2
[Kraken2](http://ccb.jhu.edu/software/kraken2/index.shtml) is available to run at `tools/kraken2_classification.cwl`.

```yaml
inputs:
  input_db: { type: File, doc: "Input TGZ containing Kraken2 database" }
  input_reads: { type: File, doc: "FA or FQ file containing sequences to be classified" }
  input_mates: { type: 'File?', doc: "Paired mates for input_reads" }
  db_path: { type: string, default: "./covid", doc: "Relative path to the folder containing the db files from input_db" }
  threads: { type: int, default: 32, doc: "Number of threads to use in parallel" }
  ram: { type: int, default: 50000, doc: "Recommended MB of RAM needed to run the job" }
  output_basename: { type: string, doc: "String to be used as the base filename of the output" }
```

```yaml
outputs:
  output: { type: File, outputBinding: { glob: "*.output" } }
  classified_reads: { type: 'File', outputBinding: { glob: "*_1.fq" } }
  classified_mates: { type: 'File?', outputBinding: { glob: "*_2.fq" } }
```
