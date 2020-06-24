cwlVersion: v1.0
class: Workflow
doc: |-
    GATK RNAseq SNV Calling Workflow

    The overall [workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) picks up from post-STAR alignment, starting at Picard mark duplicates.
    For the most part, tool parameters follow defaults from the GATK Best Practices [WDL](https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl), written in cwl with added optimatization for use on the Cavatica platform.
    The git repo serving this app and related tools can be found [here](https://github.com/d3b-center/d3b-dev-rnaseq-snv/releases/tag/0.4).
    `workflows/d3b_gatk_rnaseq_snv_wf.cwl` is the wrapper cwl used to run all tools for GATK4.

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

    ### Simulated bash calls
    An example of bash calls from each step can be found in the [git repo](https://github.com/d3b-center/d3b-dev-rnaseq-snv#gatk4-simulated-bash-calls)

id: d3b-gatk-rnaseq-snv-wf
label: "GATK RNAseq SNV Calling Workflow"
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

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

outputs:
  filtered_hc_vcf: {type: File, outputSource: gatk_filter_vcf/filtered_vcf, doc: "Haplotype called vcf with Broad-recommended FILTER values added"}
  pass_vcf: {type: File, outputSource: gatk_pass_vcf/pass_vcf, doc: "Filtered vcf selected for PASS variants"}
  anaylsis_ready_bam: {type: File, outputSource: gatk_applybqsr/recalibrated_bam, doc: "Duplicate marked, Split N trimmed CIGAR BAM, BQSR recalibratede, ready for RNAseq calling"}
  bqsr_table: {type: File, outputSource: gatk_baserecalibrator/output, doc: "BQSR table"}

steps:
  bedtools_gtf_to_bed:
    run: ../tools/bedtools_gtf_to_bed.cwl
    in:
      input_bed_gtf: call_bed_file
    out: [run_bed]
  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: bedtools_gtf_to_bed/run_bed
      reference_dict: reference_dict
      exome_flag: exome_flag
      scatter_ct: scatter_ct
      bands:
        valueFrom: ${return 80000000}
    out: [output]
  preprocess_rnaseq_bam:
    run: ../subworkflows/preprocess_rnaseq_bam.cwl
    label: "Preprocess RNAseq BAM"
    in:
      STAR_sorted_genomic_bam: STAR_sorted_genomic_bam
      reference_fasta: reference_fasta
    out:
      [sorted_md_splitn_bam]
  gatk_baserecalibrator:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge
    run: ../tools/gatk_baserecalibrator.cwl
    label: "GATK BQSR"
    in:
      input_bam: preprocess_rnaseq_bam/sorted_md_splitn_bam
      knownsites: knownsites
      reference: reference_fasta
    out: [output]
  gatk_applybqsr:
    run: ../tools/gatk_applybqsr.cwl
    label: "GATK Apply BQSR"
    in:
      reference: reference_fasta
      input_bam: preprocess_rnaseq_bam/sorted_md_splitn_bam
      bqsr_report: gatk_baserecalibrator/output
    out: [recalibrated_bam]
  gatk_haplotype_rnaseq:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.4xlarge
    run: ../tools/gatk_haplotypecaller_rnaseq.cwl
    label: "GATK Haplotype Caller"
    in:
      reference_fasta: reference_fasta
      bqsr_bam: gatk_applybqsr/recalibrated_bam
      dbsnp: dbsnp_vcf
      output_basename: output_basename
      genes_bed: gatk_intervallisttools/output
    scatter: genes_bed
    out: [hc_called_vcf]
  merge_hc_vcf:
    run: ../tools/gatk_mergevcfs.cwl
    label: "GAK Merge HC VCF"
    in:
      input_vcfs: gatk_haplotype_rnaseq/hc_called_vcf
      reference_dict: reference_dict
      output_basename: output_basename
      tool_name: tool_name
    out: [merged_vcf]
  gatk_filter_vcf:
    run: ../tools/gatk_filtervariants.cwl
    label: "GATK Hard Filter VCF"
    in:
      reference_fasta: reference_fasta
      hc_vcf: merge_hc_vcf/merged_vcf
      output_basename: output_basename
    out: [filtered_vcf]
  gatk_pass_vcf:
    run: ../tools/gatk_selectvariants.cwl
    label: "GATK Select PASS"
    in:
      input_vcf: gatk_filter_vcf/filtered_vcf
      output_basename: output_basename
      tool_name: tool_name
      mode: mode
    out: [pass_vcf]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2