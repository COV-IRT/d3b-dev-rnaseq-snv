cwlVersion: v1.0
class: Workflow
doc: |-
    VarDict Java RNAseq SNV Calling Workflow

    This [workflow](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/rnaseq/variation.py) is based on the Vardict run style of BC Bio, compatible with all releases as of 2020-Jun-17.
    `workflows/d3b_vardict_rnaseq_snv_wf.cwl` is the wrapper cwl that runs this workflow.
    Tweaking `vardict_bp_target` and `vardict_intvl_target_size` maybe be needed to improve run time in high coverage areas, by reducing their values from defaults.

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

    ### Simulated bash calls
    An example of bash calls from each step can be found in the [git repo](https://github.com/d3b-center/d3b-dev-rnaseq-snv#vardict-simulated-bash-calls)


id: d3b-vardict-rnaseq-snv-wf
label: "Vardict Java RNAseq SNV Calling Workflow"
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

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

outputs:
  vardict_prepass_vcf: {type: File, outputSource: sort_merge_vardict_vcf/merged_vcf, doc: "VarDict SNV calls"}
  vardict_pass_vcf: {type: File, outputSource: gatk_pass_vcf/pass_vcf, doc: "VarDict calls filtered on PASS"}

steps:
  bedtools_gtf_to_bed:
    run: ../tools/bedtools_gtf_to_bed.cwl
    in:
      input_bed_gtf: call_bed_file
    out: [run_bed]
  python_vardict_interval_split:
    run: ../tools/python_vardict_interval_split.cwl
    doc: "Custom interval list generation for vardict input. Briefly, ~60M bp per interval list, 20K bp intervals, lists break on chr and N regions only"
    in:
      wgs_bed_file: bedtools_gtf_to_bed/run_bed
      bp_target: vardict_bp_target
      intvl_target_size: vardict_intvl_target_size
    out: [split_intervals_bed]
  vardict:
    run: ../tools/vardict_rnaseq.cwl
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.4xlarge
    label: "VarDict Java"
    in:
      input_bam: STAR_sorted_genomic_bam
      input_name: sample_name
      padding: padding
      min_vaf: vardict_min_vaf
      cpus: vardict_cpus
      ram: vardict_ram
      reference: reference_fasta
      bed: python_vardict_interval_split/split_intervals_bed
      output_basename: output_basename
    scatter: [bed]
    out: [vardict_vcf]
  sort_merge_vardict_vcf:
    run: ../tools/gatk_sortvcf.cwl
    label: GATK Sort & merge vardict
    in:
      input_vcfs: vardict/vardict_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${return "vardict"}
    out: [merged_vcf]
  gatk_pass_vcf:
    run: ../tools/gatk_selectvariants.cwl
    label: "GATK Select PASS"
    in:
      input_vcf: sort_merge_vardict_vcf/merged_vcf
      output_basename: output_basename
      tool_name: tool_name
      mode: mode
    out: [pass_vcf]
$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2