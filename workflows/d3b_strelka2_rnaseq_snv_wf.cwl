cwlVersion: v1.0
class: Workflow
doc: |-
    Strelka2 SNV Calling Workflow
    
    This [workflow](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#rna-seq) is pretty straight forward, with a `PASS` filter step added to get `PASS` calls.
    The git repo serving this app and related tools can be found [here](https://github.com/d3b-center/d3b-dev-rnaseq-snv), compatible with all releases as of 2020-Jun-17.
    `workflows/d3b_strelka2_rnaseq_snv_wf.cwl` is the wrapper cwl that runs this workflow

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

    ### Simulated bash calls
    An example of bash calls from each step can be found in the [git repo](https://github.com/d3b-center/d3b-dev-rnaseq-snv#strelka2-simulated-bash-calls)

id: d3b-strelka2-rnaseq-snv-wf
label: "Strelka2 RNAseq SNV Calling Workflow"
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  output_basename: string
  STAR_sorted_genomic_bam: {type: File, doc: "STAR sorted alignment bam"}
  STAR_sorted_genomic_bai: {type: File, doc: "STAR sorted alignment bam index"}
  reference_fasta: {type: File, secondaryFiles: ['.fai'], doc: "Reference genome used"}
  tool_name: {type: string, doc: "description of tool that generated data, i.e. gatk_haplotypecaller"}
  mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}

outputs:
  strelka2_prepass_vcf: {type: File, outputSource: strelka2_rnaseq/output_vcf, doc: "Strelka2 SNV calls"}
  strelka2_pass_vcf: {type: File, outputSource: gatk_pass_vcf/pass_vcf, doc: "Strelka2 calls filtered on PASS"}

steps:
  strelka2_rnaseq:
    run: ../tools/strelka2_rnaseq.cwl
    label: Strelka2
    in:
      reference: reference_fasta
      input_rna_bam: STAR_sorted_genomic_bam
      input_rna_bai: STAR_sorted_genomic_bai
      # strelka2_bed: {type: File?, secondaryFiles: [.tbi], label: gzipped bed file}
      output_basename: output_basename
    out: [output_vcf]

  gatk_pass_vcf:
    run: ../tools/gatk_selectvariants.cwl
    label: "GATK Select PASS"
    in:
      input_vcf: strelka2_rnaseq/output_vcf
      output_basename: output_basename
      tool_name: tool_name
      mode: mode
      init_bam_index_style: string
    out: [pass_vcf]
