cwlVersion: v1.0
class: Workflow
doc: >-
    "Preprocessing sub workflow for GATK4 RNAseq SNV calling.
    Uses Sambamba for mark dups and sort, GATK SplitNCigarReads for final processed bam output"
id: d3b-sambamba-subwf-wf
label: "Sambamba Sort Mark Dup Sub WF"
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  STAR_sorted_genomic_bam: {type: File, doc: "STAR sorted alignment bam"}
  reference_fasta: {type: File, secondaryFiles: ['^.dict', '.fai'], doc: "Reference genome used"}
outputs:
  sorted_md_splitn_bam: {type: File, outputSource: gatk_splitntrim/cigar_n_split_bam, secondaryFiles: ['^.bai'], doc: "Dup marked, sorted, Split N trim cigar bam"}

steps:
  sambamba_md_sorted:
    run: ../tools/sambamba_markdup.cwl
    label: "Sambamba Mardup"
    in:
      input_bam: STAR_sorted_genomic_bam
    out: [markduplicates_bam]
  gatk_splitntrim:
    run: ../tools/gatk_splitncigarreads.cwl
    label: "GATK Split N Cigar"
    in:
      reference_fasta: reference_fasta
      dup_marked_bam: sambamba_md_sorted/markduplicates_bam
    out: [cigar_n_split_bam]
