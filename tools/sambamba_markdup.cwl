cwlVersion: v1.0
class: CommandLineTool
id: sambamba_sort_mark_dup

requirements:
  - class: DockerRequirement
    dockerPull: 'kfdrc/sambamba:0.7.1'
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      mkdir TMP

      sambamba markdup
      --tmpdir TMP
      -t 4
      $(inputs.input_bam.path)
      $(inputs.input_bam.nameroot).md.bam

      mv $(inputs.input_bam.nameroot).md.bam.bai $(inputs.input_bam.nameroot).md.bai

inputs:
  input_bam: File
outputs: 
  markduplicates_bam:
    type: File
    outputBinding:
      glob: '*.bam'
