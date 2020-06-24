cwlVersion: v1.0
class: CommandLineTool
id: gatk_gatherbamfiles
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.7.0R'
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 8
baseCommand: [/gatk, GatherBamFiles]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -O=$(inputs.output_bam_basename).bam
      --CREATE_INDEX=true
      --CREATE_MD5_FILE=true
inputs:
  input_bam:
    type:
      type: array
      items: File
      inputBinding:
        prefix: -I
        separate: false
    secondaryFiles: ['.bai']
  output_bam_basename: string
outputs:
  output:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [.bai, .md5]
