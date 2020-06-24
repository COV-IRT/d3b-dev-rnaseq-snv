cwlVersion: v1.0
class: CommandLineTool
id: picard_markduplicates_spark

requirements:
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.7.0R'
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 32000
    coresMin: 16
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      mkdir TMP

      /gatk MarkDuplicatesSpark
      --java-options "-Xmx30000m
      -XX:+PrintFlagsFinal
      -Xloggc:gc_log.log
      -XX:GCTimeLimit=50
      -XX:GCHeapFreeLimit=10"
      --tmp-dir TMP
      -I=$(inputs.input_bam.path)
      -O=$(inputs.input_bam.nameroot).dup_marked.bam
inputs:
  input_bam: File
outputs: 
  output_markduplicates_bam:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [.bai]
