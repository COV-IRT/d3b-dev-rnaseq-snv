cwlVersion: v1.0
class: CommandLineTool
id: gatk-splitncigarreads
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 32000
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.7.0R'
baseCommand: [/gatk, SplitNCigarReads]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx30G
      -XX:+PrintFlagsFinal
      -Xloggc:gc_log.log
      -XX:GCTimeLimit=50
      -XX:GCHeapFreeLimit=10"
      --seconds-between-progress-updates 30
      -R $(inputs.reference_fasta.path)
      -I $(inputs.dup_marked_bam.path)
      ${
        if (inputs.interval_bed != null){
          return "-L " + inputs.interval_bed.path;
        }
        else{
          return "";
        }
      }
      -OBI
      ${
          if (inputs.output_basename != null){
              return "-O " + inputs.output_basename + ".sorted.dup_marked.splitn.bam";
          }
          else{
              return "-O " + inputs.dup_marked_bam.nameroot + ".splitn.bam";
          }
      }
      
inputs:
  reference_fasta: {type: File, secondaryFiles: ['.fai', '^.dict']}
  dup_marked_bam: {type: File, secondaryFiles: ['^.bai']}
  interval_bed: {type: File?}
  output_basename: {type: string?, doc: "Can leave blank if using an inerval file"}
outputs:
  cigar_n_split_bam:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: ['^.bai']
