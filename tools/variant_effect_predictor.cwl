cwlVersion: v1.0
class: CommandLineTool
id: vep-1oo-annotate
doc: "VEP Release 100. Basic tool with no customization"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 24000
    coresMin: 14
  - class: DockerRequirement
    dockerPull: 'migbro/vep:r100'
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      tar -xzf
      $(inputs.cache.path)
      && perl /ensembl-vep/vep
      --cache --dir_cache $PWD
      --cache_version 100
      --vcf
      --symbol
      ${
        if (inputs.merged_cache){
          return " --merged ";
        }
        else{
          return " ";
        }
      }
      --canonical
      --variant_class
      --offline
      --ccds
      --uniprot
      --protein
      --numbers
      --hgvs
      --hgvsg
      --fork 14
      --sift b
      --vcf_info_field ANN
      -i $(inputs.input_vcf.path)
      -o STDOUT
      --stats_file $(inputs.output_basename)_stats.txt
      --stats_text
      --warning_file $(inputs.output_basename)_warnings.txt
      --allele_number
      --dont_skip
      --allow_non_variant
      --fasta $(inputs.reference.path) |
      /ensembl-vep/htslib/bgzip -@ 14 -c > $(inputs.output_basename).$(inputs.tool_name).vep.vcf.gz
      && /ensembl-vep/htslib/tabix $(inputs.output_basename).$(inputs.tool_name).vep.vcf.gz

inputs:
  reference: {type: File,  secondaryFiles: [.fai], label: Fasta genome assembly with index}
  input_vcf:
    type: File
    secondaryFiles: [.tbi]
  output_basename: string
  merged_cache: {type: boolean, doc: "If merged cache being used", default: true}
  tool_name: {type: string, doc: "Name of tool used to generate calls"}
  cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}

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
