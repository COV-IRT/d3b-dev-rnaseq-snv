cwlVersion: v1.0
class: CommandLineTool
id: vep-1oo-annotate
doc: "VEP Release 100. Basic annotation tool, can use cache or vcf"
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

      ${
        if (inputs.cache){
          var cmd = "tar -xzf " + inputs.cache.path + ";";
          return cmd;
          }
        else{
          return "echo No cache skipping unzip;";
          }
      }
      
      perl /ensembl-vep/vep
      --vcf
      --canonical
      --variant_class
      --ccds
      --uniprot
      --protein
      --numbers
      --hgvs
      --hgvsg
      --fork 14
      --sift b
      --symbol
      --vcf_info_field ANN
      -i $(inputs.input_vcf.path)
      -o STDOUT
      --stats_file $(inputs.output_basename)_stats.txt
      --stats_text
      --warning_file $(inputs.output_basename)_warnings.txt
      --allele_number
      --dont_skip
      --allow_non_variant
      --fasta $(inputs.reference.path)
      ${
        var args = "";
        if (inputs.cache){
          args = " --cache --dir_cache $PWD --cache_version 100 --offline ";
          if (inputs.merged_cache){
            args += " --merged ";
          }
        }
        else{
          args = " --gtf " + inputs.bgzipped_gtf.path;
        }
        return args;
      }
      |
      /ensembl-vep/htslib/bgzip -@ 2 -c > $(inputs.output_basename).$(inputs.tool_name).vep.vcf.gz
      && /ensembl-vep/htslib/tabix $(inputs.output_basename).$(inputs.tool_name).vep.vcf.gz

inputs:
  reference: {type: File,  secondaryFiles: [.fai], label: Fasta genome assembly with index}
  input_vcf:
    type: File
    secondaryFiles: [.tbi]
  output_basename: string
  merged_cache: {type: boolean?, doc: "If merged cache being used", default: true}
  tool_name: {type: string, doc: "Name of tool used to generate calls"}
  cache: {type: File?, label: tar gzipped cache from ensembl/local converted cache, doc: "Use this if not using a gtf for gene models"}
  bgzipped_gtf: {type: File?, doc: "If merged cache being used", secondaryFiles: ['.tbi'], doc: "Use this if not using a cahce, but using gtf instead for gene models"}

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
