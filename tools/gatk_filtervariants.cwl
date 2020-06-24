cwlVersion: v1.0
class: CommandLineTool
id: gatk-filter-variants
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.7.0R'
baseCommand: [/gatk, VariantFiltration]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -R $(inputs.reference_fasta.path)
      -V $(inputs.hc_vcf.path)
      --window 35
      --cluster 3
      --filter-name "FS"
      --filter "FS > 30.0"
      --filter-name "QD"
      --filter "QD < 2.0"
      ${
        if (inputs.genes_bed != null){
          return "-L " + inputs.genes_bed.path;
        }
        else{
          return "";
        }
      }
      -O $(inputs.output_basename).gatk.hc.filtered.vcf.gz

inputs:
  reference_fasta: {type: File, secondaryFiles: ['.fai', '^.dict']}
  hc_vcf: {type: File, secondaryFiles: ['.tbi']}
  genes_bed: {type: File?}
  output_basename: string
outputs:
  filtered_vcf:
    type: File
    outputBinding:
      glob: '*.gatk.hc.filtered.vcf.gz'
    secondaryFiles: ['.tbi']
