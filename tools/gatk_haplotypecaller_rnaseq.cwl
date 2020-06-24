cwlVersion: v1.0
class: CommandLineTool
id: gatk-haplotypecaller-rnaseq
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.7.0R'
baseCommand: [/gatk, HaplotypeCaller]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -R $(inputs.reference_fasta.path)
      -I $(inputs.bqsr_bam.path)
      --standard-min-confidence-threshold-for-calling 20
      -dont-use-soft-clipped-bases
      ${
        if (inputs.genes_bed != null){
          return "-L " + inputs.genes_bed.path;
        }
        else{
          return "";
        }
      }
      -O $(inputs.output_basename).gatk.hc.called.vcf.gz
      --dbsnp $(inputs.dbsnp.path)
inputs:
  reference_fasta: {type: File, secondaryFiles: ['.fai', '^.dict']}
  bqsr_bam: {type: File, secondaryFiles: ['^.bai']}
  genes_bed: {type: File?}
  dbsnp: {type: File, secondaryFiles: ['.idx'], doc: "dbSNP reference"}
  output_basename: string
outputs:
  hc_called_vcf:
    type: File
    outputBinding:
      glob: '*.gatk.hc.called.vcf.gz'
    secondaryFiles: ['.tbi']
