cwlVersion: v1.0
class: CommandLineTool
id: strelka2-rnaseq
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: ${return inputs.ram * 1000}
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'kfdrc/strelka2:2.9.10'
  - class: InitialWorkDirRequirement
    listing: [$(inputs.input_rna_bam),$(inputs.input_rna_bai)]


baseCommand: [/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --bam $(inputs.input_rna_bam.path)
      --reference $(inputs.reference.path)
      ${
        if (inputs.strelka2_bed != null){
          return "--callRegions " + inputs.strelka2_bed.path;
        }
        else{
          return "";
        }
      }
      --rna
      --runDir ./
      && ./runWorkflow.py
      -m local
      -j $(inputs.cores)
      -g $(inputs.ram)

      mv results/variants/variants.vcf.gz $(inputs.output_basename).strelka2.rnaseq.vcf.gz

      mv results/variants/variants.vcf.gz.tbi $(inputs.output_basename).strelka2.rnaseq.vcf.gz.tbi

inputs:
  reference: { type: File, secondaryFiles: [.fai] }
  input_rna_bam: {type: File}
  input_rna_bai: {type: File}
  strelka2_bed: {type: File?, secondaryFiles: [.tbi], label: gzipped bed file}
  cores: {type: ['null', int], default: 16, doc: "Num cores to use"}
  ram: {type: ['null', int], default: 30, doc: "Max mem to use in GB"}
  output_basename: string
outputs:
  output_vcf:
    type: File
    outputBinding:
      glob: '*.strelka2.rnaseq.vcf.gz'
    secondaryFiles: [.tbi]

