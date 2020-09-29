cwlVersion: v1.0
class: CommandLineTool
id: star_fusion_covirt
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'trinityctat/starfusion:1.9.0'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 64000

baseCommand: [tar]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -zxf $(inputs.genome_tar.path) &&
      /usr/local/src/STAR-Fusion/STAR-Fusion
      --genome_lib_dir ./$(inputs.genome_untar_path)
      -J $(inputs.Chimeric_junction.path)
      --output_dir STAR-Fusion_outdir
      --examine_coding_effect
      --CPU 16 &&
      mv STAR-Fusion_outdir/star-fusion.fusion_predictions.abridged.coding_effect.tsv $(inputs.SampleID).STAR.fusion_predictions.abridged.coding_effect.tsv &&
      gzip -c $(inputs.Chimeric_junction.path) > $(inputs.Chimeric_junction.basename).gz
      

inputs:
  Chimeric_junction: {type: File, doc: "Output from STAR alignment run"}
  genome_tar: {type: File, doc: "STAR-Fusion reference"}
  genome_untar_path: {type: ['null', string], doc: "This is what the path will be when genome_tar is unpackaged", default: "ctat_genome_lib_build_dir"}
  SampleID: string

outputs:
  abridged_coding:
    type: File
    outputBinding:
      glob: '*.fusion_predictions.abridged.coding_effect.tsv'
  chimeric_junction_compressed:
    type: File
    outputBinding:
      glob: "$(inputs.Chimeric_junction.basename).gz"
