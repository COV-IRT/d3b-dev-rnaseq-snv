cwlVersion: v1.0
class: CommandLineTool
id: star_fusion_ref_gen
doc: "Generate STAR Fusion reference file from custom annotation using STAR Fusion v1.9.0"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'trinityctat/starfusion:1.9.0'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 64000

baseCommand: [/usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --genome_fa $(inputs.reference_genome.path) 
      --gtf $(inputs.gtf.path)
      --pfam_db current
      --dfam_db human  
      --max_readlength $(inputs.max_readlength)
      --CPU 8

      rm ctat_genome_lib_build_dir/*.gz

      rm -rf ctat_genome_lib_build_dir/_chkpts

      tar -czf $(inputs.new_reference_name).tgz ctat_genome_lib_build_dir

inputs:
  reference_genome: File
  gtf: File
  max_readlength: {type: ['null', int], doc: "Read length of library aligned", default: 150}
  new_reference_name: {type: string, doc: "Name to use for newly created star fusion archive"}

outputs:
  star_fusion_reference:
    type: File
    outputBinding:
      glob: '*.tgz'
