# Get code and install in python env

```sh
git clone https://github.com/kmhernan/gdc-fastq-splitter
cd gdc-fastq-splitter
```
*if you don't have it yet*, get virtual env for python3

`pip3 install virtualenv # depends on your set up`

```sh
virtualenv venv --python /usr/bin/python3 # again, wherever your python3 (should be version 3.5+) is installed
./venv/bin/pip install .
venv/bin/gdc-fastq-splitter -h # if successful, you should get help menu
```

# Fun time run time
Note that processing paired end files will use two cores instead on one (likely a good thing).

```text
inputs: 
 -rw-rw-r-- 1 ubuntu ubuntu  84M May  7 04:31 COVHA-P1-D06-N.filtered.human.R1.fastq.gz
 -rw-rw-r-- 1 ubuntu ubuntu  88M May  7 04:31 COVHA-P1-D06-N.filtered.human.R2.fastq.gz
```

```sh
gdc-fastq-splitter/venv/bin/gdc-fastq-splitter COVHA-P1-D06-N.filtered.human.R1.fastq.gz COVHA-P1-D06-N.filtered.human.R2.fastq.gz -o COVHA-P1-D06-N.filtered.human_ 2> errs.log 
```


## Result outputs:

```text
 -rw-rw-r-- 1 ubuntu ubuntu  34M Jul  1 18:38 COVHA-P1-D06-N.filtered.human_H5NKGDSXY_1_R1.fq.gz
 -rw-rw-r-- 1 ubuntu ubuntu 8.9K Jul  1 18:38 COVHA-P1-D06-N.filtered.human_H5NKGDSXY_1_R1.report.json
 -rw-rw-r-- 1 ubuntu ubuntu  36M Jul  1 18:38 COVHA-P1-D06-N.filtered.human_H5NKGDSXY_1_R2.fq.gz
 -rw-rw-r-- 1 ubuntu ubuntu 8.9K Jul  1 18:38 COVHA-P1-D06-N.filtered.human_H5NKGDSXY_1_R2.report.json
 -rw-rw-r-- 1 ubuntu ubuntu  34M Jul  1 18:38 COVHA-P1-D06-N.filtered.human_H5NKGDSXY_2_R1.fq.gz
 -rw-rw-r-- 1 ubuntu ubuntu 9.6K Jul  1 18:38 COVHA-P1-D06-N.filtered.human_H5NKGDSXY_2_R1.report.json
 -rw-rw-r-- 1 ubuntu ubuntu  36M Jul  1 18:38 COVHA-P1-D06-N.filtered.human_H5NKGDSXY_2_R2.fq.gz
 -rw-rw-r-- 1 ubuntu ubuntu 9.6K Jul  1 18:38 COVHA-P1-D06-N.filtered.human_H5NKGDSXY_2_R2.report.json
```

json outputs have relevant RG info, metadata field most relevant, early barcode field is informative

from `COVHA-P1-D06-N.filtered.human_H5NKGDSXY_1_R1.report.json`

```json
   "metadata": {
     "fastq_filename": "COVHA-P1-D06-N.filtered.human_H5NKGDSXY_1_R1.fq.gz",
     "flowcell_barcode": "H5NKGDSXY",
    "lane_number": 1,
     "multiplex_barcode": "TGAACACC+CGCCTAGG",
     "record_count": 494598
   }
```

from `COVHA-P1-D06-N.filtered.human_H5NKGDSXY_2_R1.report.json`

```json
 "metadata": {
     "fastq_filename": "COVHA-P1-D06-N.filtered.human_H5NKGDSXY_2_R1.fq.gz",
     "flowcell_barcode": "H5NKGDSXY",
     "lane_number": 2,
     "multiplex_barcode": "TGAACACC+CGCCTAGG",
     "record_count": 488895
   }
```

Recommended RGs based on these two files:

```text
@RG ID:H5NKGDSXY.1	PL:ILLUMINA PU:H5NKGDSXY.1.TGAACACC+CGCCTAGG	LB:COVHA-P1-D06-N.filtered.human 	:COVHA-P1-D06-N.filtered.human
@RG	ID:H5NKGDSXY.2	PL:ILLUMINA	PU:H5NKGDSXY.1.TGAACACC+CGCCTAGG	LB:COVHA-P1-D06-N.filtered.human SM:COVHA-P1-D06-N.filtered.human
```

# Updated STAR Command, based on [Google Drive doc](https://drive.google.com/file/d/1z7euPSE77jbthKi8ty7dscq-oNx4wgF3/view?usp=sharing) follwing instruction from [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf), page 8, section 3.2:

```sh
STAR --twopassMode Basic \
    --limitBAMsortRAM 105000000000 \
--genomeDir /path/to/STAR/genome/directory \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS nM NM MD jM jI MC ch \ #same as All
    --outSAMattrRGline ID:H5NKGDSXY.1	PL:ILLUMINA PU:H5NKGDSXY.1.TGAACACC+CGCCTAGG	LB:COVHA-P1-D06-N.filtered.human 	SM:COVHA-P1-D06-N.filtered.human , ID:H5NKGDSXY.2	PL:ILLUMINA	PU:H5NKGDSXY.1.TGAACACC+CGCCTAGG	LB:COVHA-P1-D06-N.filtered.human SM:COVHA-P1-D06-N.filtered.human \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --readFilesCommand zcat \
    --runThreadN 24 \
    --chimOutType Junctions \
    --chimOutJunctionFormat 1
    --chimSegmentMin 20 \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outFileNamePrefix /path/to/output/directory/${sample} \
    --readFilesIn COVHA-P1-D06-N.filtered.human_H5NKGDSXY_1_R1.fq.gz,COVHA-P1-D06-N.filtered.human_H5NKGDSXY_2_R1.fq.gz COVHA-P1-D06-N.filtered.human_H5NKGDSXY_1_R2.fq.gz,COVHA-P1-D06-N.filtered.human_H5NKGDSXY_2_R2.fq.gz
```