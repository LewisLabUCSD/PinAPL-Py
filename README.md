# PinAPL-Py
![pinapl-py logo_small](https://cloud.githubusercontent.com/assets/23485218/21127117/6e0cac76-c0a5-11e6-8109-fce62e21c148.png)

See the **Tutorial (PinAPL-Py Tutorial.pdf)** for instructions!

Experienced users can follow the Quick Start below:

## Quick Start 

0. *Install docker* on your machine: https://docs.docker.com/engine/installation/

1. Create a *working directory* and create both a /Data and /Library folder inside it

2. Copy your read files (fastq.gz) to the /Data folder and rename them as *\<samplename>\_R\<number>.fastq.gz* 

	e.g. TreatmentA_R1.fastq.gz / TreatmentA_R2.fastq.gz / TreatmentB_R1.fastq.gz /...etc. for your treatment samples
	
	Control_R1.fastq.gz / Control_R2.fastq.gz / ... etc. for your control samples
3. Copy the *library file (.tsv)* to the /Library folder. Columns should be Gene, sgRNA_ID, Sequence. If you work with GECKO_v2 library, you can download this file from https://github.com/LewisLabUCSD/PinAPL-Py
  
4. Download *configuration.yaml* from https://github.com/LewisLabUCSD/PinAPL-Py and copy it to the working directory.

5. Edit configuration.yaml 

	NOTE: For GECKO_v2 enrichment screens, you can leave everything at default.
  
6. Start PinAPL-Py from the working directory using Terminal (Windows: Docker Quickstart Terminal)

  docker  run -v $PWD:/workingdir oharismendy/pinaplpy_docker PinAPL.py
  
7. Refer to the tutorial (PinAPL-Py Tutorial.pdf) for a full documentation on PinAPL-Py
  
