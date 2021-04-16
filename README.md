# PinAPL-Py
![pinapl-py logo_small](https://cloud.githubusercontent.com/assets/23485218/21127117/6e0cac76-c0a5-11e6-8109-fce62e21c148.png)

**!! See the Tutorial (PinAPL-Py Tutorial.pdf) for instructions !!**

Experienced users can follow the Quick Start below:

## Quick Start 

0. *Install docker* on your machine: https://docs.docker.com/engine/installation/

1. Create a *working directory* and create both a /Data and /Library subfolder inside it

2. Copy your read files (fastq.gz) to the /Data folder 

3. *Rename your read files*. To be recognized and correctly interpreted, the fastq.gz file names of control replicates need to start with *“Control_R1_...”,  “Control_R2_....”* etc. Similarly, file names of treatment replicates (e.g. treatmentX) need to start with *“TreatmentX_R1_...”, “TreatmentX_R2_...”* etc. For more information on filename requirements see Section 2 of the tutorial.

4. Copy the *library file (.tsv)* to the /Library folder. The library file is a tab delimited file of all sgRNAs in the library. Columns should be 1: gene_ID, 2: sgRNA_ID, 3: sequence. If you work with the GECKO_v21 library, you can download this file from https://github.com/LewisLabUCSD/PinAPL-Py
  
5. Download `configuration.yaml` from https://github.com/LewisLabUCSD/PinAPL-Py and copy it to the working directory.

6. Edit `configuration.yaml` 

	NOTE: For GECKO_v2 enrichment screens, you can leave everything at default.
  
7. Dockerfile belongs to v2.9 and can be found in https://hub.docker.com/r/oncogx/pinaplpy_docker 


8. Start PinAPL-Py from the working directory using Terminal (Windows: Docker Quickstart Terminal)

   `docker run -t -i --name pinaplpy_test -v $PWD:/workingdir oncogx/pinaplpy_docker`
  

9. Run PinAPL-py using python from the workingdirectory.

   `python /Scripts/PinALP-py`

  
