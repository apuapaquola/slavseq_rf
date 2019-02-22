# slavseq_rf 

Analysis pipeline for human LINE-1 targeted sequencing data obtained with the SLAV-Seq protocol.

### Publication

https://doi.org/10.1038/nn.4388


### Prerequisites

- rf 
- pyslavseq
- snakemake
- cutadapt
- bwa-mem
- tabix
- python
- perl


### How to use

This pipeline uses rf (https://github.com/apuapaquola/rf) to manage workflow.

- clone the repository
- move the pipeline directory to a suitable location
- edit _h/run scripts or set the BASE_DIR to /path/to/pipeline
- go to pipeline/genome and <code>rf run -r .</code> to download and index hg19
- add your fastq files
- go to pipeline/genome  and <code>rf run -r .</code>
