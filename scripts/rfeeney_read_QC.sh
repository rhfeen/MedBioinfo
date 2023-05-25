#!/bin/bash
echo "script start: Downloads sequencing files and run initial sequencing read quality control"
date

# List sequencing runs that I assigned myself to
sqlite3 -batch /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db /
"SELECT patient_code, username FROM sample_annot LEFT JOIN sample2bioinformatician /
USING(patient_code) WHERE username='rfeeney';"

# Save the run_accessions for these sequencing runs in a text file 
sqlite3 -batch -noheader -csv /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db \
"SELECT run_accession FROM sample_annot LEFT JOIN sample2bioinformatician \
USING(patient_code) WHERE username='rfeeney';" > /shared/projects/2314_medbioinfo/rachel/MedBioinfo/analyses/rfeeney_run_accession.txt

# Explore sra-tool with one run accession ID and first 10 reads (-x 10)
## Load sra-tools
module load sra-tools

## Use fastq-dump to download reads from NCBI's SRA 
## fastq-dump command format: fastq-dump [options] [accessions]

## Split the reads into three files (forward, reverse, singletons if any): --split-3 
## Produce gzip-ed compressed files: --gzip
## Define the definition lines as 'accession.spot.readid': --readids
## Write output files in the dedicated directory: --outdir /shared/projects/2314_medbioinfo/rachel/MedBioinfo/data/sra_fastq
## Disable multi-threading: --disable-multithreading

## Test: Download first 10 reads for accession ID "ERR6913341" using fastq-dump command with all above options included
fastq-dump --split-3 --readids --disable-multithreading -X 10 \
--gzip  --outdir /shared/projects/2314_medbioinfo/rachel/MedBioinfo/data/sra_fastq ERR6913341

## Check how the downloaded files have been named
ls -l ../data/sra_fastq/
zcat ../data/sra_fastq/ERR6913341_* 

# Download all 8 sequencing run IDs (using xargs)
cat rfeeney_run_accessions.txt | srun --cpus-per-task=2 --time=00:30:00 xargs \
fastq-dump --split-3 --readids --disable-multithreading \
--gzip --outdir ../data/sra_fastq \

## Check job usage 
sacct --format=JobID,JobName%20,ReqCPUS,ReqMem,Timelimit,State,ExitCode,Start,elapsed,MaxRSS,NodeList,Account%15 

## Check number of fastq files
ls ../data/sra_fastq | wc -l

## Count number of reads in each FASTQ file
zgrep -c "^@" ../data/sra_fastq/*fastq.gz 

## Base call quality scores 
### Found in field 4 (4th line) 
### ASCII format ranging from ! (lowest quality) to ~ (highest quality) 

# Compare number of reads obtained with grep to number of reads obtained by seqkit 
## Load seqkit
module load seqkit 

## Print statistics for each FASTQ file using seqkit 
seqkit stats --threads 2 ../data/sra_fastq/*.fastq.gz

## Compare the read count from seqkit stats with the author metadata read count
### Unsure how to do this using sqlite3 command, also unsure of location of author metadata

## Are reads un-trimmed or already quality filtered/trimmed?
# Use seqkit command to check if sequencing kit adapters have been trimmed from FASTQ files
## Search for adapters - not sure how to do this
zcat ../data/sra_fastq/*.fastq.gz | seqkit locate -i -p "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" --threads 2
zcat ../data/sra_fastq/*.fastq.gz | seqkit locate -i -p "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" --threads 2

# Quality control
## Make subdirectory for fastQC
 mkdir ./fastqc

## Load fastQC module
module load fastQC

## Test fastQC command on one pair of files (accession ERR6913341)
srun --cpus-per-task=2 --time=00:10:00 fastqc --threads 2 \
--outdir ./fastqc --noextract ../data/sra_fastq/ERR6913341_1.fastq.gz ../data/sra_fastq/ERR6913341_2.fastq.gz  

## Run fastQC command for all file pairs
srun --cpus-per-task=2 --time=00:30:00 xargs -a rfeeney_run_accessions.txt -I{} fastqc --threads 2 --outdir \ 
./fastqc --noextract ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz

# html files can be moved to local laptop using the following command (in normal console window, not while connected to server)
## scp rfeeney@core.cluster.france-bioinformatique.fr:/shared/projects/2314_medbioinfo/rachel/MedBioinfo/analyses/fastqc/*.html  /rachel/medbioinfo/fastqc_html/
## Have reads been trimmed to exclude bases with low quality scores? Yes, see high quality score for all bases in fastqc report
## Have reads been trimmed to exclude sequencing library adapters? Yes, no adapters detected in fastqc report

# Merge paired end reads
## Load flash2 module
module load flash2

## Create sub-directory for merged files 
mkdir ../data/merged_pairs

## Test merging two paired end FASTQ files (accession ERR6913341) 
srun --cpus-per-task=2 --time=00:10:00 flash2 --threads=2 \
../data/sra_fastq/ERR6913341_1.fastq.gz ../data/sra_fastq/ERR6913341_2.fastq.gz \ 
--compress --output-directory=../data/merged_pairs/ --output-prefix="ERR6913341.flash" \
2>&1 | tee -a rfeeney_flash2.log

## Check statistics with seqkit
seqkit stats --threads 2 ../data/merged_pairs/ERR6913341.flash.extendedFrags.fastq.gz
### 84.61% of reads combined successfully 

## View flash2 histogram file
head ../data/merged_pairs/ERR6913341.flash.histogram
### Suggests most library insert sizes are 39bp in length?

# Merge all paired end reads 
srun --cpus-per-task=2 --time=00:30:00 xargs -a rfeeney_run_accessions.txt -n 1 -I{} \
flash2 --threads=2 ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz \
--compress --output-directory=../data/merged_pairs/ --output-prefix={}.flash \
2>&1 | tee -a rfeeney_flash2.log

# Compare number of base pairs in unmerged reads vs remaining after merging
seqkit stats --threads 2 ../data/merged_pairs/*.flash.extendedFrags.fastq.gz
## Has information been lost or was it redundant? 
### Generally see a small decrease in sum_len in merged vs unmerged reads, but sum_len increases in some merged reads
### Dependent on overlap size? Bases lost are likely redundant and no information lost? 

# Use read mapping to check for PhiX contamination 
## Create subdirectory for PhiX genome
mkdir ../data/reference_seqs 

## Use NCBI edirect tool kit to download PhiX genome
efetch -db nuccore -id NC_001422 -format fasta > ../data/reference_seqs/PhiX_NC_001422.fnaefetch \
 -db nuccore -id NC_001422 -format fasta > ../data/reference_seqs/PhiX_NC_001422.fna

## Check PhiX file
head ../data/reference_seqs/PhiX_NC_001422.fna 

## Load bowtie2 module
module load bowtie2

## Create bowtie2 indexed database from PhiX reference seqs
mkdir ../data/bowtie2_db
srun bowtie2-build -f ../data/reference_seqs/PhiX_NC_001422.fna ../data/bowtie2_db/PhiX_bowtie2_db

## Check output index files
ls -l ../data/bowtie2_db/

## Create new subdirectory in analyses directory to collect alignment results 
mkdir bowtie

## Align merged reads against PhiX index database 
srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_db/PhiX_bowtie2_db -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz \
 -S bowtie/rachel_merged2PhiX.sam --threads 8 --no-unal 2>&1 | tee bowtie/rachel_bowtie_merged2PhiX.log 
## 0.00% alignment rate, no hits against PhiX

## Check for SARS-CoV-2 hits 
### Download SARS-CoV-2 genome from NCBI 
efetch -db nuccore -id NC_045512 -format fasta > ../data/reference_seqs/SC2_NC_045512.fnaefetch \
 -db nuccore -id NC_045512 -format fasta > ../data/reference_seqs/SC2_NC_045512.fna 

### Make reference database for SC2
srun bowtie2-build -f ../data/reference_seqs/SC2_NC_045512.fna ../data/bowtie2_db/SC2_bowtie2_db

### Align merged reads against SC2 index database
srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_db/SC2_bowtie2_db -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz \
 -S bowtie/rachel_merged2_SC2.sam --threads 8 --no-unal 2>&1 | tee bowtie/rachel_bowtie_merged2_SC2.log

### 125 reads aligned 1 time, 0.00% overall alignment rate

### View SAM file 
head bowtie/rachel_merged2_SC2.sam 

# Load multiQC module
module load mutliqc

# Create summary report with multiQC
srun multiqc --force --title "rfeeney sample sub-set" ../data/merged_pairs/ ./fastqc/ ./rfeeney_flash2.log ./bowtie/

## Save and view html on local laptop 
## scp rfeeney@core.cluster.france-bioinformatique.fr:/shared/projects/2314_medbioinfo/rachel/MedBioinfo/analyses/rfeeney-sample-sub-set_multiqc_report.html  /rachel/medbioinfo/multiqc_html/

echo "script end"

date
