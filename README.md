# IsoSeqAPA (version 0.2)
Time-stamp: <2018-02-06 Yunhao Wang, Email: yunhaowang@126.com>


## Introduction

Full-length isoform sequencing (Iso-Seq) technology originally developed by Pacific Biosciences (PacBio) has been widely applied to transcriptome study. IsoSeqAPA is a bioinformatics tool to construct full-length transcriptome (including the detection of known gene isoforms and the prediction of novel gene isoforms) and analyze APA (alternative cleavage and polyadenylation) event using PacBio Iso-Seq data. In addition, Second Generation Sequencing (SGS, e.g., Illumina platform) short reads (RNA-Seq data) can also be used to improve the construction of full-length transcriptome and the identification of polyA sites.


## Prerequisite

- Linux system

- python 2.7


## Install and Run

1. Download the package (e.g., `IsoSeqAPA-0.2.tar.gz`) to a directory (e.g., `/home/`)

2. Unpack it using the command `tar -zxvf /home/IsoSeqAPA-0.2.tar.gz`

3. Now, you can run IsoSeqAPA by the executable file `/home/IsoSeqAPA-0.2/bin/isoseqapa`. Optional, you can add IsoSeqAPA into your PATH so that you can run IsoSeqAPA without having to specify the entire path. For example, you can add one line `export PATH=/home/IsoSeqAPA-0.2/bin:$PATH` to your `~/.bashrc`


## Input

### 1. Reference genome file (FASTA format)

### 2. Gene annotation file (GTF format)

Note: considering gene duplication phenomenon, wherein multiple genic loci (from same/different chromosomes) can be transcribed into the sequence-identical transcript, thus make sure the transcript ID is unique for each derived genic loci. Suggest to use the well-annotated gene annotation libraries (e.g., RefSeq, Ensembl and Gencode).

### 3. PacBio Iso-Seq long read alignment file(s) (SAM format)

Multiple sam files can be as inputs split by "space".

Please follow the processes below to generate the Iso-Seq long read alignment files. (Here, we take SIRV (Spike-In RNA Variant by Lexogen Inc.) Iso-Seq data as an example.)

#### Raw Iso-Seq dataset

After sequencing SIRV sample by PacBio RSII platform, we get a batch of `*.bax.h5` files. Now, the BAM format (typically, `*.subreads.bam`) is the standard output of PacBio sequencer (e.g., PacBio Sequel system).

#### Step1: extract ROI

Extract ROI (Reads of Insert, also historically called CCS) by PacBio SMRT® Analysis Software or SMRT® Link Software.

In our test, we use SMRT Analysis (v2.3.0): `ConsensusTools.sh CircularConsensus --minFullPasses 0 --minPredictedAccuracy 70` to extract ROI. Now, we get the `SIRV_ROI.fasta` file. 

Alternatively, if your raw Iso-Seq data is `*.subreads.bam` file, you can use SMRT Link (v5.0.0): `ccs --minPasses 0 --minPredictedAccuracy 0.7` to get `*.ROI.bam` file.

#### Step2: classifiy ROI

Classifiy ROI by PacBio SMRT® Analysis Software or SMRT® Link Software.

In our test, we use SMRT Analysis (v2.3.0): `pbtranscript.py classify --flnc SIRV_ROI.flnc.fasta --nfl SIRV_ROI.nfl.fasta -d ./output/ -p SIRV_ROI.primer_info.csv --detect_chimera_nfl SIRV_ROI.fasta SIRV_ROI.classify.fasta`. Now, we get the FLNC (full-length non-chimera) ROI (i.e., `SIRV_ROI.flnc.fasta`) and nFLNC (non-full-length non-chimera) ROI (i.e., `./output/nflnc.fasta`)

Alternatively, if your data is `*.ROI.bam` format, you can use SMRT Link (v5.0.0): `pbtranscript.py classify --flnc *.flnc.fasta --nfl *.nfl.fasta -d ./output/ -p *.primer_info.csv --detect_chimera_nfl *.ROI.bam *.ROI.classify.fasta`. Similarly, you can get the FLNC and nFLNC.

#### Step3: separate nFLNC ROI

Separate the nFLNC ROI by the script `./IsoSeqAPA-0.2/utilities/py_isoseqcon_separate_nflnc_fasta.py -i ./output/nflnc.fasta -s ./output/SIRV_ROI_nflnc_sense.fasta -u ./output/SIRV_ROI_nflnc_unknown.fasta`. Now we have 3 fasta files: 1) `SIRV_ROI.flnc.fasta` (with strand information); 2) `./output/SIRV_ROI_nflnc_sense.fasta` (with strand information); and 3) `./output/SIRV_ROI_nflnc_unknown.fasta` (without strand information).

#### Step4: align ROI

Align 3 fasta files to reference genome by GMAP aligner (version 2016-06-09). For `SIRV_ROI.flnc.fasta` and `./output/SIRV_ROI_nflnc_sense.fasta`, we used the parameter `-z sense_force -f samse -n 0 --split-output ./SIRV_ROI_sense`. For `./output/SIRV_ROI_nflnc_unknown.fasta`, we used the parameter `-f samse -n 0 --split-output ./SIRV_ROI_unknown`. Now, we get 8 output files (SAM format) but only take 4 of them (`SIRV_ROI_sense.uniq`, `SIRV_ROI_sense.mult`, `SIRV_ROI_unknown.uniq` and `SIRV_ROI_unknown.mult` as the input of isoseqapa.

Alternatively, you can use other aligners (e.g., BLAT) to align your Iso-Seq long read to reference genome. However, some points need to be noted: (1) Iso-Seq long read sequences in `FLNC.fasta` and `nFLNC_sense.fasta` files are sense strand, and in `nFLNC_unknown.fasta` file are strand-unknown; (2) At most one alignment path is output for each Iso-Seq long read; (3) The aligned strand should be marked using the Tag:Type:Value = 'XS:A:+/-/?' in the optional fields (>=12 colunm) of SAM output file.

We strongly suggest you to use GMAP which is best aligner for Iso-Seq long read transcriptome data (see the paper "Evaluation of tools for long read RNA-seq splice-aware alignment. bioRxiv (2017)").

### 4. PacBio Iso-Seq long read primer information file (CSV format)

In the Step2 of generating the PacBio Iso-Seq long read alignment file (i.e., classifiy ROI by PacBio SMRT® Analysis Software or SMRT® Link Software), you can output the primer information (CSV format) using the parameter `-p` (e.g., we use SMRT Analysis (v2.3.0): `pbtranscript.py classify --flnc SIRV_ROI.flnc.fasta --nfl SIRV_ROI.nfl.fasta -d ./output/ -p SIRV_ROI.primer_info.csv --detect_chimera_nfl SIRV_ROI.fasta SIRV_ROI.classify.fasta`). The output file (`SIRV_ROI.primer_info.csv`) contains the primer information for all ROIs.

### 5. Second-Generation Sequencing (e.g., Illumina platform) short read alignment file(s) (SAM format).

Multiple sam files can be as inputs split by "space".

In our test, we used Hisat2 to get the `SIRV_SR.sam` file.

Alternatively, a batch of short read aligners (e.g., Tophat, Hisat2 and STAR) can be used to generate this SAM file. Note: In the CIGAR field of SAM file, only the strings 'M', 'I', 'D', 'N', 'S', and 'H' are permissible. In the optional fields (>=12 colunm) of SAM file, there is the 'Tag:Type:Value' = 'XS:A:+/-/?' showing the aligned strand.

### 6. 3End-Seq data (BED format)

3End-Seq is the abbreviation for technologies that combine polyA-enrichment method and short read sequencing technology to identify polyA sites. For easy understanding, please check the paper "Fu, Y., Sun, Y., Li, Y., Li, J., Rao, X., Chen, C., and Xu, A. (2011). Differential genome-wide profiling of tandem 3' UTRs among human breast cancer and normal cells by high-throughput sequencing. Genome Res. 21, 741-747.". For human and mouse, you can download the annotated polyA sites from some public database (such as, APADB, APASdb), then convert to BED format. Alternatively, you can generate your own 3End-Seq peak data. 3End-Seq peak data (BED format) should include 6 columns: 1) chromosome ID; 2) start position; 3) end position; 4) polyA ID; 5) score; 6) strand. Actually, only 1,2,3,6 columns are used in IsoSeqAPA. 


## Output

### 1. Constructed gene isoforms with APA analysis (modified GPD format)

(1) Gene ID

For novel singleton isoform (only including one exon) located in novel genic loci, the prefix of gene ID is "novel_sgt_loci_"

For novel multi-exon isoform (including multiple exons) located in novel genic loci, the prefix of gene ID is "novel_mlt_loci_"

(2) Isoform ID

For novel singleton isoform (only including one exon), the prefix of isoform ID is "novel_sgt_iso_"

For novel multi-exon isoform (including multiple exons), the prefix of isoform ID is "novel_mlt_iso_"

(3) Chromosome

(4) Strand

(5) Transcription start site (for "+" strand transcript)

(6) Transcription end site (for "+" strand transcript)

(7) Number of supporting Iso-Seq full-length long reads

(8) Number of supporting Iso-Seq long reads (both full-length and non-full-length)

(9) Exon number

(10) Exon start site set (for "+" strand transcript)

(11) Exon end site set (for "+" strand transcript)

(12) For novel isoform, its derived genic locus; The "-" means this novel isoform is derived novel genic locus

(13) For novel isoform, the overlap percentage between the novel isoform and its derived genic locus; The "0.0" means this novel isoform is derived novel genic locus

(14) For novel singleton isoform, if it is located at the last exon of any known isoform? If yes, show isoform ID otherwise '-'

(15) For novel singleton isoform, the overlap percentage bewteen this novel singleton isoform and the the last exon of known isoform

(16) For novel multi-exon isoform, number of its splice sites is annotated by gene annotation library and/or detected by short reads; and the total number of its splice sites. Split by '/'

(17) For novel multi-exon isoform, if the multi-exon isoform is the subset (based on splice junction combination) of any known multi-exon isoform? If yes, show isoform ID otherwise '-'

(18) For novel isoform, length of longest polyA track in defined region surrounding transcription end site

(19) For novel isoform, percentage of nucleotide 'A' in defined region surrounding transcription end site

(20) Individual polyA site information set (split by ','). The set is split by '_': polyA site postion; number of supporting long reads; number of supporting short reads

(21) Grouped polyA site information set (split by ','). The set is split by '_': polyA site postion; number of supporting long reads; number of supporting short reads

(22) Grouped polyA site information set after filtering (based on supporting long reads and short reads, split by ',').  The set is split by '_': polyA site postion; number of supporting long reads; number of supporting short reads

(23) Number of polyA sites after filtering (based on supporting long reads and short reads)

### 2. Constructed gene isoforms with APA analysis (GTF format)

In the attribute field (9th column if split by 'tab'), there are some tag-value pairs corresponding to the columns of output GPD file above:

(1) full_length_LR_count

7th column of GPD file

(2) LR_count

8th column of GPD file

(3) derived_genic_loci

12th column of GPD file

(4) derived_genic_loci_overlap_pct

13th column of GPD file

(5) sgt_last_exon_iso

14th column of GPD file

(6) sgt_last_exon_overlap_pct

15th column of GPD file

(7) mlt_splice_site

16th column of GPD file

(8) mlt_subset_iso

17th column of GPD file

(9) polya_len

18th column of GPD file

(10) polya_pct

19th column of GPD file

(11) grouped_polya_site

22th column of GPD file

(12) polya_site_count

23th column of GPD file


## Usage and Example

### 1. Have Second-Generation Sequencing (SGS) short read alignment data

`./bin/isoseqapa -g ./example/input/SIRV_reference_genome.fasta -a ./example/input/SIRV_gene_annotation.gtf -l ./example/input/SIRV_ROI_sense.uniq ./example/input/SIRV_ROI_sense.mult ./example/input/SIRV_ROI_unknown.uniq ./example/input/SIRV_ROI_unknown.mult -c ./example/input/SIRV_ROI.primer_info.csv -s ./example/input/SIRV_SR.sam --end_bed ./example/input/SIRV_3End.bed --tempdir ./example/temp_with_SR/ --output_gpd ./example/SIRV_test_with_SR.gpd --output_gtf ./example/SIRV_test_with_SR.gtf`

### 2. No SGS short read alignment data

`./bin/isoseqapa -g ./example/input/SIRV_reference_genome.fasta -a ./example/input/SIRV_gene_annotation.gtf -l ./example/input/SIRV_ROI_sense.uniq ./example/input/SIRV_ROI_sense.mult ./example/input/SIRV_ROI_unknown.uniq ./example/input/SIRV_ROI_unknown.mult -c ./example/input/SIRV_ROI.primer_info.csv --end_bed ./example/input/SIRV_3End.bed --tempdir ./example/temp_without_SR/ --output_gpd ./example/SIRV_test_without_SR.gpd --output_gtf ./example/SIRV_test_without_SR.gtf`
