# bulkRNAseq_pipeline_scripts

# SCRIPTS DESCRIPTION

---

## Pre-Alignment Quality Control (`getPreAlignmentQC.sh`)

This script performs pre-alignment quality control (QC) on raw sequencing data using **FastQC** and **MultiQC**. It processes `.fq.gz`, `.fastq.gz`, and `.fasta` files in parallel, generating individual FastQC reports for each file and a consolidated MultiQC report.

### Problems 

`getPostAlignment.QC`
- bam index must be built before the ReSEQc scripts run 
- the bed12 annotation I created might be bed, need a 12-column check and check against pre-built annotation
- my ribosomal intervals file may be smaller than the pre-built one, inspect differences 
- Picard doesn\`t build ribosomal intervals plot
- old Picard synthax
- R RSEQc scripts do not initialize properly
- Need to think of some additional tools

### Usage

```bash
./getPreAlignmentQC.sh
```

#### Options

- **FASTQ_DIR**: Directory containing FASTQ/Fasta files (default: current directory `.`).
- **FASTQC_OUTPUT_DIR**: Directory where FastQC reports will be saved (default: `./pre_align_QC/fastQC`).
- **MULTIQC_OUTPUT_DIR**: Directory where the final MultiQC report will be saved (default: `./pre_align_QC/multiQC`).
- **FASTQ_LOGS**: Directory to store logs (default: `./pre_align_QC/logs`).
- **JOB_NUMBER**: Number of parallel jobs to run (default: 6).

### What the Script Does

1. **Searches for FASTQ/Fasta files**: Finds all `.fq.gz`, `.fastq.gz`, and `.fasta` files in the specified directory.
2. **Runs FastQC**: Executes FastQC on each file in parallel, logging the results for each file.
3. **Aggregates Results with MultiQC**: Runs MultiQC to create a consolidated report from all the FastQC outputs.
4. **Logs the Process**: Keeps logs of file discovery and FastQC processing in the `logs` directory.

### Script Caveats

- **Parallel processing**: The script runs up to 6 FastQC processes in parallel (customizable via the `JOB_NUMBER` variable). Adjust this number based on system resources.
- **File types**: Only processes `.fq.gz`, `.fastq.gz`, and `.fasta` files. Other file formats need to be converted or manually included in future script updates.
- **Output paths**: Ensure the output directories have the appropriate write permissions. The script creates the directories if they don’t exist, but permission issues will cause failures.
- **Dependencies**: Requires **FastQC**, **MultiQC**, and **GNU parallel** to be installed and available in the system's `$PATH`.

---

## `runSTARalign.sh`

### Usage

```bash
./runSTARalign.sh
```

This script aligns paired-end RNA-seq FASTQ files to a reference genome using STAR, and outputs sorted BAM files along with logs.

### What Does the Script Do?

- **STAR Alignment**: The script uses the STAR aligner to align paired-end reads (R1 and R2 FASTQ files) to the specified genome directory. It processes files in parallel using GNU `parallel`.
- **Input/Output**:
  - **Input**: Paired-end FASTQ files located in the specified directory (`BATCH_1_DIR`).
  - **Output**: Sorted BAM files saved in the `STAR_alignment/bam` directory and logs in `STAR_alignment/logs`.
- **Resources**:
  - Utilizes multi-threading to speed up alignment processes.
  - Uses memory management to ensure there is always sufficient memory left for new processes.

### Caveats

1. **File Naming**: The script assumes paired-end reads follow a naming convention where the forward read is denoted with `_R1` or `_1`, and the reverse read with `_R2` or `_2`. Any deviations from this pattern may cause the script to skip files.
2. **Memory Management**: The script uses the `--memfree` option of `parallel` to manage memory allocation. Adjust this value (`GB_FREE`) based on available resources to prevent excessive memory usage.
3. **Log Output**: Each sample generates a log file in the `STAR_alignment/logs` directory. Errors encountered during alignment are logged here but also shown in the terminal for quick diagnosis.
4. **Genome Index**: Ensure that the `GENOME_DIR` variable points to a valid STAR genome index. Incorrect or missing paths will cause the alignment step to fail.
5. **Paired-end Only**: This script is designed specifically for paired-end reads. Single-end FASTQ files are not supported.

---

## `getPostAlignmentQC.sh`

### Usage

```bash
./getPostAlignmentQC.sh
```

This script performs post-alignment quality control (QC) on BAM files using Picard, RSeQC, and Samtools. It aggregates all QC metrics into a final MultiQC report.

### What Does the Script Do?

1. **Picard Tools**: Runs `CollectRnaSeqMetrics` to generate RNA-seq metrics and `MarkDuplicates` to identify and quantify duplicates without removing them.
2. **RSeQC**: Runs `infer_experiment.py` to determine strand-specificity using the BAM files.
3. **Samtools**: Runs `flagstat` to generate alignment statistics.
4. **MultiQC**: Combines the metrics from all tools into a single report, including metrics from STAR alignments located in `*_STARpass1` subdirectories.

**N.B!** The script does not inherently remove duplicates. If you need to remove duplicates, use `dedupPicard_postQC.sh`, or modify the script by adding `--REMOVE_DUPLICATES false` as specified on the [Picard manual](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard). This is done intentionally, as the main purpose of the pipeline is bulk RNAseq analysis, and without UMIs it\`s basically impossible to remove duplicates without messing up with the biology of the sample. For reference see for example [here](https://www.biorxiv.org/content/10.1101/2023.12.12.571280v1) 

### Caveats

1. **File Naming**: The script assumes all BAM files are in the specified `BAM_DIR`. Incorrectly named files or missing BAM files will result in skipped QC steps.
2. **Reference Files**: Ensure that paths to the reference genome, refFlat file, ribosomal intervals, and BED annotation are correct. Incorrect paths will lead to failures in Picard and RSeQC steps.
3. **Memory & CPU Usage**: The script runs Picard, RSeQC, and Samtools in parallel with a specified number of jobs (`JOB_NUMBER`). Ensure sufficient computational resources are available, especially for large datasets.
4. **STAR Metrics**: The script looks for STAR metrics in `*_STARpass1` subdirectories. Ensure these directories are present and correctly named; otherwise, STAR metrics won’t be included in the final MultiQC report.

---

## `dedupPicard_postQC.sh`

### Usage

```bash
./dedupPicard_postQC.sh
```

This script processes BAM files by removing duplicates using Picard's `MarkDuplicates`, runs RNA-seq QC with `CollectRnaSeqMetrics`, and generates a final report using MultiQC.

### What Does the Script Do?

1. **Picard MarkDuplicates**: Removes duplicates from BAM files, generating deduplicated BAM files and duplicate metrics.
2. **Picard CollectRnaSeqMetrics**: Runs RNA-seq QC on the deduplicated BAM files, producing various RNA-seq metrics.
3. **MultiQC**: Aggregates all the QC metrics into a unified report for easier interpretation.

### Caveats

1. **File Paths**: Ensure the paths to the reference genome and `refFlat` files are correctly specified in the script. Incorrect paths will cause QC steps to fail.
2. **Duplicate Removal**: By default, duplicates are permanently removed from the BAM files. If you need the original BAM files, make sure they are backed up before running the script.
3. **Parallelization**: The script uses GNU `parallel` to process multiple BAM files concurrently. The number of jobs is set by `NUM_JOBS` (default is 2). Adjust this value according to your system’s resources.
4. **Memory Usage**: Removing duplicates and running RNA-seq metrics can be memory-intensive. Ensure your system has enough resources, especially when working with large BAM files.

---

Here’s the README for `GTFtoRefFlat.sh`:

## `GTFtoRefFlat.sh`

### Usage

```bash
./GTFtoRefFlat.sh -i input.gtf -o output.refFlat
```

This script converts a GTF file into the RefFlat format, commonly required for certain RNA-seq QC tools like Picard.

### What Does the Script Do?

- **Input**: Takes a GTF file (using the `-i` flag) that contains gene annotation data.
- **Output**: Generates a RefFlat file (specified with the `-o` flag) by extracting gene and transcript information from the GTF file.
- The RefFlat format is essential for RNA-seq tools, as it defines gene and transcript details such as chromosomal positions and exon boundaries.

### Caveats

1. **GTF Format**: The script assumes a standard GTF format. Ensure that your GTF file contains the necessary fields (`gene_id`, `transcript_id`, `exon`) in the expected format.
2. **Transcript Names**: If no `transcript_id` is present for a record in the GTF file, that entry will be skipped in the output.
3. **Strand Information**: The script handles strand-specific data (`+` or `-`). However, if this field is not available or improperly formatted, the output might be inaccurate.
4. **Dependencies**: Ensure that `awk` is available in your environment, as this script relies heavily on it for text processing.

---

## `getRibosomalIntervals_from_gtf.sh`

### Usage

```bash
./getRibosomalIntervals_from_gtf.sh
```

This script generates a ribosomal RNA (rRNA) interval list from a GTF file, suitable for use with tools like Picard, which require interval lists for certain QC tasks.

### What Does the Script Do?

1. **Input**: It takes a GTF annotation file that contains ribosomal RNA information and a reference genome in FASTA format.
2. **Output**: It creates a Picard-compatible interval list (`ribosomal.interval_list`) that includes rRNA regions extracted from the GTF file.
3. **FASTA Index Generation**: If a `.fai` index file for the reference genome doesn't exist, the script automatically generates it using `samtools faidx`.

### Steps Performed by the Script

1. **FASTA Index Check**:
   - The script first checks if the `.fai` file (FASTA index) exists for the reference genome. If not, it creates the index using `samtools faidx`.
2. **Header Creation**:
   - The script adds a Picard-compatible header with chromosome size information extracted from the `.fai` file.
3. **rRNA Interval Extraction**:
   - The script searches the GTF file for `rRNA` features (i.e., annotations related to ribosomal RNA) and extracts their chromosome, start, end positions, and strand information.
   - The gene name is also extracted, with default fallback to `"rRNA"` if not present.
4. **Output**:
   - The extracted information is formatted and written to the output interval list in Picard's format.

### Expected Input Files

- **GTF File**: The GTF file should contain annotations for ribosomal RNA, identified by `rRNA` in the feature fields or gene names.
- **FASTA File**: The reference genome in FASTA format, along with its `.fai` index. The script will generate the `.fai` file if it’s missing.

### Output

The output is a Picard-compatible interval list file named `ribosomal.interval_list`.

### Dependencies

- **Samtools**: Required to generate the `.fai` index if it doesn't exist.
- **Grep and Awk**: Used for text processing and extracting the required fields from the GTF and `.fai` files.

### Example

```bash
./getRibosomalIntervals_from_gtf.sh
```

After running, the script will generate a file named `ribosomal.interval_list` containing ribosomal RNA intervals.

---

## `BedtoRefSeqBed.sh`

### Purpose

This script converts chromosome names in a BED file from the standard `chr1`, `chr2`, etc., format to RefSeq accession numbers (e.g., `NC_000067.7` for `chr1`), which are used in reference genomes for certain bioinformatics tools.

### Usage

```bash
./BedtoRefSeqBed.sh -I <input.bed> -O <output.bed>
```

- `-I`: Input BED file with standard chromosome names.
- `-O`: Output BED file with RefSeq accession numbers.

### Description

The script uses `sed` to search for standard chromosome names (e.g., `chr1`, `chrX`) and replaces them with the corresponding RefSeq accession numbers for the mouse genome (GRCm39). It then writes the modified content to the specified output file.

### Example

```bash
./BedtoRefSeqBed.sh -I input.bed -O output_refseq.bed
```

This will create a new BED file (`output_refseq.bed`) where the chromosome names are replaced with their corresponding RefSeq accession numbers.

### Dependencies

- **Sed**: This tool is required and is commonly available on Unix-based systems.

### RefSeq Mapping

The chromosome names are mapped as follows:

- `chr1` → `NC_000067.7`
- `chr2` → `NC_000068.8`
- `chr3` → `NC_000069.7`
- `chr4` → `NC_000070.7`
- `chr5` → `NC_000071.7`
- `chr6` → `NC_000072.7`
- `chr7` → `NC_000073.7`
- `chr8` → `NC_000074.7`
- `chr9` → `NC_000075.7`
- `chr10` → `NC_000076.7`
- `chr11` → `NC_000077.7`
- `chr12` → `NC_000078.7`
- `chr13` → `NC_000079.7`
- `chr14` → `NC_000080.7`
- `chr15` → `NC_000081.7`
- `chr16` → `NC_000082.7`
- `chr17` → `NC_000083.7`
- `chr18` → `NC_000084.7`
- `chr19` → `NC_000085.7`
- `chrX` → `NC_000086.8`
- `chrY` → `NC_000087.8`

### Output

The script generates a new BED file (`output.bed`) where the chromosome names are updated to RefSeq accession numbers. This file is suitable for workflows that require RefSeq format.

---

## `htseqCheckStrand.sh`

### Purpose

This script samples a specified number of reads from a BAM file and uses `htseq-count` to determine the most likely strandedness (i.e., whether the RNA-seq data is stranded or not) based on the provided GTF file.

### Usage

```bash
./htseqCheckStrand.sh <BAM file> <GTF file> <Number of reads to sample>
```

- `<BAM file>`: Path to the BAM file to sample and analyze.
- `<GTF file>`: Path to the GTF annotation file used by `htseq-count`.
- `<Number of reads to sample>`: Number of reads to sample from the BAM file for strandedness inference.

### Description

1. **Check for Required Tools**: Ensures that `htseq-count` and `samtools` are installed.
2. **Verify BAM File Header**: Checks if the BAM file contains `@SQ` lines in the header, which are necessary for proper processing.
3. **Sample Reads**: Creates a temporary BAM file with a specified number of reads sampled from the beginning of the original BAM file.
4. **Run `htseq-count`**: Executes `htseq-count` with different strandedness options (`-s no`, `-s yes`, `-s reverse`) to count reads.
5. **Determine Strandedness**: Compares read counts for each strandedness option and infers the most likely strandedness based on the highest read count.
6. **Output Results**: Stores the results in a text file named `<BAM file name>_strand_results.txt`.

### Example

```bash
./htseqCheckStrand.sh sample.bam annotations.gtf 1000
```

This will sample 1000 reads from `sample.bam`, run `htseq-count` for each strandedness option, and save the results to `sample_strand_results.txt`.

### Dependencies

- **htseq-count**: To count reads and infer strandedness.
- **samtools**: To handle BAM file operations.

### Notes

- Ensure the GTF file is compatible with the BAM file (i.e., it should be in the same genome assembly version).
- The sampling of reads is done from the beginning of the BAM file. If the file is large and/or the reads are not uniformly distributed, consider using a random sampling approach for better representation.

### Output

- **`<BAM file name>_strand_results.txt`**: Contains the results of the strandedness inference, including total read counts for each strandedness option and the inferred strandedness.

Feel free to modify or extend this script based on specific requirements or additional analyses.

---

## `salmon_quant.sh`

### Purpose

This script performs quantification of RNA-seq data using Salmon. It processes FASTQ files to estimate transcript abundance and generates a MultiQC report summarizing the quantification results.

### Usage

```bash
./salmon_quant.sh
```

### Description

1. **Directories and Constants**:

   - **FASTQ_DIR**: Directory where FASTQ files are located (default is the current directory).
   - **SALMON_OUTPUT_DIR**: Base directory for Salmon output.
   - **INDEX_DIR**: Directory where the Salmon index is stored or will be built.
   - **EXPERIMENT_OUTPUT_DIR**: Directory where Salmon results for different experiment types will be stored.
   - **LOG_DIR**: Directory for storing logs.
   - **TRANSCRIPTOME_FASTA**: Path to the transcriptome FASTA file used for building the Salmon index.

2. **Index Building**:

   - If the Salmon index does not exist, it is built from the provided transcriptome FASTA file.

3. **Quantification**:

   - **Function `run_salmon`**: This function runs Salmon quantification on each FASTQ file found in the specified directory. It infers strandedness automatically (`-l A` option) and stores results in a dedicated output directory for each sample.

4. **Processing FASTQ Files**:

   - Finds all FASTQ files in the `FASTQ_DIR` and processes each file using Salmon.

5. **MultiQC Report**:
   - Runs MultiQC to aggregate and visualize the Salmon quantification results.

### Arguments

This script does not require command-line arguments; it uses predefined constants for configuration.

### Example

To use this script, simply run it in the directory where your FASTQ files are located, or adjust the `FASTQ_DIR` variable in the script accordingly:

```bash
./salmon_quant.sh
```

### Dependencies

- **Salmon**: For RNA-seq quantification.
- **MultiQC**: For aggregating and visualizing results.
- **GNU Parallel**: (Optional) Can be used to run jobs in parallel if desired.

### Notes

- Ensure that the `TRANSCRIPTOME_FASTA` file matches the genome assembly used for alignment.
- The `-l A` option in Salmon allows automatic inference of the library type. If you know the specific strandedness of your samples, you can adjust this option accordingly.
- The `MultiQC` command assumes it is installed and available in your PATH.

### Output

- **Salmon Quantification Results**: Stored in `EXPERIMENT_OUTPUT_DIR`, with separate directories for each sample.
- **Logs**: Stored in `LOG_DIR`, with log files named according to the sample.
- **MultiQC Report**: Aggregates all Salmon quantification results and is stored in the `SALMON_OUTPUT_DIR`.

---

## `runQuant3p.sh`

### Purpose

This script performs RNA-seq quantification using the `quant3p` tool on BAM files. It processes BAM files to generate count matrices and logs the results.

### Usage

```bash
./runQuant3p.sh
```

### Description

1. **Directories and Constants**:

   - **JOB_NUMBER**: Number of parallel jobs for processing (default is 8).
   - **BAM_DIR**: Directory containing BAM files (default is the current directory).
   - **GTF_FILE**: Path to the GTF file for gene annotations.
   - **OUTPUT_DIR**: Directory where quant3p results will be stored.
   - **LOG_DIR**: Directory where logs will be stored.

2. **Creating Directories**:

   - The script creates the `OUTPUT_DIR` and `LOG_DIR` if they do not already exist.

3. **Quant3p Execution**:

   - **Function `run_quant3p`**: This function runs `quant3p` on each BAM file. It saves the results in `OUTPUT_DIR` and logs in `LOG_DIR`.
   - Uses `quant3p` to process each BAM file with the specified GTF file. The results and logs are saved with filenames based on the base name of each BAM file.

4. **Parallel Processing**:

   - Uses `GNU parallel` to run `quant3p` on multiple BAM files simultaneously. The number of parallel jobs is controlled by the `JOB_NUMBER` variable.

5. **Completion Logging**:
   - Logs a message indicating that the quant3p analysis has been completed.

### Arguments

This script does not take command-line arguments. It relies on predefined variables for configuration.

### Example

To use this script, run it in the directory where your BAM files are located, or adjust the `BAM_DIR` variable in the script:

```bash
./runQuant3p.sh
```

### Dependencies

- **quant3p**: For RNA-seq quantification.
- **GNU Parallel**: For running jobs in parallel.

### Output

- **Quant3p Results**: Each BAM file will generate a corresponding result file in the `OUTPUT_DIR` with the suffix `_quant3p_results.txt`.
- **Logs**: Each BAM file will generate a corresponding log file in the `LOG_DIR` with the suffix `_quant3p.log`.
- **Completion Log**: A log entry indicating the completion of the analysis will be added to `LOG_DIR/quant3p_analysis.log`.

### Notes

- Ensure that `quant3p` is installed and available in your PATH.
- Modify `JOB_NUMBER` as needed based on the number of cores available on your machine.
- Ensure the `GTF_FILE` path and `BAM_DIR` are set correctly according to your file system and experiment setup.

This script automates the process of running quant3p on multiple BAM files, facilitating efficient analysis of RNA-seq data.

---

## `get_chr_coverage.sh`

### Purpose

This script calculates the coverage of specific chromosomes (chrX and chrY, in RefSeq notation) from indexed BAM files and outputs the results to a text file.

### Usage

```bash
./get_chr_coverage.sh
```

### Description

1. **Directories and Files**:

   - **BAM_DIR**: Directory containing indexed BAM files (`.bam` files). Ensure that corresponding `.bai` files are in the same directory.
   - **OUTPUT_FILE**: File where chromosome coverage statistics will be saved.

2. **Initialize Output File**:

   - The script initializes the output file with a header line: `File\tChromosome\tCoverage`.

3. **Coverage Calculation**:

   - The script loops through each BAM file in the specified directory.
   - For each BAM file, it uses `samtools idxstats` to extract coverage statistics for chromosomes `NC_000086.8` (chrX) and `NC_000087.8` (chrY).
   - The coverage values are appended to the `OUTPUT_FILE` with the corresponding file name.

4. **Output**:
   - The results are stored in `chr_coverage.txt` in a tab-delimited format with columns for file name, chromosome, and coverage.

### Example

To use this script, simply execute it in the terminal. Make sure the `BAM_DIR` is correctly set to the directory containing your BAM files:

```bash
./get_chr_coverage.sh
```

### Dependencies

- **samtools**: For extracting coverage statistics from BAM files.

### Output

- **`chr_coverage.txt`**: Contains coverage statistics for chromosomes chrX and chrY for each BAM file in the specified directory. The file is formatted with tab-separated values:

  ```
  File\tChromosome\tCoverage
  example.bam\tNC_000086.8\t123456
  example.bam\tNC_000087.8\t654321
  ```

### Notes

- Ensure that `samtools` is installed and available in your PATH.
- Verify that the chromosome names `NC_000086.8` and `NC_000087.8` match the reference genome used for indexing the BAM files.
- Modify `BAM_DIR` and `OUTPUT_FILE` variables in the script if needed to fit your specific setup.

This script provides a quick way to assess coverage for specific chromosomes across multiple BAM files, useful for quality control and analysis in genomic studies.

---

Here's a README for `samtools_build_index.sh`:

## `samtools_build_index.sh`

### Purpose

This script creates index files (`.bai`) for BAM files located in a specified directory. Index files are essential for many downstream applications that require quick access to specific regions of BAM files.

### Usage

```bash
./samtools_build_index.sh
```

### Description

1. **Directories**:

   - **BAM_DIR**: Directory containing BAM files that need to be indexed.

2. **Indexing**:

   - The script loops through each BAM file in the specified directory.
   - For each BAM file, it runs `samtools index` to generate the corresponding index file (`.bai`).

3. **Output**:
   - Each BAM file will have an associated index file (`.bai`) created in the same directory.

### Example

To use this script, ensure that `BAM_DIR` is set to the directory containing your BAM files. Then, execute the script in the terminal:

```bash
./samtools_build_index.sh
```

### Dependencies

- **samtools**: Required for creating BAM file indexes.

### Notes

- Ensure that `samtools` is installed and available in your PATH.
- Verify that `BAM_DIR` correctly points to the directory with your BAM files.
- The script assumes that all BAM files in the directory need indexing. If you only need to index specific files, you may need to modify the script accordingly.

This script is useful for preparing BAM files for analysis by ensuring that each file has an index, which is necessary for efficient data access in many bioinformatics tools.

---

## Sorted featureCounts Output (`runFeatureCounts.sh`)
This script performs RNA-seq read quantification using **featureCounts** and processes the output to create a sorted count matrix. It handles BAM files, runs featureCounts, and then reorders the output columns lexicographically.

### Usage
```bash
./runFeatureCounts.sh
```

#### Options
- **BAM_DIR**: Directory containing BAM files (default: current directory `.`).
- **ANNOTATION_FILE**: Path to the GTF annotation file (default: `/mouse_genome/ref_genome/GCF_000001635.27_GRCm39_genomic.gtf`).
- **RAW_OUTPUT_DIR**: Directory where the raw featureCounts output will be saved (default: `./featureCounts/raw`).
- **MATRIX_OUTPUT_DIR**: Directory where the final sorted count matrix will be saved (default: `./featureCounts/count_matrices`).
- **LOG_DIR**: Directory to store logs (default: `./featureCounts/logs`).
- **THREADS**: Number of threads for featureCounts to use (default: 12).

### What the Script Does
1. **Locates and Sorts BAM Files**: Finds all `.bam` files in the specified directory and sorts them lexicographically.
2. **Runs featureCounts**: Executes featureCounts on the sorted BAM files, using the specified annotation file and parameters.
3. **Processes Output**: Uses `awk` to reorder the columns of the featureCounts output, creating a sorted count matrix. Awk logic double-checked against the script without awk sorting functionality - all correct.
4. **Generates Headers**: Creates a header row with "Geneid" as the first column, followed by simplified sample names derived from the BAM filenames.
5. **Logs the Process**: Keeps logs of the featureCounts run and the script's execution.

### Script Caveats
- **BAM file naming**: The script assumes BAM files are named descriptively (e.g., `SampleName_Aligned.sortedByCoord.out.bam`). The sample names in the output matrix are derived from these filenames.
- **Output format**: The final output is a tab-separated file with "Geneid" as the first column, followed by count data for each sample.
- **Column ordering**: Sample columns in the output are ordered lexicographically based on the BAM filenames.
- **Memory usage**: featureCounts can be memory-intensive. Ensure your system has sufficient RAM, especially for large BAM files or complex annotations.
- **Dependencies**: Requires **featureCounts** to be installed and available in the system's `$PATH`.
- **Annotation file**: Ensure the specified GTF annotation file exists and is appropriate for your organism and genome build.

### Output
- A raw featureCounts output file in `RAW_OUTPUT_DIR`.
- A sorted count matrix file in `MATRIX_OUTPUT_DIR`, with simplified column headers and lexicographically ordered sample columns.
- Log files in `LOG_DIR`, including featureCounts execution logs and script completion status.
