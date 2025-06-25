# Genomic Data Analysis Pipeline (Python Version)
Identification of somatic and germiline variants from tumor and normal samples
 Note: The bash script was added for comparism.

## Overview

This Python script is a complete recasting of a Linux-based genomic analysis pipeline. The original bash script involved downloading datasets, quality control, trimming, mapping, variant calling, merging, and annotation, all executed via command-line tools such as `wget`, `bwa`, `samtools`, and others.

The Python version aims to:

- Encapsulate these steps into functions for modular execution.
- Use Python libraries such as `subprocess` for executing shell commands.
- Use `pathlib` for file and directory management.
- Add informative print statements for progress tracking and debugging.
- Make the pipeline platform-independent; ideally, it can run on any system with the required tools installed.

---

## Why Python?

- **Portability:** Python can run across platforms, making automation scripts portable.
- **Readability:** Using functions and clear flow makes the code easier to understand and modify.
- **Control:** Better error handling, logging, and flexible control over process execution.
- **Integration:** Easier to incorporate into larger workflows or pipeline management systems.

---

## Dependencies and Environment Setup

The script relies heavily on external command-line tools. Here's what you need to prepare:

### Operating System
- **Linux (Ubuntu, Debian, CentOS, etc.)**: The script and commands are optimized for Linux environments.
- **Windows or Mac**: Might require adjustments or use of WSL (Windows Subsystem for Linux).

### Basic setup steps
1. **Install Python 3.x**:
    - Ubuntu: `sudo apt-get install python3`
2. **Install required command-line tools**:
    - Using package managers or conda:
      ```bash  
      # Example for Ubuntu:  
      sudo apt-get install bwa samtools bcftools unzip wget gunzip  
      # For conda:  
      conda create -n genomics python=3.9 fastqc trimmomatic bwa samtools bcftools  
      conda activate genomics  
Install Java Runtime Environment:
Needed for VarScan, snpEff, gemini.
Example:
bash
sudo apt-get install default-jre  
Download external tools manually if needed:
VarScan: VarScan
snpEff: snpEff
gemini: gemini
Environment variables:
Ensure commands like wget, gunzip, fastqc, bwa, samtools, etc., are available in your system's PATH.
How the Conversion Works
1. Using subprocess
Purpose: To execute shell commands exactly as they would run in a terminal.
Why: Most bioinformatics tools are CLI-based, so subprocess allows us to run them from Python.
Method: subprocess.run() is used for simplicity, with check=True to raise exceptions if commands fail.
2. Using pathlib
Purpose: To handle file and directory paths in an OS-independent way.
Why: More readable and robust compared to using string concatenation.
For example, Path("raw_data") creates a directory path object, and / operator joins paths seamlessly.
3. Using Print Statements
Why: To provide real-time feedback on script progress.
Helps with debugging—for example, knowing which step the script is executing or where it failed.
Also improves transparency during long-running processes.
4. Managing Downloads and Execution
For reproducibility, the script downloads all necessary datasets and references if they are not already present.
Checks if files exist before downloading to avoid redundancy.
Commands like wget, gunzip, are executed via subprocess.
5. Why not just run bash?
Pure bash scripts are less flexible and harder to debug programmatically.
Wrapping commands in Python allows integration with error handling, logging, and potential extensions.

## Step-by-Step Explanation of the Script


a. Directory Setup
Creates a structured workspace with raw_data, Fastqc_Reports, Mapping, and Variants.
Ensures directories are present before file operations, preventing runtime errors.
b. Downloading Data & References
Downloads raw FASTQ files and the reference genome from public repositories
c. Downloading and Preparing Reference Genome
Checks if the reference genome (hg19.chr5_12_17.fa.gz) exists.
If absent, downloads using wget.
Unzips the reference file with gunzip, preparing it for downstream processing.
The reference unzipped path is used consistently across tools like BWA, samtools, etc.
d. Quality Control with FastQC
Runs fastqc on each FASTQ file to generate quality reports.
Reports are saved in the Fastqc_Reports/ directory.
Providing feedback via print statements helps monitor progress and troubleshoot issues if quality issues are detected.
e. Trimming with Trimmomatic
Reads are processed to remove adapters and low-quality bases.
The outputs are new paired FASTQ files saved in the trimmed_reads/ directory.
Post-trimming quality is checked again with FastQC.
f. Mapping Reads
Uses bwa mem to align reads to the reference genome.
The read group (@RG) info is added for downstream analysis.
Output SAM files are saved in the Mapping/ directory.
g. Conversion, Sorting, and Indexing
Converts SAM to BAM (samtools view).
Sorts BAM files for efficient access (samtools sort).
Indexes BAM files (samtools index) to enable quick random access.
This step readies files for variant calling.
h. BAM Filtering & Processing
Filters BAM files to retain high-quality, properly paired reads (samtools view with -q, -f, -F).
Optional: Remove duplicates using samtools markdup or similar tools.
Uses additional tools like bamleftalign and bamtools for realignment and filtering—integrated via subprocess.
i. Variant Calling
Generates pileup files with samtools mpileup.
Calls somatic variants with VarScan (java -jar VarScan.v2.3.9.jar somatic).
Outputs VCF files for snp and indel variants.
j. Merging Variants
Compresses VCFs with bgzip.
Indexes with tabix.
Merges VCFs into a combined VCF with bcftools merge.
k. Annotation with snpEff
Downloads snpEff latest core database.
Annotates the merged VCF file to add functional annotations.
Output is a .ann.vcf file located in Variants/.
l. Clinical Annotation with Gemini
Downloads gemini install script.
Loads annotated variants into gemini database (gemini load).
Enables further clinical or genetic interpretation.
Notable Design Choices
## Use of Print Statements
Added throughout for tracking execution flow and debugging.
Useful for long pipelines where steps may take hours; it helps identify which step is in progress or stuck.
## Why use subprocess?
Provides close control over command execution.
Enables capturing output, errors, and checking success (via check=True).
Better than os.system() for robustness and error handling.
## Use of pathlib
Cross-platform file path management.
Simplifies path manipulations (/ operator).
Improves code readability and reduces bugs.
Downloading External Data
Downloads datasets, reference genome, annotation tools.
Checks if files exist to avoid redundant downloads.
Ensures reproducibility and automation.
## Why combine these tools?
Most bioinformatics workflows are command-line driven. Using Python allows automating, error checking, and integrating multiple commands seamlessly, making the entire pipeline reproducible and manageable.

## Final Remarks
Adjust paths, URLs, and parameters to your environment and datasets.
Install all external dependencies as specified.
Run the script from a Linux environment with appropriate permissions.
For Windows users, consider using WSL or adjusting commands accordingly.

## Conclusion
This Python script offers an automated, modular, and transparent way to perform complex genomic analyses. It provides clear instructions, progress updates, and control over each step—all while leveraging the powerful suite of command-line bioinformatics tools.