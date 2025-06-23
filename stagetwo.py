import os  
import subprocess  
from pathlib import Path  

# Define directories  
raw_data_dir = Path("raw_data")  
fastqc_reports_dir = Path("Fastqc_Reports")  
mapping_dir = Path("Mapping")  
reference_fa = "hg19.chr5_12_17.fa"  
reference_gz = "hg19.chr5_12_17.fa.gz"  
reference_unzipped = "hg19.chr5_12_17.fa"  

def create_dirs():  
    print("Creating directories...")
    raw_data_dir.mkdir(parents=True, exist_ok=True)  
    fastqc_reports_dir.mkdir(parents=True, exist_ok=True)  
    mapping_dir.mkdir(parents=True, exist_ok=True)  
    print("Directories created.\n")
def download_files(): 
    print("Starting file downloads...") 
    datasets = [  
        # URLs  
        "https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz",  
        "https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz",  
        "https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz",  
        "https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz"  
    ]  
    for url in datasets:  
        filename = url.split('/')[-1]  
        filepath = raw_data_dir / filename  
        if not filepath.exists():  
            print(f"Downloading {filename}...")  
            subprocess.run(["wget", "-O", str(filepath), url], check=True)  
        else:  
            print(f"{filename} already exists. Skipping download.")  
        print("File downloads completed.\n")  
     

def download_reference():  
    print("Checking reference genome...") 
    ref_url = "https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz"  
    ref_path = Path(reference_gz)  
    if not ref_path.exists():  
        print("Downloading reference genome...")  
        subprocess.run(["wget", "-O", str(ref_path), ref_url], check=True)  
    else:  
        print("Reference genome already exists. Skipping download.")  
    print("Reference genome check complete.\n")

def unzip_reference(): 
    print("Unzipping reference genome...") 
    ref_path = Path(reference_gz)  
    if ref_path.exists():  
        subprocess.run(["gunzip", "-f", str(ref_path)], check=True)  
        print("Unzipping completed.\n")  
    else:  
        print("Reference gz file not found, skipping unzipping.\n")  

def run_fastqc():  
    print("Running FastQC on samples...")  
    samples = ["SLGFSK-N_231335", "SLGFSK-T_231336"]  
    for sample in samples:  
        r1 = raw_data_dir / f"{sample}_r1_chr5_12_17.fastq.gz"  
        r2 = raw_data_dir / f"{sample}_r2_chr5_12_17.fastq.gz"  
        subprocess.run(["fastqc", str(r1), "-o", str(fastqc_reports_dir)], check=True)  
        subprocess.run(["fastqc", str(r2), "-o", str(fastqc_reports_dir)], check=True) 
    print("FastQC completed.\n")  
        
def index_reference():  
    print("Indexing reference genome...")
    subprocess.run(["bwa", "index", reference_unzipped], check=True) 
    print("Indexing completed.\n")

def map_samples(): 
    print("Mapping samples to reference genome...")   
    samples = ["SLGFSK-N_231335", "SLGFSK-T_231336"]  
    for sample in samples:  
        r1 = raw_data_dir / f"{sample}_r1_chr5_12_17.fastq.gz"  
        r2 = raw_data_dir / f"{sample}_r2_chr5_12_17.fastq.gz"  
        sam_file = mapping_dir / f"{sample}.sam"  
        subprocess.run([  
            "bwa", "mem", "-R", f"@RG\\tID:{sample}\\tSM:{sample}",  
            reference_unzipped, str(r1), str(r2)  
        ], stdout=open(sam_file, 'w'), check=True)  
    print("Mapping completed.\n")  

def convert_sort_index_bam():
    samples = ["SLGFSK-N_231335", "SLGFSK-T_231336"]
    for sample in samples:
        print(f"Processing sample: {sample}")
        sam_file = mapping_dir / f"{sample}.sam"
        bam_file = mapping_dir / f"{sample}.sorted.bam"
        
        print(f"Converting SAM to BAM for {sample}...")
        subprocess.run(["samtools", "view", "-@20", "-S", "-b", str(sam_file)],
                       stdout=open(str(mapping_dir / f"{sample}.bam"), 'wb'), check=True)
        print(f"Conversion complete for {sample}.")

        print(f"Sorting BAM for {sample}...")
        subprocess.run(["samtools", "sort", "-@32", "-o", str(bam_file), str(mapping_dir / f"{sample}.bam")], check=True)
        print(f"Sorting complete for {sample}.")

        print(f"Indexing BAM for {sample}...")
        subprocess.run(["samtools", "index", str(bam_file)], check=True)
        print(f"Indexing complete for {sample}.\n")


def filter_bam():
    samples = ["SLGFSK-N_231335", "SLGFSK-T_231336"]
    for sample in samples:
        print(f"Filtering BAM for {sample}...")
        sorted_bam = mapping_dir / f"{sample}.sorted.bam"
        filtered_bam = mapping_dir / f"{sample}.filtered1.bam"
        subprocess.run(["samtools", "view", "-q", "1", "-f", "0x2", "-F", "0x8", "-b", str(sorted_bam)],
                       stdout=open(str(filtered_bam), 'wb'))
        print(f"Filtering complete for {sample}.\n")

def run_variant_calling_and_merging():
    print("Creating 'Variants' directory...")
    Path("Variants").mkdir(parents=True, exist_ok=True)

    print("Generating pileup files...")
    for sample in ["SLGFSK-N_231335", "SLGFSK-T_231336"]:
        print(f"Processing pileup for {sample}...")
        subprocess.run([
            "samtools", "mpileup", "-f", reference_unzipped,
            f"Mapping/{sample}.refilter.bam",
            "--min-MQ", "1", "--min-BQ", "28",
            "-o", f"Variants/{sample}.pileup"
        ], check=True)

    print("Running variant calling with VarScan...")
    subprocess.run([
        "java", "-jar", "VarScan.v2.3.9.jar",
        "somatic",
        "Variants/SLGFSK-N_231335.pileup",
        "Variants/SLGFSK-T_231336.pileup",
        "Variants/SLGFSK",
        "--normal-purity", "1",
        "--tumor-purity", "0.5",
        "--output-vcf", "1"
    ], check=True)

    print("Compressing VCF files...")
    for vcf in ["Variants/SLGFSK.snp.vcf", "Variants/SLGFSK.indel.vcf"]:
        print(f"Compressing {vcf}...")
        subprocess.run(["bgzip", vcf], check=True)

    print("Indexing compressed VCF files...")
    for gz in ["Variants/SLGFSK.snp.vcf.gz", "Variants/SLGFSK.indel.vcf.gz"]:
        print(f"Indexing {gz}...")
        subprocess.run(["tabix", gz], check=True)

    print("Merging VCF files into final VCF...")
    with open("Variants/SLGFSK.vcf", 'w') as outfile:
        subprocess.run([
            "bcftools", "merge",
            "Variants/SLGFSK.snp.vcf.gz",
            "Variants/SLGFSK.indel.vcf.gz"
        ], stdout=outfile, check=True)

def download_and_prepare_snpeff():
    print("Downloading snpEff...")
    subprocess.run([
        "wget", "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip",
        "-O", "snpEff_latest_core.zip"
    ], check=True)

    print("Unzipping snpEff...")
    subprocess.run(["unzip", "-o", "snpEff_latest_core.zip"], check=True)

    print("Downloading snpEff database...")
    subprocess.run([
        "java", "-jar", "snpEff.jar", "download", "hg19"
    ], check=True)

    print("Annotating variants with snpEff...")
    with open("Variants/SLGFSK.ann.vcf", 'w') as outfile:
        subprocess.run([
            "java", "-Xmx8g", "-jar", "snpEff/snpEff.jar", "hg19", "Variants/SLGFSK.vcf"
        ], stdout=outfile, check=True)

def run_gemini_annotation():
    print("Downloading gemini install script...")
    subprocess.run([
        "wget", "https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py",
        "-O", "gemini_install.py"
    ], check=True)

    # Uncomment below line if you want to run the install automatically
    # print("Installing gemini...")
    # subprocess.run(["python", "gemini_install.py", "/usr/local", "/usr/local/share/gemini"], check=True)

    print("Loading gemini database...")
    subprocess.run([
        "gemini", "load",
        "-v", "Variants/SLGFSK.ann.vcf",
        "-t", "snpEff",
        "Annotation/gemini.db"
    ], check=True)

create_dirs()
download_files()
download_reference()
unzip_reference()
run_fastqc()
index_reference()
map_samples()
convert_sort_index_bam()
filter_bam()
run_variant_calling_and_merging()
download_and_prepare_snpeff()
run_gemini_annotation()