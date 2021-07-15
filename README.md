# InfoGenomeR
- InfoGenomeR is the Integrative Framework for Genome Reconstruction that uses a breakpoint graph to model the connectivity among genomic segments at the genome-wide scale. InfoGenomeR integrates cancer purity and ploidy, total CNAs, allele-specific CNAs, and haplotype information to identify the optimal breakpoint graph representing cancer genomes.

<p align="center">
    <img height="700" src="https://github.com/YeonghunL/InfoGenomeR/blob/master/doc/overview.png">
  </a>
</p>

# Requirements
- Executables in your PATH
    - bwa (version 0.7.15)
    - samtools (version 1.3)
    - bedtools (version 1.3)
    - blat (version 36)
    - blastn (version 2.2.30+)
- R (version 3.4.3) and libraries
    - lpSolveAPI (version 5.5.2.0.17)
    - ABSOLUTE (version 1.0.6)
    - fitdistrplus (version 1.1.11)
    - plyr (version 1.8.4)
- BIC-seq2 (version 0.7.2)
# Environment settings
- download InfoGenomeR and BIC-seq2 in your working directory. 
```
## For example, the working directory is /home/dmcblab.
cd /home/dmcblab

#### Download InfoGenomeR
git clone https://github.com/dmcblab/InfoGenomeR.git

#### Download BIC-seq2
wget http://compbio.med.harvard.edu/BIC-seq/NBICseq-seg_v0.7.2.tar.gz
tar -xvf NBICseq-seg_v0.7.2.tar.gz

```
- download haplotype DAGs, and bwa-indexed reference genome if it doesn't exist (GRCh37 or GRCh38).
```
## GRCh37
wget https://zenodo.org/record/5105505/files/GRCh37.tar.xz
tar Jxvf GRCh37.tar.xz
wget https://zenodo.org/record/5105505/files/haplotype_1000G.tar.xz
tar Jxvf haplotype_1000G.tar.xz

## GRCh38
wget https://zenodo.org/record/5105505/files/GRCh38.tar.xz
tar Jxvf GRCh38.tar.xz
wget https://zenodo.org/record/5105505/files/haplotype_1000G_GRCh38.tar.xz?download=1
tar Jxvf haplotype_1000G_GRCh38.tar.xz
```
- set the BIC-seq2 path, InfoGenomeR_lib, and Haplotype path.
```
export BICseq2_path=/home/dmcblab/NBICseq-seg_v0.7.2
export InfoGenomeR_lib=/home/dmcblab/InfoGenomeR
export Haplotype_path=/home/dmcblab/haplotype_1000G ## for GRCh37. If GRCh38 was used, export Haplotype_path=/home/dmcblab/haplotype_1000G_38
```
- set executables (bwa, samtools ...) in your PATH. You may use precompiled binaries in the InfoGenomeR/ext folder if they work in your computer. Otherwise, please install them.
```
export PATH=$InfoGenomeR_lib/ext:$PATH
```
- set the PATH environment.
```
export PATH=$InfoGenomeR_lib/breakpoint_graph:$InfoGenomeR_lib/allele_graph:$InfoGenomeR_lib/haplotype_graph:$PATH
```

# Inputs
- genome-binning read depths (cn_norm, cn_norm_germ)
- Initial SV calls (delly.format, manta.format, novobreak.format)
- Non-properly paired reads (NPE.fq1, NPE.fq2)
- Reference genome (GRCh37.fa, .fa.fai, .2bit, .bwt2) or (GRCh38.fa,.fa.fai, .2bit, .bwt2)
- SNP calls (het_snps.format, hom_snps.format)
- Haplotype DAGs (haplotype_1000G) or (haplotype_1000G_GRCh38)

# Outputs
- Haplotype-resolved SVs and CNAs (SVs.CN_opt.phased, copy_number.CN_opt.phased)
- Haplotypes (haplotype)
- Purity and ploidy (purity_ploidy)
- Haplotype graph (node_keys, edge_information.txt)
- Karyotypes (Eulerian_path.0)
- SV cluster and topology (cluster_sv)

# Running InfoGenomeR
- Breakpoint graph construction
```
Usage: breakpoint_graph <SVs> <cn_norm> [options]

Options:
        -m, --mode (required)
                 Select the mode (germline, total, somatic, simplification)
        -i, --lambda_ini (required)
                 Initial lambda for the first round iterations (default: 1)
        -f, --lambda_fi (required)
                 Final lambda for the second round iterations (default: 4)
        -t, --cancer_type (optional)
                 Cancer type used for ABSOLUTE estimation (BRCA, GBM, OV, ...). If unknown, write null.
        -n, --min_ploidy (optional)
                 minimum cancer ploidy used for ABSOLUTE estimation (default: 1.5)
        -x, --max_ploidy (optional)
                 maximum cnacer ploidy used for ABSOLUTE estimation (default: 5)
        -d, --npe_dir (optional)
                 Directory that contains NPE.fq1 and NPE.fq2 (non-properly paired reads). If it is not assigned, InfoGenomeR runs without NP reads mapping.
        -g, --ref_genome (required for NP reads mapping)
                 Fasta prefix (hg19 or hg38). Enter the prefix without .2bit and .fa extension.
        -c, cn_norm_germ (required for somatic mode)
                 Directory that contains copy number bins from a control genome.
        -s, --germ_LocSeq_result (required for somatic mode)
                 Local segmentation results from a control genome.
        -o, --out_dir
                 If it already exists, results are overlaid (default: InfoGenomeR_job)
        -h, --help
```
- Allele-specific graph construction
```
Usage: allele_graph <hom_snps.format> <het_snps.format> [options]

Options:
        -m, --mode (required)
                 Select the mode (germline, total, somatic)
        -s, --germ_LocSeq_result (required for somatic mode)
                 Local segmentation results from a control genome.
        -g, --ref_genome (required for NP reads mapping)
                 Fasta prefix (hg19 or hg38). Enter the prefix without .2bit and .fa extension.
        -o, --breakpoint_graph_dir
                 The output directory of breakpoint graph construciton
        -t, --threads
                 The number of threads
        -h, --help
```
- Haplotype graph construction
```
Usage: haplotype_graph [options]

        -o, --allele_graph_dir
                 The output directory of allele graph construciton
        -t, --threads
                 The number of threads
        -h, --help
```
# How to generate inputs from BAM
- We provide scripts for generating inputs using BIC-seq2, DELLY2, Manta, and novobreak, but we recommend that a user follow the guideline of each tool.
- Generate cn_norm and cn_norm_germ.
    - Install BICseq2 and modifiedSamtools.
    - Set configs in the script `bicseq_preprocess.sh` and run the script for cn_norm.
    - Set configs in the script `bicseq_preprocess_germ.sh` and run the script for cn_norm_germ.
```
reference=hg19.fa
norm_script=./NBICseq-norm_v0.2.4/NBICseq-norm.pl
map_file=./NBICseq-norm_v0.2.4/hg19.CRG.50bp/
read_length=100
fragment_size=350
tumor_bam=tumor.bam
```
 `./preprocessing/bicseq_preprocess.sh`
        
 ```
reference=hg19.fa
norm_script=./NBICseq-norm_v0.2.4/NBICseq-norm.pl
map_file=./NBICseq-norm_v0.2.4/hg19.CRG.50bp/
read_length=100
fragment_size=350
normal_bam=normal.bam
```
`./preprocessing/bicseq_preprocess_germ.sh`
- Generate NPE.fq1 and NPE.fq2
```
samtools view -h -F 2 -b simulated_sorted.bam >simulated_sorted_NPE.bam
samtools sort -n simulated_sorted_NPE.bam > simulated_sorted_NPE_sorted.bam
bamToFastq  -i simulated_sorted_NPE_sorted.bam  -fq NPE.fq1 -fq2 NPE.fq2
```
- Generate delly SV calls (delly.format).
    - Install DELLY2, set configs, and run the script.
```
ref=hg19.fa
tsv=samples.tsv
tumor_bam=tumor.bam
normal_bam=normal.bam
```
`./preprocessing/delly_somatic.sh 3 20`
- Generate manta SV calls (manta.format).
    - Install Manta, set configs, and run the script.
```
reference=hg19.fa
script1=./manta-1.1.0.centos5_x86_64/bin/configManta.py
script2=./runWorkflow.py
normal_bam=normal.bam
tumor_bam=tumor.bam
```
`./preprocessing/manta_somatic.sh`
- Generate novobreak SV calls (novobreak.format).
    - Install novobreak and copy the scripts `./preprocessing/run_novoBreak_orient_partial_corrected.sh` and `./preprocessing/infer_sv_orientation.pl` to the novobreak-installed folder.
    - Set configs, and run the script.
```
novo_dir=./novoBreak_distribution_v1.1.3rc
reference=hg19.fa
tumor_bam=tumor.bam
normal_bam=normal.bam
```
`./preprocessing/novobreak_somatic.sh`
- Generate SNP calls (het_snps.format and hom_snps.format).
    - Set configs and run the script.
```
reference=hg19.fa
normal_bam=normal.bam
tumor_bam=tumor.bam
```
`./preprocessing/het_SNP_detection_somatic.sh`
 
# Tutorials
- Download demo files. Demo contains input files for InfoGenomeR. 
- Tutorial 1:Germline and somatic mode (GRCh37). a simiulated cancer genome (haplotype coverage 5X, triploidy, purity 75%) that has 162 somatic SVs (true_SV_sets_somatic). [Tutorial_1](https://zenodo.org/record/5105505/files/tutorial1.tar.xz).
- Tutorial 2:Somatic mode (GRCh38). a simiulated cancer genome (haplotype coverage 10X, triploidy, purity 75%). [Tutorial_2](https://zenodo.org/record/4545666/files/GRCh38.triploidy.f10.p0.75.tar.xz)

# Tutorial 1 (GRCh37)
```
wget https://zenodo.org/record/5105505/files/tutorial1.tar.xz
tar Jxvf tutorial1.tar.xz
cd tutorial1
cp $InfoGenomeR_lib/etc/SV_performance.R ./
```
- Check baselines for SVs.
```
## delly
Rscript SV_performance.R delly.format true_SV_sets_somatic
precision:0.690217 recall:0.798742 fmeasure: 0.740525
## manta
Rscript SV_performance.R manta.format true_SV_sets_somatic
precision:0.949580 recall:0.710692 fmeasure: 0.812950
## novobreak
Rscript SV_performance.R novobreak.format true_SV_sets_somatic
precision:0.902174 recall:0.496855 fmeasure: 0.640801
```
- Running InfoGenomeR.
```
## Generate a germline copy number profile
breakpoint_graph -m germline germline_delly.format -o germline_job cp_norm_germ
cp germline_job/copy_numbers ./copy_numbers.control

## Merge SV calls.
cat delly.format manta.format novobreak.format > SVs
## breakpoint graph construction
breakpoint_graph -m somatic -d ./ SVs ./cp_norm -g /home/dmcblab/GRCh37/GRCh37 -c cp_norm_germ -s copy_numbers.control -o somatic_job
## allele graph construction
allele_graph -m somatic -s copy_numbers.control hom_snps.format het_snps.format -o somatic_job -g /home/dmcblab/GRCh37/GRCh37 -t 23
## haplotype graph construction
haplotype_graph -o somatic_job -t 6 ## 40GB memory per one thread.
```
It takes a few hours during five iterations and outputs SVs, copy numbers and a breakpoint graph at the haplotype level.

- Check performance for SV calls from InfoGenomeR.
```
## The output directory is InfoGenomeR_output
Rscript SV_performance.R somatic_job/InfoGenomeR_output/SVs.CN_opt.phased true_SV_sets_somatic
precision:0.955224 recall:0.805031 fmeasure: 0.873720
```

# Tutorial 2 (GRCh38)
```
wget://zenodo.org/record/4545666/files/GRCh38.triploidy.f10.p0.75.tar.xz
tar Jxvf GRCh38.triploidy.f10.p0.75.tar.xz
cd GRCh38.triploidy.f10.p0.75 
cp $InfoGenomeR_lib/etc/SV_performance.R ./
### Change the environment to GRCh38.
export Haplotype_path=/home/dmcblab/haplotype_1000G_GRCh38
```
- Running InfoGenomeR.
```
## Merge SV calls.
cat delly.format manta.format novobreak.format > SVs
## breakpoint graph construction
breakpoint_graph -m somatic -d ./ SVs ./cp_norm -g /home/dmcblab/GRCh38/GRCh38 -c cp_norm_germ -s copy_numbers.control -o somatic_job
## allele graph construction
allele_graph -m somatic -s copy_numbers.control hom_snps.format het_snps.format -o somatic_job -g /home/dmcblab/GRCh38/GRCh38 -t 23
## haplotype graph construction
haplotype_graph -o somatic_job -t 6 ## 40GB memory per one thread.
```
- Check performance for SV calls from InfoGenomeR.
```
## The output directory is InfoGenomeR_output
Rscript SV_performance.R somatic_job/InfoGenomeR_output/SVs.CN_opt.phased true_SV_sets_somatic

```
