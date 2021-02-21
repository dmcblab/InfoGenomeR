
# InfoGenomeR
- InfoGenomeR is the Integrative Framework for Genome Reconstruction that uses a breakpoint graph to model the connectivity among genomic segments at the genome-wide scale. InfoGenomeR integrates cancer purity and ploidy, total CNAs, allele-specific CNAs, and haplotype information to identify the optimal breakpoint graph representing cancer genomes.

<p align="center">
    <img height="150" src="https://github.com/YeonghunL/InfoGenomeR/blob/master/doc/overview.png">
  </a>
</p>

# Requirements
- bwa
- samtools
- bedtools
- blat
- R
- lpSolveApi
- ABSOLUTE
- BICseq2

# Inputs
- genome-binning read depths (cn_norm, cn_norm_germ)
- Initial SV calls (delly.format, manta.format, novobreak.format)
- Non-properly paired reads (NPE.fq1, NPE.fq2)
- SNP calls (het_snps.format, hom_snps.format)
# Outputs
- Haplotype-resolved SVs and CNAs (SVs.CN_opt.phased, copy_number.CN_opt.phased)
- Haplotype graph (node_keys, edge_information.txt)
- Karyotypes (Eulerian_path.0)
- SV cluster and topology (cluster_sv)
# How to generate inputs from BAM
- We provide scripts for generating inputs using BIC-seq2, DELLY2, Manta, and novobreak, but we recommend that a user follow the guideline of each tool.
- Generate cn_norm and cn_norm_germ.
    - Install BICseq2 and modifiedSamtools.
    - Set configs in the script `bicseq_preprocess.sh` and run the script for cn_norm.
    - Set configs in the script `bicseq_preprocess_germ.sh` and run the script for cn_norm_germ.
```reference=hg19.fa
norm_script=./NBICseq-norm_v0.2.4/NBICseq-norm.pl
map_file=./NBICseq-norm_v0.2.4/hg19.CRG.50bp/
read_length=100
fragment_size=350
tumor_bam=tumor.bam
```
 `./preprocessing/bicseq_preprocess.sh`
        
 ```reference=hg19.fa
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
 
# Running InfoGenomeR
- Breakpoint graph construction.\
`./breakpoint_graph/breakpoint_graph.sh <mode> <sample_name> <cancer_type> <min_ploidy> <max_ploidy> <bicseq_norm> <bicseq_norm_germ> <copy_numbers.control> <fasta_prefix> <haplotype_coverage> <tumor_bam> <normal_bam> <chr_prefix>`
- Allele-specific graph construction.\
` ./allele_graph/allele_graph.sh <mode> <copy_numbers.control> <hom_snps.format> <het_snps.format> <fasta>`
- Haplotype graph construction.\
`./haplotype_graph/haplotype_graph.sh`

# Demos
- Download demo files. Demo contains input files for InfoGenomeR. 
- Demo 1: a simiulated cancer genome (haplotype coverage 5X, triploidy, purity 75%) that has 162 somatic SVs (true_SV_sets_somatic). [Demo 1](http://www.gcancer.org/InfoGenomeR/demo1.tar.gz).
- Demo 2: A549 cancer cell line that is triploidy with der(11)t(8;11) and der(19)t(8;19)x2. It has chromothripsis on chromosome 15. [Demo 2](http://www.gcancer.org/InfoGenomeR/demo2.tar.gz).
From initial 4054 SV calls, InfoGenomeR reconstructs der(11)t(8;11) and der(19)t(8;19)x2.

# Demo 1
- Check baselines for SVs.\
`Rscript SV_performance.R delly.format true_SV_sets_somatic 0 0`\
precision: 0.6902174 recall: 0.7987421 fmeasure: 0.7405248\
`Rscript SV_performance.R manta.format true_SV_sets_somatic 0 0`\
precision: 0.9495798 recall: 0.7106918 fmeasure: 0.8129496\
`Rscript SV_performance.R novobreak.format true_SV_sets_somatic 0 0`\
precision: 0.9021739 recall: 0.4968553 fmeasure: 0.6408014

- Running InfoGenomeR.\
Merge SV calls.\
`cat delly.format manta.format novobreak.format > SVs`\
Run scripts for breakpoint graph construction.\
`./breapkpoint_graph/breakpoint_graph.sh somatic sample1 sample1 2 4 bicseq_norm bicseq_norm_germ copy_numbers.control hg19 5 null null 0`\
`./allele_graph/allele_graph.sh somatic copy_numbers.control hom_snps.format het_snps.format hg19.fa`\
`./haplotype_graph/haplotype_graph.sh`\
It takes a few hours during five iterations and outputs SVs, copy numbers and a breakpoint graph at the haplotype level.

- Running JaBbA for comparison study.\
`jba junctions.vcf coverage.txt -p 3.4 -q 0.75`

- Check performance for SV calls from InfoGenomeR.\
`Rscript SV_performance.R SVs true_SV_sets_somatic 0 0`\
precision: 0.9552239 recall: 0.8050314 fmeasure: 0.8737201
- Check performance for SV calls from JaBbA.\
`Rscript SV_performance.R jabba.simple.vcf.format true_SV_sets_somatic 0 0`\
precision: 0.9444444 recall: 0.7484277 fmeasure: 0.8350877
# Demo 2
- Run scripts for breakpoint graph construction.\
`./breakpoint_graph/breakpoint_graph.sh total A549 LUSC 2 4 bicseq_norm bicseq_norm_germ null hg19 5 null null 0`\
`./breakpoint_graph/breakpoint_graph_simplifying.sh total A549 LUSC 2 4 bicseq_norm bicseq_norm_germ null hg19 5 null null 0`\
`./allele_graph/allele_graph.sh total null hom_snps.format het_snps.format hg19.fa`\
`./haplotype_graph/haplotype_graph.sh`\
It takes a day during 25 iterations (maximum 24 threads and 256Gb memory).
- Run the script for plotting the haplotype graph.\
`Rscript ACN_GRAPH.R 3,8,11,15,19`
<p align="center">
    <img height="300" src="https://github.com/YeonghunL/InfoGenomeR/blob/master/doc/haplotype_graph.png">
  </a>
</p>

- Enter a directory for a subset of chromosomes, run the Eulerian path finding script, and run the plotting script for the candidate karyotype.\
`cd euler.15.19`\
`./DAG.sh F`\
`Rscript KAR_GRAPH.R`

<p align="center">
    <img height="150" src="https://github.com/YeonghunL/InfoGenomeR/blob/master/doc/karyotype.png">
  </a>
</p>

