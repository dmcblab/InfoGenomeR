# InfoGenomeR
- InfoGenomeR is the Integrative Framework for Genome Reconstruction that uses a breakpoint graph to model the connectivity among genomic segments at the genome-wide scale. InfoGenomeR integrates cancer purity and ploidy, total CNAs, allele-specific CNAs, and haplotype information to identify the optimal breakpoint graph representing cancer genomes.

# Requirements
- bwa
- samtools
- blat
- R
- lpSolveApi
- ABSOLUTE
- BICseq2

# Inputs
- BICseq norm files (bicseq_norm, bicseq_norm_germ)
- Initial SV calls (delly.format, manta.format, novobreak.format)
- Non-properly paired reads (NPE.fq1, NPE.fq2)
- Segmentation file for a control genome (copy_numbers.control)
- SNP calls (het_snps.format, hom_snps.format)
# How to generate inputs from BAM
Please follow the guideline.
# Usages
- Breakpoint graph construction.\
`./breakpoint_graph/breakpoint_graph.sh <mode> <sample_name> <cancer_type> <min_ploidy> <max_ploidy> <bicseq_norm> <bicseq_norm_germ> <copy_numbers.control> <fasta_prefix> <haplotype_coverage> <tumor_bam> <normal_bam> <chr_prefix>`
- Allele-specific graph construction.\
` ./allele_graph/allele_graph.sh <mode> <copy_numbers.control> <hom_snps.format> <het_snps.format> <fasta>`
- Haplotype graph construction.\
`./haplotype_graph/haplotype_graph.sh`

# Tutorials
- Download demo files. Demo contains input files for InfoGenomeR. 
- Tutorial 1: a simiulated cancer genome (haplotype coverage 5X, purity 75%) that has 162 somatic SVs (true_SV_sets_somatic).
- Tutorial 2: A549 cancer cell line that is triploidy with der(11)t(8;11) and der(19)t(8;19)x2. It has chromothripsis on chromosome 15.
From initial 4054 SV calls (2009 translocations), InfoGenomeR reconstructs der(11)t(8;11) and der(19)t(8;19)x2.

# Tutorial 1
# Check baselines for SVs
`Rscript SV_performance.R delly.format true_SV_sets_somatic 0 0`\
precision: 0.6902174 recall: 0.7987421 fmeasure: 0.7405248\
`Rscript SV_performance.R manta.format true_SV_sets_somatic 0 0`\
precision: 0.9495798 recall: 0.7106918 fmeasure: 0.8129496\
`Rscript SV_performance.R novobreak.format true_SV_sets_somatic 0 0`\
precision: 0.9021739 recall: 0.4968553 fmeasure: 0.6408014

# Running InfoGenomeR
- Merge SV calls.\
`cat delly.format manta.format novobreak.format > SVs`
- Run scripts for breakpoint graph construction.\
`./breakpoint_graph.sh somatic sample1 sample1 1.5 5 bicseq_norm bicseq_norm_germ copy_numbers.control hg19 5 null null 0`\
`./allele_graph/allele_graph.sh somatic copy_numbers.control hom_snps.format het_snps.format hg19.fa`\
`./haplotype_graph/haplotype_graph.sh`\
It takes a few hours during five iterations and outputs SVs, copy numbers and a breakpoint graph at the haplotype level.\
# Check performance for SV calls.
`Rscript SV_performance.R SVs true_SV_sets_somatic 0 0`\
precision: 0.9552239 recall: 0.8050314 fmeasure: 0.8737201

# Tutorial 2
