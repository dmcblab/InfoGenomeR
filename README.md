
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
Please follow the guideline.
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
    <img height="300" src="https://github.com/YeonghunL/InfoGenomeR/blob/master/haplotype_graph.png">
  </a>
</p>

- Enter a directory for a subset of chromosomes, run the Eulerian path finding script, and run the plotting script for the candidate karyotype.\
`cd euler.15.19`\
`./DAG.sh F`\
`Rscript KAR_GRAPH.R`

<p align="center">
    <img height="150" src="https://github.com/YeonghunL/InfoGenomeR/blob/master/karyotype.png">
  </a>
</p>

