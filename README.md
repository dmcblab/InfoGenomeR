# InfoGenomeR
# Requirements
- bwa
- samtools
- blat
- R
- lpSolveApi
- ABSOLUTE
- BICseq2
# Tutorials
- Download demo files.
A simiulated cancer genome (haplotype coverage 5X, purity 75%) has 162 somatic SVs (true_SV_sets_somatic).\
Demo contains input files for InfoGenomeR. 

# Inputs
- BICseq norm files (bicseq_norm, bicseq_norm_germ)
- Initial SV calls (delly.format, manta.format, novobreak.format)
- Non-properly paired reads (NPE.fq1, NPE.fq2)
- Segmentation file for a control genome (copy_numbers.control)
- SNP calls (het_snps.format, hom_snps.format)
# How to generate inputs from BAM
Please follow the guideline.
# Check baselines for SVs
`Rscript SV_performance.R delly.format true_SV_sets_somatic 0 0 fmeasure`\
precision:0.6902174 recall:0.7987421 fmeasure:0.7405248\
`Rscript SV_performance.R manta.format true_SV_sets_somatic 0 0 fmeasure`\
precision:0.9495798 recall:0.7106918 fmeasure:0.8129496\
`Rscript SV_performance.R novobreak.format true_SV_sets_somatic 0 0 fmeasure`\
precision:0.9021739 recall:0.4968553 fmeasure:0.6408014\

0.6408014
# Running InfoGenomeR
Merge SV calls.\
`cat delly.format manta.format novobreak.format > SVs`\
Run a script for breakpoint graph construction.\
`./breakpoint_graph.sh somatic sample1 sample1 1.5 5 bicseq_norm bicseq_norm_germ copy_numbers.control hg19 5 null null 0`\
It takes an hour during five iterations and outputs SVs, copy numbers and a breakpoint graph.
# Check performance for SV calls.
`Rscript SV_performance.R SVs 0 0 fmeasure 1000`\
precision:0.9552239 recall:0.8050314 fmeasure:0.8737201
