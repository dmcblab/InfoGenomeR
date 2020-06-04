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
- Download demo files
A simiulated cancer genome has 162 somatic SVs (true_SV_sets_somatic)
Demo contains input files for InfoGenomeR. 

# Inputs
- BICseq norm files (bicseq_norm, bicseq_norm_germ)
- Initial SV calls (delly.format, manta.format, novobreak.format)
- Non-properly paired reads (NPE.fq1, NPE.fq2)
- Segmentation file for a control genome (copy_numbers.control)
# How to generate inputs from BAM
# Check baselines for SVs
`Rscript SV_performance.R delly.format true_SV_sets_somatic 0 0 fmeasure 1000`
`Rscript SV_performance.R manta.format true_SV_sets_somatic 0 0 fmeasure 1000`
`Rscript SV_performance.R novobreak.format true_SV_sets_somatic 0 0 fmeasure 1000`

# Running InfoGenomeR
Merge SV calls.
`cat delly.format manta.format novobreak.format > SVs`
