samtools sort data/bt2_pair_end.bam data/bt2_pair_end.sorted
samtools index data/bt2_pair_end.sorted.bam
samtools faidx data/genome_NCHROMS_5_CL_100000_ALUFRAC_0.10.fa
