bowtie  bowtie_index/genome_NCHROMS_5_CL_100000_ALUFRAC_0.10 -1 data/genome_NCHROMS_5_CL_100000_ALUFRAC_0.10.fa.person.fa.pair1 -2 data/genome_NCHROMS_5_CL_100000_ALUFRAC_0.10.fa.person.fa.pair2 -S -f --all -X 600 -v 0 > data/bt1_pair_end.sam
samtools view -bS  data/bt1_pair_end.sam > data/bt1_pair_end.bam
samtools sort data/bt1_pair_end.bam data/bt1_pair_end.sorted
samtools index data/bt1_pair_end.sorted.bam
