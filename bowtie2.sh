if [ ! -n "$1" ]
then
  echo "usage: `basename $0` data-dir"
  exit $E_BADARGS
fi

cd ${@: -1}

#bowtie2 --quiet -f --local -x bowtie2_index/genome -1 pair1.fa -2 pair2.fa -S bt2_pair_end.sam
#samtools view -bS bt2_pair_end.sam > bt2_pair_end.bam
#samtools sort bt2_pair_end.bam bt2_pair_end.sorted
#samtools index bt2_pair_end.sorted.bam


echo "$(date): bowtie2 pairs in ALUs vs genome"

#bowtie2 --rdg 100,100 --rfg 100,100 -p 4 --quiet --local -x genome/bowtie2_index/genome -1 pair1.in_alu.fastq -2 pair2.in_alu.fastq | samtools view -bS - | samtools sort - bt2_pair.in_alu.sorted
bowtie2 -p 4 --quiet --local -x genome/bowtie2_index/genome -1 pair1.in_alu.fastq -2 pair2.in_alu.fastq | samtools view -bS - | samtools sort - bt2_pair.in_alu.sorted

samtools index bt2_pair.in_alu.sorted.bam


echo "$(date): bowtie2 all pairs vs genome"

#bowtie2 --rdg 100,100 --rfg 100,100 -p 4 --quiet --local -x genome/bowtie2_index/genome -1 pair1.fastq -2 pair2.fastq | samtools view -bS - | samtools sort - bt2_pair_end.sorted
bowtie2 -p 4 --quiet --local -x genome/bowtie2_index/genome -1 pair1.fastq -2 pair2.fastq | samtools view -bS - | samtools sort - bt2_pair_end.sorted

samtools index bt2_pair_end.sorted.bam


cd -