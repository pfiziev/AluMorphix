if [ ! -n "$1" ]
then
  echo "usage: `basename $0` data-dir"
  exit $E_BADARGS
fi

cd ${@: -1}
echo "$(date) concatenating alu pairs"
cat pair1.in_alu.fa pair2.in_alu.fa > alu_reads.fa

echo "$(date) building index for alu reads"
bowtie2-build --quiet alu_reads.fa bowtie2_index/alu_reads

echo "$(date): bt2 pair1 vs alu_reads"
bowtie2  -p 4 --quiet -f -a --local -x bowtie2_index/alu_reads -U pair1.fa | samtools view -bS -F 4 - | samtools sort - pair1_vs_alu_reads.sorted
samtools index pair1_vs_alu_reads.sorted.bam

echo "$(date): bt2 pair2 vs alu_reads"
bowtie2  -p 4 --quiet -f -a --local -x bowtie2_index/alu_reads -U pair2.fa | samtools view -bS -F 4 - | samtools sort - pair2_vs_alu_reads.sorted
samtools index pair2_vs_alu_reads.sorted.bam

cd -