if [ ! -n "$1" ]
then
  echo "usage: `basename $0` data-dir"
  exit $E_BADARGS
fi

cd ${@: -1}
echo "$(date): bowtie2-build"
bowtie2-build --quiet alu/alu.fa alu/bowtie2_index/alu_bt2

echo "$(date): bt2 pair1 vs alu"
bowtie2  -p 4 --quiet --local -x alu/bowtie2_index/alu_bt2 -U pair1.fastq | samtools view -bS - | samtools sort - pair1_vs_alu.sorted
samtools index pair1_vs_alu.sorted.bam


echo "$(date): bt2 pair2 vs alu"
bowtie2  -p 4 --quiet --local -x alu/bowtie2_index/alu_bt2 -U pair2.fastq | samtools view -bS - | samtools sort - pair2_vs_alu.sorted
samtools index pair2_vs_alu.sorted.bam

cd -

echo "$(date): bt2_vs_alu.sh finished"