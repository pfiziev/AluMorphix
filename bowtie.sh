if [ ! -n "$1" ]
then
  echo "usage: `basename $0` data-dir"
  exit $E_BADARGS
fi

cd ${@: -1}

bowtie --quiet -a -f bowtie_index/genome -1 pair1.fa -2 pair2.fa -S bt_pair_end.sam
samtools view -bS bt_pair_end.sam > bt_pair_end.bam
samtools sort bt_pair_end.bam bt_pair_end.sorted
samtools index bt_pair_end.sorted.bam

cd -