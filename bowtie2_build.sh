if [ ! -n "$1" ]
then
  echo "usage: `basename $0` data-dir"
  exit $E_BADARGS
fi

cd ${@: -1}
rm -rf bowtie2_index
mkdir bowtie2_index
bowtie2-build --quiet genome.fa bowtie2_index/genome
samtools faidx genome.fa
samtools faidx alu.fa

cd -