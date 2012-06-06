if [ ! -n "$1" ]
then
  echo "usage: `basename $0` data-dir"
  exit $E_BADARGS
fi

cd ${@: -1}
rm -rf bowtie_index
mkdir bowtie_index
bowtie-build --quiet genome.fa bowtie_index/genome
samtools faidx genome.fa

cd -