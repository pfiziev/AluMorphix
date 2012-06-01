if [ ! -n "$1" ]
then
  echo "usage: `basename $0` data-dir"
  exit $E_BADARGS
fi

rm -rf $1
mkdir $1

echo "python generate_genome.py $1"
python generate_genome.py $1

./bowtie2_build.sh $1

echo "python generate_random_reads.py $1"
python generate_random_reads.py $1

./bt2_vs_alu.sh $1

echo "python mark_ALU_reads.py $1"
python mark_ALU_reads.py $1


./bowtie2.sh $1

#./bowtie_build.sh $1
#./bowtie.sh $1



