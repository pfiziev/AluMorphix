if [ ! -n "$1" ]
then
  echo "usage: `basename $0` data-dir"
  exit $E_BADARGS
fi

echo "$(date): python generate_subject_genome.py $1"
python generate_subject_genome.py $1

echo "$(date): python generate_reads.py $1"
python generate_reads.py $1

echo "$(date): ./bt2_vs_alu.sh $1"
sh ./bt2_vs_alu.sh $1


echo "$(date): python mark_ALU_reads.py $1"
python mark_ALU_reads.py $1

echo "$(date): ./bowtie2.sh $1"
sh ./bowtie2.sh $1


echo "$(date): python alumorphix.py $1"
python alumorphix.py $1


echo "$(date): finished"



