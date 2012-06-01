import json
import random
import sys
import pysam
import utils
from utils import *

if __name__ == '__main__':
    if len(sys.argv) == 1:
        error('usage: %s data-dir' % sys.argv[0])

    utils.DATA_DIR = sys.argv[1]

    alu_info = json.load(open(data_d('alu.json')))

    reads_per_alu_bp = 0

    for fname in [data_d('pair1_vs_alu.sorted.bam'), data_d('pair2_vs_alu.sorted.bam')]:
        with pysam.Samfile(fname) as pair:
            for col in pair.pileup():
                reads_per_alu_bp += len(col.pileups)


    reads_per_alu_bp /= float(len(alu_info['seq']))


    inputf = pysam.Samfile(data_d('bt2_pair_end.sorted.bam'))
    reads_per_base = 0
    total_bases = 0
    for col in inputf.pileup():
        reads_per_base += len(col.pileups)
        total_bases += 1

    reads_per_base /= float(total_bases)

    print reads_per_alu_bp, reads_per_base
    print 'Estimated number of ALUs:', reads_per_alu_bp/reads_per_base
