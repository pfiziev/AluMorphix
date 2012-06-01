import json
import random
import sys

from utils import *
import utils

NCHROMS = 1
CHROM_LEN = 100000
ALU_LEN = 300
ALU_FRAC = 0.1 # percentage of genome that is Alu
ALU_PER_CHROM = int(CHROM_LEN*ALU_FRAC/ALU_LEN)
NON_ALU_BASES = int(CHROM_LEN*(1-ALU_FRAC))

if __name__ == '__main__':
    if len(sys.argv) == 1:
        error("usage: %s data-dir" % __file__)

    utils.DATA_DIR = sys.argv[1]



#    alu = ''.join(['A','C','T','G'][random.randint(0,3)] for j in xrange(ALU_LEN))

    alu = read_alu()
    print alu

    positions = {}

    f = open(data_d('genome.fa'), 'w')

    total_alus = 0
    for i in xrange(NCHROMS):
        chrom = []
        alu_positions = []
        c_pos = 0
        for j in xrange(NON_ALU_BASES + ALU_PER_CHROM):
            if random.random() <= float(ALU_PER_CHROM)/(NON_ALU_BASES + ALU_PER_CHROM):

                to_insert, strand = orient(mutate(alu))

                alu_positions.append({'pos': c_pos, 'strand' : strand})
                c_pos += len(to_insert)
                chrom.append(to_insert)
                total_alus += 1

            else:
                chrom.append(['A','C','T','G'][random.randint(0,3)])
                c_pos += 1


        f.write('>chr%d\n%s\n' % (i + 1, ''.join(chrom)))
        positions['chr%d'%(i+1)] = alu_positions

    f.close()
    print 'Total number of ALUs:', total_alus
    with open(data_d('alu.json'), 'w') as out:
        json.dump({'seq': alu, 'positions' : positions, 'TOTAL' : total_alus}, out, indent = 1)

    with open(data_d('alu.fa'), 'w') as out:
        out.write('>ALU\n%s\n' % alu)