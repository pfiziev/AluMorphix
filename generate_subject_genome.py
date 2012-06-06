__author__ = 'pf'
import sys, random, string, json
from utils import *
import utils

NEW_ALUS = 1000
ALU_MUTATION_RATE = 0.05


N_HAPS = 2

SNP_FREQ = 10**-3
INDELS = 10**-4

RANDOM_INSERTS = 10**-5
RANDOM_INSERT_SIZE = 300

MAX_FRAGLEN = 500

def add_alus(seq, alu, n_alus):

    positions = set()
    #    positions = sorted(random.sample(xrange(len(seq)), n_alus))
    while len(positions) != n_alus:

        while True:
            pos = random.randint(3*MAX_FRAGLEN, len(seq) - 3*MAX_FRAGLEN)
            pos_is_ok = True
            for i in xrange(pos - MAX_FRAGLEN, pos + MAX_FRAGLEN):
                if seq[i] == 'N':
                    pos_is_ok = False
                    break
            if pos_is_ok:
                break
        positions.add(pos)

    positions = sorted(positions)

    ins_len = 0
    for i in xrange(n_alus):
        pos = positions[i] + ins_len

        to_insert, strand = orient(mutate(alu, mut_prob = ALU_MUTATION_RATE))

        ins_len += len(to_insert)

        positions[i] = {'ref_pos' : positions[i], 'pers_pos' : pos, 'strand' : strand}
        seq = seq[:pos] + to_insert + seq[pos:]

    return seq, positions

nucs = ['A', 'C', 'T', 'G']
def add_structural_variations(seq):
    parts = []
    prev_part = 0
    for i in xrange(1, len(seq)-1):
        if random.random() <= INDELS:
            if random.random() <= 0.5:
                parts.append(seq[prev_part:i-1])
            else:
                parts.append(seq[prev_part:i])
                parts.append(random.choice(nucs))
            prev_part = i
#        elif random.random() <= RANDOM_INSERTS:
#            parts.append(seq[prev_part:i])
#            parts.append(''.join(random.choice(nucs) for j in xrange(RANDOM_INSERT_SIZE)))
#            prev_part = i
#        print parts
    parts.append(seq[prev_part:])
#    print parts
    return ''.join(parts)

def generate_diploid_chromosomes(chrom_id, chrom_parts, alu, pers):

    chrom_seq = ''.join(chrom_parts)

    alu_positions = []

    for haplotype in xrange(N_HAPS):
        logm('generating haplotype %d from %s' % (haplotype, chrom_id))
        hap_id = chrom_id +'_' + str(haplotype)
        seq, positions = add_alus(add_structural_variations(add_snps(chrom_seq, snp_prob = SNP_FREQ)), alu, NEW_ALUS)
        alu_positions.append(positions)

        pers.write('>%s\n%s\n' % (hap_id, seq))

        logm('haplotype length: %d' % len(seq))


    return alu_positions


if __name__ == '__main__':
    if len(sys.argv) == 1:
        error('usage: %s data-dir' % sys.argv[0])

    utils.DATA_DIR = sys.argv[1]


    pers_genome = data_d('subject_genome.fa')

    pers = open(pers_genome, 'w')

    alu = read_alu()
    alu_positions = {}

    chrom_parts = []
    chrom_id  = None

    for line in open(data_d('genome/genome.fa')):
        if line[0] == '>':
            if chrom_id is not None:
                alu_positions[chrom_id] = generate_diploid_chromosomes(chrom_id, chrom_parts, alu, pers)

            chrom_id = line[1:].split()[0]
            chrom_parts = []
            logm('reading %s' % chrom_id)
        else:
            chrom_parts.append(line.strip().upper())

    alu_positions[chrom_id] = generate_diploid_chromosomes(chrom_id, chrom_parts, alu, pers)

    pers.close()

    json.dump(alu_positions, open(pers_genome + '.alu_positions.json', 'w'), indent = 1)

    logm('Output stored in: %s'% pers_genome)
    logm('inserted ALUs: %d' %  sum(len(hap_positions) for c in alu_positions for hap_positions in alu_positions[c]))
