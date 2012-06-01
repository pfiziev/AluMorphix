import sys, random, string, json
import math
from utils import *
import utils

COVERAGE = 10
RLEN = 100 # read length
MAX_FRAGLEN = 500
MIN_FRAGLEN = 150

NEW_ALUS = 1000

SEQ_ERROR_RATE = 0.05

N_HAPS = 2
SNP_FREQ = 0.001

def add_seq_errors(seq, error_rate = 0):
    return mutate(seq, mut_prob = error_rate)


def gen_reads(seq, seq_id, nreads, pair1, pair2):
    for n in xrange(nreads):

        while True:
            pos = random.randint(0, len(seq) - MAX_FRAGLEN)
            strand = random.choice('+-')
            endp = pos + random.randint(MIN_FRAGLEN, MAX_FRAGLEN)
            read = seq[pos : endp]

            if strand == '-': 
                read = reverse_compl_seq(read)

            if 'N' not in read:
                break



        p1_st = pos + 1
        p1_en = pos + RLEN
        
        p2_st = endp - RLEN + 1
        p2_en = endp

        if strand == '-':
            p1_st, p1_en, p2_st, p2_en = p2_st, p2_en, p1_st, p1_en
        
        seq1 = add_seq_errors(read[:RLEN], error_rate = SEQ_ERROR_RATE)
        seq2 = add_seq_errors(reverse_compl_seq(read[len(read)- RLEN:]), error_rate = SEQ_ERROR_RATE)

#        p1_st, p1_en, seq1 = mutate(p1_st, p1_en, seq1)
#        p2_st, p2_en, seq2 = mutate(p2_st, p2_en, seq2)

        pair1.write('>%s_%d_%s_%d_%d_%d_%d_p1\n%s\n' % (seq_id, n, strand, p1_st, p1_en, p2_st, p2_en, seq1))

        pair2.write('>%s_%d_%s_%d_%d_%d_%d_p2\n%s\n' % (seq_id, n, strand, p1_st, p1_en, p2_st, p2_en, seq2))


def add_alus(seq, alu, n_alus):

    positions = []
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
        positions.append(pos)

    positions.sort()

    ins_len = 0
    for i in xrange(n_alus):
        pos = positions[i] + ins_len

        to_insert, strand = orient(mutate(alu, mut_prob = 0.05))

        ins_len += len(to_insert)

        positions[i] = {'ref_pos' : positions[i], 'pers_pos' : pos, 'strand' : strand}
        seq = seq[:pos] + to_insert + seq[pos:]

    return seq, positions



def generate_haplotype_reads(chrom_id, chrom_parts, alu, pair1, pair2, pers):

    chrom_seq = ''.join(chrom_parts)

    alu_positions = []

    for haplotype in xrange(N_HAPS):
        logm('generating haplotype %d from %s' % (haplotype, chrom_id))
        hap_id = chrom_id +'_' + str(haplotype)
        seq, positions = add_alus(add_snps(chrom_seq, snp_prob = SNP_FREQ), alu, NEW_ALUS)
        alu_positions.append(positions)

        pers.write('>%s\n%s\n' % (hap_id, seq))

        reads_per_chrom = COVERAGE * len(seq)/(2*RLEN)

        logm('haplotype length: %d' % len(seq))
        logm('generating %d pair-end reads' % reads_per_chrom)

        gen_reads(seq, hap_id, reads_per_chrom, pair1, pair2)


    return alu_positions


if __name__ == '__main__':
    if len(sys.argv) == 1:
        error('usage: %s data-dir' % sys.argv[0])

    utils.DATA_DIR = sys.argv[1]


    pers_genome = data_d('subject_genome.fa')
    
    pair1_fname = data_d('pair1.fa') #+'.M%dI%fMR%dCT%d'%(MUTATIONS, INDELS, MUT_RANGE, CELL_TYPES)+'.pair1'
    pair2_fname = data_d('pair2.fa') #+'.M%dI%fMR%dCT%d'%(MUTATIONS, INDELS, MUT_RANGE, CELL_TYPES)+'.pair2'

    pers = open(pers_genome, 'w')
    pair1 = open(pair1_fname, 'w')
    pair2 = open(pair2_fname, 'w')

    alu = read_alu()
    alu_positions = {}

    chrom_parts = []
    chrom_id  = None

    for line in open(data_d('genome.fa')):
        if line[0] == '>':
            if chrom_id is not None:
                alu_positions[chrom_id] = generate_haplotype_reads(chrom_id, chrom_parts, alu, pair1, pair2, pers)

            chrom_id = line[1:].split()[0]
            chrom_parts = []
            logm('reading %s' % chrom_id)
        else:
            chrom_parts.append(line.strip().upper())

    alu_positions[chrom_id] = generate_haplotype_reads(chrom_id, chrom_parts, alu, pair1, pair2, pers)

    pair1.close()
    pair2.close()
    pers.close()


    json.dump(alu_positions, open(pers_genome + '.alu_positions.json', 'w'), indent = 1)

    logm('Output stored in: %s and %s'% (pair1_fname, pair2_fname))
    logm('inserted ALUs: %d' %  sum(len(hap_positions) for c in alu_positions for hap_positions in alu_positions[c]))