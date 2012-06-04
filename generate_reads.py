import sys, random, string, json
import math
from utils import *
import utils

COVERAGE = 30
RLEN = 100 # read length
MAX_FRAGLEN = 500
MIN_FRAGLEN = 150


SEQ_ERROR_MEAN = 0.05
SEQ_ERROR_VARIANCE = 0.01
SEQ_ERROR_BINS = 5


def add_seq_errors(seq):
    bases = list(seq)
    error_rate = min(max(random.gauss(SEQ_ERROR_MEAN, SEQ_ERROR_VARIANCE), 0.0001 ), 1)
    error_probs = [(i/(RLEN/SEQ_ERROR_BINS) + 1)*error_rate/((SEQ_ERROR_BINS+1)/2.) for i in xrange(RLEN)]
    for i in xrange(len(seq)):
        if random.random() <= error_probs[i]:
            bases[i] = random.choice([n for n in NUCS if n != bases[i]])

    return ''.join(bases), ''.join(chr(int(33 - 10*math.log(p, 10))) for p in error_probs)


def gen_reads(haplotypes, nreads, pair1, pair2):
    for n in xrange(nreads):
        seq_id, seq = random.choice(haplotypes)
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

        seq1, qual1 = add_seq_errors(read[:RLEN])
        seq2, qual2 = add_seq_errors(reverse_compl_seq(read[len(read) - RLEN:]))

        seq1_id = '%s_%d_%s_%d_%d_%d_%d_p1' % (seq_id, n, strand, p1_st, p1_en, p2_st, p2_en)
        seq2_id = '%s_%d_%s_%d_%d_%d_%d_p2' % (seq_id, n, strand, p1_st, p1_en, p2_st, p2_en)

        pair1.write('@%s\n%s\n+%s\n%s\n' % (seq1_id, seq1, seq1_id, qual1))
        pair2.write('@%s\n%s\n+%s\n%s\n' % (seq2_id, seq2, seq2_id, qual2))



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



def generate_reads(haplotypes, pair1, pair2):

    reads_per_chrom = COVERAGE * ((len(haplotypes[0][1])+len(haplotypes[1][1]))/2)/(2*RLEN)

    logm('haplotype lengths: %d %d' % (len(haplotypes[0][1]), len(haplotypes[1][1])))
    logm('generating %d pair-end reads' % reads_per_chrom)

    gen_reads(haplotypes, reads_per_chrom, pair1, pair2)




if __name__ == '__main__':
    if len(sys.argv) == 1:
        error('usage: %s data-dir' % sys.argv[0])

    utils.DATA_DIR = sys.argv[1]


    pers_genome = data_d('subject_genome.fa')

    pair1_fname = data_d('pair1.fastq') #+'.M%dI%fMR%dCT%d'%(MUTATIONS, INDELS, MUT_RANGE, CELL_TYPES)+'.pair1'
    pair2_fname = data_d('pair2.fastq') #+'.M%dI%fMR%dCT%d'%(MUTATIONS, INDELS, MUT_RANGE, CELL_TYPES)+'.pair2'

    pair1 = open(pair1_fname, 'w')
    pair2 = open(pair2_fname, 'w')


    chrom_parts = []
    chrom_id  = None
    haplotypes = []
    for line in open(pers_genome):
        if line[0] == '>':
            if chrom_id is not None:
                haplotypes.append([chrom_id, ''.join(chrom_parts)])

                if len(haplotypes) == 2:
                    generate_reads(haplotypes, pair1, pair2)
                    haplotypes = []

            chrom_id = line[1:].split()[0]
            chrom_parts = []
            logm('reading %s' % chrom_id)
        else:
            chrom_parts.append(line.strip().upper())

    haplotypes.append([chrom_id, ''.join(chrom_parts)])
    generate_reads(haplotypes, pair1, pair2)

    pair1.close()
    pair2.close()



    logm('Output stored in: %s and %s'% (pair1_fname, pair2_fname))
