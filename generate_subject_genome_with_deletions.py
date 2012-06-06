import json
import utils
from utils import *
from generate_subject_genome import *

ALUS_TO_DELETE = 0.5


def delete_alus(seq, to_delete):
    alu_len = len(read_alu())
    parts = []
    prev_part = 0

    positions = []
    for i, alu in enumerate(sorted(to_delete, key = lambda a: a['pos'])):
        parts.append(seq[prev_part: alu['pos']])
        prev_part = alu['pos'] + alu_len + 1
        positions.append({'ref_pos' : alu['pos'], 'pers_pos' : alu['pos'] - i * alu_len, 'strand' : alu['strand']})
    parts.append(seq[prev_part : ])
    return ''.join(parts), positions

def generate_diploid_chromosomes(chrom_id, chrom_parts, known_alus, pers):
    chrom_seq = ''.join(chrom_parts)

    alu_positions = []

    for haplotype in xrange(N_HAPS):

        logm('generating haplotype %d from %s' % (haplotype, chrom_id))

        hap_id = chrom_id +'_' + str(haplotype)

        seq, positions = delete_alus(add_structural_variations(add_snps(chrom_seq, snp_prob = SNP_FREQ)),
                                     random.sample(known_alus, int(ALUS_TO_DELETE*len(known_alus))))

        alu_positions.append(positions)

        pers.write('>%s\n%s\n' % (hap_id, seq))

        logm('haplotype length: %d' % len(seq))
    return alu_positions


if __name__ == '__main__':
    if len(sys.argv) == 1:
        error('usage: %s data-dir' % sys.argv[0])

    utils.DATA_DIR = sys.argv[1]
    known_alus = json.load(open(data_d('genome/alu_positions.json')))

    deleted_alus = {}
    pers_genome = data_d('subject_genome.fa')

    pers = open(pers_genome, 'w')


    chrom_parts = []
    chrom_id  = None

    for line in open(data_d('genome/genome.fa')):
        if line[0] == '>':
            if chrom_id is not None:
                deleted_alus[chrom_id] = generate_diploid_chromosomes(chrom_id, chrom_parts, known_alus[chrom_id], pers)

            chrom_id = line[1:].split()[0]
            chrom_parts = []
            logm('reading %s' % chrom_id)
        else:
            chrom_parts.append(line.strip().upper())

    deleted_alus[chrom_id] = generate_diploid_chromosomes(chrom_id, chrom_parts, known_alus[chrom_id], pers)

    pers.close()

    json.dump(deleted_alus, open(pers_genome + '.deleted_alus.json', 'w'), indent = 1)

    logm('Output stored in: %s'% pers_genome)
    logm('deleted ALUs: %d' %  sum(len(hap_positions) for c in deleted_alus for hap_positions in deleted_alus[c]))

