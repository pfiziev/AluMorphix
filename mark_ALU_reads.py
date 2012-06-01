import cPickle as pickle

__author__ = 'pf'

import json
import random
import sys
import pysam
import utils
from utils import *



def clipped_bases(cigar):
    return sum(bases for op, bases in cigar if op == BAM_SOFTCLIP)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        error('usage: %s data-dir' % sys.argv[0])

    utils.DATA_DIR = sys.argv[1]

    filtered_reads = {}
    p1_fname = data_d('pair1_vs_alu.sorted.bam')
    p2_fname = data_d('pair2_vs_alu.sorted.bam')

    reads = {}

    full = set()
    partial = set()
    almost = set()

    for fname in [p1_fname, p2_fname]:
        filtered_reads[fname] = set()
        reads[fname] = {}
        logm('reading '+fname)
        with pysam.Samfile(fname, 'rb') as pair:
            for aln in pair:
                reads[fname][aln.qname] = (aln.seq, aln.qual)

                if not aln.is_unmapped:
                    clipped_ratio = float(clipped_bases(aln.cigar))/len(aln.seq)

                    if clipped_ratio < 0.1:
                        full.add(aln.qname)
                    elif clipped_ratio <= 0.5:
                        almost.add(aln.qname)
                    else:
                        partial.add(aln.qname)
                    filtered_reads[fname].add(aln.qname)

    both_full = 0
    pair1_outf = open(data_d('pair1.in_alu.fastq'), 'w')
    pair2_outf = open(data_d('pair2.in_alu.fastq'), 'w')
    alu_reads = set()
    logm('processing reads from ' + p1_fname)
    for rname in filtered_reads[p1_fname]:
        mate_rname = rname[:-1]+'2'

        if not (rname in full and mate_rname in full):
            pair1_outf.write('@%s\n%s\n+%s\n%s\n' % (rname, reads[p1_fname][rname][0], rname, reads[p1_fname][rname][1]))
            pair2_outf.write('@%s\n%s\n+%s\n%s\n' % (mate_rname, reads[p2_fname][mate_rname][0], mate_rname, reads[p2_fname][mate_rname][1]))
#            pair2_outf.write('>%s\n%s\n' % (mate_rname, reads[p2_fname][mate_rname]))
        else:
            both_full += 1

        alu_reads.add(rname)
        alu_reads.add(mate_rname)
    logm('processing reads from ' + p2_fname)
    for rname in filtered_reads[p2_fname]:
        mate_rname = rname[:-1]+'1'
        if rname not in alu_reads and mate_rname not in alu_reads:

            pair2_outf.write('@%s\n%s\n+%s\n%s\n' % (rname, reads[p2_fname][rname][0], rname, reads[p2_fname][rname][1]))
            pair1_outf.write('@%s\n%s\n+%s\n%s\n' % (mate_rname, reads[p1_fname][mate_rname][0], mate_rname, reads[p1_fname][mate_rname][1]))

#            pair2_outf.write('>%s\n%s\n' % (rname, reads[p2_fname][rname]))
#            pair1_outf.write('>%s\n%s\n' % (mate_rname, reads[p1_fname][mate_rname]))
            alu_reads.add(rname)
            alu_reads.add(mate_rname)

    pair1_outf.close()
    pair2_outf.close()

    pickle.dump(alu_reads, open(data_d('alu_reads.pickle'), 'w'), protocol= pickle.HIGHEST_PROTOCOL)

    logm('Full: %d\tAlmost: %d\tPartial: %d' % (len(full), len(almost), len(partial)))

    logm('Both mates are full ALUs: %d' % both_full)

