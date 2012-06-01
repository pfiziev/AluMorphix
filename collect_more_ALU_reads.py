from pyfasta.fasta import Fasta
import pysam
from generate_subject_genome_and_reads import MAX_FRAGLEN
from utils import *
import utils
import sys


if __name__ == '__main__':
    if len(sys.argv) == 1:
        error('usage: %s data-dir' % sys.argv[0])

    utils.DATA_DIR = sys.argv[1]

    genomef = Fasta(data_d('genome.fa'), flatten_inplace=True)
    inputf = pysam.Samfile(data_d('bt2_pair.in_alu.sorted.bam'))
    outf = open(data_d('alu_regions.fa'), 'w')

    for chrom in inputf.references:
        logm('processing ' + chrom)
        regions = []
        max_reg_start = -1
        max_reg_end = -1

        for col in inputf.pileup(reference = chrom):
            if len(col.pileups) < 4:
                continue
            c_reg_start = col.pos - MAX_FRAGLEN
            c_reg_end = col.pos + MAX_FRAGLEN

            if max_reg_start <= c_reg_start <= max_reg_end:
                max_reg_end = c_reg_end
            else:
                regions.append((max_reg_start, max_reg_end))

                max_reg_start = c_reg_start
                max_reg_end   = c_reg_end

        regions.append((max_reg_start, max_reg_end))

        logm('cutting %d regions from %s' % ( len(regions)-1, chrom))
        chrom_seq = genomef[chrom]
        for start, end in regions:
            if start == -1: continue
            outf.write('>%s_%d_%d\n%s\n' % (chrom, start, end, chrom_seq[start:end]))




    outf.close()
    inputf.close()
    logm('done')