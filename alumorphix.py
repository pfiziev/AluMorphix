import copy
import json
from pprint import pformat
import sys
import math
from numpy.ma.core import mean, var
import pysam
import numpy as np
from generate_subject_genome_and_reads import RLEN
from utils import *
import utils
import cPickle as pickle

__author__ = 'pf'

WIN_LENGTH = 30
MIN_MAPQ = 10
MIN_PEAK = 5
MIN_ALU_READS = 2

rnames = lambda col: set([r.alignment.qname for r in col.pileups if r.alignment.mapq >= MIN_MAPQ])


def window_stats(window):
    return float(len(window[0]['reads'] & window[-1]['reads']))/(len(window[0]['reads']))
#    return float(len(window[0]['reads'] & window[-1]['reads']))/(len(window[0]['reads']) + len(window[-1]['reads']))

def potential_ALU_insert(window, mean_heads, mean_tails):
    max_heads = 0
    max_tails = 0
    max_heads_pos = -1
    max_tails_pos = -1
    w_alu_reads = []
    reason = ''

    for i in xrange(WIN_LENGTH - 10): #site in window[: WIN_LENGTH - 10]:

        if window[i]['insert_deviations']:
            reason = 'insert_deviation'

        w_alu_reads.extend(window[i]['alu_reads'])
        if window[i]['heads'] > max_heads:
            max_heads = window[i]['heads']
            max_heads_pos = i

        if window[i]['tails'] > max_tails:
            max_tails = window[i]['tails']
            max_tails_pos = i


    if (max_heads >= MIN_PEAK and
                        #            window[max_heads_pos]['alu_heads'] >= MIN_ALU_READS_IN_PEAK and
                        max_tails >= MIN_PEAK and
                        #            window[max_tails_pos]['alu_tails'] >= MIN_ALU_READS_IN_PEAK and
                        len(set(w_alu_reads)) >= MIN_ALU_READS and
                        -2 <= (max_tails_pos - max_heads_pos) <= 10):
        reason += ' peak'

    return reason




boundary = lambda w: WIN_LENGTH <= (w[-1]['pos'] - w[0]['pos'] + 1) <= WIN_LENGTH + 3

#enough_coverage = lambda w, mean: abs(len(w[0]['reads']) - mean) <= 2*math.sqrt(mean) and abs(len(w[-1]['reads']) - mean) <= 2*math.sqrt(mean)
enough_coverage = lambda w, mean: len(w[0]['reads']) >= mean/3 and len(w[-1]['reads']) >= mean/3

alu_reads = None

def site(col):
    heads = 0
    alu_heads = 0
    tails = 0
    alu_tails = 0
    reads = []
    site_alu_reads = []
    insert_deviations = False
    for r in col.pileups:
        aln = r.alignment

        if aln.mapq >= MIN_MAPQ:
            is_alu_read = False
            qname = aln.qname
            reads.append(qname)

            if qname in alu_reads:
                site_alu_reads.append(qname)
                is_alu_read = True

            if r.is_tail:
                tails += 1
                if qname in alu_reads:
                    alu_tails += 1

            elif r.is_head:
                heads += 1
                if qname in alu_reads:
                    alu_heads += 1

            if is_alu_read and not aln.mate_is_unmapped and aln.mate_is_reverse != aln.is_reverse and aln.tid == aln.rnext and abs(aln.pos - aln.pnext) < 50:
                insert_deviations = True

    return {'pos'   : col.pos,
            'reads' : set(reads),
            'heads' : heads,
            'alu_heads' : alu_heads,
            'tails'     : tails,
            'alu_tails' : alu_tails,
            'alu_reads' : site_alu_reads,
            'insert_deviations'   : False } #insert_deviations }






if __name__ == '__main__':
    if len(sys.argv) == 1:
        error( "usage %s data-dir" % __file__)

    utils.DATA_DIR = sys.argv[1]

    not_found = json.load(open(data_d('subject_genome.fa.alu_positions.json')))
    pers_alu_info = copy.deepcopy(not_found)
    count_alus = lambda nf: sum(len(pos) for c in nf for pos in nf[c])
    total_alus = count_alus(not_found)
    false_positives = []

    alu_reads = pickle.load(open(data_d('alu_reads.pickle')))


    print data_d('bt2_pair_end.sorted.bam')
    inputf = pysam.Samfile(data_d('bt2_pair_end.sorted.bam'), 'rb')
#    inputf = pysam.Samfile(data_d('bt2_pair.in_alu.sorted.bam'))

    window = []


    total_wins = 0
    bgr_spanning = 0
    reads_per_base = 0
    heads_per_base = 0
    tails_per_base = 0
    total_bases = 0
    p_heads = 0
    p_tails = 0
    col_no = 0

    for col in inputf.pileup():
        col_no += 1
        if col_no > 100000:
            break
        s = site(col)
        window.append(s)

        reads_per_base += len(s['reads'])
        total_bases += 1

        heads_per_base += s['heads']
        tails_per_base += s['tails']

        p_heads += float(s['heads'])/ len(s['reads'])
        p_tails += float(s['tails'])/ len(s['reads'])

        if len(window) > WIN_LENGTH:
            window = window[1:]

            if boundary(window) and enough_coverage(window, 12): # coverage should be at least 4 reads
                bgr_spanning += window_stats(window)
                total_wins += 1


    window = []

    bgr_spanning = float(bgr_spanning)/total_wins
    reads_per_base = float(reads_per_base)/total_bases
    heads_per_base = float(heads_per_base)/total_bases
    tails_per_base = float(tails_per_base )/total_bases

    p_heads /= total_bases
    p_tails /= total_bases

    logm('%f\t%f\t%f\t%f' % (bgr_spanning, reads_per_base, heads_per_base, tails_per_base))
    logm('p_heads: %.4lf\tp_tails: %.4lf' % (p_heads, p_tails))

#    exit(0)

    h0_mean = bgr_spanning
    h1_mean = bgr_spanning/2
    h2_mean = 0

    fp_scores, tp_scores = [],[]

    skip_until = -1
    for col in inputf.pileup():
        if col.pos < skip_until:
            continue


        if len(window) == WIN_LENGTH:
            window = window[1:WIN_LENGTH]

        window.append(site(col))

#        if col.pos < 36265890:
#            continue
        if boundary(window) and enough_coverage(window, reads_per_base):
            reason = potential_ALU_insert(window, heads_per_base, tails_per_base)

            if reason:
                spanning = window_stats(window)

                if spanning >= h0_mean:
                    continue
    #            hyp, _ = min((None, h0_mean), ('heterozygous', h1_mean), ('homozygous', h2_mean),
    #                             key = lambda (hyp, mean): abs(spanning - mean))
                hyp, _ = min( ('heterozygous', h1_mean), ('homozygous', h2_mean),
                                 key = lambda (hyp, mean): abs(spanning - mean))

                if hyp:
                    chrom = inputf.getrname(col.tid)
                    window = []
                    skip_until = col.pos + RLEN

                    fp = True
                    for hap_no, haplotype_positions in enumerate(pers_alu_info[chrom]):

                        for inserted in haplotype_positions:
                            if abs(inserted['ref_pos'] - col.pos) < 300:
                                if inserted in not_found[chrom][hap_no]: not_found[chrom][hap_no].remove(inserted)
                                fp = False


                    if fp:
                        false_positives.append([hyp, spanning, chrom, col.pos, reason])
                        fp_scores.append(spanning)
                    else:
                        tp_scores.append(spanning)

                    logm('\t'.join(map(str, [hyp, not fp, spanning, chrom, col.pos, reason])))


    print 'False positives:\n', pformat(sorted(false_positives, key = lambda fp: (fp[-1], fp[-2])))
    print 'False negatives:\n', pformat(not_found)
    print

    print bgr_spanning, reads_per_base
    print 'True positives: %d (%.2lf%%)\t' % (len(tp_scores), float(100*len(tp_scores))/(len(tp_scores) + len(fp_scores))), mean(tp_scores), var(tp_scores)
    if fp_scores:
        print 'False positives: %d (%.2lf%%)\t' % (len(fp_scores), float(100*len(fp_scores))/(len(tp_scores) + len(fp_scores))), mean(fp_scores), var(fp_scores)

    print 'False negatives: %d(%.2f%%)\n' % (count_alus(not_found), float(100*count_alus(not_found))/total_alus)
    logm('done')
    inputf.close()


