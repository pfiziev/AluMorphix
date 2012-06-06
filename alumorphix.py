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
MIN_PEAK = 2
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





count_alus = lambda nf: sum(len(pos) for c in nf for pos in nf[c])


def detect_insertions( h0_mean, h1_mean, h2_mean, heads_per_base, inputf, reads_per_base, tails_per_base):
    fp_scores, tp_scores = [], []
    not_found_insertions = json.load(open(data_d('subject_genome.fa.alu_positions.json')))

    if not not_found_insertions:
        return

    pers_alu_info = copy.deepcopy(not_found_insertions)
    total_alus = count_alus(not_found_insertions)
    false_positives = []
    skip_until = -1
    window = []

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
                hyp, _ = min(('heterozygous', h1_mean), ('homozygous', h2_mean),
                                key = lambda (hyp, mean): abs(spanning - mean))

                if hyp:
                    chrom = inputf.getrname(col.tid)
                    window = []
                    skip_until = col.pos + RLEN

                    fp = True
                    for hap_no, haplotype_positions in enumerate(pers_alu_info[chrom]):
                        for inserted in haplotype_positions:
                            if abs(inserted['ref_pos'] - col.pos) < 300:
                                if inserted in not_found_insertions[chrom][hap_no]: not_found_insertions[chrom][
                                                                                    hap_no].remove(inserted)
                                fp = False

                    if fp:
                        false_positives.append([hyp, spanning, chrom, col.pos, reason])
                        fp_scores.append(spanning)
                    else:
                        tp_scores.append(spanning)

                    logm('\t'.join(map(str, [hyp, not fp, spanning, chrom, col.pos, reason])))

    print 'False positives:\n', pformat(sorted(false_positives, key=lambda fp: (fp[-1], fp[-2])))
    print 'False negatives:\n', pformat(not_found_insertions)
    print
    #print bgr_spanning, reads_per_base
    print 'True positives: %d (%.2lf%%)\t' % (len(tp_scores), float(100 * len(tp_scores)) / (len(tp_scores) + len(fp_scores))), mean(tp_scores), var(tp_scores)

    if fp_scores:
        print 'False positives: %d (%.2lf%%)\t' % (

        len(fp_scores), float(100 * len(fp_scores)) / (len(tp_scores) + len(fp_scores))), mean(fp_scores), var(fp_scores)

    print 'False negatives: %d(%.2f%%)\n' % (count_alus(not_found_insertions), float(100 * count_alus(not_found_insertions)) / total_alus)



def detect_deletions(h0_mean, h1_mean, h2_mean, heads_per_base, inputf, reads_per_base, tails_per_base):
#    known_alus = json.load(open(data_d('genome/alu_positions.json')))
    known_alus = json.load(open(data_d('genome/all_known_alus.json')))

    deleted_alus = json.load(open(data_d('subject_genome.fa.deleted_alus.json')))
    not_found = set(['%s %d' % (chr_id, alu['ref_pos']) for chr_id in deleted_alus for hap in deleted_alus[chr_id] for alu in hap])

    alu_len = len(read_alu())

    total_alus = count_alus(deleted_alus)

    false_positives = []

    def detect_tails_or_heads(pileup, detect_tails = True):
        window = []
        hyp = None
        start = None
        max_heads_or_tails = 0
        for col in pileup:
            start = start or col.pos
            if col.pos > start + 200:
#                print start, col.pos, detect_tails, max_heads_or_tails
                break

            if len(window) == WIN_LENGTH:
                window = window[1:WIN_LENGTH]

            window.append(site(col))
#            print col.pos

            if boundary(window) and enough_coverage(window, reads_per_base):

                for s in window[:WIN_LENGTH - 10]:
                    max_heads_or_tails = max(s['tails' if detect_tails else 'heads'], max_heads_or_tails)

                if max_heads_or_tails >= MIN_PEAK+1:
                    spanning = window_stats(window)

                    if spanning >= h0_mean:
                        continue
                        #            hyp, _ = min((None, h0_mean), ('heterozygous', h1_mean), ('homozygous', h2_mean),
                        #                             key = lambda (hyp, mean): abs(spanning - mean))
                    hyp, _ = min(('heterozygous', h1_mean), ('homozygous', h2_mean),
                                     key = lambda (hyp, mean): abs(spanning - mean))
        return hyp

    def set_pileup_position(pos, pileup):
        col = pileup.next()
        if col.pos > pos:
            print 'cpos', col.pos, 'rewinding to', pos
            pileup = inputf.pileup(reference = inputf.gettid(chr_id), start = pos - 5)
            col = pileup.next()

        while col.pos < pos:
            col = pileup.next()
        return pileup

    for chr_id in known_alus:
        print chr_id
        to_check = sorted(known_alus[chr_id], key = lambda a: a['pos'])

        # detect homozygous
        homozygous = set()
        pileup = inputf.pileup(reference = inputf.gettid(chr_id))
        fp_keys = set()
        while to_check:
            print 'to_check:', len(to_check)
            prev_pos = 0

            for alu in list(to_check):

                if alu['pos'] - prev_pos < 200:
                    continue

                to_check.remove(alu)

                prev_pos = alu['pos'] + alu_len + 200


                pileup = set_pileup_position(alu['pos'], pileup)
                window_reads_per_base = 0
                total_bases = 0
                for col in pileup:
                    window_reads_per_base += len(site(col)['reads'])
                    total_bases += 1
                    if col.pos > alu['pos'] + alu_len:
                        break

                window_reads_per_base /= float(total_bases)
                if window_reads_per_base < reads_per_base/5:
    #                print window_reads_per_base
                    key = '%s %d' % (chr_id, alu['pos'])
                    if key in not_found:
                        print alu, 'Homozygous\tTrue Postive'
                        not_found.remove(key)
                        homozygous.add(key)
                    else:
                        print alu, 'Homozygous\tFalse Positve'
                        fp_keys.add(key)
                        false_positives.append(alu)

        elapsed('homozygous')
        print 'TP homozygous:', len(homozygous)
        print 'FP homozygous:', len(false_positives)

        pileup = inputf.pileup(reference = inputf.gettid(chr_id))

        to_check = sorted(known_alus[chr_id], key = lambda a: a['pos'])


        while to_check:
            print 'to_check:', len(to_check)
            prev_pos = 0
            for alu in list(to_check):

                if alu['pos'] - prev_pos < 200:
                    continue

                to_check.remove(alu)
                key = '%s %d' % (chr_id, alu['pos'])
                if key in homozygous or key in fp_keys:
                    continue

                prev_pos = alu['pos'] + alu_len + 200

#            print 'testing ', alu

                pileup = set_pileup_position(alu['pos'] - 100, pileup)
                hyp1 = detect_tails_or_heads(pileup, detect_tails=True)
    #            print hyp1

                pileup = set_pileup_position(alu['pos'] + alu_len - 100, pileup)
                hyp2 = detect_tails_or_heads(pileup, detect_tails=False)
    #            print hyp2

                if hyp1 and hyp2:
                    if key in not_found:
                        print alu, 'True Postive'
                        not_found.remove(key)
                    else:
                        print alu, 'False Positve'
                        false_positives.append(alu)
                else:
                    if key in not_found:
                        print alu, 'False Negative'




    print 'True positives: %d' % (total_alus - len(not_found))
    print 'False positives: %d' % len(false_positives)
    print 'False negatives: %d' % len(not_found)




if __name__ == '__main__':
    if len(sys.argv) == 1:
        error( "usage %s data-dir" % __file__)

    utils.DATA_DIR = sys.argv[1]



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

    detect_insertions(h0_mean, h1_mean, h2_mean, heads_per_base, inputf, reads_per_base, tails_per_base)
#    detect_deletions(h0_mean, h1_mean, h2_mean, heads_per_base, inputf, reads_per_base, tails_per_base)
    logm('done')
    elapsed('done')
    inputf.close()


