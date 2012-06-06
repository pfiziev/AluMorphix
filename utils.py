import os
import random
import string
import sys
import datetime

__author__ = 'pf'


DATA_DIR = './data'

data_d = lambda path: os.path.join(DATA_DIR, path)


def error(msg):
    print >>sys.stderr, msg
    exit(1)



global_stime = datetime.datetime.now()
def elapsed(msg = None):
    print "[%s]" % msg if msg is not None else "+", "Last:" , datetime.datetime.now() - elapsed.stime, '\tTotal:', datetime.datetime.now() - global_stime

    elapsed.stime = datetime.datetime.now()

elapsed.stime = datetime.datetime.now()

BAM_MATCH = 0
BAM_INS = 1
BAM_DEL = 2
BAM_SOFTCLIP = 4



def reverse_compl_seq(strseq):
    return strseq.translate(string.maketrans("ATCG", "TAGC"))[::-1]



NUCS = ['A','C','G','T']
def mutate(seq, mut_prob = 0.0, mut_range = 0):

    for m in xrange(int(mut_prob*len(seq))):
        pos = random.choice(xrange(mut_range, len(seq) - mut_range))
        n = seq[pos]
        seq = seq[:pos] + random.choice([q for q in NUCS if q != n])+ seq[pos+1:]

    return seq

def add_snps(chrom_seq, snp_prob = 0.0001):
    parts = []
    prev_pos = 0
    for i in xrange(len(chrom_seq)):
        if random.random() < snp_prob:
            parts.append(chrom_seq[prev_pos:i])
            parts.append(random.choice([q for q in NUCS if q != chrom_seq[i]]))
            prev_pos = i + 1
    parts.append(chrom_seq[prev_pos:])
    return ''.join(parts)


def add_indels(chrom_seq, indels = 0.0001):
    error('sucks! reimplement this!')
    for m in xrange(int(indels*len(chrom_seq))):
        n = ''
        if random.random() > .5: # insertion
            n = random.choice(NUCS)

        pos = random.randint(0, len(chrom_seq) - 1)
        chrom_seq = chrom_seq[:pos] + n + chrom_seq[pos+1:]
    return chrom_seq


def orient(seq):
    strand = random.choice('+-')
    return seq if strand == '+' else reverse_compl_seq(seq), strand


def read_alu():
    with open(data_d('alu/alu.fa')) as aluf:
        alu = ''.join(l.strip().upper() for l in aluf if l[0] != '>')
    return str(alu) # for some reason it may come as unicode, which breaks other code..

def logm(message):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') , message)


from math import factorial
def binom(n,k):
    return factorial(n) / (factorial(k)*factorial(n-k))