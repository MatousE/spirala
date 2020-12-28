"""
This module contains all the functions for the bioloop website
"""
import random
import os
from Bio.Seq import Seq
from Bio import motifs


def calculatecg(seq):
    """"This function calculates the GC content of num_a user inputted sequence
    This is achieved by counting the number of occurances for the bases
    A, C, G & T and then applies the following formula
    Count(G + C)/Count(A + T + G + C) * 100%.
    """

    # Counts all occurrences of the letters A, C, G & T
    num_a = seq.count('A')
    num_c = seq.count('C')
    num_g = seq.count('G')
    num_t = seq.count('T')

    # Returns the gc content after applying the  following formula
    # Count(G + C)/Count(A + T + G + C) * 100%
    return str((num_g + num_c) / (num_a + num_c + num_g + num_t) * 100)


def sequencegen(length):
    """"This function generates a random DNA sequence of a length specified by
    the user.
    """
    bases = 'ACGT'
    return ''.join([random.choice(bases) for i in range(length)])


def dnamotifsearch(motif, file):
    """" This functions generates the locations of the motifs in the
    file
    """
    motif = motif.upper()

    all_motifs = []

    if not 'X' in motif:
        all_motifs.append(motif)
    else:
        no_x = motif.count('X')
        for i in range(no_x):
            for base in 'ACGT':
                all_motifs.append(motif.replace('X', base))

    instances = [Seq(i.strip()) for i in all_motifs]

    m = motifs.create(instances)

    file.save(file.filename)
    open_file = open(file.filename, 'r')

    seq = open_file.read()

    line_no = 0
    for line in seq.splitlines():
        line_no += 1
        for pos, seq in m.instances.search(line):
            print("%i %s %i" % (pos, seq, line_no))

    open_file.close()
    os.remove(file.filename)

    return 1
