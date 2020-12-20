"""
This module contains all the functions for the bioloop website
"""

import random


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
    seq = [random.choice(bases) for i in range(length)]
    return ''.join(seq)
