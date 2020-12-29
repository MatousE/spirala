"""
This module contains all the functions for the bioloop website
"""
import random
import os
from Bio.Seq import Seq
from Bio import motifs
import re

def calculatecg(seq):
    """This function calculates the GC content of num_a user inputted sequence
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
    """Returns a string of random bases
    This function generates a random DNA sequence of a length specified by
    the user.
    """

    return ''.join([random.choice('ACGT') for i in range(length)])


def dnamotifsearch(motif, file):
    """Returns list of strings
    This functions generates the locations of the motifs in the
    file
    """
    # Motif is turned to caps
    motif = motif.upper()

    all_motifs = []
    # If the motif has no x's it will by default append the all_motif
    if not 'X' in motif:
        all_motifs.append(motif)
    # If the motif has a x will create the possible motifs where x can be
    # A C G T
    else:
        for base in 'A C G T'.split():
            all_motifs.append(motif.replace('X', base, 1))

        for motiv in all_motifs:
            if 'X' in motiv:
                for base in 'A C G T'.split():
                    all_motifs.append(motiv.replace('X', base, 1))

        to_remove = []

        for mot in all_motifs:
            if 'X' in mot:
                to_remove.append(mot)

        for remove in to_remove:
            all_motifs.remove(remove)

    # Removes spaces at start of the motifs and converts then to a Seq
    instances = [Seq(i.strip()) for i in all_motifs]

    # Creates instance of Bio.motifs with the given instances
    m = motifs.create(instances)

    # Temporarily saves the file to the os
    file.save(file.filename)
    # Opens the file for read format
    open_file = open(file.filename, 'r')

    # Reads in the contents of the file into seq
    seq = open_file.read()
    motif_occurrences = []
    line_no = 0
    # Checking each line individually for the motifs in m
    for line in seq.splitlines():
        # Increments the line
        line_no += 1
        # Searching for the instances in the line
        for pos, seq in m.instances.search(line):
            # Stores the line, position and sequence in motif_occurrences
            motif_occurrences.append(f"Line: {line_no} Position: {pos} Seq: {seq}")

    # Closes the openfile as it is no longer needed
    open_file.close()
    # Removes the file from the os as it is no longer needed
    os.remove(file.filename)

    return motif_occurrences


def proteinmotifsearch(imp, file):
    imp = imp.upper()
    instances = []
    if not 'X' in imp:
        instances.append(imp.strip())
    else:
        for base in 'G A L M F W K Q E S P V I C Y H R N D T'.split():
            instances.append(imp.replace('X', base, 1))

        for motiv in instances:
            if 'X' in motiv:
                for base in 'G A L M F W K Q E S P V I C Y H R N D T'.split():
                    instances.append(motiv.replace('X', base, 1))

        to_remove = []

        [to_remove.append(mot) for mot in instances if 'X' in mot]

        [instances.remove(rem) for rem in to_remove]

        instances = [mo.lstrip() for mo in instances]

    # Temporarily saves the file to the os
    file.save(file.filename)
    # Opens the file for read format
    open_file = open(file.filename, 'r')

    # Reads in the contents of the file into seq
    seq = open_file.read()

    motif_occurrences = []
    motif_locations = []
    return_vals = []

    for instance in instances:
        for line in seq.splitlines():
            if instance in line:
                motif_locations.append([m.start() for m in re.finditer(instance, seq)])
                motif_occurrences.append(instance)

    for seq in motif_occurrences:
        for loc in motif_locations:
            return_vals.append(f"Pos: {str(loc)} Motif: {seq}")

    open_file.close()

    os.remove(file.filename)

    return return_vals
