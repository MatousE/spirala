"""
This module contains all the functions for the bioloop website
"""
import random
import os
import re
from Bio.Seq import Seq
from Bio import motifs
from Bio.SeqUtils import six_frame_translations


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

    return ''.join([random.choice('ACGT') for _ in range(length)])


def dnamotifsearch(motif, file):
    """Returns list of strings
    This functions generates the locations of the motifs in the
    file
    """
    # Motif is turned to caps
    motif = motif.upper()

    all_motifs = []
    # If the motif has no x's it will by default append the all_motif
    if 'X' not in motif:
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
    motifs_instance = motifs.create(instances)

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
        for pos, seq in motifs_instance.instances.search(line):
            # Stores the line, position and sequence in motif_occurrences
            motif_occurrences.append(f"Line: {line_no} Position: {pos} Seq: {seq}")

    # Closes the openfile as it is no longer needed
    open_file.close()
    # Removes the file from the os as it is no longer needed
    os.remove(file.filename)

    return motif_occurrences


def proteinmotifsearch(imp, file):
    """Returns a list of positions within the FASTA file
        where a user inputted motif has occurred
    """
    imp = imp.upper()
    instances = []
    if 'X' not in imp:
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


def codonusage(file):
    """ Counts the codons and its associative
        amino acid and gives a percentage for
        that codon
    """
    codon_dict = {
        "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0,
        "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0,
        "ATT": 0, "ATC": 0, "ATA": 0, "ATG": 0,
        "GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0,
        "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
        "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0,
        "AAT": 0, "AAC": 0, "AAA": 0, "AAG": 0,
        "GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
        "TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0,
        "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
        "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0,
        "GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0,
        "TGT": 0, "TGC": 0, "TGA": 0, "TGG": 0,
        "CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0,
        "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
        "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}

    synonymous_codons = {
        "CYS": ["TGT", "TGC"],
        "ASP": ["GAT", "GAC"],
        "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
        "GLN": ["CAA", "CAG"],
        "MET": ["ATG"],
        "ASN": ["AAC", "AAT"],
        "PRO": ["CCT", "CCG", "CCA", "CCC"],
        "LYS": ["AAG", "AAA"],
        "STOP": ["TAG", "TGA", "TAA"],
        "THR": ["ACC", "ACA", "ACG", "ACT"],
        "PHE": ["TTT", "TTC"],
        "ALA": ["GCA", "GCC", "GCG", "GCT"],
        "GLY": ["GGT", "GGG", "GGA", "GGC"],
        "ILE": ["ATC", "ATA", "ATT"],
        "LEU": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
        "HIS": ["CAT", "CAC"],
        "ARG": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
        "TRP": ["TGG"],
        "VAL": ["GTA", "GTC", "GTG", "GTT"],
        "GLU": ["GAG", "GAA"],
        "TYR": ["TAT", "TAC"],
    }

    file.save(file.filename)

    codon_count = codon_dict.copy()

    with open(file.filename) as handle:
        dna_sequence = handle.read().splitlines()
        for seq in dna_sequence:
            if not seq[0] == '>':
                for i in range(0, len(seq), 3):
                    codon = seq[i: i + 3]
                    if codon in codon_count:
                        codon_count[codon] += 1

    total = 0.0
    for codon in codon_count:
        total += codon_count[codon]

    return_vals = ["AmAcid   Codon    Number   Fraction"]

    for amino_acid in synonymous_codons:

        codons = synonymous_codons[amino_acid]

        for codon in codons:
            return_vals.append(f"{amino_acid}        {codon}        {float(codon_count[codon])}           {float(codon_count[codon] / total * 100)}")
        return_vals.append("")

    os.remove(file.filename)

    return return_vals


def sixframetranslation(seq):
    """ Returns a string showing the translation map
        of a user inputted RNA sequence
    """
    return six_frame_translations(seq)


def cpgisland(file, search_window, search_frame):
    """ Returns a list of location and gc content of a
        frame of the given window size
    """
    file.save(file.filename)

    cpg_island = []

    with open(file.filename) as handle:
        data = handle.read()

        start = 0
        end = search_frame
        seq = data[:search_window]

        while start <= len(seq)-end:
            frame = seq[start:end]
            c_count = frame.count('C')
            g_count = frame.count('G')
            cg_count = frame.count('CG')

            cpg_ratio = 0
            gc_content = float((c_count + g_count) / search_frame)
            if c_count != 0 and g_count != 0:
                cpg_ratio = cg_count / ((c_count * g_count) / search_frame)
            if gc_content >= 0.5 and cpg_ratio >= 0.6:
                cpg_island.append(f"CpG island detected in region {start} to {end} (Obs/Exp = {round(cpg_ratio * 100, 2)} and %GC = {round(gc_content*100, 2)})")
            start += 1
            end += 1
            
    os.remove(file.filename)

    return cpg_island
