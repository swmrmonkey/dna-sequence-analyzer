# DNA Sequence Analyzer
# Analyzes a user-provided DNA sequence, including nucleotide composition,
# complementary strands, transcription, translation, codon usage, and motif searching.

CODON_TABLE = {
    # Phenylalanine (F) / Leucine (L)
    "UUU": "F", "UUC": "F",
    "UUA": "L", "UUG": "L",

    # Serine (S)
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",

    # Tyrosine (Y) / Stop
    "UAU": "Y", "UAC": "Y",
    "UAA": "STOP", "UAG": "STOP",

    # Cysteine (C) / Tryptophan (W) / Stop
    "UGU": "C", "UGC": "C",
    "UGA": "STOP", "UGG": "W",

    # Leucine (L)
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",

    # Proline (P)
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",

    # Histidine (H) / Glutamine (Q)
    "CAU": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",

    # Arginine (R)
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    # Isoleucine (I) / Methionine (M)
    "AUU": "I", "AUC": "I", "AUA": "I",
    "AUG": "M",

    # Threonine (T)
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",

    # Asparagine (N) / Lysine (K)
    "AAU": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",

    # Serine (S) / Arginine (R)
    "AGU": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",

    # Valine (V)
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",

    # Alanine (A)
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",

    # Aspartic acid (D) / Glutamic acid (E)
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",

    # Glycine (G)
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


def get_dna_sequence():
    """Prompt the user until a valid DNA sequence containing only A, C, G, and T is entered."""
    not_valid = True

    while not_valid:
        dna = input("Enter a DNA sequence: ").upper()

        if not dna:
            print("No DNA sequence entered. Please enter a valid DNA sequence.")
            continue

        not_valid = False

        for nucleotide in dna:
            if nucleotide not in "ACGT":
                print("Invalid DNA sequence. Please enter a sequence containing only A, C, G, and T.")
                not_valid = True
                break

    return dna


def nucleotide_count(dna):
    """Count the occurrences of each nucleotide in a DNA sequence."""
    counts = {nucleotide: dna.count(nucleotide) for nucleotide in "ACGT"}
    return counts


def calculate_nucleotide_percentages(counts, sequence_length):
    """Calculate the percentage of the sequence represented by each nucleotide."""
    percentages = {nucleotide: (count/sequence_length) * 100 for nucleotide, count in counts.items()}
    return percentages


def gc_and_at_content(percentages):
    """Calculate GC and AT content from nucleotide percentages."""
    gc_content = percentages['G'] + percentages['C']
    at_content = percentages['A'] + percentages['T']
    return gc_content, at_content


def complementary_strand(dna):
    """Generate the complementary strand of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complementary_dna = "".join(complement[nucleotide] for nucleotide in dna)
    return complementary_dna


def reverse_complementary_strand(dna):
    """Generate the reverse complement of a DNA sequence."""
    complementary_dna = complementary_strand(dna)
    reverse_complementary_dna = complementary_dna[::-1]
    return reverse_complementary_dna


def transcribe_dna_to_rna(dna):
    """Transcribe a coding DNA sequence into RNA."""
    rna = dna.replace('T', 'U')
    return rna


def translate_rna_to_protein(rna):
    """Translate an RNA sequence into a protein sequence using the standard genetic code.
    
    Translation begins at the first nucleotide of the RNA sequence and does not
    search for a start codon. Translation ends at the first stop codon or at the
    end of the sequence. Incomplete final codons are ignored.
    """
    protein = ""

    # Process the RNA sequence one codon at a time.
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]

        # Ignore an incomplete codon at the end of the sequence.
        if len(codon) < 3:
            break

        if codon in CODON_TABLE:
            amino_acid = CODON_TABLE[codon]

            # Translation terminates at the first stop codon.
            if amino_acid == 'STOP':
                break

            protein += amino_acid

    return protein


def codon_count(rna):
    """Count complete codons in an RNA sequence according to its current reading frame."""
    codon_counts = {}

    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]

        # Do not count an incomplete final codon.
        if len(codon) < 3:
            break

        codon_counts[codon] = codon_counts.get(codon, 0) + 1

    return codon_counts


def get_valid_motif():
    """Prompt the user until a valid DNA motif containing only A, C, G, and T is entered."""
    not_valid = True

    while not_valid:
        motif = input("Enter a motif to search for: ").upper()

        if not motif:
            print("No motif entered. Please enter a valid motif.")
            continue

        not_valid = False

        for nucleotide in motif:
            if nucleotide not in "ACGT":
                print("Invalid motif. Please enter a motif containing only A, C, G, and T.")
                not_valid = True
                break

    return motif


def motif_search(dna, motif):
    """Return all motif locations using 1-based biological sequence positions."""
    positions = []

    for i in range(len(dna) - len(motif) + 1):
        if dna[i:i + len(motif)] == motif:
            positions.append(i + 1)

    return positions


def main():
    # Perform sequence analyses.
    dna = get_dna_sequence()
    sequence_length = len(dna)
    counts = nucleotide_count(dna)
    percentages = calculate_nucleotide_percentages(counts, sequence_length)
    gc_content, at_content = gc_and_at_content(percentages)

    complementary_dna = complementary_strand(dna)
    reverse_complementary_dna = reverse_complementary_strand(dna)

    rna = transcribe_dna_to_rna(dna)
    protein = translate_rna_to_protein(rna)
    codon_counts = codon_count(rna)

    motif = get_valid_motif()
    motif_positions = motif_search(dna, motif)
    rev_complementary_motif_pos = motif_search(reverse_complementary_dna, motif)

    # Display general sequence information.
    print("DNA sequence:", dna)
    print(f"Length of the DNA sequence: {sequence_length} nt")

    print("----Nucleotide counts----")
    for nucleotide, count in counts.items():
        print(f"Count of {nucleotide}: {count}")

    print("----Nucleotide percentages----")
    for nucleotide, percentage in percentages.items():
        print(f"Percentage of {nucleotide}: {percentage:.2f}%")

    print("----GC and AT content----")
    print(f"GC Content: {gc_content:.2f}%")
    print(f"AT Content: {at_content:.2f}%")

    # Display derived DNA, RNA, and protein sequences.
    print("----Derived sequences----")
    print("Complementary DNA strand:", complementary_dna)
    print("Reverse Complementary DNA strand:", reverse_complementary_dna)
    print("Transcribed RNA sequence:", rna)
    print("Translated Protein sequence:", protein)

    print("----Codon counts in the RNA sequence----")
    for codon, count in codon_counts.items():
        print(f"Count of {codon}: {count}")

    # Report motif matches on the forward and reverse-complementary strands.
    print("----Motif search results----")
    if motif_positions and rev_complementary_motif_pos:
        print(f"Motif found at positions: {', '.join(map(str, motif_positions))}")
        print(f"Motif found in reverse complementary strand at positions: {', '.join(map(str, rev_complementary_motif_pos))}")

    elif motif_positions:
        print(f"Motif found at positions: {', '.join(map(str, motif_positions))}")

    elif rev_complementary_motif_pos:
        print(f"Motif found in reverse complementary strand at positions: {', '.join(map(str, rev_complementary_motif_pos))}")

    else:
        print("Motif not found in either the DNA sequence or its reverse complementary strand.")


main()