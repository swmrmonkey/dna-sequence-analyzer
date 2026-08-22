# DNA Sequence Analyzer

A Python program that performs basic analysis of a user-provided DNA sequence.

## Features

- Counts nucleotides and calculates nucleotide percentages
- Calculates GC and AT content
- Generates complementary and reverse-complementary DNA strands
- Transcribes DNA into RNA
- Translates RNA into a protein sequence
- Counts codon occurrences
- Searches for motifs in both the original and reverse-complementary strands
- Detects overlapping motif matches

## Requirements

- Python 3.11 or later
- No external libraries required

## Usage

Run the program from the project directory:

```bash
python3 main.py
```

Enter a DNA sequence containing only `A`, `C`, `G`, and `T` when prompted. The program will analyze the sequence and then prompt for a DNA motif to search for.

Example:

```text
Enter a DNA sequence: ATGGCTTTTGGACCGAAATGCTAA
Enter a motif to search for: ATG
```

## Notes

- DNA input is treated as a coding strand in the 5' to 3' direction.
- Translation begins at the first nucleotide rather than searching for a start codon.
- Translation stops at the first stop codon.
- Incomplete codons at the end of a sequence are ignored.
- Motif positions are reported using 1-based indexing.
- Reverse-complementary motif positions refer to positions within the reverse-complementary strand.