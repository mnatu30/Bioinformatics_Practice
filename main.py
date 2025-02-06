#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 10:55:40 2025

@author: mnat
"""

from DNAToolkit import *
from Utilities import color
import random
import Structures

# Create a random DNA sequence for testing 



randDNAstr = ''.join([random.choice(Nucleotides)
                        for nuc in range(50)])

DNAStr = validateSeq(randDNAstr)

print(f'\n Sequence: {color(DNAStr)}\n')
print(f'[1] + Sequence Length: {len(DNAStr)}\n')
print(color(f'[2] + Nucleotide Frequency: {countNucFrequency(DNAStr)}\n'))
print(f'[3] + DNA/RNA Transcription: {transcription(color(DNAStr))}\n')
print(f"[4] + DNA String + Complement + Reverse Complement:\n 5' {color(DNAStr)} 3' ")
print(f"    {''.join (['|' for c in range(len(DNAStr))])}")
print(f"3' {reverse_complement(DNAStr)[::-1]} 5' [Complement]")
print(f"5' {color(reverse_complement(DNAStr))} 3' [Reverse Complement]\n")
print(f"[5] + GC Content: {gc_content(DNAStr)}%\n")
print(
      f'[6] + GC Content in Subsection k=5: {gc_content_subsec(DNAStr, k=5)}\n')
print(
      f'[7] + Amino Acids Sequence from DNA: {translate_seq(DNAStr, 0)}\n')
print(
      f'[8] + Codon Frequency (L): {codon_usage(DNAStr, "L")}\n')

print('[9] + Reading_frames:')
for frame in gen_reading_frames(DNAStr):
    print(frame)
    
print('\n[10] + All proteins from 6 ORFs:')
for prot in all_proteins_from_orfs(NM_000207_3, 0, 0, True):
    print(f'{prot}')
    



    

