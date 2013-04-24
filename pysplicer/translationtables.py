#!/usr/bin/env python3
'''A set of codon table objects providing a bidirectional way to translate or reverse-translate from DNA to Amino sequences.'''
table1 = {
    "starts": ["TTG", "CTG", "ATG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAG", "TAA", "TGA"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "R", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "R", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "*", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table2 = {
    "starts": ["ATT", "ATC", "ATA", "ATG", "GTG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["AGG", "AGA", "TAG", "TAA"], "M": ["ATA", "ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG", "TGA"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "M", "AGG": "*", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "*", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "W", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table3 = {
    "starts": ["ATA", "ATG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAG", "TAA"], "M": ["ATA", "ATG"], "L": ["TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"], "T": ["CTT", "ACC", "ACA", "ACT", "CTG", "CTA", "CTC", "ACG"], "W": ["TGG", "TGA"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "T", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "M", "AGG": "R", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "R", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "T", "CTA": "T", "CTC": "T", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "W", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table4 = {
    "starts": ["TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAG", "TAA"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG", "TGA"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "R", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "R", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "W", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table5 = {
    "starts": ["TTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAG", "TAA"], "M": ["ATA", "ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGG", "AGC", "AGA", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG", "TGA"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "M", "AGG": "S", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "S", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "W", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table6 = {
    "starts": ["ATG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TGA"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG", "TAG", "TAA"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "R", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "R", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "Q", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "Q", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "*", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table9 = {
    "starts": ["ATG", "GTG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAG"], "*": ["TAG", "TAA"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAA", "AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGG", "AGC", "AGA", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG", "TGA"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "N", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "S", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "S", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "W", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table10 = {
    "starts": ["ATG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC", "TGA"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAG", "TAA"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "R", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "R", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "C", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table11 = {
    "starts": ["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAG", "TAA", "TGA"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "R", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "R", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "*", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table12 = {
    "starts": ["CTG", "ATG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAG", "TAA", "TGA"], "M": ["ATG"], "L": ["CTT", "CTA", "CTC", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "CTG", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "R", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "R", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "S", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "*", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table13 = {
    "starts": ["TTG", "ATA", "ATG", "GTG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["AGG", "AGA", "GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAG", "TAA"], "M": ["ATA", "ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG", "TGA"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "M", "AGG": "G", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "G", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "W", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table14 = {
    "starts": ["ATG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAG"], "*": ["TAG"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAA", "AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGG", "AGC", "AGA", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG", "TGA"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC", "TAA"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "N", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "S", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "S", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "Y", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "W", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table15 = {
    "starts": ["ATG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAA", "TGA"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG", "TAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "R", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "R", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "Q", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "*", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table16 = {
    "starts": ["ATG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAA", "TGA"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TAG", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "R", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "R", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "L", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "*", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table21 = {
    "starts": ["ATG", "GTG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATT"], "H": ["CAT", "CAC"], "K": ["AAG"], "*": ["TAG", "TAA"], "M": ["ATA", "ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAA", "AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGG", "AGC", "AGA", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG", "TGA"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "N", "ATC": "I", "AAC": "N", "ATA": "M", "AGG": "S", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "S", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "W", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table22 = {
    "starts": ["ATG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAA", "TGA", "TCA"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TAG", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC"], "R": ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "R", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "R", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "L", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "*", "GAC": "D", "TCC": "S", "TCA": "*", "ATG": "M", "CGC": "R"}}
table23 = {
    "starts": ["ATT", "ATG", "GTG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AAG"], "*": ["TAG", "TTA", "TAA", "TGA"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "R", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "R", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "*", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "*", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}
table24 = {
    "starts": ["TTG", "CTG", "ATG", "GTG"],
    "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"], "C": ["TGT", "TGC"], "E": ["GAA", "GAG"], "D": ["GAT", "GAC"], "G": ["GGT", "GGG", "GGA", "GGC"], "F": ["TTT", "TTC"], "I": ["ATC", "ATA", "ATT"], "H": ["CAT", "CAC"], "K": ["AAA", "AGG", "AAG"], "*": ["TAG", "TAA"], "M": ["ATG"], "L": ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N": ["AAC", "AAT"], "Q": ["CAA", "CAG"], "P": ["CCT", "CCA", "CCG", "CCC"], "S": ["AGC", "AGA", "AGT", "TCT", "TCG", "TCC", "TCA"], "R": ["CGA", "CGG", "CGT", "CGC"], "T": ["ACC", "ACA", "ACT", "ACG"], "W": ["TGG", "TGA"], "V": ["GTA", "GTC", "GTG", "GTT"], "Y": ["TAT", "TAC"]},
    "codons": {"CTT": "L", "ACC": "T", "ACA": "T", "AAA": "K", "ATC": "I", "AAC": "N", "ATA": "I", "AGG": "K", "CCT": "P", "ACT": "T", "AGC": "S", "AAG": "K", "AGA": "S", "CAT": "H", "AAT": "N", "ATT": "I", "CTG": "L", "CTA": "L", "CTC": "L", "CAC": "H", "ACG": "T", "CAA": "Q", "AGT": "S", "CCA": "P", "CCG": "P", "CCC": "P", "TAT": "Y", "GGT": "G", "TGT": "C", "CGA": "R", "CAG": "Q", "TCT": "S", "GAT": "D", "CGG": "R", "TTT": "F", "TGC": "C", "GGG": "G", "TAG": "*", "GGA": "G", "TGG": "W", "GGC": "G", "TAC": "Y", "TTC": "F", "TCG": "S", "TTA": "L", "TTG": "L", "CGT": "R", "GAA": "E", "TAA": "*", "GCA": "A", "GTA": "V", "GCC": "A", "GTC": "V", "GCG": "A", "GTG": "V", "GAG": "E", "GTT": "V", "GCT": "A", "TGA": "W", "GAC": "D", "TCC": "S", "TCA": "S", "ATG": "M", "CGC": "R"}}

# Aliases
# 1   The Standard Code
standard = table1
# 2   The Vertebrate Mitochondrial Code
vertebrate_mitochondrial = table2
# 3   The Yeast Mitochondrial Code
yeast_mitochondrial = table3
# 4   The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
mold_mitochondrial = table4
protozoan_mitochondrial = table4
coelenterate_mitochondrial = table4
# 5   The Invertebrate Mitochondrial Code
invertebrate_mitochondrial = table5
# 6   The Ciliate, Dasycladacean and Hexamita Nuclear Code
ciliate = table6
dasycladacean = table6
hexamita = table6
# 9   The Echinoderm and Flatworm Mitochondrial Code
echinoderm_mitochondrial = table9
flatworm_mitochondrial = table9
# 10   The Euplotid Nuclear Code
euplotid = table10
# 11   The Bacterial, Archaeal and Plant Plastid Code
bacterial = table11
archaeal = table11
plant_plastid = table11
chloroplast = table11
# 12   The Alternative Yeast Nuclear Code
alternative_yeast = table12
alt_yeast = table12
# 13   The Ascidian Mitochondrial Code
ascidian_mitochondrial = table13
# 14   The Alternative Flatworm Mitochondrial Code
alternative_flatworm_mitochondrial = table14
alt_flatworm_mitochondrial = table14
# 15   Blepharisma Nuclear Code
blepharisma = table15
# 16   Chlorophycean Mitochondrial Code
chlorophycean_mitochondrial = table16
# 21   Trematode Mitochondrial Code
trematode_mitochondrial = table21
# 22   Scenedesmus Obliquus Mitochondrial Code
scenedesmus_obliquus_mitochondrial = table22
# 23   Thraustochytrium Mitochondrial Code
thraustochytrium_mitochondrial = table23
# 24   Pterobranchia Mitochondrial Code
pterobranchia_mitochondrial = table24