'''A set of weighted reverse translation tables for use with PySplicer.
User-defined tables can be added to this list in JSON format, though it may
be more convenient to simply specify such tables as separate files with the -T
switch in the pysplicer script.'''

# Codon frequency/bias tables from the Creative Commons licensed paper:
# "Design Parameters to Control Synthetic Gene Expression in Escherichia coli", by
# Mark Welch, Sridhar Govindarajan, Jon E. Ness, Alan Villalobos, Austin Gurney,
# Jeremy Minshull, Claes Gustafsson, PLoS, 2009.
# From figure: "Two codon frequency tables were constructed as averages of the
# frequencies found in the best variants:
# FreqA from the 4 genes comprised of the 2 most highly expressing variants of
# each set (A1, A21, P19, and P20);
# FreqB from 10 of the most highly expressing variants (P19, P20, A1, A21,
# A1_14, A17_17_1, A17_1_1, A1_17_1, A1_11_11, and A_FreqA) and
# FreqC approximates the bias used to create polymerase variant P19.
# A fourth set of frequencies, “HiCAI” used codons that are most common in
# highly expressed native E. coli genes."

welch_freqa = {'A': {'GCA': 0.15, 'GCC': 0.18, 'GCG': 0.51, 'GCT': 0.15},
 'C': {'TGC': 0.65, 'TGT': 0.35},
 'D': {'GAC': 0.59, 'GAT': 0.4},
 'E': {'GAA': 0.41, 'GAG': 0.59},
 'F': {'TTC': 0.58, 'TTT': 0.42},
 'G': {'GGA': 0.01, 'GGC': 0.41, 'GGG': 0.01, 'GGT': 0.57},
 'H': {'CAC': 0.58, 'CAT': 0.42},
 'I': {'ATA': 0.0, 'ATC': 0.46, 'ATT': 0.53},
 'K': {'AAA': 0.44, 'AAG': 0.56},
 'L': {'CTA': 0.0,
  'CTC': 0.05,
  'CTG': 0.61,
  'CTT': 0.04,
  'TTA': 0.06,
  'TTG': 0.24},
 'M': {'ATG': 1.0},
 'N': {'AAC': 0.56, 'AAT': 0.44},
 'P': {'CCA': 0.15, 'CCC': 0.0, 'CCG': 0.72, 'CCT': 0.13},
 'Q': {'CAA': 0.51, 'CAG': 0.48},
 'R': {'AGA': 0.09,
  'AGG': 0.01,
  'CGA': 0.0,
  'CGC': 0.4,
  'CGG': 0.0,
  'CGT': 0.49},
 'S': {'AGC': 0.38,
  'AGT': 0.01,
  'TCA': 0.02,
  'TCC': 0.25,
  'TCG': 0.21,
  'TCT': 0.12},
 'T': {'ACA': 0.0, 'ACC': 0.63, 'ACG': 0.25, 'ACT': 0.11},
 'V': {'GTA': 0.08, 'GTC': 0.22, 'GTG': 0.38, 'GTT': 0.31},
 'W': {'TGG': 1.0},
 'Y': {'TAC': 0.57, 'TAT': 0.42},
 '*': {'TAA':0.64,'TAG':0.0,'TGA':0.36}}

welch_freqb = {'A': {'GCA': 0.24, 'GCC': 0.19, 'GCG': 0.44, 'GCT': 0.12},
 'C': {'TGC': 0.58, 'TGT': 0.42},
 'D': {'GAC': 0.54, 'GAT': 0.46},
 'E': {'GAA': 0.43, 'GAG': 0.57},
 'F': {'TTC': 0.55, 'TTT': 0.45},
 'G': {'GGA': 0.0, 'GGC': 0.39, 'GGG': 0.0, 'GGT': 0.6},
 'H': {'CAC': 0.62, 'CAT': 0.38},
 'I': {'ATA': 0.0, 'ATC': 0.51, 'ATT': 0.49},
 'K': {'AAA': 0.51, 'AAG': 0.49},
 'L': {'CTA': 0.0,
  'CTC': 0.03,
  'CTG': 0.78,
  'CTT': 0.02,
  'TTA': 0.03,
  'TTG': 0.14},
 'M': {'ATG': 1.0},
 'N': {'AAC': 0.53, 'AAT': 0.47},
 'P': {'CCA': 0.1, 'CCC': 0.0, 'CCG': 0.81, 'CCT': 0.09},
 'Q': {'CAA': 0.45, 'CAG': 0.55},
 'R': {'AGA': 0.03,
  'AGG': 0.0,
  'CGA': 0.0,
  'CGC': 0.35,
  'CGG': 0.0,
  'CGT': 0.62},
 'S': {'AGC': 0.68,
  'AGT': 0.0,
  'TCA': 0.02,
  'TCC': 0.13,
  'TCG': 0.05,
  'TCT': 0.1},
 'T': {'ACA': 0.0, 'ACC': 0.57, 'ACG': 0.33, 'ACT': 0.1},
 'V': {'GTA': 0.03, 'GTC': 0.28, 'GTG': 0.34, 'GTT': 0.35},
 'W': {'TGG': 1.0},
 'Y': {'TAC': 0.58, 'TAT': 0.42},
 '*': {'TAA':0.64,'TAG':0.0,'TGA':0.36}}

welch_freqc = {'A': {'GCA': 0.0, 'GCC': 0.0, 'GCG': 0.73, 'GCT': 0.26},
 'C': {'TGC': 1.0, 'TGT': 0.0},
 'D': {'GAC': 0.61, 'GAT': 0.38},
 'E': {'GAA': 0.36, 'GAG': 0.63},
 'F': {'TTC': 0.5, 'TTT': 0.5},
 'G': {'GGA': 0.0, 'GGC': 0.43, 'GGG': 0.0, 'GGT': 0.56},
 'H': {'CAC': 0.5, 'CAT': 0.5},
 'I': {'ATA': 0.0, 'ATC': 0.77, 'ATT': 0.22},
 'K': {'AAA': 0.35, 'AAG': 0.64},
 'L': {'CTA': 0.0,
  'CTC': 0.18,
  'CTG': 0.27,
  'CTT': 0.0,
  'TTA': 0.11,
  'TTG': 0.43},
 'M': {'ATG': 1.0},
 'N': {'AAC': 0.64, 'AAT': 0.35},
 'P': {'CCA': 0.13, 'CCC': 0.0, 'CCG': 0.86, 'CCT': 0.0},
 'Q': {'CAA': 0.53, 'CAG': 0.46},
 'R': {'AGA': 0.0, 'AGG': 0.0, 'CGA': 0.0, 'CGC': 0.2, 'CGG': 0.0, 'CGT': 0.8},
 'S': {'AGC': 0.16,
  'AGT': 0.0,
  'TCA': 0.0,
  'TCC': 0.58,
  'TCG': 0.04,
  'TCT': 0.2},
 'T': {'ACA': 0.0, 'ACC': 0.91, 'ACG': 0.0, 'ACT': 0.08},
 'V': {'GTA': 0.0, 'GTC': 0.0, 'GTG': 0.4, 'GTT': 0.59},
 'W': {'TGG': 1.0},
 'Y': {'TAC': 0.57, 'TAT': 0.42},
 '*': {'TAA':0.64,'TAG':0.0,'TGA':0.36}}

welch_hicai = {'A': {'GCA': 0.11, 'GCC': 0.0, 'GCG': 0.44, 'GCT': 0.44},
 'C': {'TGC': 0.75, 'TGT': 0.25},
 'D': {'GAC': 0.58, 'GAT': 0.41},
 'E': {'GAA': 1.0, 'GAG': 0.0},
 'F': {'TTC': 0.55, 'TTT': 0.44},
 'G': {'GGA': 0.0, 'GGC': 0.45, 'GGG': 0.0, 'GGT': 0.54},
 'H': {'CAC': 0.66, 'CAT': 0.33},
 'I': {'ATA': 0.0, 'ATC': 0.58, 'ATT': 0.41},
 'K': {'AAA': 1.0, 'AAG': 0.0},
 'L': {'CTA': 0.0, 'CTC': 0.0, 'CTG': 1.0, 'CTT': 0.0, 'TTA': 0.0, 'TTG': 0.0},
 'M': {'ATG': 1.0},
 'N': {'AAC': 1.0, 'AAT': 0.0},
 'P': {'CCA': 0.0, 'CCC': 0.0, 'CCG': 1.0, 'CCT': 0.0},
 'Q': {'CAA': 0.07, 'CAG': 0.92},
 'R': {'AGA': 0.0,
  'AGG': 0.0,
  'CGA': 0.0,
  'CGC': 0.38,
  'CGG': 0.0,
  'CGT': 0.61},
 'S': {'AGC': 0.02,
  'AGT': 0.0,
  'TCA': 0.0,
  'TCC': 0.45,
  'TCG': 0.0,
  'TCT': 0.51},
 'T': {'ACA': 0.0, 'ACC': 0.65, 'ACG': 0.0, 'ACT': 0.34},
 'V': {'GTA': 0.0, 'GTC': 0.0, 'GTG': 0.41, 'GTT': 0.58},
 'W': {'TGG': 1.0},
 'Y': {'TAC': 0.61, 'TAT': 0.38},
 '*': {'TAA':0.64,'TAG':0.0,'TGA':0.36}}

ecoli_optimal = welch_freqb