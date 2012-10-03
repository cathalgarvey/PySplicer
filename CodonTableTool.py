#!/usr/bin/env python3
'''CodonTableTool - A PySplicer tool by Cathal Garvey
   Find on GitHub: github.com/cathalgarvey/pysplicer
   Reads a flat text file containing copy/pasted data
   from the codon usage database, and converts the contents
   into PySplicer-style JSON. Exports either to stdout or
   to a specified save-file. Warning: this clobbers savefile
   if it already exists!
'''

from pysp_classes import CodonTable
from sys import exit
import argparse
import json

desc = \
'''CodonTableTool: A codon usage table import tool for pysplicer.\r\n
Imports a codon usage table from the extended format used by the http://www.kazusa.or.jp/codon/ codon usage database (read from a flat text file where it should be copy pasted by the user) into the extensible dictionary format used by pysplicer.'''

epil = \
'''CodonTableTool is part of the PySplicer toolkit by Cathal Garvey.\r\n
PySplicer is licensed under the GNU General Public License: https://www.gnu.org/licenses/gpl.html.\r\n
Code is hosted primarily on Github: https://www.github.com/cathalgarvey/pysplicer.\r\n
Huge thanks to the maintainers of the Codon Usage Database!'''

Parser = argparse.ArgumentParser(description=desc, epilog=epil, argument_default=None)
Parser.add_argument("input", help="Input file, which should contain codon usage table. If in Codon Usage Database format, it will be imported.")
Parser.add_argument("-s", "--splice", help="A pre-formatted PySplicer-style codon table file to splice with the input table during import.")
Parser.add_argument("-c", "--compromise-frequency", help="Only valid when splicing tables. Sets the minimum common frequency to accept. Default is 0.1, or 10% usage.", type=float)
Parser.add_argument("-m", "--minimum-frequency", help="Set codons from input table with a frequency lower than this threshold to 0.0 and rebalance frequencies to accomodate. When splicing, this occurs to the final table.", type=float)
Parser.add_argument("-o", "--output", help="Output file, which will be created or overwritten. Leave blank for stdout.")

Args = Parser.parse_args()
# Argparser delivers a "namespace", which is like a dictionary but more esoteric. vars gives a dictionary.
Args = vars(Args)

# ===== Make sure arguments aren't being mangled/misused =====
if Args['compromise_frequency'] and not Args['splice']:
    print("Error: -c/--compromise-frequency is only valid when splicing, but no table to splice input with was specified.")
    exit(1)
if Args['minimum_frequency']:
    if not 1.0 > Args['minimum_frequency'] > 0.0:
        print("Error: --minimum-frequency must be a float between 0 and 1.")
        exit(1)
if Args['compromise_frequency']:
    if not 1.0 > Args['compromise_frequency'] > 0.0:
        print("Error: --compromise-frequency must be a float between 0 and 1.")
        exit(1)

# ===== Create Codon Table, Read Input File =====
Table = CodonTable()

with open(Args['input']) as ImportFile:
    ImportTable = ImportFile.read()

# ===== Determine whether this is a new CUD-CUT or a pre-formatted JSON file =====
try:
    ImportTable = json.loads(ImportTable)
    Table.check_cut(ImportTable)
    Table.set_codon_table(ImportTable)
except ValueError:
    if not isinstance(ImportTable, dict):
        # Input not JSON; failed JSON loads, import CUD CUT
        Table.import_cud_cut(ImportTable)
    else:
        # Input is JSON-imported dict but fails check_cut
        print("Input file specified appears to be a JSON-formatted dict but is not a valid PySplicer Codon Usage Table.")
        exit(1)

# ===== If splicing with another table, make said other table and splice, already =====
if Args['splice']:
    Table2 = CodonTable()
    try:
        Table2.set_codon_table_from_file(Args['splice'])
        if Args['compromise_frequency']:
            compfreq = Args['compromise_frequency']
        else:
            compfreq = 0.1
        Table.splice_table(Table2.codontable, lower_threshold=compfreq)
    except Exception as e:
        print(e, "Error probably fatal, quitting.")
        exit(1)

if Args['minimum_frequency']:
    Table.purge_by_frequency(Args['minimum_frequency'])

if Args['output']:
    Table.export_codon_table_to_file(Args['output'])
else:
    Table.printout_cut()
