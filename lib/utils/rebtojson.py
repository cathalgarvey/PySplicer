#!/usr/bin/env python3
# Imports and reformats REBASE enzyme lists to JSON.
# --------------------------------------------------
# Reference:
# <1>: Enzyme Name
# <2>: Prototypical enzyme (first discovered with same specificity)
# <3>: Microbe
# <4>: Source (Institution or Individual)
# <5>: Recognition Sequence (formatting information below)
# <6>: Methylation Site (if known)
# <7>: Commercial Availability (single letter code)
# <8>: Reference Numbers
# 
# Enzyme Formatting: A carat is used to indicate top-strand cleavage:
#           "ACGCG^T"
#  Alternatively, Parentheses indicate (top/bottom) cleavage, numbered
#   from the last letter of the defined recognition sequence. So:
# For example HgaI GACGC (5/10) indicates cleavage as follows:
#          5' GACGCNNNNN^      3'
#          3' CTGCGNNNNNNNNNN^ 5'
# All DNA sequences use the IUPAC Ambiguous DNA alphabet.
# Enzymes should be listed in the REBASE style as shown above, with a blank line
# between each entry: the script cleaves enzyme entries by double-line-breaks.
# References should be presented as a separate file with refs enumerated like so:
# 11. A reference, Garvey et al, etc..

import sys, json

if len(sys.argv) != 4:
    print("Usage: python3 rebtojson.py enzymefile referencesfile savefile")
    exit(1)
JSONData = {}
with open(sys.argv[1], encoding='utf-8', mode='r') as RebaseFile:
    REBASEData = RebaseFile.read()
REBASEData = REBASEData.split("\n\n") #Split into 8-entry enzyme strings.

SuppliersDict = {
            'B':        'Invitrogen Corporation (1/11)',
            'C':        'Minotech Biotechnology (1/11)',
            'E':        'Stratagene (3/10)',
            'F':        'Fermentas International Inc. (1/11)',
            'H':        'American Allied Biochemical, Inc. (4/10)',
            'I':        'SibEnzyme Ltd. (1/11)',
            'J':        'Nippon Gene Co., Ltd. (1/11)',
            'K':        'Takara Bio Inc. (1/11)',
            'M':        'Roche Applied Science (1/11)',
            'N':        'New England Biolabs (1/11)',
            'O':        'Toyobo Biochemicals (9/09)',
            'Q':        'Molecular Biology Resources - CHIMERx (1/11)',
            'R':        'Promega Corporation (7/10)',
            'S':        'Sigma Chemical Corporation (1/11)',
            'U':        'Bangalore Genei (1/11)',
            'V':        'Vivantis Technologies (11/10)',
            'X':        'EURx Ltd. (1/11)',
            'Y':        'CinnaGen Inc. (1/11)'
}

with open(sys.argv[2], encoding='utf-8', mode='r') as ReferencesFile:
    References = ReferencesFile.read()

References = References.split('\n')
RefDict = {}
for Ref in References:
    Ref = Ref.strip() # Kill whitespace
    # Fetch number by grabbing characters *until* first occurrance of period character
    FirstPeriodIndex = Ref.find('.')
    RefNum = Ref[:FirstPeriodIndex]
    RefContent = Ref[FirstPeriodIndex+1:].strip() #Get everything *after* period, minus whitespace.
    RefDict[RefNum] = RefContent

for RebEntry in REBASEData:
    Components = RebEntry.split('\n')
    ThisEnzyme = {'name':		Components[0][3:],
                  'prototype':		Components[1][3:],
                  'organism':		Components[2][3:],
                  'source':		Components[3][3:],
                  'target_site':	Components[4][3:],
                  'methylation_site':	Components[5][3:],
                  'suppliers':		[],
                  'references':		[]
                  }
    RefLine = Components[7][3:].strip()
    Refs = RefLine.split(',')
    for Ref in Refs:
        ThisEnzyme['references'].append(RefDict[Ref])

    Suppliers = Components[6][3:].strip().upper()
    for char in range(0,len(Suppliers)): # For # of times = # of chars in "suppliers"
        if Suppliers[char] in SuppliersDict.keys():
            ThisEnzyme['suppliers'].append(SuppliersDict[Suppliers[char]])

    print("Adding enzyme "+ThisEnzyme['name']+" to JSON Data.")
    JSONData[ThisEnzyme['name'].lower()] = ThisEnzyme

print("Saving JSON Data to file '"+sys.argv[3]+"'.")
with open(sys.argv[3], encoding='utf-8', mode='w') as JSONFile:
    JSONFile.write(json.dumps(JSONData))
print("Done!")
