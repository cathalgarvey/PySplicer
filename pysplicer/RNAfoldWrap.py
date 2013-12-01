#!/usr/bin/env python3
'''A dumb wrapper for RNAfold from the Vienna RNA Package for Python 3.
Uses subprocess.Popen to open a stdin/stdout/stderr stream to RNAfold, with
the --noPS option, allowing inline use of RNAfold on RNA strings with output
wrapped in a simple object to separate sequence, energy and structure.
'''
import collections
import subprocess
import random
# Warn if Vienna not found.
from shutil import which
if not which("RNAfold"):
    from sys import stderr
    print("WARNING: PySplicer could not locate the RNAfold binary from the ViennaRNA package. Without ViennaRNA installed, structural issues cannot be diagnosed or resolved, and DNA/RNA instability may result. To put it bluntly, your designs may suck. Install ViennaRNA for best results (or any results at all)!", file=stderr)
del(which)

# Used for testing:
randseq = lambda n: ''.join([random.choice("ACGU") for x in range(0,n)])

RNAStructure = collections.namedtuple("RNAStructure",["structure","energy"])

class RNAFoldError(BaseException):
    'Used to wrap and raise messages from stderr when calling RNAfold.'
    pass

class RNAFoldOutput:
    'Wraps the two-line output from RNAfold and extracts sequence, structure and energy.'
    def __init__(self, rnafold_output):
        output_lines = rnafold_output.strip().splitlines()
        self.sequence = output_lines[0]
        structure = output_lines[1].split(None,1)[0].strip()
        energy = float(output_lines[1].rsplit("(",1)[1].strip("()").strip())
        self.folding = RNAStructure(structure, energy)

def RNAfold(sequence, *args):
    # Note that RNAfold auto-converts "T" to "U" so this is unnecessary in Python
    # This behaviour can be overridden if calculation of DNA is desired.
    rnaf = subprocess.Popen(["RNAfold","--noPS"] + list(args),
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            # Universal Newlines effectively allows string IO.
                            universal_newlines=True)
    foldout, folderr = rnaf.communicate(sequence)
    if folderr:
        raise RNAFoldErr(folderr)
    output = RNAFoldOutput(foldout)
    #print("Debug: Energy was", output.folding.energy,"- RNAfold output was:\n",foldout) # Debug
    return output

class RNASuboptError(BaseException):
    'Used to wrap and raise messages from stderr when calling RNAsubopt.'
    pass

class RNASuboptOutput:
    'Wraps the muti-line output from RNAsubopt and extracts sequence, structures and energies.'
    def __init__(self, rnafold_output):
        output_lines = rnafold_output.strip().splitlines()
        # Top line of RNAsubopt has sequence and two numbers; split, take first.
        self.sequence = output_lines.pop(0).strip().split()[0]
        self.foldings = []
        for structure in output_lines:
            # Create a namedtuple of structure/energy by splitting line.
            structure, energy = structure.strip().split()
            newstructure = RNAStructure(structure, float(energy))
            self.foldings.append(newstructure)

def RNAsubopt(sequence, *args):
    rnaf = subprocess.Popen(["RNAsubopt"] + list(args),
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            # Universal Newlines effectively allows string IO.
                            universal_newlines=True)        
    foldout, folderr = rnaf.communicate(sequence)
    if folderr:
        raise RNASuboptError(folderr)
    return RNASuboptOutput(foldout)

# Base function to filter by free energy of folding (sign irrelevant)
filter_dg = lambda R,dg: [x for x in R.foldings if abs(x.energy) >= abs(dg)]

# Function to find "key" location of a secondary structure to return to mapper.
# Should probably default to finding biggest contiguous section of parentheses.
# Parse through and tabulate every "run" of parens with start index, then sort/pop?
def find_key_structure(folding, within=None):
    if not within:
        within = len(folding.structure)
    keys = []
    # As in, "on a run right now"
    run = False
    thisrun = 0
    thisrun_start = 0
    str_ind = 0
    for l in folding.structure:
        if str_ind > within: break
        if l in "()":
            if not run:
            #    print("Found new run at",str_ind)
                thisrun_start = str_ind
                run = True
            #print("Extending run through",str_ind)
            thisrun += 1
        else:
            # Append length first, then start index, to ease sorting.
            if run:
                #print("Run ends at index",str_ind,"giving:",
                #    thisrun_start,"->",thisrun_start+thisrun,
                #    "- or: ",folding.structure[thisrun_start:thisrun_start+thisrun])
                keys.append((thisrun, thisrun_start))
                thisrun = 0
                run = False
        str_ind += 1
    keys = sorted(keys,reverse=True)
    # Return start and span.
    if not keys: return 0, 0
    return keys[0][1], keys[0][0]

def structure_occludes(structure, start, end, occlusion_ratio = 0.5):
    '''Returns true if a structure occludes the area covered by start/end in structure.
    This is to be used mainly to test for RBS occlusion by secondary structures.
    start, end are slice indices for structure to check for occlusion at.
    occlusion_ratio is the threshold of occlusion which will return True.'''
    num_parens = structure[start:end].count("(") + structure[start:end].count(")")
    if num_parens / abs(end - start) >= occlusion_ratio:
        return True
    return False

# Deprecate and replace: this isn't really a very sensible system, as structures
# from multiple suboptimal structures are considered which may have no likeness
# or comparability.
def map_structures(some_seq, min_energy, prioritise_5prime = True):
    '''Designed to return useful information to PySplicer so it can attempt to
    modify key parts of stable structures to reduce stability. Returns a dictionary
    of the format PySplicer is tailored to work with.
    Min_energy should be a minimum energy value PER N that is considered worrisome.'''
    subopt_foldings = RNAsubopt(some_seq.upper())
    # Need to adjust min_energy to be sequence-length based.
    valid_foldings = filter_dg(subopt_foldings, min_energy*len(some_seq))
    results = {}
    for folding in valid_foldings:
        # If the 5' end of the sequence is occluded, later structures are less
        # relevant, so kill structures here first; assume RBS is here.
        if structure_occludes(folding.structure, 0, 10):
            key_n, span = find_key_structure(folding, within=10)
        else:
            key_n, span = find_key_structure(folding)
        if key_n == 0 and span == 0:
            continue
        results[key_n] = {'span':span}
    return results

def calc_fe(some_seq):
    if isinstance(some_seq, list):
        return RNAfold(''.join(some_seq)).folding.energy
    else:
        return RNAfold(some_seq).folding.energy

def sort_by_fe(seq_list):
    return sorted(seq_list, key = calc_fe)
    
def map_structure(some_seq, min_energy_to_consider = -9):
    '''Returns the apparent key point of a structure, if that structure's free-
    energy is greater (in absolute terms) than abs(min_energy_to_consider).'''
    seq_struct = RNAfold(some_seq)
    if abs(seq_struct.folding.energy) >= abs(min_energy_to_consider):
        # Map "key" part of structure, returning a dict indexed by "main" structure
        # and containing a "span" key.
        start, span = find_key_structure(seq_struct.folding)
        return {start:{"span":span}}
    else:
        return {}
