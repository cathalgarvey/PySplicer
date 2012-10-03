#!/usr/bin/env python3
'''A set of classes that comprise core functionality in PySplicer.'''
from copy import deepcopy # Needed for codon table imports..

class IUPACTool:
    '''IUPACTool provides methods for dealing with DNA using the IUPAC notation,
    offer translation methods to regular expressions, expansion of redundant
    sequences to possible permutations, conversion between RNA and DNA, and
    generation of redundant base-pairing strings or regular expressions.'''
    # B=Not A, V=Not T, D=Not C, H=Not G
    # Strong and Weak: GC vs. AT
    # "Keto" and "aMino": GT vs. AC
    # "puRine" and "pYrimidine"; AG vs CT
    # "-", "N", "." : Wildcards; any N.
    iupac_re = {
           'A': 'A',       'B': '[CGTU]',
           'C': 'C',       'D': '[AGTU]',
           'G': 'G',       'H': '[ACTU]',
           'K': '[GTU]',   'M': '[AC]',
           'N': '[ACGTU]', 'R': '[AG]',
           'S': '[CG]',    'T': 'T', 'U': 'U',
           'V': '[ACG]',   'W': '[ATU]',
           'Y': '[CTU]',   '.': '[ACGTU]',
           '-': '[ACGTU]'
           }
    iupac_complement = {
           "A":	"T", "T": "A",
           "C":	"G", "G": "C",
           "W":	"W", "S": "S",
           "B":	"V", "V": "B",
           "H":	"D", "D": "H",
           "M":	"K", "K": "M",
           "R":	"Y", "Y": "R",
           "V":	"B", "B": "V",
           "N":	"N", ".": ".",
           "-":	"-"
           }

    def __init__(self, strict_type=''):
        '''Arguments for strict_type are "dna" and "rna". In either case, the class is forced to deal only in
        the appropriate type of nucleotide, raising a ValueError if given the wrong type.'''
        self.strictmode = False
        if strict_type:
            self.strictmode = True
            if strict_type.lower() == 'dna':
                # Remove all mentions of Uracil.
                del(self.iupac_re['U'])
                for key in self.iupac_re.keys():
                    self.iupac_re[key].replace("U", "")
            elif strict_type.lower() == 'rna':
                # Remove all mentions of Thymine.
                del(self.iupac_re['T'])
                for key in self.iupac_re.keys():
                    self.iupac_re[key].replace("T", "")
            elif strict_type.lower() == 'strict':
                # This is only here to allow a mode that is strict on inputs but agnostic to type.
                pass
            else:
                raise ValueError("Incorrect argument for strict_type in IUPACTool: must be either 'dna', 'rna', 'strict', or left unspecified.")
        # Below list of characters is used to determine which letters are permissible.
        self.iupac_characters = list(self.iupac_re.keys())

    def nucleic_iupac_purify(self, iupac_sequence, strict=False, strip_flanking_wildcards=False):
        '''Prunes only legal IUPAC DNA/RNA characters from a given string and returns those.

        If strict is set to True, this will raise a ValueError if any illegal characters are encountered.
        If strip_flanking_wildcards is set to True, this will strip "-", "." or "N" from either
        end of the given nucleotide sequence. This may be useful if the permutations of the resulting
        string are to be computed, in which case external wildcards are inefficient.'''
        try:
            assert isinstance(iupac_sequence, str) and isinstance(strict, bool) and isinstance(strip_flanking_wildcards, bool)
        except AssertionError:
            raise TypeError("iupac_sequence must be a string, and strict/strip_flanking_wildcards must be booleans.")
        iupac_sequence = iupac_sequence.upper()
        output_char_list = []
        for char in iupac_sequence:
            if char in self.iupac_characters:
                output_char_list.append(char)
            else:
                if strict:
                    raise ValueError("Illegal character passed to nucleic_iupac_purify in strict mode.")
        output_sequence = ''.join(output_char_list)
        if strip_flanking_wildcards:
            # Not default behaviour, as nucleotides passed to this arg are not always user determined, and size/position
            # of nucleotides in downstream assemblies may depend heavily on the flanking wildcards.
            for char in ["N", "-", "."]:
                output_sequence = output_sequence.replace(char, '')
        return output_sequence

    def iupac_to_regex(self, sequence, return_pattern=False):
        '''Uses a dictionary to assemble a (poorly formatted) python regex function from an IUPAC string.

        For example, given the DNA sequence "AWWSSTGGC", this method returns a regular expression object with
        this pattern: "A[ATU][ATU][GC][GC]TGGC". If the IUPACTools instance is set to strict_type dna/rna, the
        dictionary that powers this method will only output DNA or RNA regexes, and will raise a ValueError if
        an illegal character is detected.'''
        import re
        # First off, get rid of the crap
        sequence = self.nucleic_iupac_purify(sequence, self.strictmode)
        sequence_re_list = []
        for base in sequence:
            sequence_re_list.append(iupac_re[base])
        # Flatten the list into a string
        sequence_re_str = ''.join(sequence_re_list)
        # Compile a regex function from the string
        sequence_re = re.compile(sequence_re_str)
        if return_pattern:
            return sequence_re, sequence_re_str # Former is the regex, latter is just for verbose readout.
        return sequence_re

    def iupac_complement(self, iupac_string):
        'Returns the complement of a redundant iupac notation string in iupac notation.'
        assert isinstance(sequence,str)
        iupac_string = self.nucleic_iupac_purify(iupac_string, self.strictmode)
        complement_sequence_list = []
        for base in iupac_string:
            complement_sequence_list.append(self.iupac_complement[base])
        complement_sequence = ''.join(complement_sequence_list)
        return complement_sequence

    def iupac_rev_complement(self, iupac_string):
        complement_string = self.iupac_complement(iupac_string)
        return complement_string[::-1]

    def all_iupac_permutations(self, iupac_string):
        '''Returns a list of possible permutations for a given IUPAC nucleotide string.

        It is suggested that this only be used when in one of the "strict modes": "dna" or "rna",
        both to avoid confusion and to reduce the number of permutations returned.
        Warning: no size/sanity limits; if you pass large or highly redundant sequences to
        this method it may choke.'''
        # Will hold a list item of possible characters for each character in iupac_string.
        snp_list = []
        # Will hold the permutations for export.
        perm_list = []
        iupac_string = self.nucleic_iupac_purify(iupac_string, self.strictmode)
        for char in iupac_string:
            # Call up the regex function for a character, minus the brackets.
            snp_list.append(self.iupac_re[char].strip("[]"))
        perm_generator = itertools.product(*snp_list) # "*" character dumps list contents as args.
        for item in perm_generator:
            perm_list.append(''.join(item))
        return perm_list

class NucleotideMapper:
    '''DNAMapper is a class that can be fed and updated with DNA/RNA sites
    that are undesirable in a sequence, and which provides search functions
    for mapping any matches in a target sequence by string or list index
    (depending on which method is called), so that other class methods or
    functions can resolve these issues.

    Supported features that this can discover: consensus sites or redundant
    DNA/RNA sequences.
    Planned features: Detection of secondary structure by fold-and-walk of
    the sequence with a minimum structure-seeking threshold. Detection of
    potential G-quadroplexes. Detection of conspicuous sites of DNA torsion?'''
    pass

class CodonTable:
    '''Contains and provides data for codon usage tables. Methods can be used to return encoded amino acids,
    to return frequencies for codons, or to return whole sets of codons with frequencies for given aminos.
    Can also be used to convert codon tables in XXX or YYY format to internal JSON-based structure.
    Methods can be used to remove codons entirely, automatically adjusting the frequencies of remaining
    synonymous codons to compensate.
    Methods employing the preceding function can cleave all codons from the table whose frequency is below
    a certain threshold.
    The reverse operation is also possible; injecting an (unused!) codon '''
    def __init__(self, init_table=None):
        self.validaminos = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', '*',
                            'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
        self.validcodons = ['GCA', 'GCC', 'GCU', 'GCG', 'UGC', 'UGU', 'GAG', 'GAA', 'GAU', 'GAC',
                            'GGU', 'GGG', 'GGA', 'GGC', 'UUU', 'UUC', 'AUA', 'AUC', 'AUU', 'CAC',
                            'CAU', 'AAG', 'AAA', 'UAA', 'UGA', 'UAG', 'AUG', 'CUU', 'CUG', 'CUC',
                            'CUA', 'UUG', 'UUA', 'AAU', 'AAC', 'CAA', 'CAG', 'CCU', 'CCG', 'CCA',
                            'CCC', 'UCU', 'AGC', 'UCG', 'UCC', 'UCA', 'AGU', 'CGA', 'CGC', 'AGA',
                            'AGG', 'CGG', 'CGU', 'ACC', 'ACA', 'ACG', 'ACU', 'UGG', 'GUC', 'GUA',
                            'GUG', 'GUU', 'UAC', 'UAU']
        if init_table: self.set_codon_table(init_table)

    def frequency_of_codon(self, codon):
        'Take a three-letter string identifying a codon, return the frequency of incidence.'
        codon = self._validate_codon(codon)
        return self.codonkeyed[codon]['frequency']

    def encodes(self, codon):
        'Take a three-letter string identifying a codon, and return the one-letter Amino/Stop encoded by it.'
        codon = self._validate_codon(codon)
        return self.codonkeyed[codon]['amino']

    def codons_for(self, amino):
        'Take a single-letter amino acid or stop codon identifier, and return a dictionary of codons with frequencies.'
        try:
            assert isinstance(amino, str)
            amino = amino.upper()
            assert amino in self.validaminos
        except AssertionError: raise ValueError("Invalid or improperly formatted amino. Single letter code only.")
        return self.codontable[amino]

    def remove_codon(self, codon):
        self.codontable = self._remove_codon(codon, self.codontable)
        self.codonkeyed = self.invert_cut(self.codontable)

    def _remove_codon(self, codon, codon_table):
        'Removes a codon from the specified codon table and returns the table after rebalancing frequencies.'
        codon = self._validate_codon(codon)
        try:
            self.check_cut(codon_table)
        except ValueError:
            print("Error testing codon table before removal.")
        #if self._is_valid_codon(codon):
        codons_in_table = self.invert_cut(codon_table)
        if codon in codons_in_table.keys():
            codon_encodes = codons_in_table[codon]["amino"]
            del(codon_table[codon_encodes][codon])
        codon_table = self.balance_frequencies(codon_table)
        return codon_table

    def add_codon(self, target_codon, target_amino, frequency):
        self.codontable = self._add_codon(target_codon, target_amino, frequency, self.codontable)
        self.codonkeyed = self.invert_cut(self.codontable)

    def _add_codon(self, target_codon, target_amino, frequency, codon_table):
        '''Adds a codon to a specified amino and rebalances frequencies to accomodate specified relative frequency.
        Clobbers over codon if it already exists.'''
        target_codon = self._validate_codon(target_codon)
        self.check_cut(codon_table)
        # Get inverted table to check that codon isn't in use
        inverted_table = self.invert_cut(codon_table)
        if target_codon in inverted_table:
            # Remove and rebalance
            del(codon_table[inverted_table[target_codon]['amino']][target_codon])
            codon_table = self.balance_frequencies(codon_table)
        # Get amino dict so it can be edited wholesale and re-inserted
        amino_dict = codon_table[target_amino]
        for existing_codon in amino_dict:
            # Reduce frequencies of all the other codons to match the desired insert frequency.
            amino_dict[existing_codon]['frequency'] = amino_dict[existing_codon]['frequency'] * (1 - frequency)
        # Now add new codon at desired frequency
        amino_dict[target_codon] = {}
        amino_dict[target_codon]['frequency'] = frequency
        # Re-insert amino dictionary, replacing old one.
        codon_table[target_amino] = amino_dict
        # Rebalance codons after all the above trauma.
        codon_table = self.balance_frequencies(codon_table)
        return codon_table

    def _validate_codon(self, codon):
        'Just to verify that a valid three-letter code has been given.'
        try:
            assert isinstance(codon, str)
            codon = codon.upper()
            assert codon in self.validcodons
        except AssertionError:
            raise ValueError("Specified codon is invalid or badly formatted.")
        return codon

    def invert_cut(self, cut_dict):
        'Returns a dict keyed by codon, where value is a dict containing the "amino" and "frequency".'
        codon_keyed = {}
        for Amino in cut_dict:
            for Codon in cut_dict[Amino]:
                codon_keyed[Codon] = {'amino':Amino, 'frequency':cut_dict[Amino][Codon]['frequency']}
        return codon_keyed

    def purge_by_frequency(self, threshold):
        'If any codons have a usage frequency below threshold, they are set to 0 and other codons compensated.'
        if not isinstance(threshold, float) and 1 > threshold > 0:
            raise ValueError("Threshold frequency must be a float between 0 and 1.")
        for Codon in self.validcodons:
            try:
                if self.frequency_of_codon(Codon) < threshold:
                    CodonAmino = self.codonkeyed[Codon]['amino']
                    self.add_codon(Codon, CodonAmino, 0.0)
            except KeyError:
                # If this is a spliced table or whatever, then there may be missing codons, and that's OK.
                pass
    
    def set_codon_table(self, cut_dict, copy=True):
        'Import a codon table in the format used internally, probably previously saved from this software.'
        if copy:
            cut_dict = deepcopy(cut_dict) # So external table, once imported, is not edited by internal methods.
        try: self.check_cut(cut_dict)
        except:
            raise ValueError("Error checking codon table in import_cut_dict; aminos/codons may be missing or frequencies may be incorrect.")
        self.codontable = self.balance_frequencies(cut_dict)
        self.codonkeyed = self.invert_cut(self.codontable)

    def set_codon_table_from_file(self, cut_file):
        import json
        try:
            with open(cut_file) as cut_infile:
                cut = json.loads(cut_infile.read())
            self.set_codon_table(cut, copy=False)
            return True
        # "Patch" likely errors with more informative ones.
        except IOError as e:
            raise IOError("Error: Cannot open specified file; does it exist?", e)
        except ValueError as v:
            raise ValueError("Error: File contents probably not valid JSON.", v)

    def printout_cut(self):
        'Print codon table, JSON-formatted, to stdout for piping etc.'
        import json
        print(json.dumps(self.codontable))

    def export_codon_table_to_file(self, cut_file):
        'Create or overwrite specified file with codon usage table in json format.'
        import json
        with open(cut_file, mode='w') as SaveFile:
            SaveFile.write(json.dumps(self.codontable),'',sep='\r\n')

    def import_cud_cut(self, cut_string):
        '''Parses, corrects and verifies validity of a codon usage table from Codon Usage Database (which must be presented with
        Aminos and relative frequencies in the format obtained by selecting "Format" with a genetic code and selecting
        "Codon Usage Table with Amino Acids".
        '''
        candidate_cut = self.parse_cud_cut(cut_string)
        candidate_cut = self.balance_frequencies(candidate_cut)
        try:
            self.check_cut(candidate_cut)
        except:
            raise ValueError("Error checking codon table in import_cut; aminos/codons may be missing or frequencies may be incorrect.")
        self.codontable = candidate_cut
        self.codonkeyed = self.invert_cut(self.codontable)

    def parse_cud_cut(self, cut_string):
        '''Parse a codon usage table (CUT) as provided by the Codon Usage Database (http://www.kazusa.or.jp/codon/) when asked to present
        a table with amino acids, where each codon is presented as: [triplet] [frequency: per thousand] ([number]).
        The whole table should be provided to this method as a string. This method discards frequency per thousand and
        codon numbers and only uses the first two fields of the CUT; the codon and the relative frequency vs. synonymous codons.
        Output is a dictionary containing the same values; keys are capitalised codons, value is relative frequency.'''
        # String methods magic: split the given codon-usage-table (cut_string) by right-brackets to get one field per codon,
        #  then for each thing in that list, strip flanking whitespace and remove the other bracket/redundant space. Result is a
        #  string for each codon of format 'triplet amino relative_frequency frequency_per_1000 number_of_codons_in_dataset',
        #  space-delimited. This is then sliced [0:3] to get just [triplet amino relative_frequency]
        # Note that "strip(")") is used first on cut_string before splitting; this removes trailing bracket, which otherwise
        #  would lead to an empty entry in the resulting list after splitting by this character.
        # Beginners to python: please forgive me. Regard this as a learning experience rather than a justification for murder.
        cut_as_list = cut_string.strip().strip(")").split(")")
        cut_as_list = [x.strip().replace("(", " ").split()[0:3] for x in cut_as_list]
        codon_frequency_dict = {}
        for codon in cut_as_list:
            # format for each codon is: [triplet, amino, relative_frequency], all as strings.
            # Can use first two values as is, third can be cast as a float using float()
            # Below: format this list into an amino-keyed dictionary of codon-choice subdictionaries,
            #  where each codon key yields its relative frequency (as a float) as a value.
            if codon[1] in codon_frequency_dict:
                # Add this codon to the amino entry in the table, with the float value of its relative frequency as value.
                codon_frequency_dict[codon[1]][codon[0]] = {'frequency':float(codon[2])}
            else:
                # Create entry for this amino/stop, starting with this codon in there as value.
                codon_frequency_dict[codon[1]] = {codon[0]:{'frequency':float(codon[2])}}
        return codon_frequency_dict

    def balance_frequencies(self, cut_dict):
        '''Takes a codon frequency in internal object format and rebalances codon frequencies for each amino acid proportionally.

        This can be used after deletion of a codon to re-normalise remaining codons, or to resolve minor statistical errors in
        unmodified tables.'''
        for Amino in cut_dict:
            frequencysum = 0.0
            for Codon in cut_dict[Amino]:
                frequencysum += cut_dict[Amino][Codon]['frequency']
            for Codon in cut_dict[Amino]:
                cut_dict[Amino][Codon]['frequency'] = cut_dict[Amino][Codon]['frequency'] / frequencysum
        return cut_dict

    def check_cut(self, cut_dict, strict=False):
        '''Expects a dictionary of format {Single-letter-amino:{Codon1:float_frequency, Codon2:float_frequency}..} -
        confirms that all necessary Aminos/Stops and Codons are present as keys, and that the sums of all listed codons sum to 1.0.
        If a table that does not assign certain codons is desired, it is recommended to assign unused codons to a null key such as "$".'''
        try:
            presentcodons = []
            freq = 0.0
            for X in self.validaminos:
                freq = 0.0
                assert X in cut_dict.keys()
                # Strict: Aminos must have at least one codon encoding them.
                if strict: assert len(cut_dict[X]) > 0
                presentcodons.extend(list(cut_dict[X].keys()))
                for Codon in cut_dict[X]:
                    freq += cut_dict[X][Codon]['frequency']
                #print(freq)
                assert 0.99 <= freq <= 1.01
            # Nonstrict: only care that all keys for aminos are valid codons.
            for X in presentcodons:
                assert X in self.validcodons
            # Strictmode: All codons must be in use.
            if strict:
                for X in self.validcodons:
                    assert X in presentcodons
        except AssertionError: raise ValueError("Codon Table is not valid.")

    def splice_table(self, other_table, lower_threshold=0.1):
        'Attempts to find a compromise table between internal table and a given table. Uses self._splice_table method.'
        table_attempt = self._splice_table(self.codontable, other_table, lower_threshold)
        if table_attempt:
            # Strict mode false: this is one of a few legitimate cases where codons may be entirely absent:
            self.check_cut(table_attempt, strict=False)
            self.codontable = table_attempt
            self.codonkeyed = self.invert_cut(self.codontable)
            return True
        else:
            return False

    def _splice_table(self, table1, table2, threshold, verbose=False):
        '''Create a compromise table between table1 and table2, attempting to omit codons that fall below lower_threshold frequency.

        Compromise is reached using the following strategy:
          - First, check that codons exist for all amino acids.
          - Second, check that codons are used for same aminos; remove any that conflict.
          - Third, remove codons with a frequency below lower_threshold
          - All removals should be removed from *both* tables
          - For remaining codons in both tables..average their relative frequencies?
        If compromise is impossible (polar opposite usage of codons for a given amino), this returns False.'''
        def vprint(*args, **vargs):
            if verbose: print(*args, **vargs)
        CompromiseTable = {}
        for A in table1.keys():
            vprint("==========","\nAmino:", A)
            CompromiseTable[A] = {}
            for C in table1[A].keys():
                Freq1 = table1[A][C]['frequency']
                Freq2 = table2[A][C]['frequency']
                vprint('\t---', C, "---\n\t Freq1:", Freq1, "\n\t Freq2", Freq2)
                try:
                    if Freq1 > threshold and Freq2 > threshold:
                        vprint("\t  Freq in both greater than threshold: Compromise frequency set to Average")
                        CFreq = Freq1 + Freq2 / 2
                        CompromiseTable[A][C] = {'frequency':CFreq}
                    else:
                        vprint("\t  Freq in both tables below threshold: Compromise frequency set to Zero")
                        CompromiseTable[A][C] = {'frequency':0.0}
                except KeyError:
                    vprint("Codon not common to both tables at this Amino; skipping entirely.")
            if len(CompromiseTable[A]) == 0:
                # Would this ever happen?
                raise ValueError("Compromise Table at Amino "+A+" has no codons with given parameters and tables.")
            AtLeastOneValidCodon = False
            for Codon in CompromiseTable[A]:
                if CompromiseTable[A][Codon]['frequency'] > 0.0:
                    AtLeastOneValidCodon = True
                    break
            if not AtLeastOneValidCodon:
                raise ValueError('''No valid codons (f>0.0) remain in Compromise Table at Amino {0}:
this probably indicates an overambitious threshold frequency or totally incompatible tables.'''.format(A))
        CompromiseTable = self.balance_frequencies(CompromiseTable)
        return CompromiseTable

class CodonJuggler:
    '''Can translate codons or reverse-translate them, using an embedded CodonTable object.

    Provides methods to:
    * Analyse codon strings for relative frequency, 
    * Select, according to relative frequency, a codon for a given amino
    * Select "best" amino for a given amino
    * Extend the above to lists of aminos, returning lists of codons
    * Return X permutations of RNA codons for a given peptide or codon list

    Methods are provided to import and parse codon adaptation tables, to use these to provide information on
    codon lists or sequence strings, and to translate or reverse translate, with regard to codon frequency.

    Therefore, removal of entire codons is possible, and ratio of remaining codons to one another will
    be retained. This might be useful if a codon-replacement strategy is desirable for later re-use of codons
    for alternative amino acids.
    Methods can also return permutations of given sequences that retain amino identity, which can help when
    resolving problems detected by NucleotideMapper instances for example.'''
    def __init__(self, codontable):
        self.codontable = codontable

    def analyse_cai_str(self, codonstring, frameoffset=0):
        pass

    def analyse_cai_list(self, codonlist):
        pass

class bioseq(str):
    'Inherits string methods in addition to some biological sequence-specific methods.'
    iupac_re_chars = {
           'A': 'A',       'B': '[CGTU]',
           'C': 'C',       'D': '[AGTU]',
           'G': 'G',       'H': '[ACTU]',
           'K': '[GTU]',   'M': '[AC]',
           'N': '[ACGTU]', 'R': '[AG]',
           'S': '[CG]',    'T': 'T', 'U': 'U',
           'V': '[ACG]',   'W': '[ATU]',
           'Y': '[CTU]',   '.': '[ACGTU]',
           '-': '[ACGTU]' }
    iupac_complement_chars = {
           "A":	"T", "T": "A",
           "C":	"G", "G": "C",
           "W":	"W", "S": "S",
           "B":	"V", "V": "B",
           "H":	"D", "D": "H",
           "M":	"K", "K": "M",
           "R":	"Y", "Y": "R",
           "V":	"B", "B": "V",
           "N":	"N", ".": ".",
           "-":	"-" }
    def __new__(self, value):
        'Overrides default string behaviour to make all characters uppercase.'
        # Strings, being immutable, don't use "__init__", they use "__new__".
        strself = str.__new__(self, value.upper())
        strself.rnachars = ["A","U","C","G"]
        strself.dnachars = ["A","T","C","G"]
        strself.aminochars = ['A','C','E','D','G','F','I','H','K','*',
                              'M','L','N','Q','P','S','R','T','W','V','Y']
        strself.iupacchars = ['A','C','B','D','G','H','K','-','M','N',
                               'S','R','U','T','W','V','Y','.']
        return strself
    def as_codon_list(self, offset=0, trailing=False):
        '''Slices the string into three-letter partitions with an optional starting offset.
        Optional "trailing" argument tells the method to return additional letters after the last
        set of three; this can cause bugs in parsers that expect only codons, so default is False.'''
        if trailing:
            from itertools import zip_longest
            codonlist = [''.join(x) for x in zip_longest(self.upper()[offset::3],
                                                         self.upper()[offset+1::3],
                                                         self.upper()[offset+2::3],
                                                         fillvalue='') ]
        else:
            codonlist = [''.join(x) for x in zip(self.upper()[offset::3],
                                                 self.upper()[offset+1::3],
                                                 self.upper()[offset+2::3]) ]
        return codonlist
    def _is_valid_x(self,alphabet):
        'Helper method for type-validation methods; DNA/RNA/Aminos/IUPAC etc.'
        for char in self:
            if char not in alphabet:
                return False
        return True
    def is_valid_dna(self):
        '''Checks that string contains only DNA chars.'''
        return self._is_valid_x(self.dnachars)
    def is_valid_rna(self):
        '''Checks that string contains only RNA chars.'''
        return self._is_valid_x(self.rnachars)
    def is_valid_amino(self):
        '''Checks that string contains only Amino Acid chars.'''
        return self._is_valid_x(self.aminochars)
    def is_valid_iupac(self):
        '''Checks that string contains only IUPAC DNA/RNA chars.'''
        return self._is_valid_x(self.iupacchars)
    def iupac_only(self, strip_flanking_wildcards=False, strict=False):
        '''Prunes only legal IUPAC DNA/RNA characters from a given string and returns those.

        If strict is set to True, this will raise a ValueError if any illegal characters are encountered.
        If strip_flanking_wildcards is set to True, this will strip "-", "." or "N" from either
        end of the given nucleotide sequence. This may be useful if the permutations of the resulting
        string are to be computed, in which case external wildcards are inefficient.'''
        if strict and not self.is_valid_iupac():
            raise ValueError("Value of bioseq does not appear to be valid IUPAC.")
        output_char_list = []
        for char in self:
            if char in self.iupacchars:
                output_char_list.append(char)
        output_sequence = ''.join(output_char_list)
        if strip_flanking_wildcards:
            output_sequence = output_sequence.strip("N-.")
        return output_sequence
    def iupac_to_regex(self, strip_wildcards = False, return_pattern=False):
        '''Uses a dictionary to assemble a (poorly formatted) python regex function from an IUPAC string.

        For example, given the DNA sequence "AWWSSTGGC", this method returns a regular expression object with
        this pattern: "A[ATU][ATU][GC][GC]TGGC". If the IUPACTools instance is set to strict_type dna/rna, the
        dictionary that powers this method will only output DNA or RNA regexes, and will raise a ValueError if
        an illegal character is detected.'''
        if not self.is_valid_iupac():
            raise ValueError("Cannot convert a non-iupac bioseq to an IUPAC regex.")
        import re
        # First off, get rid of the crap
        sequence_re_list = []
        for base in self.iupac_only(strip_wildcards):
            sequence_re_list.append(self.iupac_re_chars[base])
        # Flatten the list into a string
        sequence_re_str = ''.join(sequence_re_list)
        # Compile a regex function from the string
        sequence_re = re.compile(sequence_re_str)
        if return_pattern:
            return sequence_re, sequence_re_str # Former is the regex, latter is just for verbose readout.
        return sequence_re
    def iupac_complement(self):
        'Returns the complement of a redundant iupac notation string in iupac notation.'
        if not self.is_valid_iupac():
            raise ValueError("Cannot get IUPAC Complement of a non-IUPAC string.")
        complement_sequence_list = []
        for base in self:
            complement_sequence_list.append(self.iupac_complement_chars[base])
        return ''.join(complement_sequence_list)
    def iupac_rev_complement(self):
        return self.iupac_complement()[::-1]
    def to_rna(self):
        if not self.is_valid_iupac():
            raise ValueError("Cannot convert invalid IUPAC code to RNA.")
        return self.replace("U","T")
    def to_dna(self):
        if not self.is_valid

Tim Buckley - 0214385956
0862748776


mytestdna = bioseq("agcctagtctggtgtcggcagcataggctattaataattaaataa")
