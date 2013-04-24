'''nucutils - Nucleotide-wrangling functions and constants.
by Cathal Garvey
Part of the DNAmespace project. License accessible as nucutils.license.
'''
# A set of dictionaries with keys "codons", "aminos" and "starts", containing
# respectively a dictionary mapping of codons to aminos, aminos to lists of
# corresponding codons, and start codons.
#from pysplicer2 import translationtables
from pysplicer import translationtables
import random
import re
import itertools

def _chunks(l, n):
    "Yield successive n-sized chunks from l."
    for i in range(0, len(l), n): yield l[i:i+n]
iupac_dna_chars =    "ATCGBVDHSWKMRYN"
iupac_rna_chars =    "AUCGBVDHSWKMRYN"
iupac_nucleotides = "ATUCGBVDHSWKMRYN"
#iupac_nucleotides = [ 'A', 'T', 'C', 'G', 'U',    # Canonical bases
#                     'B', 'V', 'D', 'H',     # B=Not A, V=Not T, D=Not C, H=Not G
#                     'S', 'W',               # Strong and Weak: GC vs. AT
#                     'K', 'M',               # "Keto" and "aMino": GT vs. AC
#                     'R', 'Y',               # "puRine" and "pYrimidine": AG vs CT
#                     'N' ]         # Wildcards: any N.

aminoiupac = {'A': 'Alanine',              'B': 'Aspartic Acid or Asparagine',
              'C': 'Cysteine',             'D': 'Aspartic Acid',
              'E': 'Glutamic Acid',        'F': 'Phenylalanine',
              'G': 'Glycine',              'H': 'Histidine',
              'I': 'Isoleucine',           'K': 'Lysine',
              'L': 'Leucine',              'M': 'Methionine',
              'N': 'Asparagine',           'P': 'Proline',
              'Q': 'Glutamine',            'R': 'Arginine',
              'S': 'Serine',               'T': 'Threonine',
              'V': 'Valine',               'W': 'Tryptophan',
              'X': 'Any',                  'Y': 'Tyrosine',
              'Z': 'Glutamine or Glutamic Acid', '*': 'Stop'}
amino_non_iupac = {'A': 'Alanine',   'C': 'Cysteine',      'D': 'Aspartic Acid',
              'E': 'Glutamic Acid',  'F': 'Phenylalanine', 'G': 'Glycine',
              'H': 'Histidine',      'I': 'Isoleucine',    'K': 'Lysine',
              'L': 'Leucine',        'M': 'Methionine',    'N': 'Asparagine',
              'P': 'Proline',        'Q': 'Glutamine',     'R': 'Arginine',
              'S': 'Serine',         'T': 'Threonine',     'V': 'Valine',
              'W': 'Tryptophan',     'Y': 'Tyrosine' }

# All complements are in lowercase so that string substitution is simple:
# just set string to uppercase (if not already), replace all characters
# with lowercase complements, then if desired set resulting string to
# uppercase again.
dnacomplement = {"A":"t","T":"a","G":"c","C":"g"}
rnacomplement = {"A":"u","U":"a","G":"c","C":"g"}
dnaiupaccomplement = {   # Translates to an IUPAC sequence of potential complements.
              "A":   "t", "T": "a",
              "C":   "g", "G": "c",
              "W":   "w", "S": "s",
              "B":   "v", "V": "b",
              "H":   "d", "D": "h",
              "M":   "k", "K": "m",
              "R":   "y", "Y": "r",
              "V":   "b", "B": "v",
              "N":   "n"}

rnaiupaccomplement = {   # Translates to an IUPAC sequence of potential complements.
              "A":   "u", "U": "a", "C":   "g", "G": "c",
              "W":   "w", "S": "s", "B":   "v", "V": "b",
              "H":   "d", "D": "h", "M":   "k", "K": "m",
              "R":   "y", "Y": "r", "V":   "b", "B": "v",
              "N":   "n"}

bpdict = {'A':'UT','C':'G','G':'UC','U':'GA','T':'A'}
def bp(l,m):
    'Returns true if the two (explicit, non-extended-IUPAC) nucleotides can base-pair.'
    if m in bpdict.get(l,''): return True
    else: return False

def _uniquify(string):
    '''Reduces a string down to its component characters.
    This is a fast, order-preserving function for removing duplicates
    from a sequence/list, by "Dave Kirby"
    Found here: http://www.peterbe.com/plog/uniqifiers-benchmark'''
    seen = set()
    return [x for x in string if x not in seen and not seen.add(x)]

def random_dna(length):
    return ''.join([random.choice("GCAT") for x in range(0,length)])

def random_rna(length):
    return ''.join([random.choice("GCAU") for x in range(0,length)])

def random_iupac_dna(length):
    return ''.join([random.choice(iupac_dna_chars) for x in range(0,length)])

def random_iupac_rna(length):
    return ''.join([random.choice(iupac_rna_chars) for x in range(0,length)])

def random_iupac_chars(length):
    return ''.join([random.choice(iupac_nucleotides) for x in range(0,length)])

def deduce_alphabet(string):
    string = string.upper()
    charset = _uniquify(string)
    couldbenuc, couldbeaminos = True, True
    for character in charset:
        if character not in iupac_nucleotides:
            couldbenuc = False
        if character not in aminoiupac:
            couldbeaminos = False
    if couldbenuc:
        if "T" in charset and "U" in charset:
            if not couldbeaminos:
                return "hybrid_nuc"
        elif "T" in charset:
            return "dna"
        elif "U" in charset:
            return "rna"
        else:
            # Assume DNA if ambiguous.
            return "dna"
    elif couldbeaminos:
        return "aminos"
    else:
        raise ValueError("Could not deduce an IUPAC alphabet from the input "
                         "string. Charset is '{}'".format(charset))

def get_complement_alphabet(string):
    string = string.upper()
    alphabet = deduce_alphabet(string)
    if alphabet == "aminos":
        raise ValueError("Cannot get the complement of an amino sequence.")
    elif alphabet == "rna":
        return rnaiupaccomplement
    elif alphabet == "dna":
        return dnaiupaccomplement
    else:
        raise ValueError("Could not get complement for string: "+string+"\nAlphabet deduced was: "+str(alphabet))

def get_complement(nucleotides):
    'Given a string of nucleotides (RNA *or* DNA), return reverse complement.'
    # Reverse the string:
    nucleotides = nucleotides.upper()[::-1]
    # Determine molecule type:
    basedict = get_complement_alphabet(nucleotides)
    for base in basedict.keys():
        # Replace each base with its complementary base in lowercase:
        # (lowercase means previously substituted bases aren't replaced)
        nucleotides = nucleotides.replace(base, basedict[base])
    # Return sequence in uppercase:
    return nucleotides.upper()

def translate(sequence, table, frame=1):
    sequence = sequence.upper()
    frame -= 1
    translation_table = translationtables.__dict__[table]
    aminos = []
    for codon in _chunks(sequence[frame:], 3):
        if len(codon) < 3: break
        encoded = translation_table['codons'][codon]
        aminos.append(encoded)
        if encoded == "*": break
    return ''.join(aminos)

def dumb_backtranslate(sequence, table):
    'Using the chosen table, return a back-translation of an amino sequence without codon weighting.'
    sequence = sequence.upper()
    translation_table = translationtables.__dict__[table]
    codons = []
    for amino in sequence:
        candidate_codons = translation_table['aminos'][amino]
        codons.append(random.choice(candidate_codons))
    return ''.join(codons)

class bioseq(str):
    'Inherits string methods in addition to some biological sequence-specific methods.'
    # Can be used to assemble regular expressions from IUPAC notation
    iupac_re_chars = {
           'A': 'A',       'B': '[CGTU]',    'C': 'C',       'D': '[AGTU]',
           'G': 'G',       'H': '[ACTU]',    'K': '[GTU]',   'M': '[AC]',
           'N': '[ACGTU]', 'R': '[AG]',      'S': '[CG]',    'T': 'T',
           'U': 'U',       'V': '[ACG]',   'W': '[ATU]',     'Y': '[CTU]',
           '.': '[ACGTU]', '-': '[ACGTU]'}
    # Same as above but with only T or U
    iupac_re_dna_chars = {'A':'A','B':'[CGT]','C':'C','D':'[AGT]','G':'G',
        'H':'[ACT]','K':'[GT]','M':'[AC]','N':'[ACGT]','R':'[AG]','S':'[CG]',
        'T':'T','V':'[ACG]','W':'[AT]','Y':'[CT]','.':'[ACGT]','-':'[ACGT]'}
    iupac_re_rna_chars = {'A':'A','B':'[CGU]','C':'C','D':'[AGU]','G':'G',
        'H':'[ACU]','K':'[GU]','M':'[AC]','N':'[ACGU]','R':'[AG]','S':'[CG]',
        'U':'U','V':'[ACG]','W':'[AU]','Y':'[CU]','.':'[ACGU]','-':'[ACGU]'}
    dna_iupac_complement = dnaiupaccomplement
    rna_iupac_complement = rnaiupaccomplement
    rnachars = ["A","U","C","G"]
    dnachars = ["A","T","C","G"]
    aminochars = ['A','C','E','D','G','F','I','H','K','*',
                    'M','L','N','Q','P','S','R','T','W','V','Y']
    iupacchars = ['A','C','B','D','G','H','K','-','M','N',
                    'S','R','U','T','W','V','Y','.']

    @classmethod
    def _iupac_nuc_only(cls, sequence, strip_wildcards=False):
        '''Prunes only legal IUPAC DNA/RNA characters from a given string and returns those.
        If strip_wildcards is set to True, this will strip "-", "." or "N" from ends.'''
        # This is called BY other classmethods to give an object, it's not supposed to
        # wrap returned output in cls!
        output_char_list = []
        for char in sequence.upper():
            if char in cls.iupacchars:
                output_char_list.append(char)
        output_seq = ''.join(output_char_list)
        if strip_wildcards:
            output_seq = output_seq.strip("N-.")
        return output_seq

    @classmethod
    def fromstring(cls, sequence, strict=False):
        sequence = sequence.upper()
        if strict: return cls(cls._iupac_nuc_only(sequence))
        return cls(sequence)

    def iupac_nuc_only(self, *args):
        '''Prunes only legal IUPAC DNA/RNA characters from a given string and returns those.
        If strip_wildcards is set to True, this will strip "-", "." or "N" from ends.'''
        return bioseq(self._iupac_nuc_only(self, *args))

    def infer_type(self):
        'Calls deduce_alphabet on self, returns same values: dna/rna/aminos/hybrid_nuc.'
        return deduce_alphabet(self)

    def _is_valid_iupac(self):
        'Simple look-before-you-leap method for when try/except is too much work.'
        for char in self:
            if char not in self.iupacchars: return False
            return True

    def as_codon_list(self, offs=0, trailing=False):
        '''Slices the string into three-letter partitions with an optional starting offset.
        Optional "trailing" argument tells the method to return additional letters after the last
        set of three; this can cause bugs in parsers that expect only codons, so default is False.'''
        strV = self.upper()
        # Have to do this likely-inefficient workaround to get this to work with
        # negative offsets, where normally this would kill zips due to one of the
        # fields being empty.
        if offs < 0:
            offs, strV = 0, strV[abs(offs):]
        frames = (strV[offs+0::3],
                  strV[offs+1::3],
                  strV[offs+2::3])
        if trailing:
            return [''.join(x) for x in itertools.zip_longest(*frames,fillvalue='')]
        else:
            return [''.join(x) for x in zip(*frames)]

    def codons(self, frame=0, leading=False, trailing=True):
        # To deal with negative frames, just convert to positive frames, and
        # do or do not return leading codons depending on that argument.
        if frame < 0 and trailing: frame, seq = abs(frame), self
        elif frame < 0: seq, frame = self[abs(frame):], 0
        else: seq = self
        # Give an empty list if not "leading" or if there are no prefix codons.
        # Otherwise would get a list with an empty string, which is bug-bait.
        prefix = [seq[:frame]] if leading and seq[:frame] else []
        # Concatenate prefix list to a generated list, keeping or discarding
        # bits smaller than 3n long depending on the "trailing" argument.
        if trailing: return prefix + [x for x in _chunks(seq[frame:],3)]
        else: return prefix + [x for x in _chunks(seq[frame:],3) if len(x)==3]

    def iupac_rev_complement(self):
        return bioseq(get_complement(self))

    def to_rna(self):
        return bioseq(self.upper().replace("T","U"))

    def to_dna(self):
        return bioseq(self.upper().replace("U","T"))

    def iupac_to_regex(self, include_complement=False, strip_wildcards = False,
                             alphabet='', return_type="regex"):
        '''Returns a (poorly formatted) regex object or pattern for Sequence.
        E.g.: The sequence "AWWSSTGGC" becomes "A[ATU][ATU][GC][GC]TGGC".
        return_type: "regex" (default), "string", "list" (of patterns for each char).'''
        if not self._is_valid_iupac(): raise ValueError("Can't convert non-nucleotides to regex: "+self)
        # If alphabet is not given, infer from character content.
        alphabet = alphabet.lower() or self.infer_type()
        if   alphabet == "dna":        alphabet = self.iupac_re_dna_chars
        elif alphabet == "rna":        alphabet = self.iupac_re_rna_chars
        elif alphabet == "hybrid_nuc": alphabet = self.iupac_re_chars
        elif alphabet == "aminos": raise NotImplementedError("Regex for Aminos not done yet.")
        else: raise NotImplementedError("Alphabet "+str(alphabet)+" is not recognised.")
        # Get RE characters from a crap-stripped sequence:
        sequence_re_list = []
        iupac_pattern = self.iupac_nuc_only(strip_wildcards)
        for base in iupac_pattern:
            sequence_re_list.append(alphabet[base])
        if include_complement:
            # Make an "or" type regular expression of form "(this|that)"
            sequence_re_list.insert(0,"(")
            sequence_re_list.append("|")
            for base in get_complement(iupac_pattern):
                sequence_re_list.append(alphabet[base])
            sequence_re_list.append(")")
        # Return favourite representation:
        if return_type.lower() == "list": return sequence_re_list
        elif return_type.lower() == "string": return ''.join(sequence_re_list)
        elif return_type.lower() == "regex": return re.compile(''.join(sequence_re_list))
        else: raise ValueError('"'+str(return_type)+'" is not a valid return type.')

    def all_iupac_permus(self, alphabet=''):
        '''Returns a generator for all possible literal RNA/DNA permutations of an IUPAC string.
        The generator is itertools.product and provides list output.
        "alphabet" is passed directly to self.iupac_re(), which is used to generate
        the permutation list. If not specified, it will be inferred.
        Warning: this generator may go on for a very long time for highly redundant sequences;
        it would be unwise to unpack this generator into a list or similar.'''
        # Will hold a list item of possible characters for each character in iupac_string.
        snp_list = []
        iupac_string = self.iupac_nuc_only()
        for regex_char in self.iupac_to_regex(alphabet=alphabet, return_type="list"):
            # Regex function will return list of all possible characters at each
            # locus, where each set of 2+ size is flanked in square braces.
            snp_list.append(regex_char.strip("[]"))
        # Return an interator for snp_list's possible combinations. This iterator
        # is *not* calculated here, it is calculated on-the-fly as it is used.
        return (''.join(x) for x in itertools.product(*snp_list))
