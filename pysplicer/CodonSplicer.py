from pysplicer import sequtils

import sys

class Splicer:
    def __init__(self, table, mapper, verbose=False):
        '''Expects a translators.ReverseTranslator table to farm out codon operations,
        and a DNAMapper.DNAMapper object to find objectionable codons/areas to avoid.'''
        self.mapper = mapper
        self.rtable = table
        # Alias for convenience.
        self.ftable = self.rtable.for_table
        self.verbose = verbose

    def verbose_msg(self, *args, **kwargs):
        if self.verbose:
            print(*args,**kwargs)

    def assert_synonymity(self, seq_list):
        'This function checks that seq_list is a non-zero list and that contents are valid and equivalent codon lists.'
        if not isinstance(seq_list, list) or not seq_list:
            raise TypeError("List that should contain synonymous codon lists is either not a list or empty:\n\n"+ str(seq_list))
        for i in seq_list:
            if not isinstance(i, list) or not i:
                raise TypeError("Entry in synonymous codon-list set is either not a list or is empty:\n\n"+str(i))
        # Use a method of the embedded forward table of the embedded reverse table.. D:
        First = self.ftable.translate_codon_list_to_str(seq_list[0])
        for i in seq_list:
            if not self.ftable.translate_codon_list_to_str(i) == First:
                raise ValueError("Not all sequences in synonymous codon-sequence set are identical when translated.")
        return True

    def find_best_span(self, work_list, codon_offset):
        'Working from codon_offset, compare each sequence in work_list to find one with longest-stretch-til-next-issue.'
        # This will be the index of which sequence has best stretch of clear codons.
        at_end = False
        best_so_far = {'Seq':-1, 'codon':-1}
        for i in work_list:
            # Filter out indices smaller than the current offset.
            i_indices = [x for x in filter(lambda x: x >= codon_offset, i['mapIndices'])]
#            print("For sequence",i,"undesired site indexes are:\n",str(i_indices))
            if not i_indices:
                # If, in other words, there are no remaining issues with this
                # codon list: declare the "codon" entry to be the last one in
                # the list. For some reason, keep iterating over other lists
                best_so_far = {'Seq':work_list.index(i), "codon":len(i['fullcodonlist'])-1}
                at_end = True
                break
            if i_indices[0] > best_so_far['codon']:
                best_so_far = {'Seq':work_list.index(i), 'codon':i_indices[0]}
        if (best_so_far['codon'] == codon_offset) or (best_so_far['Seq'] < 0):
            # This means we couldn't get past even the start
            # so we instead return next codon anyway.
            if not at_end:
                print("Warning:\tCould not avoid inclusion of mapped issue at position",
                  codon_offset,"during splicing.\n\t\tThis error may be fixed at a "+\
                  "later stage in the optimisation process.",file=sys.stderr)
            self.verbose_msg("Index of best list:",best_so_far['Seq'],
                    "; sequence:",work_list[best_so_far["Seq"]]['fullcodonlist'][abs(codon_offset):abs(codon_offset)+1])
            return work_list[best_so_far['Seq']]['fullcodonlist'][codon_offset:codon_offset+1]
        # Return the longest found span this round..
        self.verbose_msg("Index of best list:",best_so_far['Seq'],
                "; sequence:",work_list[best_so_far["Seq"]]['fullcodonlist'][abs(codon_offset):best_so_far['codon']])
        best_full_length = work_list[ best_so_far['Seq'] ][ 'fullcodonlist' ]
        next_codons = best_full_length[ codon_offset : best_so_far['codon'] ]
        return next_codons

    def splice_sequences(self, sequence_list):
        '''Takes a list of candidate sequences (which must encode same amino sequence)
        and attempts to splice a site-free form.

        This method uses the given DNAMapper to identify codon sites to be avoided
        although it actually doesn't use the content of the returned dictionary,
        just the dictionary indices.

        It does this by repeatedly scanning from the "current position" forward
        through every given synonymous codon list, returning the next longest
        span of codons without problems.

        In doing so, three problems may arise:
            A) The method may encounter an issue/site that is universal to all
                sequences, and it must ignore said site and continue.
            B) The method may create new issues/sites at the splice junctions,
                but resolving these issues can be left to another function/method.
        '''
        self.assert_synonymity(sequence_list)
        seq_len = len(sequence_list[0]) # Used later for a while loop.
        # Get an ordered list of codons to work around for each sequence "i" of sequence_list.
        working_list = []
        for i in sequence_list:
            working_list.append( { "mapIndices": sorted(self.mapper.codon_indices(i).keys()),
                                   "fullcodonlist": i } )
        new_codons = []
        # While new codons aren't as long as one of the (already known to be synonymous) sequences,
        # keep appending "next best" codons as returned by self.find_best_span.
        while len(new_codons) < len(sequence_list[0]):
            # Use len(new_codons) because we're starting at the next index after the
            # end of the previous one; the increment of 1 implied by len vs. list indices
            # is a feature, not a bug. :)
            # NB: using %timeit in ipython, use of += operator appears faster than list.extend().
            new_codons += self.find_best_span(working_list, len(new_codons))

        self.verbose_msg("Splicer: Sanity checking that new_codons matches amino sequence of input sequence set..")
        try:
            self.assert_synonymity([new_codons, sequence_list[0]])
        except ValueError as E:
            # To make it clear this bug happened *after* everything else, and not during
            # initial sanity-check call to this function.
            raise ValueError("Error occurred while verifying finished spliced-sequence: "+str(E))
        return new_codons

def splice_aminos_to_codons(aminostring, table, mapper, attempts_number=50, suggestions = [], verbose=False):
    '''Creates attempts_number versions of aminostring using the given table,
    then splices them to find one with minimal issues as determined by mapper.'''
    splicer = Splicer(table, mapper, verbose)
    query_sequences = [table.get_codons(aminostring) for x in range(0,attempts_number)]
    if suggestions: query_sequences + suggestions
    return splicer.splice_sequences(query_sequences)

def tests(aminos):
    ecoli_table='{"A": {"GCA": 0.24242424242424243, "GCC": 0.1919191919191919, "GCT": 0.12121212121212122, "GCG": 0.4444444444444444}, "C": {"TGC": 0.5772005772005773, "TGT": 0.42279942279942284}, "E": {"GAG": 0.29558634020618557, "GAA": 0.7044136597938144}, "D": {"GAT": 0.46, "GAC": 0.54}, "G": {"GGT": 0.6060606060606061, "GGG": 0.0, "GGA": 0.0, "GGC": 0.3939393939393939}, "F": {"TTC": 0.55, "TTT": 0.45}, "I": {"ATT": 0.5820752914198357, "ATC": 0.34702847315115615, "ATA": 0.07089623542900822}, "H": {"CAC": 0.45275181723779856, "CAT": 0.5472481827622014}, "K": {"AAG": 0.49, "AAA": 0.51}, "M": {"ATG": 1.0}, "L": {"CTT": 0.02, "CTG": 0.78, "CTA": 0.0, "CTC": 0.03, "TTA": 0.03, "TTG": 0.14}, "N": {"AAT": 0.4726604711476119, "AAC": 0.5273395288523882}, "Q": {"CAA": 0.45, "CAG": 0.55}, "P": {"CCT": 0.09, "CCG": 0.81, "CCA": 0.1, "CCC": 0.0}, "S": {"TCT": 0.10204081632653061, "AGC": 0.6938775510204082, "TCG": 0.05102040816326531, "AGT": 0.0, "TCC": 0.1326530612244898, "TCA": 0.02040816326530612}, "R": {"AGG": 0.02671690357938003, "CGC": 0.4447679397157047, "CGG": 0.07021750299708854, "CGA": 0.073642747045727, "AGA": 0.023462921733173492, "CGT": 0.3611919849289262}, "T": {"ACC": 0.57, "ACA": 0.0, "ACG": 0.33, "ACT": 0.09999999999999999}, "W": {"TGG": 1.0}, "V": {"GTA": 0.1735462488701416, "GTC": 0.17640855679421516, "GTT": 0.2529376318168123, "GTG": 0.397107562518831}, "Y": {"TAT": 0.5342029907731467, "TAC": 0.46579700922685335}, "*": {"TAG": 0.0, "TGA": 0.3576642335766423, "TAA": 0.6423357664233577}}'
    import DNAMapper
    import translators
    import random
    print("\nRunning tests for Splicer object and splice_amino_to_codons function.")
    testtable = translators.ReverseTranslator(ecoli_table)
    testtable.set_
    maptargets = [sequtils.random_iupac_dna(random.randint(6,10)) for x in range(2,7)]
    print("Using the optimised E.coli table as default.")
    print("Targeted sites to splice around are:",'; '.join(maptargets))
    testmapper = DNAMapper.DNAMapper(maptargets)
    print(''.join(splice_aminos_to_codons(aminos, testtable, testmapper, 100, verbose=True)),"\n")

if __name__ == "__main__":
    tests(sys.argv[1])
