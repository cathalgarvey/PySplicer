from pysplicer import sequtils
from pysplicer import translators
from pysplicer import DNAMapper
import textwrap # For errmsgs.
#import copy

class Substitutor:
    def __init__(self, table, mapper, verbose=False):
        self.rtable = table
        self.mapper = mapper
        self.verbose = verbose
        # Aliases
        self.ftable = self.rtable.for_table
        translate = self.ftable.translate_codon_list_to_str
        reverse_translate = self.rtable.get_codons

    def verbose_msg(self, *args, **kwargs):
        if self.verbose: print(*args, **kwargs)

    def shuffle_codons(self, codon_list, avoids=[]):
        '''Given a codon list and a set of synonymous codon lists to avoid, return first non-avoided one.'''
        # This method inverts codon_list so that iterator cycles primarily on the
        # *first* codon and not the last, which may prevent wasting lots of cycles
        # on codon permutations far from the most affected codons at the start.
        # It then re-assigns "codon_list" to the first permutation found that doesn't
        # occur in "avoids", a list of (presumably) previous attempts with the present
        # attempt added automatically.
        lavoids = [x[::-1] for x in avoids]# copy.deepcopy(avoids)
        lcodon_list = codon_list[::-1] # Copy, just in case changes propagate back to caller.
        lavoids.append(lcodon_list)
        aminos = self.ftable.translate_codon_list_to_str(lcodon_list)
        # Go through permutations until a non-borky permutation is found. If none, cry.
        for perm in self.rtable.get_codon_permutations(aminos):
            perm = list(perm)
            if perm in lavoids: continue
            lcodon_list = perm
            break
        else:
            raise ValueError(("No permutation of codon_list could be found that"
                              " was not already in avoids list!"))
        # Reverse back to correct order!
        lcodon_list.reverse()
        return lcodon_list

    def remove_mapped(self, codon_list):
        '''Attempts to substitute codons at where mapper detects sites of interest.
        Uses sequtils.bioseq methods, embedded mapper and translation tables.'''
        # 1. Map codon_list.
        # 2. For each codon in the mapped region, iterate over valid permutations
        #     (according to the permutation generator in self.rtable), remapping
        #     the results and ending iteration only when a match at the *current*
        #     locus is not made.
        # 3. Remap codon list for each change made.
        # 4. If a codon has no permutations for which a site mapping can be excluded,
        #     don't change, print to stderr with a Warning and continue.
        # =====
        # Keep a copy of the original for the final sanity-check.
        original_codons = codon_list[:]
        # Mapping (1):
        map_results = self.mapper.codon_indices(codon_list)
        self.verbose_msg("Attempting to resolve mapped problems at these loci:",*map_results.keys())
        # Iterating (2):
        # A list of codons we couldn't fix, so we'll ignore them when fixing future ones.
        ignored_codons = []
        # Can't iterate over map_results as we'll be editing it as we go. However,
        # as map_results is index-keyed, using codon_index works just as well.
        for codon_index in range(0,len(codon_list)):
            # Start: skip OK codons, break at end of sequence.
            if codon_index not in map_results: continue
            self.verbose_msg("At codon",str(codon_index)+".")
            # Middle: Iterate over permutations, avoiding those that map an issue
            #     at the same or previous loci (but not next, we don't care about that yet).
            attempted_perms = []
            codons_to_shuffle = codon_list[codon_index:codon_index+3]
            self.verbose_msg(" Codons to be shuffled:",codons_to_shuffle)
            number_of_permutations = self.rtable.count_codon_permus(codons_to_shuffle)
            self.verbose_msg(" There are",number_of_permutations,"synonymous",
                             "permutations for this and the next two codons.")
            for i in range(0,number_of_permutations-1):
                # Shuffle 3 codons, balance odds of killing sites vs. computer time.
                this_perm = self.shuffle_codons(codons_to_shuffle,attempted_perms)
                perm_in_context = codon_list[:codon_index]+this_perm+codon_list[codon_index+3:]
                perm_map = self.mapper.codon_indices(perm_in_context)
                perm_map_keys = [x for x in sorted(perm_map.keys()) if x not in ignored_codons]
                # If this didn't solve the issue, or introduced new/different issues:
                if codon_index in perm_map_keys:
                   # self.verbose_msg(" Permutation",this_perm,"introduced same"
                   #     "or new issues at current locus. Skipping.")
                    attempted_perms.append(this_perm)
                # OR if this introduced a new issue to a prior codon not previously
                # written off as intractable (filtered those out above):
                elif [x for x in perm_map_keys if x < codon_index]:
                   # self.verbose_msg(" Permutation",this_perm,"introduced new"
                   #     "issues at prior sites. Skipping.")
                    attempted_perms.append(this_perm)
                # Otherwise, we win! Break and continue.
                else:
                    # Replace whole codon_list with the "fix".
                    self.verbose_msg(" Found OK permutation:",this_perm)
                    codon_list = perm_in_context
                    break
            else:
                # With for-loops, else-clauses only get executed if nothing breaks
                # the for-loop prior to end-of-iteration. So, if that happens, it
                # means none of the attempted permutations worked.
                self.verbose_msg(" Have exhausted possibilities without avail."
                    "Skipping and ignoring this locus.")
                ignored_codons.append(codon_index)

            # End: remap in case future codons have new issues due to changes.
            map_results = self.mapper.codon_indices(codon_list)

        # Finished the iteration, should be OK but we'll just sanity-check:
        original_aminos = self.ftable.translate_codon_list_to_str(original_codons)
        new_aminos = self.ftable.translate_codon_list_to_str(codon_list)
        if original_aminos != new_aminos:
            errMsg = textwrap.dedent("""\
                After codon substitutions, codon set is no longer synonymous
                with input:
                \tOriginal Aminos:\t{0}
                \tNew Aminos:\t\t{1}
                \tOriginal Codons:\t{2}
                \tNew codons:\t\t{3}""".format(original_aminos,new_aminos,
                                '-'.join(original_codons),'-'.join(codon_list)))
            raise ValueError(errMsg)

        self.verbose_msg("Ignored codons with remaining issues:",ignored_codons)
        self.verbose_msg("Remaining issues:",*self.mapper.codon_indices(codon_list).items(),sep="\n\t")
        return codon_list

def tests(aminos):
    ecoli_table='{"A": {"GCA": 0.24242424242424243, "GCC": 0.1919191919191919, "GCT": 0.12121212121212122, "GCG": 0.4444444444444444}, "C": {"TGC": 0.5772005772005773, "TGT": 0.42279942279942284}, "E": {"GAG": 0.29558634020618557, "GAA": 0.7044136597938144}, "D": {"GAT": 0.46, "GAC": 0.54}, "G": {"GGT": 0.6060606060606061, "GGG": 0.0, "GGA": 0.0, "GGC": 0.3939393939393939}, "F": {"TTC": 0.55, "TTT": 0.45}, "I": {"ATT": 0.5820752914198357, "ATC": 0.34702847315115615, "ATA": 0.07089623542900822}, "H": {"CAC": 0.45275181723779856, "CAT": 0.5472481827622014}, "K": {"AAG": 0.49, "AAA": 0.51}, "M": {"ATG": 1.0}, "L": {"CTT": 0.02, "CTG": 0.78, "CTA": 0.0, "CTC": 0.03, "TTA": 0.03, "TTG": 0.14}, "N": {"AAT": 0.4726604711476119, "AAC": 0.5273395288523882}, "Q": {"CAA": 0.45, "CAG": 0.55}, "P": {"CCT": 0.09, "CCG": 0.81, "CCA": 0.1, "CCC": 0.0}, "S": {"TCT": 0.10204081632653061, "AGC": 0.6938775510204082, "TCG": 0.05102040816326531, "AGT": 0.0, "TCC": 0.1326530612244898, "TCA": 0.02040816326530612}, "R": {"AGG": 0.02671690357938003, "CGC": 0.4447679397157047, "CGG": 0.07021750299708854, "CGA": 0.073642747045727, "AGA": 0.023462921733173492, "CGT": 0.3611919849289262}, "T": {"ACC": 0.57, "ACA": 0.0, "ACG": 0.33, "ACT": 0.09999999999999999}, "W": {"TGG": 1.0}, "V": {"GTA": 0.1735462488701416, "GTC": 0.17640855679421516, "GTT": 0.2529376318168123, "GTG": 0.397107562518831}, "Y": {"TAT": 0.5342029907731467, "TAC": 0.46579700922685335}, "*": {"TAG": 0.0, "TGA": 0.3576642335766423, "TAA": 0.6423357664233577}}'
    import random
    print("\nRunning tests for Substitutor.")
    testtable = translators.ReverseTranslator(ecoli_table)
    maptargets = [sequtils.random_iupac_dna(random.randint(5,7)) for x in range(2,10)]
    print("Using the optimised E.coli table as default.")
    print("Targeted sites to edit away are:",'; '.join(maptargets))
    testmapper = DNAMapper.DNAMapper(maptargets)
    Subber = Substitutor(table=testtable, mapper=testmapper, verbose=True)
    test_codon_list = testtable.get_codons(aminos)
    scrubbed_codons = Subber.remove_mapped(test_codon_list)
    print("Final codon list:",scrubbed_codons)
    print("Translates to:",testtable.for_table.translate_codon_list_to_str(scrubbed_codons))
    print("Original Heatmap:\t",'-'.join([str(int(x*100)).zfill(2) for x in testtable.heat_map(test_codon_list)]))
    print("Heatmap:\t\t",'-'.join([str(int(x*100)).zfill(2) for x in testtable.heat_map(scrubbed_codons)]))

if __name__ == "__main__":
    import sys
    tests(sys.argv[1])
