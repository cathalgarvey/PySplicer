import math
import itertools
from pysplicer import sequtils

class DNAMapper:
    '''Created with a list of IUPAC Sequence Patterns, this provides a set of
    methods for analysing sequence objects or codon lists to find matched sites.
    It can be used for sequence analysis, or to help identify problem areas
    during automated DNA design for synthesis, for example by identifying sites
    that may be recognised by endogenous factors in the target host so they can
    be replaced if possible by synonymous sequences.'''

    @staticmethod
    def uniquify(seq):
        'Fast, order-preserving uniquifier, only works on sequences of immuntable stuff (strings).'
        seen = set() # From http://www.peterbe.com/plog/uniqifiers-benchmark
        return [x for x in seq if x not in seen and not seen.add(x)]

    @staticmethod
    def allindices(string, sub, offset=0):
        'Returns a list of ALL indices of a substring, not just first.'
        listindex = []
        i = string.find(sub, offset)
        while i >= 0: # Because find returns -1 when it finds nothing.
            listindex.append(i)
            i = string.find(sub, i + 1)
        return listindex

    def verbose_msg(self, *args, **kwargs):
        if self.verbose: print(*args, **kwargs)

    def __init__(self, map_sites, search_complement=True, structure_mapper=None, verbose=False):
        '''Create a DNAMapper using the IUPAC pattern list in map_sites.
        map_sites should be a list of strings in IUPAC format.
        search_complement directs the object to build regexes that will also match reverse complement.
        structure_mapper, if provided, will be called with sequence as argument
            to get a locus-indexed dict of secondary structures found in the sequence.
        verbose instructs the object to print debug info and return redundant information when mapping codons.
        '''
        self.verbose = verbose
        self.map_sites = map_sites
        self.search_complement = search_complement
        self.structure_mapper = structure_mapper or self._dummy_sequence_mapper
        self.regexes = self.build_regex_list(self.map_sites)
        self.verbose_msg("Created DNAMapper object with patterns:",[x.pattern for x in self.regexes])

    def _dummy_sequence_mapper(self, *args, **kwargs):
        'If no real mapper is given, this takes its place and just returns an empty dict.'
        # Mappers which *do* return something should return
        return {}

    def build_regex_list(self, iupac_list):
        regex_list = []
        for iupac_pattern in iupac_list:
            pattern_object = sequtils.bioseq(iupac_pattern.upper())
            regex_list.append(pattern_object.iupac_to_regex(
                                include_complement=self.search_complement,
                                strip_wildcards=True))
        return regex_list

    def map_string(self, query):
        '''Given a query string, this method returns a dictionary of matches, where
        each key is an integer corresponding to the string index of the match, and
        the value is a sub-dict containing the following:
        - patterns: [regex_patterns..] # Patterns matching at this locus
        - span: int # Length of the largest match at this locus, to aid in recovering match.'''
        query = query.upper()
        results = {}
        for pattern in self.regexes:
            matches = pattern.findall(query)
            for match in matches:
                # Regexes return the sequence they match, which doesn't always
                # match the pattern itself if it's not a simple string, so find
                # each substring matching the regex match.
                matchindices = self.allindices(query, match)
                for sequence_locus in matchindices:
                    # setdefault allows us to pretend that the key already exists
                    # as the second argument if it hasn't already been set.
                    # First: Get (create?) the sub-dict for this sequence locus:
                    dict_locus = results.setdefault(sequence_locus, {})
                    # Then: Append to (after creating?) the list of pattern matches.
                    # Duplicates are removed below using the fast uniquify static method.
                    dict_locus.setdefault("patterns",[]).append(pattern.pattern)
                    # If the span of this match is greater than existing one, replace prior.
                    dict_locus['span'] = max(dict_locus.get("span",0),len(match))
        # Without a mapper given at instantiation, structure_mapper is a dummy method
        # that returns an empty dict, so the below codeblock will be skipped.
        structural_results = self.structure_mapper(query)
        for structure_locus in structural_results:
            # Make sure it doesn't overwrite existing entries
            structure_d = structural_results[structure_locus]
            dict_locus = results.setdefault(structure_locus, {})
            dict_locus['span'] = max(dict_locus.get("span",0), structure_d['span'])
            hairpin_r_nucleotides = []
            for hairpin_span in structure_d.get('hairpin',[]):
                # each subunit of "hairpin" key is a tuple of (span, r_start, f_start, score, [(basepair),(basepair)])
                if not hairpin_span: continue
                    # hairpin base pairs are ordered outward from hairpin loop, so
                    # the proximal "r" bases are ordered "backwards". So, copy backwards.
                this_span = [x[0] for x in hairpin_span[4][::-1]]
                hairpin_r_nucleotides.extend(this_span)
                hairpin_r_nucleotides.append("=")
            dict_locus.setdefault('patterns',[]).append(''.join(hairpin_r_nucleotides)+"=0")

        # Remove duplicate pattern matches:
        for l in results:
            results[l]['patterns'] = self.uniquify(results[l]['patterns'])
        return results

    @staticmethod
    def codon_index(index, frame=0):
        'Returns which codon a nucleotide falls into, from the start of an arbitrary sequence at the given frameshift.'
        # Currently truncates start of sequence, so "codons" start counting at this start.
        # This can return negative codons where index is less than frameshift.
        return math.floor((index-frame)/3)

    def codon_indices(self, sequence, frame=0):
        '''Wraps map_string to provide codon-indexed, and codon-spanning, data.
        If given a string and optional frame, then converts map_string output to
        codon-mapped data. If given a list of codons, assumes frame 0 and does
        likewise. Can return negative codon indices for positive frames if map_string
        returns mapped sites in the frame-truncated leader sequence.
        Wheareas map_string returns indices for matches at the first nucleotide,
        this function may group together matches within the same codon. That is,
        the returned set of patterns may have matched any of the nucleotides in
        the first codon, not just the first.
        Returned dict is of form (int):{"span":(int) codons, "patterns":(list) regexes},
        for example {2:{'spans':3,'patterns':['GC[AT][AT][AT]GCGG']}'''
        # Accept either str or list for sequence, but convert to str.
        if isinstance(sequence, list):
            sequence = ''.join(sequence)
            frame = 0 # If codons are given then it is assumed they are in frame
        if not isinstance(sequence, str): # Is true of str subclasses like bioseq.
            raise TypeError("codon_indices accepts either string or list-of-codons input.")
        self.verbose_msg("codon_indices called on sequence",
            ' '.join([str(x[1])+'-'+x[0] for x in zip(sequtils.bioseq(sequence).codons(),range(0,100))]),
            ", frame",frame)
        # First get string indices, patterns and spans.
        stringindices = self.map_string(sequence)
        # This is int-indexed dicts containing keys "span" and "patterns".
        # Need to convert the indexes and "span" value to codons instead.
        # Do this by converting index and index+span into codon indices,
        # then converting span into codon_index(index+span)-codon_index(index).
        codonindices = {}
        for index in stringindices:
            if index < frame: continue # Ignore string indices "before" the current frame.
            self.verbose_msg("Mapping hit at string index",index,
                 "({0}) to codons.".format(sequence[abs(index-2):abs(index+3)]))
            nuc_span = stringindices[index]['span']
            nuc_patterns = stringindices[index].get('patterns',[])
            nuc_structures = stringindices[index].get('structures',[])
            codon_n = self.codon_index(index,frame)
            span_nuc_n = index+nuc_span + 1
            span_codon_n = self.codon_index(span_nuc_n,frame)
            spans_codons = 1 + span_codon_n - codon_n
            # From here, have to be careful not to overwrite entries from prior
            # nucleotides in the same codon!
            codon_entry = codonindices.setdefault(codon_n, {})
            # Dicts being mutable, it should be OK to edit codon_entry and have
            # changes map back to its value in codonindices without further reference.
            codon_entry.setdefault("patterns",[]).extend(nuc_patterns)
            codon_entry['span'] = max(codon_entry.get("span",0),spans_codons)
            if self.verbose:
                # In verbose mode, spit out an additional dict in each codon index
                # containing subdicts for each nucleotide index processed:
                nucdebug = codon_entry.setdefault('debug',{}).setdefault(index,{})
                nucdebug['sequence'] = sequence
                nucdebug['string_vs_codon'] = (index, codon_n)
                nucdebug['original'] = stringindices[index]
                nucdebug['span_nuc_n'] = span_nuc_n
                nucdebug['span_codon_n'] = span_codon_n
        return codonindices

class HairpinMapper:
    bpdict = {'A':'UT','C':'G','G':'UC','U':'GA','T':'A'}
    bp_scores = {
                 ("A","U"):2,
                 ("U","A"):2,
                 ("A","T"):2,
                 ("T","A"):2,
                 ("G","C"):3,
                 ("C","G"):3,
                 ("G","U"):1,
                 ("U","G"):1,
                 }

    def __init__(self, min_hairpin=5, endnum=3, min_score=10, max_loop=6, verbose=False):
        self.min_hairpin = min_hairpin
        self.endnum = endnum
        self.min_score = min_score
        self.max_loop = max_loop
        self.verbose = verbose

    def verbose_msg(self, *args, **kwargs):
        if self.verbose: print(*args, **kwargs)

    def bp(self, l,m):
        if m in self.bpdict.get(l,''): return True
        else: return False

    def _score_contig(self, bp_tuples):
        'Returns a rough measure of stem stability based on base-pair type and position relative to other BPs.'
        scores = []
        for bp in bp_tuples:
            try:
                scores.append(self.bp_scores[bp])
            except KeyError as E:
                raise KeyError("Base pair that isn't in bp_scores: "+E)
        return sum(scores)

    def _extend_contig(self, seq, r_init, f_init):
        'Given start-bases r and f, extends a base-paired contig along seq, returns number of uninterrupted base pairs encountered.'
        span = 0
        bps_in_contig = []
        try:
            while True:
                if self.bp(seq[r_init-span], seq[f_init+span]):
                    #self.verbose_msg("Extending contig..")
                    bps_in_contig.append((seq[r_init-span],seq[f_init+span]))
                    span += 1
                else:
                    break
        except IndexError:
            score = self._score_contig(bps_in_contig)
#            self.verbose_msg("Broken by indexerror, score is:",score)
            return span, r_init+1, f_init, score, bps_in_contig
        score = self._score_contig(bps_in_contig)
        # Returns span and init points so method can be used with anonymous inputs
        # and extract the one that worked with max().
#        self.verbose_msg("Final score is:",score)
        # Increment r_init by one because usual string usage after this method
        # will involve slicing rather than specific indexing, and not doing so
        # gives the incorrect index for dictionary indexing and manual lookup.
        return span, r_init+1, f_init, score, bps_in_contig

    def _map_hairpin(self, seq, index, boffs=0, foffs=0):
        'Walks along a potential hairpin using contig-extension, favouring longest contigs found at each step.'
        empty_span = 0
        initial_tolerance = self.max_loop
        r_current, f_current = index - boffs, index + foffs
        # Begin as you mean to continue: test the three "normal" "frames" of a hairpin
        # and pick the one with the best contig or inc/dec/rement the location "pins"
        found_hairpin = []
        # While we have seen fewer nucleotides in a gap than endnum:
        while empty_span < self.endnum or empty_span < initial_tolerance:
            # Assemble a set of queries where each represents a different possible path forward.
            # Query spans is a list of r/f startpoints at every possible offset within 0,self.endnum,
            # which effectively defines the size of bulges or kinks allowed within a hairpin.
            query_spans =  [(r_current-z,    f_current+z)    for z in range(0,self.endnum)]
            query_spans += [(r_current-z[0], f_current+z[1]) for z in itertools.permutations(range(0,self.endnum),2)]
            next_span = max(self._extend_contig(seq, *x) for x in query_spans)

            if next_span[0] >= self.min_hairpin:
                # If this is the first match, yay! No further need for initial tolerance.
                # (%timeit in ipython reveals that dumb assignment is faster than
                # conditional assignment like if initial_tolerance > 0: initial_tolerance = 0.
                initial_tolerance = 0
                # Add to growing hairpin set.
#                self.verbose_msg("Found",next_span[0],"bp hairpin during extension:",next_span)
                found_hairpin.append(next_span)
                r_current = next_span[1] - next_span[0]
                f_current = next_span[2] + next_span[0]
                empty_span = 0
                continue
            else:
                empty_span += 1
                r_current -= 1
                f_current += 1
#        self.verbose_msg("Final hairpin:",found_hairpin)
        return found_hairpin

    def map_hairpins(self, seq_to_map, *args, **kwargs):
        'Builds and returns a list of hairpins found in query sequence. Intensive; char-by-char.'
        found_hairpins = {}
        for charN in range(0,len(seq_to_map)):
            this_hp = self._map_hairpin(seq_to_map, charN, *args, **kwargs)
            if not this_hp: continue
            last_span = this_hp[len(this_hp)-1]
            hp_begins = last_span[1] - last_span[0]
            hp_ends = last_span[2] + last_span[0]
            hp_span = hp_ends - hp_begins
            score = sum([x[3] for x in this_hp])
            if score > self.min_score:
                if hp_begins not in found_hairpins:
                    # This returns a "patterns" key so as to be drop-in compatible with DNAMapper.
                    self.verbose_msg("Found hairpin at",hp_begins,"spanning",hp_span,"nucleotides, with a score of",score)
                    found_hairpins[hp_begins] = {"span":hp_span, "score":score, "patterns":[], "hairpin":this_hp}
                else:
                    self.verbose_msg("Hairpin at index",hp_begins,
                        "already in dict, but should not be a duplicate?\nDict contains:",
                        found_hairpins[hp_begins], "\nNew Entry:", this_hp)
            else:
                self.verbose_msg("Hairpin at index",hp_begins,"rejected as its score of",score,
                    "is below the minimum threshold of",self.min_score,":",this_hp)

        for key in sorted(found_hairpins.keys()):
            # Discard keys that are so early they can only have been accepted by
            # rolling back through the sequence, not allowed unless sequence is
            # circular, and even if it is.. that's not properly handled by the code
            # yet. Better not to even accomodate a thing if it's not fully correct!
            if key < self.min_hairpin:
                found_hairpins.pop(key)
            else:
                # sorted keys, so break at first acceptable key.
                break
        return found_hairpins

    def print_hairpin(self, hairpin_list):
        'If given a hairpin, returns the base-paired regions of the hairpin for pretty-printing.'
        out_string_lines = [[" /¯"],["{  "],[" \\_"]]
        # Want: measurements of loop-size, and bulge-size.
        # Loop size should be forward_start - backward_start of first stem:
        loop_size = hairpin_list[0][2] - hairpin_list[0][1]
        # Used to contain easily accessed data on when last contig ended,
        # for comparison with start of new contig, and therefore useful output.
        prior_contig_ends = (0,0)
        for hp in hairpin_list:
            # If this is the second or later but not last stem of the hairpin, then add the
            # number of bases on either strand prior to this stem:
            if 0 < hairpin_list.index(hp):
                prior_hp = hairpin_list[hairpin_list.index(hp) - 1]
                # On the "reverse" end, string indexes are getting lower, while
                # on the "forward" end they're getting larger. Thankfully, we
                # can ignore this and use abs() to get a positive value of any difference.
                r_span_since_contig = abs(prior_hp[1] - prior_hp[0]) - hp[1]
                f_span_since_contig = abs(prior_hp[2] + prior_hp[0] - hp[2])
                out_string_lines[0].append("-"+str(r_span_since_contig)+"-")
                out_string_lines[1].append("   ")
                out_string_lines[2].append("-"+str(f_span_since_contig)+"-")
            # hp index 4 is the list of base-pair tuples, where for list x:
            #   x30 x20 x10 x00      x01 x11 x21 x31
            #   (   (   (   (        )   )   )   )
            # ======================================
            for x in hp[4]:
                out_string_lines[0].append(x[0])
                out_string_lines[1].append("|")
                out_string_lines[2].append(x[1])
            prior_contig_ends = (hp[1],hp[2])
        return '\n'.join([''.join(x).rstrip(" .") for x in out_string_lines])

def tests():
    print("Running tests on DNAMapper.")
    import re, random
    testsequence = sequtils.bioseq(sequtils.random_dna(100))
    maptargets = [sequtils.random_iupac_dna(random.randint(2,10)) for x in range(2,7)]
    mapper = DNAMapper(maptargets,verbose=True)

    def verify_stringmap():
        result = True
        string_map = mapper.map_string(testsequence)
        print("\t ----- Checking Hits for Consistency ----- ")
        for x in string_map:
            x_dict = string_map[x]
            for pattern in x_dict['patterns']:
                thisre = re.compile(pattern)
                if not thisre.match(testsequence[x:x+x_dict['span']]):
                    print("\tFailed to match",pattern,"with substring",testsequence[x:])
                    result = False
                else:
                    print("\tMatched",pattern,"successfully at target site.")
        return result

    def verify_framemap(frame):
        result = True
        #mycodons = testsequence.as_codon_list(trailing=True)
        mycodons = testsequence.codons(leading=False,trailing=True) # Leading?
        codon_map = mapper.codon_indices(testsequence,frame)
        print("\t ----- Checking Hits for Consistency ----- ")
        for indexC in sorted(codon_map.keys()):
            region = ''.join(mycodons[indexC:indexC+codon_map[indexC]['span']])
            for pattern in codon_map[indexC]["patterns"]:
                thisre = re.compile(pattern)
                if not thisre.findall(region):
                    print("\tNon-match with pattern",thisre.pattern,"in string",region)
                    print("\t\tDebug dict:",codon_map[indexC]['debug'])
                    result = False
                else:
                    pass
                    #print("\tMatch with pattern",thisre.pattern,"in string",region)
        return result

    allpassing = True
    if verify_stringmap():
        print("Passed stringmap tests.\n","==============================\n")
    else:
        print("Failed stringmap tests.\n","==============================\n")
        allpassing = False
    for x in range(0,3):
        if verify_framemap(x):
            print("Passed framemap test at frame",str(x)+".\n","==============================\n")
        else:
            print("Failed framemap test at frame",str(x)+".\n","==============================\n")
            allpassing = False
    if allpassing: print("\n >>> All Tests Passed! <<<")
    else: print("\n >>> Tests Failed <<<")

if __name__ == "__main__": tests()
