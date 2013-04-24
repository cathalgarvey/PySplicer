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

    def __init__(self, min_hairpin=5, endnum=3, min_score=10, max_loop=9, debug=False):
        self.min_hairpin = min_hairpin
        self.endnum = endnum
        self.min_score = min_score
        self.max_loop = max_loop
        self.debug = debug

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
                    if self.debug: print("Extending contig..")
                    bps_in_contig.append((seq[r_init-span],seq[f_init+span]))
                    span += 1
                else:
                    break
        except IndexError:
            score = self._score_contig(bps_in_contig)
            if self.debug: print("Broken by indexerror, score is:",score)
            return span, r_init, f_init, score, bps_in_contig
        score = self._score_contig(bps_in_contig)
        # Returns span and init points so method can be used with anonymous inputs
        # and extract the one that worked with max().
        if self.debug: print("Final score is:",score)
        return span, r_init, f_init, score, bps_in_contig

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
            # What I want: a set of (x+z,y+w) where z/w are incremented between 1 and self.endnum

            if next_span[0] >= self.min_hairpin:
                # If this is the first match, yay! No further need for initial tolerance.
                # (%timeit in ipython reveals that dumb assignment is faster than
                # conditional assignment like if initial_tolerance > 0: initial_tolerance = 0.
                initial_tolerance = 0
                # Add to growing hairpin set.
                if self.debug: print("Found",next_span[0],"bp hairpin during extension:",next_span)
                found_hairpin.append(next_span)
                r_current = next_span[1] - next_span[0]
                f_current = next_span[2] + next_span[0]
                empty_span = 0
                continue
            else:
                empty_span += 1
                r_current -= 1
                f_current += 1
        if self.debug: print("Final hairpin:",found_hairpin)
        return found_hairpin

    def map_hairpins(self, seq_to_map, *args, **kwargs):
        'Builds and returns a list of hairpins found in query sequence. Intensive; char-by-char.'
        found_hairpins = {}
        for charN in range(0,len(seq_to_map)):
            this_hp = self._map_hairpin(seq_to_map, charN, *args, **kwargs)
            if not this_hp: continue
            hp_begins = this_hp[len(this_hp)-1][1]
            hp_ends = this_hp[len(this_hp)-1][2]
            hp_span = hp_ends - hp_begins
            score = sum([x[3] for x in this_hp])
            if score > self.min_score:
                if hp_begins not in found_hairpins:
                    # This returns a "patterns" key so as to be drop-in compatible with DNAMapper.
                    found_hairpins[hp_begins] = {"span":hp_span, "score":score, "patterns":[], "hairpin":this_hp}
                else:
                    if self.debug: print("Hairpin at index",hp_begins,
                        "already in dict, but should not be a duplicate?\nDict contains:",
                        found_hairpins[hp_begins], "\nNew Entry:", this_hp)
            elif self.debug:
                print("Hairpin at index",hp_begins,"rejected as its score of",score,
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
