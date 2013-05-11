import random
import collections
import json
import itertools
from pysplicer import sequtils
#from sequtils.translationtables import *
# Provides a set of tables of format:
# table24 = {"starts": ["TTG", "CTG", "ATG", "GTG"],
#            "aminos": {"A": ["GCA", "GCC", "GCG", "GCT"]...},
#            "codons": {"CTT": "L", "ACC": "T"...}}
# "starts" and "codons" subkeys can be directly used by Translator classes.
# Tables are presented by table number, and also by these aliases:
#  standard; vertebrate_mitochondrial; yeast_mitochondrial; mold_mitochondrial;
#  protozoan_mitochondrial; coelenterate_mitochondrial; invertebrate_mitochondrial;
#  ciliate; dasycladacean; hexamita; echinoderm_mitochondrial; flatworm_mitochondrial;
#  euplotid; bacterial; archaeal; plant_plastid; chloroplast; alternative_yeast;
#  alt_yeast; ascidian_mitochondrial; alternative_flatworm_mitochondrial;
#  alt_flatworm_mitochondrial; blepharisma; chlorophycean_mitochondrial;
#  trematode_mitochondrial; scenedesmus_obliquus_mitochondrial;
#  thraustochytrium_mitochondrial; pterobranchia_mitochondrial
default_table = sequtils.translationtables.bacterial

class Translator:
    @classmethod
    def from_reverse_table(cls, reverse_table):
        'Create a translator from the type of table used by ReverseTranslator.'
        table = {}
        for amino in reverse_table:
            for codon in reverse_table[amino]:
                table[codon] = amino
        return cls(table)

    @classmethod
    def from_translation_table(cls, translation_table):
        return cls(translation_table["codons"], translation_table["starts"])

    def __init__(self, translation_table, start_codons=['ATG']):
        'Give a table and a set of potential start codons (default ["ATG"]).'
        if isinstance(translation_table, str):
            translation_table = json.loads(translation_table)
        self.table = translation_table
        self.start = start_codons

    def translate_codon(self, codon):
        'Returns amino encoded by DNA/RNA codon.'
        return self.table[codon.upper().replace("U","T")]

    def translate_codon_list(self, codon_list):
        for codon in codon_list:
            yield self.translate_codon(codon)

    def translate_codon_list_to_str(self, codon_list):
        return ''.join(self.translate_codon_list(codon_list))

    def translate_string(self, string, frame=0):
        'Yields aminos for codons in the incoming string.'
        string = sequtils.bioseq.fromstring(string.upper().replace("U","T"))
        for codon in string.codons(frame, leading=False, trailing=False):
            yield self.translate_codon(codon)

    def get_orfs(self, string):
        'Uses "dumb" string comparison and bioseq.codon() to return ORFs in 3 frames.'
        string = sequtils.bioseq.fromstring(string.upper().replace("U","T"))
        buffering, amino_buffer = False, []
        orfs = [[],[],[]]
        for frame in range(0,3):
            frame_codons = string.codons(frame, leading=False, trailing=False)
            for codon in frame_codons:
                if codon in self.start: buffering = True
                # Yea, "elif": I'm sure there's at least one weird species out
                # there that uses a stop codon as an alternative start codon.
                elif self.translate_codon(codon) == "*":
                    buffering = False
                    amino_buffer.append("*")
                    orfs[frame].append(''.join(amino_buffer))
                    amino_buffer = []
                if buffering:
                    amino_buffer.append(self.translate_codon(codon))
        return orfs

    def sixframe(self, string):
        'Returns frames 0,1,2,c0,c1,c2, covering all six frames of "string".'
        string = sequtils.bioseq.fromstring(string.upper().replace("U","T"))
        orfs = self.get_orfs(string)
        orfs.extend(self.get_orfs(string.iupac_rev_complement()))
        return orfs

class ReverseTranslator:
    '''An object representing a codon usage frequency table, providing means to
    do weighted selection of individual or many aminos, with an embedded forward
    table generated from the codon frequency table. Can be instantiated from JSON,
    copy/pasted tables in standard format from the codon usage database (though
    species codon alphabet must be specified), or even from frequency-less forward
    tables in which all codons are assumed perfectly proportionate.'''
    def __init__(self, frequency_table, normalise=True, verbose=False):
        '''frequency_table should be an amino-indexed set of codon:frequency mappings, eg:
         {"A": {"GCA": 0.24, "GCC": 0.19, "GCU": 0.12, "GCG": 0.44}...}'''
        if isinstance(frequency_table, str):
            frequency_table = json.loads(frequency_table)
        if normalise: self.import_table(frequency_table)
        else: self.table = frequency_table
        # For convenience, embed a forward translator in the reverse transcriptor.
        self.for_table = Translator.from_reverse_table(self.table)
        self.verbose = verbose

    def verbose_msg(self, *args, **kwargs):
        if self.verbose: print(*args, **kwargs)


    @classmethod
    def from_forward_table(cls, forward_table):
        'Returns a "flat" reverse-translation table by reversing a forward table.'
        if isinstance(forward_table, str):
            forward_table = json.loads(forward_table)
        reversed_table = {}
        for codon in forward_table:
            encoded_amino = forward_table[codon]
            reversed_table.setdefault(encoded_amino,{})[codon] = 1.0
        return cls(reversed_table)

    @classmethod
    def from_rich_table(cls, richtable):
        'Parses a richer table format where codons point to dicts containing (at least) "frequency" subkeys.'
        # Earlier versions of pysplicer imported and exported this more
        # complex format to leave wiggle-room for additional data.
        if isinstance(richtable, str):
            richtable = json.loads(richtable)
        freq_dict = {}
        for amino in richtable:
            amino_subdict = {}
            for codon in richtable[amino]:
                codon_freq = richtable[amino][codon]['frequency']
                amino_subdict[codon] = codon_freq
            freq_dict[amino] = amino_subdict
        return cls(freq_dict, normalise=True)

    @classmethod
    def from_cud_table(cls, cud_table, translation_table=default_table['codons']):
        '''Imports frequency table for codons from a codon-usage-database formatted table.
        Needs a mapping of codons to aminos to use in calculating amino-set membership
        for each codon in the input table, which is needed for frequency normalisation.
        By default, it uses sequtils.translationtables.bacterial['codons'] for this.'''
        # Define a convenience class for attribute access of codon data.
        CodonInfo = collections.namedtuple("CodonInfo",["triplet","freq_per_1000","sample_number"])
        # Crack table by right-brackets, but only after stripping trailing right-bracket.
        table_split = cud_table.strip(")").split(")")
        table_data = []
        for x in table_split:
            # Crack by whitespace after killing remaining bracket
            codon_details_list = x.strip().replace("(", " ").split()
            if not len(codon_details_list)==3: raise ValueError("Error parsing CUD table: "+str(codon_details_list))
            # Amend RNA to DNA to match translation tables in sequtils
            codon_details_list[0] = codon_details_list[0].replace("U","T")
            # Import to appropriate datatypes
            codon_details_list[1] = float(codon_details_list[1])
            codon_details_list[2] = int(codon_details_list[2])
            # Create convenience named tuple
            Codon = CodonInfo(*codon_details_list)
            table_data.append(Codon)
        # Now parse codons, and use translation table to sort into amino-keyed dictionaries.
        reversed_table = {}
        for codon in table_data:
            amino_encoded = translation_table[codon.triplet]
            reversed_table.setdefault(amino_encoded,{})[codon.triplet] = codon.freq_per_1000
        # Leave it to the class initialiser to normalise the f/1000 values for each amino.
        return cls(reversed_table, normalise=True)

    def get_codon(self, amino, exclude_codons=[], exclude_pattern='', force_excludes=False):
        '''Performs weighted-selection of codons for target amino, optionally excluding
        specified codons or an IUPAC codon pattern. Exclusion may be useful for (e.g.)
        adding only non-NGG codons to the first portion of a CDS, where NGG can encourage
        translation abortion.'''
        # Copy because if exclude_codons is used, it'll be modifying the table!
        codons = self.table[amino].copy()
        # Expand exclude_pattern and extend exclude_codons with the results
        if exclude_pattern:
            exclude_codons.extend(sequtils.bioseq(exclude_pattern).all_iupac_permus("dna"))
        # Remove all codons we want to exclude
        if exclude_codons:
            for key in exclude_codons: codons.pop(key, None)
            if not sum([codons[x] for x in codons]):
                # That is, if there are no codons remaining, or none with frequency > 0!
                if force_excludes:
                    raise IndexError("No valid codons left for amino "+amino+" after exclusions: "+str(codons))
                else:
                    # Print a log of the problem and carry on. Can't win 'em all.
                    import sys
                    print("Warning: Could not exclude codons",codons,
                        ", and still select a codon for",amino+".",
                        "Have reverted to unmodified codon table in this case.",
                        file=sys.stderr)
                    codons = self.table[amino].copy()
            codons = self.normalise_codon_frequencies(codons)
        return self.random_category(codons)

    def get_codons(self, amino_string, *args, **kwargs):
        'Wraps get_codon to provide a list of results. Accepts same args and forwards.'
        return [self.get_codon(x, *args, **kwargs) for x in amino_string]

    def get_codon_permutations(self, aminos):
        'Returns an iterator that provides every synonymous permutation of codons for aminos.'
        if isinstance(aminos,list):
            aminos = self.for_table.translate_codon_list_to_str(aminos)
        perm_list = []
        for A in aminos:
            # This filters out codons that have a frequency of zero, either due
            # to natural bias or the use of remove_codon or set_min_frequency methods.
            Valid_Codons = filter(lambda x:bool(self.table[A][x]),self.table[A].keys())
            perm_list.append([x for x in Valid_Codons]) #Unpack iterator.
        return itertools.product(*perm_list)

    def count_codon_permus(self, aminos):
        'Counts available permutations of valid codons for a string of aminos.'
        if isinstance(aminos,list):
            aminos = self.for_table.translate_codon_list_to_str(aminos)
        p, perm_list = 1, []
        for A in aminos:
            # This filters out codons that have a frequency of zero, either due
            # to natural bias or the use of remove_codon or set_min_frequency methods.
            Valid_Codons = filter(lambda x:bool(self.table[A][x]),self.table[A].keys())
            perm_list.append(len([x for x in Valid_Codons])) #Unpack iterator.
        for n in perm_list:
            p *= n
        return p

    @staticmethod
    def random_category(prob_dict):
        '''Accepts a dictionary of {opt:wgt}, e.g.: random_category({'a':.15, 'b':.35, 'c':.5})
        Returns a selection based on weight. Weights must be normalised, should sum to 1!'''
        import random
        r = random.random() # range: 0,1
        total = 0
        for value,prob in prob_dict.items():
            if prob <= 0: continue # ignore items with a 0 weight.
            total += prob
            if total>r: return value
        raise ValueError('Distribution not normalized: {probs}'.format(probs=str(prob_dict)))

    def import_table(self, table):
        'Normalises every amino in a table.'
        self.table = table.copy()
        for amino in self.table:
            self.table[amino] = self.normalise_codon_frequencies(self.table[amino])

    def export_table(self, *args, **kwargs):
        return json.dumps(self.table, *args, **kwargs)

    @staticmethod
    def normalise_codon_frequencies(codon_table):
        'Given a dict of format {"XXX":float, "YYY":float}, returns with floats balanced to sum to 1.'
        renamed_table = {}
        for codon in codon_table:
            newcodon = codon.upper().replace("U","T")
            renamed_table[newcodon] = codon_table[codon]
        codon_table = renamed_table
        # Get current total:
        total_input_frequencies = sum(codon_table.values())
#        total_input_frequencies = 0.0
#        for codon in codon_table.keys():
#            total_input_frequencies += codon_table[codon]
        # Map new codon table with normalised values:
        total_output_frequencies = 0.0
        output_table = {}
        for codon in codon_table.keys():
            # Should give a proportional value between 0, 1, but with some error
            # due to C's handling of recurring floats, so give wiggle-room in
            # comparisons before raising errors.
            output_table[codon] = codon_table[codon] / total_input_frequencies
            if output_table[codon] > 1.0000001:
                raise ValueError("While attempting to re-normalise frequency "+\
                    "table, encountered error at codon "+codon+". Dumping dicts: "+\
                    str(codon_table)+" - Output so far: "+str(output_table))
            total_output_frequencies += output_table[codon]
        if not 0.99999 < total_output_frequencies < 1.00001:
            raise ValueError("Adjusted frequency total is outside accepted error: "+\
                str(total_output_frequencies)+", with table: "+str(output_table))
        return output_table
        
    def edit_codon(self, codon, frequency):
        'Manually set a codon frequency and automatically readjust synonymous codon frequencies.'
        if not 0 < frequency < 1:
            raise ValueError("Frequency to set manually must be between 0 and 1.")
        codon = codon.upper().replace("U","T")
        amino_id = self.for_table.translate_codon(codon)
        if amino_id == None:
            raise ValueError("Codon '"+codon+"' not found.")
        amino_set = self.table[amino_id]
        new_set = {codon:frequency}
        # First, normalise all the codons besides the one to be set:
        partial_set = {}
        for sub_codon in amino_set:
            if codon == sub_codon: continue
            partial_set[sub_codon] = amino_set[sub_codon]
        partial_set = self.normalise_codon_frequencies(partial_set)
        # Then, add each element of partial_set to new_set, adjusting the frequency
        # to fit into the remaining frequency space (1 - desired_codon_frequency).
        for sub_codon in partial_set:
            new_set[sub_codon] = partial_set[sub_codon] * (1 - frequency)
        self.table[amino_id] = new_set

    def remove_codon(self, codon):
        'Given a codon, find it, remove it, and re-normalise the remaining set.'
        codon = codon.upper().replace("U","T")
        amino_id = self.for_table.translate_codon(codon)
        if amino_id == None:
            raise ValueError("Codon '"+codon+"' not found.")
        amino_set = self.table[amino_id]
        if len(amino_set) == 1:
            self.verbose_msg("Cannot remove codon",codon,"from set, as it is only codon for amino",amino_id+".")
            return None
        amino_set[codon] = 0.0
        amino_set = self.normalise_codon_frequencies(amino_set)
        self.table[amino_id] = amino_set

    def set_min_frequency(self, freq):
        '''Sets all codons whose frequency is below "freq" to 0.0, and re-normalises.
        It is advised not to set freq above a very low value intended only to prune
        exceptionally rare codons, as the "best pick" codon optimisation strategy is
        largely discredited.'''
        # Do this two-step to avoid editing-while-iterating bugs.
        # First, id the codons below freq
        hitlist = []
        for amino in self.table.values():
            for codon in amino:
                if amino[codon] < freq:
                    hitlist.append(codon)
        # then, remove each in turn
        for codon in hitlist:
            self.remove_codon(codon)

    def favour_base(self, base, bonus = 0.15, minimum_f = 0.15):
        '''Allows a table to be re-written to favour codons rich in a specified base.
        bonus is the frequency to boost chosen codons by.
        minimum_f is the minimum frequency below which a codon will not be
        considered viable for boosting.
        This operation was written to be used in Adenine-enrichment of early
        mRNA, which is one simple strategy used to avoid secondary structures.'''
        base = base.upper().replace("U","T")
        for amino in self.table:
            if amino == "*": continue
            self.verbose_msg("Favouring base",base,"for amino",amino)
            # Get list of codons and sort by occurrance of favoured base.
            codons = list(self.table[amino].keys())
            sorted(codons, key=lambda s: s.count(base))
            # Boost the first codon in this ordered list that is above minimum_f.
            if len(codons) == 1: continue
            for codon in codons:
                if base not in codon: continue
                current_frequency = self.table[amino][codon]
                if current_frequency < minimum_f:
                    self.verbose_msg("Skipping",codon,
                        "as frequency is below threshold:",
                        "{0} < {1}".format(current_frequency, minimum_f))
                    continue
                else:
                    self.edit_codon(codon, min(0.99, current_frequency + bonus))
                    self.verbose_msg("Boosted codon",codon,"by {}.".format(bonus))
                    break

    def heat_map(self, sequence, frame=0):
        '''Returns an iterator for codons in sequence[frame:] that yields (codon, frequency).
        This can be used to analyse the codon usage bias in a sequence.'''
        if isinstance(sequence, list):
            frame = 0 # Accept in frame as given.
            sequence = ''.join(sequence)
        sequence = sequtils.bioseq(sequence.upper().replace("U","T"))
        codons = sequence.codons(frame, leading=False, trailing=False)
        results = []
        for codon in codons:
            referred_amino = self.for_table.translate_codon(codon)
            # Get the float value of codon's frequency relative to synonymous codons.
            results.append(self.table[referred_amino][codon])
        return results

    def splice_tables(self, other_table, min_frequency):
        '''Given another object of this type and a minimum frequency, first create
        an averaged table between the two, then for each value which is used in either
        table at a frequency below min_frequency, or which is used for different
        amino acids between the two tables, remove it from the compromise table
        and re-normalise frequencies. The result is a table that should function at
        some level in both represented species.'''
        raise NotImplementedError("Have not yet done this.")
