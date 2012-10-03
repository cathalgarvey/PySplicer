#!/usr/bin/env python3
# pysplicer.py - A codon-optimising script with support for enzyme site exclusion.
# Uses a weighted-random codon selection system with optional minimum frequency threshold to exclude very rare codons.
# Comes with presets including empirically determined E.coli optimal codon tables from Welch et al, 2009.
# Applies miscellaneous rules from "Designing Genes for Successful Protein Expression" by Welch et al, 2011:
#  - TODO: Minimum frequency threshold is increased for first 25 codons.
#  - TODO: "NGG" codons are penalised in first 25 codons.

import random
import sys
import os
import argparse
import json
import itertools
import re
from math import ceil, floor

# Below: Some files and directories are required for the script to correctly import codon tables etc.,
# so first the expected directories are set to variables, and then it's checked whether or not they exist.
script_dir = sys.path[0] # Gives script directory. Default/imported codon tables will be here.
lib_dir = os.path.join(script_dir, 'lib')
codon_table_dir = os.path.join(lib_dir, 'codontables')
working_dir = os.getcwd() # Gives current working directory.
for needed_dir in [script_dir, lib_dir, codon_table_dir]:
    try:
        assert os.path.exists(needed_dir)
    except AssertionError:
        print("Could not find necessary directory: "+needed_dir+". Aborting.")
        sys.exit(1)

# Below: Initialise ArgumentParser object, which provides the nice CLI interface.
arg_parser = argparse.ArgumentParser(description=('Performs reverse-translation and weighted-random '
                                         'codon optimisation on an amino acid sequence. '
                                         'Default setting is an E.coli table derived from '
                                         'empirical synthetic DNA research by Welch et al, '
                                         '2009.'))
arg_parser.add_argument('infile', type=str,
                       help='Protein sequence, in FASTA file format, to be reverse translated and codon optimised.')
arg_parser.add_argument('-d', '--output-dna', action='store_true',
                       help='Output codon-optimised sequence as DNA, not RNA.')
arg_parser.add_argument('-v', '--verbose', action='store_true',
                       help='Print far more information during processing.')
arg_parser.add_argument('-s', '--species', type=str, default='default_ecoli',
                       help='Species to optimise DNA for. Default is E.coli.')
arg_parser.add_argument('-x', '--exclude-sites',
                       help=('A json-formatted list of DNA or RNA substrings to forbid, '
                             'usually restriction enzyme target sites.'))
arg_parser.add_argument('-m', '--min-codon-frequency', type=float, default=0.0,
                       help=('A minimum codon frequency to accept; frequencies below this are discarded and'
                            ' remaining codon frequencies re-normalised.'))
arg_parser.add_argument('-n', '--splice-candidates', type=int, default=0,
                       help=('Number of random candidate codon sets to generate before attempting to splice '
                             'together a suitably enzyme-free set. Larger numbers of candidates may be useful '
                             'when attempting to exclude many sites, or highly redundant sites.'))
arg_parser.add_argument('--usage', action='store_true',
                       help='Ignore other flags/input and provide detailed usage information.')

args = vars(arg_parser.parse_args()) # vars makes a dict of the results.
# ============================================================================================== #
# Below: The class that does most of the heavy lifting here:

class rev_translate_rna:

    def __init__(self, aminos, codontable, excludesites, min_codon_frequency=0.0, outputdna=False, candidates=0, verbosity=False):
        '''aminos: A string containing just the peptide sequence of a desired protein.
           codontable: Should be a dict mapping amino acid letters to codons with frequencies. Example entry:
              "C": { "TGC": { "frequency": 0.5772005772005773 },
                     "TGT": { "frequency": 0.42279942279942284}  }
           exclude_sites: Should be a list of sites to exclude; list is iterated over and compared to output.
           output_dna: Boolean: output is converted to DNA before export.
           candidates: Number of permutations of raw codon list to create by weighted-random selection before splicing.
        '''
        self.verbose = verbosity
        with open("pysplicer_report", mode='w') as Reportfile:
            Reportfile.write("New job.\r\n")

        def err_then_exit(String):
            print(String)
            sys.exit(1)

        try:
            self.peptidesequence = self.load_aminos(aminos)
        except AssertionError:
            err_then_exit("Error: Provided 'aminos' not a string!")
        self.report("Peptide sequence:\r\n"+self.peptidesequence)

        try:
            self.excludesites, self.largestexclude = self.parse_excludes(excludesites)
            self.excluderegexes, self.largestexclude = self.parse_excludes_to_regex(excludesites)
        except AssertionError:
            err_then_exit("Error: Provided 'excludesites' not a list of strings!")
        self.report("Forbidden sequences:\r\n"+json.dumps(self.excludesites, indent=1))

        try:
            self.codontable = self.import_codon_table(codontable, min_codon_frequency)
        except AssertionError:
            err_then_exit("Error: Provided 'codontable' not a dictionary!")

        self.translationtable = self.build_translation_table(codontable)
        self.buildwindowsize = ceil(self.largestexclude/3) # i.e. round up.
        self.outputdna = outputdna
        self.rnacodonlist = []
        self.ignoredcodons = [] # Used to ignore intransigent sites.
        self.result = ''
        candidatesnum = candidates or 20
        self.generate_candidate_list(candidatesnum)
        self.verbose_msg("Initialised rev_translate_rna object.")

    def verbose_msg(self, msg):
        'For messages that are unimportant unless user wants verbose output.'
        if self.verbose:
            print(msg)

    def report(self, msg):
        'For messages that are important enough to show the user.'
        # If it's important enough to report it's definitely important enough for verbose output.
        assert isinstance(msg, str)
        self.verbose_msg(msg)
        with open("pysplicer_report", mode='a') as Reportfile:
            Reportfile.write("\r\n"+msg)

    def translate(self, codonlist):
        'Translate a list of codons into a string of amino acids. No regard for start/stop: translates all.'
        outstring = ''
        for codon in codonlist:
            encodes = self.translationtable[codon.upper()]
            outstring += encodes
        return outstring


    def reverse_translate(self, amino_string):
        'Returns a list of codons encoding the amino_string, chosen by weighted-random from self.codontable.'
        codon_list = []
        stringindex = 0
        for amino in amino_string:
            codon = self.select_codon(amino, index=stringindex)
            codon_list.append(codon)
            stringindex += 1
        return codon_list

    def generate_codon_list(self):
        return self.reverse_translate(self.peptidesequence)

    def generate_codon_sequence(self):
        'Shortcut method that just sets self.rnacodonlist to a raw reverse-translated form of self.peptidesequence.'
        self.rnacodonlist = self.generate_codon_list()

    def generate_candidate_list(self, nums=20):
        assert nums > 0 # Better an assertionerror than an infinite loop on negative nums..
        self.candidatecodonlists = []
        while nums > 0:
            nums -= 1 # Count down.
            self.candidatecodonlists.append(self.generate_codon_list())

    def pretty_codons_list(self, codons):
        'An improvement to the codon output format; tabulates codons in rows of ten for ease of manual error checking.'
        num_codons = len(codons)
        format_string = '{}\t{} {} {} {} {} {} {} {} {} {}\r\n'
        print_string = ''
        for x in range(0, num_codons, 10): # x counts from 0 in increments of 10
            sublist = codons[x:x+10]
            while len(sublist) < 10: # Prevent format errors on last line.
                sublist.append('')
            print_string += format_string.format(x,*sublist)
        return print_string

    def splice_sequences(self, sequence_list=[]):
        '''Takes a list of candidate sequences (which must encode same amino sequence) and attempts to splice a site-free form.

        This method uses self.issueMap to do the job of mapping each candidate sequence for issues, and compares sequences
        to find which one has the longest stretch of codons without any issues to address. It appends all the codons prior to
        the first discovered issue to an internal candidate sequence, and then starts again from this point, attempting to find
        another candidate with a useful stretch of site-free nucleotides.

        In doing so, two problems may arise:
            A) The method may encounter an issue/site that is universal to all sequences, and it must ignore said site and continue.
            B) The method may create new issues/sites at the splice junctions, but resolving these issues can be left to the
               self.optimise_sequence() method.
        '''
        def assert_equivalence(seq_list):
            'This function checks that seq_list is a non-zero list and that contents are valid and equivalent codon lists.'
            assert isinstance(seq_list, list)
            assert len(seq_list) > 0
            for i in seq_list:
                assert isinstance(i, list)
            First = self.translate(seq_list[0])
            #nonlocal seq_len
            #seq_len = len(seq_list[0]) # Used later for a while loop.
            for i in seq_list:
                assert self.translate(i) == First
            return True

        if not sequence_list:
            sequence_list = self.candidatecodonlists

        try:
            assert_equivalence(sequence_list)
            seq_len = len(sequence_list[0]) # Used later for a while loop.
        except AssertionError:
            self.verbose_msg("Splicer: ERROR: sequence_list must be a list of codon lists, which all encode the same amino sequence.")
            self.verbose_msg("Splicer: ERROR: Can't work under these conditions, returning None.")
            return None

        # Below: We're going to iterate over sequence_list, making a new working list with dictionary values instead of plain codon lists.
        # Each dictionary is simply the issueMap of the codon list, with the codon list embedded into it under key "fullcodonlist".
        working_list = []
        for i in sequence_list:
            i_issuemap = {}
            i_issuemap['issueMap'] = self.map_excluded_sites(i, self.excluderegexes)
            i_issuemap['fullcodonlist'] = i
            i_issuemap['mapIndices'] = list(i_issuemap['issueMap'].keys())
            i_issuemap['mapIndices'].sort() # Probably unnecessary, as dict-keys are already sorted..
            working_list.append(i_issuemap)
            # i_issuemap['issueMap'] is a dict with numerical keys.
            # For each key (which corresponds to a site with issue(s)), the values returned are:
            #   {'spanning_codons', 'encodes', 'contains_sites', 'ensuing_codons', 'preceding_codons'}
            # Splicer will (probably) ignore most of these keys, instead comparing the 'mapIndices' lists simply to check which
            # Sequence has the longest span of uninterrupted codons.

        new_codons = []
        self.verbose_msg("Splicer: About to attempt splicing of the following sequences:")
        for i in working_list:
            i_number = working_list.index(i)
            self.verbose_msg("Splicer: Seq #"+str(i_number)+" issues: "+str(i['mapIndices']))

        def find_best_span(work_list, codon_offset):
            'Working from codon_offset, compare each sequence in work_list to find one with longest-stretch-til-next-issue.'

            self.verbose_msg("Splicer: Working from Offset "+str(codon_offset)+"...")
            # This will be the index of which sequence has best stretch of clear codons.
            best_so_far = {'Seq':-1, 'codon':-1}

            for i in range(0,len(work_list)):
                i_indices = work_list[i]['mapIndices']
                # Below: scan over a slice-copy of the list so we can edit the list as we iterate.
                for x in i_indices[:]:
                    if x < codon_offset: del(i_indices[i_indices.index(x)]) # Delete if X less than codon_offset.
                if len(i_indices) < 1:
                    best_so_far = {'Seq':i, 'codon':seq_len-1}
                    continue
                if i_indices[0] > best_so_far['codon']:
                    best_so_far = {'Seq':i, 'codon':i_indices[0]}
            self.verbose_msg("Splicer: Best so far is Seq #"+str(best_so_far['Seq'])+\
                             ", which spans up until codon "+str(best_so_far['codon'])+".")
            if (best_so_far['codon'] == codon_offset) or (best_so_far['Seq'] < 0):
                self.verbose_msg("Splicer: had to skip codon #"+str(codon_offset)+\
                                 ": this is either a recalcitrant enzyme site (bad) or the end of "+\
                                 "the sequence (good).")
                return work_list[best_so_far['Seq']]['fullcodonlist'][codon_offset:codon_offset+1] #Return next codon of "best" anyway.

            best_full_length = work_list[best_so_far['Seq']]['fullcodonlist']
            next_codons = best_full_length[codon_offset:best_so_far['codon']]
            return next_codons

        while len(new_codons) < seq_len:
            new_codons += find_best_span(working_list, len(new_codons))

        self.verbose_msg("Splicer: Sanity checking that new_codons matches amino sequence of input sequence set..")
        assert self.translate(new_codons) == self.translate(sequence_list[0])
        self.rnacodonlist = new_codons

    def create_codon_permutations(self, amino_sequence, sanity_ceiling = 300):
        '''Return a list of codon-lists that encode the desired amino_sequence. Warning: for long
        sequences, this might be a really dumb thing to do, so this implements a sanity check,
        comparing the predicted number of permutations to the "sanity_ceiling" and aborting if
        greater.'''
        # Will be populated with sub-lists of codon choices at every codon.
        codon_choices = []
        for amino in amino_sequence:
            possible_codons = list(self.codontable[amino.upper()].keys())
            codon_choices.append(possible_codons)

        # Sanity check warning:
        predicted_permutations = 1
        for codon_set in codon_choices:
            # Calculate # of possible permutations by multiplying choices
            predicted_permutations *= len(codon_set)
        # Warn user if we're shortcutting through massive permutation-space.
        # Later, this could be used to set a flag that changes behaviour of generator, below.
        if predicted_permutations > sanity_ceiling:
            self.verbose_msg("Predicted permutations: "+str(predicted_permutations)+\
                             ".. only "+str(sanity_ceiling)+" will be returned.")
        self.verbose_msg("Reverse Translator: Generating "+str(predicted_permutations)+\
                         " permutations for amino sequence "+amino_sequence+"..")
        permutation_maker = itertools.product(*codon_choices) # Yields generator.
        perms = []

        # Extract vals from generator, stopping short at the sanity-ceiling value.
        for permutation in permutation_maker:
            perms.append(permutation)
            # Break the loop once we have sanity_ceiling permutations.
            # TODO: Make this skip through the generator to fetch more diverse permutations
            #     if it knows it'll be hitting the ceiling; otherwise options may be pretty limited
            #     in large permutation sets with >1000s of possibilities; first 300 will suck.
            if len(perms) >= sanity_ceiling:
                self.verbose_msg("Generated more than 500 permutations for "+amino_sequence+\
                                 "..generating no further.")
                break
        return perms

    def choose_best_codon_permutation(self, permutations, codon_position=0):
        '''This is called to choose a good codon out of a given set. It returns the chosen codon.

        In PySplicer, this is called after removing all the unsuitable codons, freeing this Method to
        pick the "best" of the remainder where a preference makes any sense.

        This should pick with some regard to the desired frequency table, or perhaps in such a way as to
        "mend" deviations from the desired frequencies.

        If a codon locational index is given this can be used to bias the choice; for example, NGG codons
        are apparently a bad thing for expression when found within the first ~25 codons in E.coli.'''
        import random
        # If a codon position is provided, check against further rules relating to position.
        if codon_position:
            # Below rule: some empirical work suggests that first 25 codons in E.coli are sensitive in particular
            # to rare codons or to NGG codons, and this can cause ribosomal dropoff. So, check whether the codon
            # choices include non-NGG codons and pick from those if available.
            if codon_position <= 25:
                self.verbose_msg("Choosing best codon for one of first 25 codons; trying to avoid NGG codons.")
                print("Choosing best codon for one of first 25 codons; trying to avoid NGG codons.")
                non_NGG_codons = []
                for codon in permutations:
                    if codon not in ["AGG", "UGG", "GGG", "CGG"]:
                        non_NGG_codons.append(codon)
                if non_NGG_codons:
                    choice = random.choice(non_NGG_codons)
                    # If one or more of the codon choices are non-NGG, choice from them:
                    self.verbose_msg("Non-NGG codons exist for this locus: selecting "+choice)
                    print("Non-NGG codons exist for this locus: selecting "+choice)
                    return choice
                else:
                    self.verbose_msg("No Non-NGG codons existed for this locus after other filters were applied, sorry.")
                    print("No Non-NGG codons existed for this locus after other filters were applied, sorry.")
            # If rules cannot apply without eliminating all choices, or if there were no rules to apply, just proceed with normal choice.
            return random.choice(permutations)
        # If no codon position is provided, then just proceed with normal choice.
        else:
            return random.choice(permutations)

    def optimise_sequence(self):
        def shuffle_codons(codon_list):
            'Translates and then reverse-translates the target region.'
            aminos = self.translate(codon_list)
            new_codons = self.reverse_translate(aminos)
            return new_codons

        def attempt_permutations(issue_dict, index):
            # issue_context_output is used for pretty codon display in verbose mode.
            issue_context_output = ' '.join(issue_dict['preceding_codons']) + '{0}' + ' '.join(issue_dict['ensuing_codons'])
            # issue_context, with no spaces, is used for formatting permutations into the local context, checking.
            issue_context = ''.join(issue_dict['preceding_codons']) + '{0}' + ''.join(issue_dict['ensuing_codons'])

            # Brute-force create all available permutations for the amino sequence encoded by the area, as a list.
            # This will stop short of a "sanity_ceiling", in case you give it too much possibility-stuff
            permutations = self.create_codon_permutations(issue_dict['encodes'])

            self.verbose_msg("Reverse Translator: Trimming permutations to remove those containing excluded sites..")
            self.verbose_msg("Reverse Translator: permutations are: "+json.dumps(permutations))
            # Copy list: can parse over and delete from original as desired without interrupting iteration.
            for perm in permutations[:]:
                # Create a string representation of the local area with this permutation in it.
                perm_sequence = issue_context.format(''.join(perm)).upper()
                # Now scan through all the excluded sites and see if they're in the area when this permutation is added.
                for exclude in self.excludesites:
                    if exclude in perm_sequence:
                        # Remove this permutation from the list of possibilities.
                        # NB: This is removing perm from the original list, not the copy we've made above for this loop.
                        self.verbose_msg("Permutation "+''.join(perm)+" was found to harbour forbidden sequences, and will be "+\
                                         "removed from the set of possibilities for this locus.")
                        try: permutations.remove(permutations.index(perm))
                        except ValueError: self.verbose_msg("permutation "+''.join(perm)+" not found in set of permutations from "+\
                                                            "which it was originally copied. How odd.")

            self.verbose_msg("Reverse Translator: permutations are: "+json.dumps(permutations))
            # If no viable options remain, then flag this locus as an intractable site to ignore in future.
            if len(permutations) < 1:
                print("Reverse Translator: Could not find an enzyme-free permutation for issue at index "+str(index)+", "+\
                      "corresponding to nucleotide "+str(index*3)+"! You must address this manually. Sorry!")
                self.ignoredcodons.append(index)
            # If there are suitable permutations:
            else:
                # Select a suitable codon (currently random) and replace problem area with it.
                suitable_codons = self.choose_best_codon_permutation(permutations, index)
                num_codons = len(suitable_codons)
                self.verbose_msg("Reverse Translator: Choosing: "+\
                                 ''.join(suitable_codons) + " for codons "+\
                                 str(index)+" to "+str(index+num_codons-1)+".")
                self.verbose_msg("New Context is:\r\n"+issue_context_output.format('|'+' '.join(suitable_codons)+'|').upper())
                self.rnacodonlist[index:index+num_codons] = suitable_codons

        self.verbose_msg("Reverse Translator: Generating first-pass codon list for amino acid sequence.")
        if not self.rnacodonlist: # If there's nothing to work with, generate a candidate:
            self.rnacodonlist = self.reverse_translate(self.peptidesequence)
            self.verbose_msg("Reverse Translator: First pass codon list generation complete.")

        # Now, try three cycles of issue-finding and issue-resolving.
        i = 3
        while i > 0:
            # Count down:
            i -= 1
            # Get issueMap for this cycle. If it's empty, we're done, break the loop and continue.
            issueMap = self.map_excluded_sites(self.rnacodonlist, self.excluderegexes)
            if not issueMap:
                break

            # issueMap is a dict with numerical keys. The values returned are:
            #   {'spanning_codons', 'encodes', 'contains_sites', 'ensuing_codons', 'preceding_codons'}
            
            # Create an ordered list of indices:
            index_list = list(issueMap.keys())
            self.verbose_msg("Reverse Translator: index_list is: "+str(index_list))
            # Iterate over list of troublesome codon indices. Recover the details of the problem from issueMap,
            # then pass these details to the attempt_permutations function, above.
            for index in index_list:
                self.verbose_msg("Reverse Translator: Addressing issue at codon "+str(index)+" (nucleotide "+str(index*3)+"):")
                ProblemDict = issueMap[index]
                attempt_permutations(ProblemDict, index)

        # Finished:
        self.verbose_msg("Reverse Translator: Finished addressing issues, concatenating codons to result!")
        self.result = ''.join(self.rnacodonlist)
        if self.outputdna:
            self.result = self.result.replace("U", "T")

    def map_excluded_sites(self, codon_list, regex_sites):
        '''Returns a dictionary containing integer keys, where each key corresponds to a list index for
        codon_list where an excluded site was found. The dictionary contains the number of codons the
        excluded site spans, the amino acids encoded by those codons, and any other sub-sites found at
        the same index.'''
        def uniquify(seq): # Dave Kirby
            '''This is a fast, order-preserving function for removing duplicates from a sequence/list.
            Found here: http://www.peterbe.com/plog/uniqifiers-benchmark'''
            seen = set()
            return [x for x in seq if x not in seen and not seen.add(x)]

        def allindices(string, sub, offset=0):
            'Returns a list of ALL indices of a substring, not just first.'
            listindex = []
            i = string.find(sub, offset)
            while i >= 0: # Because find returns -1 when it finds nothing.
                listindex.append(i)
                i = string.find(sub, i + 1)
            return listindex

        def find_problem_codons(sequence, excludedsite):
            'Maps string indices to codon-list indices.'
            stringindices = allindices(sequence, excludedsite)
            codonlistindices = []
            for stringindex in stringindices:
                codonindex = floor(stringindex/3) # floor = imported from math.
                if codonindex not in self.ignoredcodons:
                    if codonindex not in codonlistindices:
                        codonlistindices.append(codonindex)
                else:
                    self.verbose_msg("Mapper: Ignoring issue at index "+str(codonindex)+" because it can't be resolved.")
            return codonlistindices

        self.verbose_msg("Mapping issues in RNA codon map using regex-equivalents of IUPAC excluded sequences.")
        self.verbose_msg(self.pretty_codons_list(codon_list))

        found_issues = {}
        sequence_string = ''.join(codon_list)
        regex_identified_sites = []
        for regex in regex_sites:
            # This appends a list of all identified permutations of the regex to master list.
            identified_sites = regex.findall(sequence_string)
            identified_sites = uniquify(identified_sites)
            regex_identified_sites.extend(identified_sites)

        # Now that we have a list of identified permutations of the regexes, we can be more specific.
        for exclude in regex_identified_sites:
            exclude_codon_length = ceil(len(exclude) / 3) # Round up.
            problem_codons = find_problem_codons(sequence_string, exclude)
            for problem_codon in problem_codons: # problem_codon is a list index for the codon in question.
                self.verbose_msg("Mapper: Found site: "+exclude+" @ codon "+str(problem_codon)+".")

                # Below: sites don't always align perfectly with boundary of "problemcodon", so must determine if overhang exists.
                #  So, join codons, search for problem site, and if it's not in joined codons, probably need extra codon.
                if exclude in ''.join(codon_list[problem_codon : problem_codon + exclude_codon_length]):
                    spanned_codons = codon_list[problem_codon : problem_codon + exclude_codon_length]
                    ensuing = problem_codon + exclude_codon_length
                else:
                    spanned_codons = codon_list[problem_codon : problem_codon + exclude_codon_length + 1]
                    ensuing = problem_codon + exclude_codon_length + 1
                    assert exclude in ''.join(spanned_codons) # Throw a fit if this fails.

                preceding_codons = codon_list[problem_codon-self.largestexclude:problem_codon]
                ensuing_codons = codon_list[ensuing:ensuing+self.largestexclude]

                problem_string = '\tContext: {0}|{1}|{2}'.format(' '.join(preceding_codons), '-'.join(spanned_codons), ' '.join(ensuing_codons))
                self.verbose_msg(problem_string)

                if problem_codon not in found_issues.keys():

                    found_issues[problem_codon] = { 'spanning_codons': exclude_codon_length,
                                                    'encodes': self.translate(spanned_codons),
                                                    'preceding_codons': preceding_codons,
                                                    'ensuing_codons': ensuing_codons,
                                                    'contains_sites': [exclude] }

                elif found_issues[problem_codon]['spanning_codons'] >= exclude_codon_length:
                    found_issues[problem_codon]['contains_sites'].append(exclude)
                elif found_issues[problem_codon]['spanning_codons'] < exclude_codon_length:
                    found_issues[problem_codon]['spanning_codons'] = exclude_codon_length
                    found_issues[problem_codon]['contains_sites'].append(exclude)
                    found_issues[problem_codon]['encodes'] = self.translate(spanned_codons)
                else:
                    print("Error in map_excluded_sites.")
        self.verbose_msg("Mapper: Finished! Returning to reverse translator.")
        self.verbose_msg("Mapper: Issues are: "+json.dumps(found_issues,indent=1))
        return found_issues

    def select_codon(self, amino, index=None, floor=0.0):
        '''Returns a weighted-random-selected codon with frequency greater than "floor", applying rules if index is not None.'''
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
            raise Exception('Distribution not normalized: {probs}'.format(probs=str(prob_dict)))
        # Above function is re-usable. Below wraps it for specific codon-selection purposes.
        assert isinstance(amino, str)
        candidate_codons = self.codontable[amino]
        if index is not None:
            trimmed_candidate_codons = {}
            if index < 25:
            # Try to exclude NGG codons from first 25 codons; can cause ribosomal dropoff in E.coli at least.
                for codon in candidate_codons.keys():
                    if codon not in ["AGG", "UGG", "CGG", "GGG"]:
                        trimmed_candidate_codons[codon] = candidate_codons[codon]
                if trimmed_candidate_codons:
                    if len(trimmed_candidate_codons) != len(candidate_codons):
                        self.verbose_msg("Trimmed away some NGG codons in first 25 codons.")
                    trimmed_candidate_codons = self.normalise_codon_frequencies(trimmed_candidate_codons)
                    chosen_codon = random_category(trimmed_candidate_codons)
                    return chosen_codon
                else:
                    self.report("Cannot choose non-NGG codon at position "+str(index)+".")
        chosen_codon = random_category(candidate_codons)
        return chosen_codon

    def iupac_reverse_complement(self, sequence):
        'Accepts a DNA string in IUPAC format and returns the reverse complement.'
        assert isinstance(sequence,str)
        inverse_iupac = { # Translates to an IUPAC sequence of potential complements.
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
        complement_sequence_list = []
        for base in sequence:
            try:
                complement_sequence_list.append(inverse_iupac[base])
            except KeyError:
                # This shouldn't happen, as sequences have already been passed through self.nucleic_iupac_purify by now..
                ErrMsg = "Error: non-IUPAC character detected in sequence"+sequence+"..ignoring."
                print(ErrMsg, file=sys.stderr)
                self.verbose_msg(ErrMsg)
                continue
        complement_sequence = ''.join(complement_sequence_list)
        reverse_complement_sequence = complement_sequence[::-1]
        return reverse_complement_sequence

    def nucleic_iupac_purify(self, iupac_sequence, strip_flanking_wildcards=False):
        '''Prunes only legal IUPAC DNA/RNA characters from a given string and returns those.'''
        assert isinstance(iupac_sequence, str)
        iupac_sequence = iupac_sequence.upper()
        iupac_characters = ['A', 'T', 'C', 'G', 'U',	# Canonical bases
                            'B', 'V', 'D', 'H',		# B=Not A, V=Not T, D=Not C, H=Not G
                            'S', 'W',			# Strong and Weak: GC vs. AT
                            'K', 'M',			# "Keto" and "aMino": GT vs. AC
                            'R', 'Y',			# "puRine" and "pYrimidine": AG vs CT
                            'N', '.', '-']		# Wildcards: any N.
        output_char_list = []
        for char in iupac_sequence:
            if char in iupac_characters:
                output_char_list.append(char)
        output_sequence = ''.join(output_char_list)
        if strip_flanking_wildcards:
            # Having "NNNN" before or after a sequence is not only useless,
            # it can cause downstream hiccups where the program may generate huge sets
            # of redundant permutations, or unnecessarily large regular expressions.
            for char in ["N", "-", "."]:
                output_sequence = output_sequence.replace(char, '')
        return output_sequence

    def parse_excludes(self, exclude_list):
        '''Generates a list of permutations for all excluded/forbidden sequences in exclude_list.

        Returns this list and the size of the largest encountered sequence.'''
        self.verbose_msg("Parsing through excluded enzyme list.")
        assert isinstance(exclude_list, list)
        output_list = []
        largest = 0
        for entry in exclude_list:
            assert isinstance(entry, str)
            entry = self.nucleic_iupac_purify(entry, strip_flanking_wildcards=True)
            entry_complement = self.iupac_reverse_complement(entry)
            entry_len = len(entry)
            output_list.extend(self.create_rna_permutations(entry))
            output_list.extend(self.create_rna_permutations(entry_complement))
            if entry_len > largest:
                largest = entry_len
        self.verbose_msg("excluded sites are: \r\n"+str(output_list))
        return output_list, largest

    def dna_iupac_to_regex(self, sequence, return_pattern=False):
        '''Uses a dictionary to assemble a (poorly formatted) python regex function from an IUPAC DNA string.

        For example, given the DNA sequence "AWWSSTGGC", this method returns a regular expression object with
        this pattern: "A[AT][AT][GC][GC]TGGC"'''
        # For those not familiar with regex, square brackets mean "any one of the contents".
        # Therefore, [AGC] means either A, G or C.
        iupac_re = { 'A': 'A',    'B': '[CGT]',
                     'C': 'C',    'D': '[AGT]',
                     'G': 'G',    'H': '[ACT]',
                     'K': '[GT]',   'M': '[AC]',
                     'N': '[ACGT]', 'R': '[AG]',
                     'S': '[GC]',   'T': 'T',
                     'V': '[ACG]',  'W': '[AT]',
                     'Y': '[CT]',   '.': '[ACGT]',
                     '-': '[ACGT]' }
        # First off, get rid of the crap
        sequence = self.nucleic_iupac_purify(sequence, strip_flanking_wildcards=True)
        # Then, create an empty list, and append it with regex substitutes for each IUPAC letter
        sequence_re_list = []
        for base in sequence:
            sequence_re_list.append(iupac_re[base])
        # Flatten the list into a string
        sequence_re_str = ''.join(sequence_re_list)
        # Compile a regex function from the string
        # Again for those not familiar with regex in python, this is an object with lots of methods
        # for parsing arbitrary strings of data and finding anything matching its "pattern".
        sequence_re = re.compile(sequence_re_str)
        if return_pattern:
            return sequence_re, sequence_re_str # Former is the regex, latter is just for verbose readout.
        return sequence_re

    def parse_excludes_to_regex(self, exclude_list):
        '''Translate IUPAC DNA exclusions into regex functions, return a list of compiled regexes for forbidden sites.'''
        assert isinstance(exclude_list, list)
        # output_list is the regular expressions to return, patterns_list is for verbose output.
        output_list = []
        patterns_list = []
        # largest is to tally the largest encountered sequence in exclude_list, which is returned.
        largest = 0
        for entry in exclude_list:
            assert isinstance(entry, str)
            entry = self.nucleic_iupac_purify(entry, strip_flanking_wildcards=True)
            reverse_complement = self.iupac_reverse_complement(entry)
            entry_len = len(entry)
            entry_re, entry_pattern = self.dna_iupac_to_regex(entry, return_pattern=True)
            rev_comp_re, rev_comp_pattern = self.dna_iupac_to_regex(reverse_complement, return_pattern=True)
            output_list.extend([entry_re, rev_comp_re])
            patterns_list.extend([entry_pattern, rev_comp_pattern])
            if entry_len > largest:
                largest = entry_len
        self.verbose_msg("excluded sites are: \r\n"+str(patterns_list))
        return output_list, largest

    def build_translation_table(self, codontable):
        'Inverts a codon table, omitting frequencies, to infer the translation table of the target species.'
        trans_table = {}
        self.verbose_msg("Building forward translation table from codon table.")
        for amino in codontable.keys():
            for codon in codontable[amino]:
                trans_table[codon] = amino
        return trans_table

    def load_aminos(self, aminos):
        assert isinstance(aminos, str)
        # NOTE: This function merely strips illegal (non-IUPAC) characters.
        #  It expects JUST the amino sequence, not markup e.g. FASTA titles! Sanitise input!
        permittedaminos = 'ACDEFGHIKLMNPQRSTVWY' # TODO: implement B, X and Z?
        aminos = aminos.upper()
        for char in aminos:
            if char not in permittedaminos:
                print("Found illegal character (not an IUPAC amino acid): "+char+"..deleting all instances.")
                aminos = aminos.replace(char, "")
        return aminos

    def normalise_codon_frequencies(self, codon_table):
        'Given a dict of format {"XXX":float, "YYY":float}, returns with floats balanced to sum to 1.' 
        assert isinstance(codon_table, dict)
        total_input_frequencies = 0.0
        for codon in codon_table.keys():
            # Sum the floats to get the current total.
            assert isinstance(codon_table[codon], float)
            total_input_frequencies += codon_table[codon]
        total_output_frequencies = 0.0
        for codon in codon_table.keys():
            # Re-balance floats to give a reasonable frequency spread summing to 1.
            codon_table[codon] = codon_table[codon] / total_input_frequencies
            try: assert codon_table[codon] < 1.0000001
            except AssertionError: print("Codons not normalised:", codon, codon_table[codon], file=sys.stderr)
            total_output_frequencies += codon_table[codon]
        return codon_table            

    def import_codon_table(self,codon_table, floor = 0.0):
        '''Imports codon tables that may have excess metadata into simple codon:frequency dicts.
           Also removes codons whose frequency is below floor, and re-normalises remaining frequencies.
           i.e. convert this:
              "C": { "TGC": { "frequency": 0.5772005772005773, 'info':'More highly charged under starvation' },
                     "TGT": { "frequency": 0.42279942279942284, 'synthetase':{'seq':'mabcd','protID=foo1234'}
                    }
           ...into this:
              "C":{"TGC":0.5772005772005773,
                   "TGT":0.42279942279942284}'''
        self.verbose_msg("Importing specified codon table, discarding frequencies below "+str(floor))
        assert isinstance(codon_table, dict)
        frequency_dict = {}
        for amino in codon_table.keys():
            frequency_dict[amino] = {}
            for codon in codon_table[amino].keys():
                # This For-loop gets the frequency for each codon, and copies it if it's greater than floor.
                try:
                    Freq = codon_table[amino][codon]['frequency']
                    if Freq > floor:
                        frequency_dict[amino][codon] = Freq
                except:
                    print("Error while importing codon table; entry {0} in {1} \
                           has no frequency value. Skipping {0}.".format(codon,amino))
            frequency_dict[amino] = self.normalise_codon_frequencies(frequency_dict[amino])
        return frequency_dict

    def create_dna_permutations(self, site):
        'Returns a list of possible permutations for a given IUPAC DNA string.'
        self.verbose_msg("Creating a set of DNA permutations for "+site)

        IUPACDNA = { 'A': 'A',    'B': 'CGT',
                     'C': 'C',    'D': 'AGT',
                     'G': 'G',    'H': 'ACT',
                     'K': 'GT',   'M': 'AC',
                     'N': 'ACGT', 'R': 'AG',
                     'S': 'GC',   'T': 'T',
                     'V': 'ACG',  'W': 'AT',
                     'Y': 'CT',   '.': 'ACGT',
                     '-': 'ACGT' }
        snp_list = []
        perm_list = []
        for char in site:
            snp_list.append(IUPACDNA[char])
        perm_generator = itertools.product(*snp_list) # "*" character dumps list contents as args.
        for item in perm_generator:
            perm_list.append(''.join(item))
        return perm_list

    def create_rna_permutations(self, site):
        'Returns a list of possible RNA permutations for a given IUPAC DNA string.'
        DNAsites = self.create_dna_permutations(site)
        self.verbose_msg("Converting DNA permutations to RNA.")
        RNAsites = []
        for DNAsite in DNAsites:
            RNAsites.append(DNAsite.replace("T", "U"))
        return RNAsites

    def convert_rna_to_dna(self,RNAinput):
        assert isinstance(RNAinput, str)
        DNAoutput = RNAinput.upper().replace("T","U")
        return DNAoutput

    def give_output(self):
        'Runs quick sanity checks and returns result.'
        newaminosequence = self.translate(self.rnacodonlist)
        #print(newaminosequence)
        if not self.peptidesequence == newaminosequence:
            print(("ERROR: Input peptide sequence and translated peptide sequence of output RNA/DNA do not match!\n")
                  ("Please report this to the developer on Github: https://github.com/cathalgarvey\n")
                  ("--Input--\r\n{0}\n".format(self.peptidesequence))
                  ("--Output--\r\n{0}".format(newaminosequence)))
            sys.exit(1)
        if self.outputdna:
            self.result += 'TAA'
        else:
            self.result += 'UAA'
        return self.result

# ======================================================================== #
# Below: Some functions that might be better off as methods for the sake of portability:
from fasta_utils import read_fasta_file

def import_codon_table(codon_table_file):
    '''Expects a file containing a codontable in json format.

    codontable: Should be a dict mapping amino acid letters to codons with frequencies. Example entry:
    "C": { "TGC": { "frequency": 0.5772005772005773 },
           "TGT": { "frequency": 0.42279942279942284}  }
    The file can actually contain other metadata for each codon, but should not contain keys other
    than codons for each amino acid; keys are assumed to correspond to codons only. Likewise, keys
    in the "master" dictionary are assumed to be amino acids and will be imported as such, although
    there is less risk of errors here as keys will only become relevant if called by a function/method.'''
    import json
    codon_table_fileDir = os.path.join(codon_table_dir, codon_table_file + '.json')
    with open(codon_table_fileDir) as TableFile:
        codon_table = json.loads(TableFile.read())
    assert isinstance(codon_table, dict)
    ReturnTable = {}
    for amino in codon_table.keys(): # Assumes all keys are amino acids.
        ReturnTable[amino] = {}
        for codon in codon_table[amino].keys(): # Assumes all keys are codons.
            # What does this achieve, exactly? It just prunes codontable and verifies that each
            #  codon has a frequency value..useful?
            try: ReturnTable[amino][codon] = {'frequency':codon_table[amino][codon]['frequency']}
            except: print("codon "+codon+" of amino "+amino+" has no 'frequency' key! Skipped.")
    return ReturnTable

def import_exclusions(exclude_file):
    '''Expects a JSON formatted file containing a simple list of IUPAC DNA sequences to avoid.

    sites are expanded into permutations and stored in the object for comparison to draft sequences.
    Providing too many sites with ambiguous characters from the IUPAC set may result in huge exclusion
    lists, potential program slowdown or crashing. For the love of all that's holy, strip unnecessary "N"
    characters from TypeIIs or Homing Endonucleases, where sites may often be written with lots of Ns outside
    the actual recognition area; this program will assume these are relevant and will generate 4^(number of Ns)
    permutations.'''

    if not isinstance(exclude_file, str):
        return []
    with open(exclude_file) as excludes_file:
        excludesList = json.loads(excludes_file.read())
    assert isinstance(excludesList, list)
    IUPACDNA = 'ABCDGHKMNRSTVWY.-'
    for exclude in excludesList:
        for char in exclude:
            assert char in IUPACDNA
    return excludesList

# ======================================================================== #
# Below: Parse all those args from above.

# args contains: 'infile', 'output_dna', 'species', 'exclude_sites', 'min_codon_frequency', 'usage', 'splice_candidates'
# Usage: RevTranslatedRNA(aminos, codontable, excludesites, output_dna=False):
fasta_dict = read_fasta_file(args['infile'])
amino_sequence = fasta_dict['sequence']

if not amino_sequence:
    print("No valid amino Acid sequence found in target file. Quitting.")
    sys.exit(1)
try:
    species_codon_table = import_codon_table(args['species'])
except AssertionError:
    print("codon table file specified is not a dict. Quitting.")
    sys.exit(1)
except KeyError:
    print("Encountered a missing dictionary key (amino? codon?) while parsing and pruning codon Table. Quitting.")
    sys.exit(1)

try:
    site_exclusion_list = import_exclusions(args['exclude_sites'])
except AssertionError:
    print(("Either the specified file doesn't encode a list,"
           "or the list entries aren't all valid IUPAC DNA strings.\r\n"
           "NB: permitted IUPAC DNA codes are 'ABCDGHKMNRSTVWY.-'"
           "\r\nProgram will quit."))
    sys.exit(1)

# Below: Argparser does funny things with some optional arguments if unset. For e.g., "store true" arguments, if unset,
# are stored as "None", not as "False". To ward against obscure bugs due to such behaviour, the below just sets args
# in stone before applying to the object.
def make_none_false(argparse_output):
    if argparse_output == None:
        return False
    else:
        return argparse_output

dna_arg = make_none_false(args['output_dna'])
verbose_arg = make_none_false(args['verbose'])
if args['min_codon_frequency']: freq_arg = args['min_codon_frequency']
else: freq_arg = 0.0
if args['splice_candidates']: splices_arg = args['splice_candidates']
else: splices_arg = 0

# Below: Make the magic happen. Instantiate reverse_translator and perform usual operations to get desired output.
reverse_translator = rev_translate_rna(amino_sequence, species_codon_table, site_exclusion_list, freq_arg, dna_arg, splices_arg, verbose_arg)
reverse_translator.splice_sequences()
reverse_translator.optimise_sequence()

# Deliver results of algorithmic goodness:
reverse_translator.verbose_msg("codons with remaining issues: "+json.dumps(reverse_translator.ignoredcodons))
reverse_translator.verbose_msg("\r\nOutput:\r\n")
print(reverse_translator.give_output())
