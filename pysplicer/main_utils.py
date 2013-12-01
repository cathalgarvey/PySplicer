'''Contains the main stuff PySplicer *does*, which is called by scripts or by
client libraries to actually make use of PySplicer.'''

# Needs a "main" function which replicates the behaviour of the pysplicer script
# -> Might just copy the code from there and wrap in "main()" for starters.

from pysplicer import DNAMapper
from pysplicer import translators
from pysplicer import CodonSplicer
from pysplicer import CodonJuggler
from pysplicer import sequtils
from pysplicer import builtin_frequency_tables
from pysplicer import enzyme_lib

# Won't work unless Vienna package is installed to path but importing does no harm:
from pysplicer import RNAfoldWrap

import json
import sys

def table_by_name(tablen):
    'For getting species tables by name at command prompt or when otherwise lazy'
    tablen = tablen.lower()
    try:
        table_to_use = builtin_frequency_tables.__dict__[tablen]
    except KeyError:
        raise KeyError("Could not find species table '"+tablen+\
            "' - Available tables are: "+\
            str(sorted(builtin_frequency_tables.__dict__.keys())))
    return table_to_use

def splice_compile(input_sequence,
                     excludes,                      # IUPAC DNA patterns to (attempt to) exclude from output
                     table_to_use,                  # Codon frequency table to use (default is E.coli empirical/optimal)
                     output_rna=False,              # Output RNA rather than DNA after processing
                     max_early_codon_frequency=0.7, # Frequency cap above which early codons should be excluded (on-ramp hypothesis!)
                     min_codon_frequency=0.025,     # Frequency value (between 0 and 1) below which codons should be excluded
                     verbose=False,                 

                     avoid_ngg_span=15,             # Number of leading codons in which to avoid NGG codons, default 15
                     enrich_adenine=False,          # Crude way to avoid early secondary structures, favour adenine-rich early codons
                     ignore_structure=False,        # Skip searching for and correcting minor secondary structures
                     use_vienna=True,               # If Vienna is installed this should be True. Perform check?
                     max_early_free_energy=6,       # If using ViennaRNA, this is the maximum free energy allowable for the
                                                    #     first several codons (as specified in #avoid-ngg-span option)
                     heatmap = False,               # Output a frequency heatmap after compilation
                     candidates=200,                # Number of splice candidates to generate during each round of splicing
                     ):
    
    def vp(*args,**kwargs):
        if verbose: print(*args, **kwargs)
    
    # Tell the world!
    vp('''splice_compile called with extra args:
    output_rna = {0}
    max_early_codon_frequency = {1}
    min_codon_frequency = {2}
    avoid_ngg_span = {3}
    enrich_adenine = {4}
    ignore_structure = {5}
    use_vienna = {6}
    max_early_free_energy = {7}
    heatmap = {8}
    candidates = {9}
    '''.format(output_rna, max_early_codon_frequency, min_codon_frequency, avoid_ngg_span, enrich_adenine, ignore_structure, use_vienna, max_early_free_energy, heatmap, candidates))
    
    # Leader table for early NGG avoidance.
    leader_table = translators.ReverseTranslator(table_to_use, verbose=verbose)
    # These won't all be removed, usually; TGG is only codon for "W" in some tables,
    # for example. If verbose is on, you'll see messages saying so..
    leader_table.remove_codon("TGG")
    leader_table.remove_codon("AGG")
    leader_table.remove_codon("GGG")
    leader_table.remove_codon("CGG")
    leader_table.set_max_frequency(max_early_codon_frequency)
    if enrich_adenine:
        leader_table.favour_base("A")
    # leader_table.set_min_frequency(min_codon_frequency)
    
    translationtable = translators.ReverseTranslator(table_to_use, verbose=verbose)
    translationtable.set_min_frequency(min_codon_frequency)
    
    # A structure-mapper for identifying and removing unwanted secondary structures.
    # At present, this is super-inefficient, and *only* removes fairly simple hairpins.
    # HairpinMapper init args: (self, min_hairpin=5, endnum=3, min_score=10, max_loop=6, verbose=False):
    if ignore_structure:
        vp("Ignoring structure.")
        Structure_Mapper = lambda seq: {}
    elif use_vienna:
        Structure_Mapper = lambda seq: RNAfoldWrap.map_structure(seq, max_early_free_energy)
    else:
        _Structure_Mapper   = DNAMapper.HairpinMapper(verbose=verbose)
        Structure_Mapper = _Structure_Mapper.map_hairpins
    
    # Now, either use or don't use structure-mapper as part of overall mapper!
    struct_mapper = DNAMapper.DNAMapper(excludes, structure_mapper=Structure_Mapper, verbose=verbose)
    mapper = DNAMapper.DNAMapper(excludes, verbose=verbose)
    
    Struct_Subber = CodonJuggler.Substitutor(table=translationtable, mapper=struct_mapper, verbose=verbose)
    Subber = CodonJuggler.Substitutor(table=translationtable, mapper=mapper, verbose=verbose)
    
    # Iterate/Alternate between splicing large numbers of weighted-randomly generated
    # candidates and trying to remove troublesome codons by substitution. If this doesn't
    # work out first-time around, submit the best effort of prior round(s) to the next
    # set of splice-candidates. Do this maximum 3 times.
    def optimise_chunk(seq, localmapper, localtable, localsubber, local_use_vienna=False):
        prior_efforts = []
        for i in range(0,3):
            vp("Attempt ",i,"/3 - Generating splice candidates (plus any",
                      " prior good results) and optimising..",sep='')
            order_function = RNAfoldWrap.sort_by_fe if local_use_vienna else None
            if ignore_structure: order_function = lambda x:x
            spliced_codon_list = CodonSplicer.splice_aminos_to_codons(seq, localtable,
                                    localmapper, candidates, prior_efforts,
                                    verbose, order_function)
            scrubbed_codons = localsubber.remove_mapped(spliced_codon_list)
            if not localmapper.codon_indices(scrubbed_codons):
                break
            prior_efforts.append(scrubbed_codons)
        return scrubbed_codons
    
    # get the sequence to convert:
    try:
        with input_sequence as InputSequence:
            seq_stream = []
            hitcontent = False
            for line in InputSequence:
                line = line.strip()
                if not line:
                    if not hitcontent:
                        # Don't quit if there is leading whitespace in the file.
                        continue
                    else:
                        # If content has been reached and finished with, break
                        # at first empty line (i.e. don't concatenate a multifasta
                        # file into one long sequence, just use first block.)
                        # TODO: Alternative batch mode? Redundant with Fastac2 I guess
                        break
                else: # i.e. line has contents
                    if line[0] not in ">;#":
                        # Ignore comments and titles.
                        # When meaningful non-comment/title/whitespace content is
                        # hit, mark hitcontent True. This allows starting comment
                        # matter and clarifying whitespace prior to initial block.
                        seq_stream.append(line)
                        hitcontent = True
            seq = ''.join(seq_stream).upper()
    except KeyboardInterrupt:
        print("\nKeyboardInterrupt caught while parsing input. Aborting pysplicer.")
        sys.exit(1)
    
    # First do leader codons:
    lead_seq = optimise_chunk(seq[:avoid_ngg_span], struct_mapper, leader_table, Struct_Subber, use_vienna)
    
    # And latter:
    if len(seq) > avoid_ngg_span:
        latter_seq = optimise_chunk(seq[avoid_ngg_span:], mapper, translationtable, Subber)
    
        # Now join and fix, without reverting to NGG codons: use leader-table minus
        # structure mapping.
        # No changes outside the splice junction should occur, as all options are already
        # exhausted for this codon table, unless structure conflicted with raw sequence,
        # in which case at this point certain substitutions may take place within the
        # leader sequence now that structure is no longer a constraint.
        lead_junction_subber = CodonJuggler.Substitutor(table=leader_table, mapper=mapper, verbose=verbose)
        final = lead_junction_subber.remove_mapped(lead_seq + latter_seq)
    else:
        final = lead_seq
    
    finalseq = ''.join(final)
    if output_rna: finalseq = finalseq.replace("T","U")
    
    vp("Final codon list:",final)
    vp("Translates to:",translationtable.for_table.translate_codon_list_to_str(final))
    if verbose or heatmap:
        print("Heatmap:\t\t",'-'.join([str(int(x*100)).zfill(2) for x in translationtable.heat_map(final)]))
        if use_vienna:
            print("Free energy of avoid-ngg-span region:", RNAfoldWrap.calc_fe(finalseq[:avoid_ngg_span*3]))
    
    return finalseq

def parse_excludes_file(f):
    'Expects *open file object*, parses it for either IUPAC sequences or quoted enzyme names, ignoring "#" comments.'
    excludes = []
    for line in f:
        line = line.strip()
        # Allow empty lines and "#" delimited comments.
        if not line:            continue
        elif line[0] == "#":    continue
        elif line[0] in "\"'":  # Allow named enzymes if quoted with either " or '
            try:
                excludes.append(enzyme_lib.get_enzyme(line.strip("\"'")).target_site)
            except KeyError:
                print("Could not find enzyme {0} from excludes file in enzyme library: Skipping.".format(line), file=sys.stderr)
        else: 
            excludes.append(line.upper().replace("U","T")) 
    return excludes
    
def main(args):
    'Args is an argparse namespace object; this is the function called by pysplicer main script.'

    # Get excludes from invocation and from file, if any.
    excludes = [x.strip().upper().replace("U","T") for x in args.exclude]
    if args.exclude_file:
        excludes.extend(parse_excludes_file(args.exclude_file))
        # All of the below now functionalised; test before deleting.
        #for line in args.exclude_file:
        #    line = line.strip()
        #    # Allow empty lines and "#" delimited comments.
        #    if not line:            continue
        #    elif line[0] == "#":    continue
        #    elif line[0] in "\"'":  # Allow named enzymes if quoted with either " or '
        #        try:
        #            excludes.append(enzyme_lib.get_enzyme(line.strip("\"'")))
        #        except KeyError:
        #            print("Could not find enzyme {0} from excludes file in enzyme library: Skipping.".format(line), file=sys.stderr)
        #    else: 
        #        excludes.append(line.upper().replace("U","T"))
        args.exclude_file.close()
    
    if args.exclude_enzymes:
        # For each name, converts to lowercase, seeks in builtin enzyme list.
        # Prints to stderr if not found and seeks next enzyme.
        excludes.extend(enzyme_lib.get_target_sites(args.exclude_enzymes))

    if args.table:
        with open(args.table) as InFile:
            table_to_use = json.load(InFile)
    elif args.species:
        table_to_use = table_by_name(args.species)
    else:
        raise Exception("Either a species name or a manually created Codon Frequency/Usage table must be provided for PySplicer to begin!")

    seq = splice_compile(args.input_sequence, excludes, table_to_use, # Positional
                            output_rna = args.output_rna,
                            max_early_codon_frequency = args.max_early_codon_frequency,
                            min_codon_frequency = args.min_codon_frequency,
                            verbose = args.verbose,
                            avoid_ngg_span = args.avoid_ngg_span,
                            enrich_adenine = args.enrich_adenine,
                            ignore_structure = args.ignore_structure,
                            use_vienna = args.use_vienna,
                            max_early_free_energy = args.max_early_free_energy,
                            heatmap = args.heatmap,
                            candidates = args.candidates
                            )

    if args.fasta_title: print(">",args.fasta_title)
    if args.wrap:
        seq = "\n".join(sequtils._chunks(seq, args.wrap))
    print(seq)
