#!/usr/bin/env python3
'''
A simple enzyme profile object. Should be supplied with a master table in JSON format (typically a large file kept with the script), and can parse through the table and supply matching search keys as desired.
Search keys can be given individually or as a list. If as a list, matching keys are output first, and errors and suggestions second.
If a search key is incomplete, a list of potential matches is returned as a prompt. If the search key ends in *, all matching items are returned.
'''

import json
import sys

class EnzymeDirectory:
    def __init__(self, enzyme_dict, verbose=True, logfile=None):
        self.enzymes = enzyme_dict
        self.verbose = verbose
        self.logfile = logfile

    def uniquify(self, seq): # Dave Kirby
        '''This is a fast, order-preserving function for removing duplicates 
        from a sequence/list.
        Found here: http://www.peterbe.com/plog/uniqifiers-benchmark'''
        seen = set()
        return [x for x in seq if x not in seen and not seen.add(x)]

    def report(self,string):
        if not self.verbose: return None
        if not isinstance(string, str):
            string = str(string)
        print(string)
        if self.logfile:
            try:
                with open(self.logfile, mode='a') as LogOut:
                    LogOut.write("\r\n"+string)
            except IOError:
                with open(self.logfile, mode='w') as LogOut:
                    LogOut.write(string)

    def sites_from_file(self, enzyme_list_file, iupac_only=False, strip_wildcards=False,
                         partial_matches=False):
        '''Reads from a flat text file with a single enzyme on each line,
        returns list of sites. Ignores lines beginning with "#" as comments.''' 
        enzymes_to_search = []
        with open(enzyme_list_file) as EnzymeFile:
            for line in EnzymeFile:
                line = line.strip()
                if not line:
                    continue
                elif line[0] == "#":
                    continue
                else:
                    enzymes_to_search.append(line.lower())
        return self.enzyme_sites(enzymes_to_search, iupac_only, strip_wildcards)

    def enzyme_sites(self, enzyme_list, iupac_only=True, strip_wildcards=False,
                        partial_matches=False):
        '''Fetches all enzymes in the target list and returns list of target sites.
        If iupac_only is True (default), then non-iupac characters are stripped;
        this will remove, for example, cleavage-site indicators.'''
        # Keys: 'target_site', 'name', 'suppliers', 'source', 'references',
        #       'prototype', 'organism', 'methylation_site'
        enzymes = self.parse_list(enzyme_list, partial=partial_matches)
        sites = []
        for enzyme in enzymes:
            target_site = enzymes[enzyme]['target_site']
            target_site = target_site.strip().upper()
            if not target_site or target_site == "?":
                self.report("Restriction site for enzyme "+enzyme+" is unknown.")
                continue
            if iupac_only:
                for x in target_site:
                    if x not in 'ABCDGHKMNRSTUVWY.-':
                        target_site = target_site.replace(x, "")
            if strip_wildcards:
                target_site = target_site.strip("N-.")
            if not target_site: continue
            sites.append(target_site)
        # Remove duplicates
        sites = self.uniquify(sites)
        self.report("Found sites: "+str(sites))
        return sites

    def print_specified_keys(self,enzyme_list,wantedkeys=['name','organism','target_site'],
                              partial=False, indentnum=1):
        'Wraps return_specified_keys and prints the output with indentation.'
        PrunedDicts = self.return_specified_keys(enzyme_list, partial)
        print(json.dumps(PrunedDicts, indent=indentnum))

    def return_specified_keys(self,enzyme_list, wantedkeys=['name','organism','target_site'],
                              partial=False):
        '''Acts as a filter for parse_list, returning pruned dictionaries.
        Keys may be any number of:
             'target_site': Target site of the enzyme.
             'name': Name of the enzyme.
             'suppliers': Potential suppliers of the enzyme.
             'source': Source strain culture collection (eg ATCC) identifier.
             'references': References relating to discovery, characterisation.
             'prototype': Prototypical enzyme name, if any.
             'organism': Source organism/strain.
             'methylation_site': Recognition site for cognate methylase.'''
        FoundEnzymes = self.parse_list(enzyme_list, partial)
        for Enzyme in FoundEnzymes.keys():
            # Get dict
            EnzymeDict = FoundEnzymes[Enzyme]
            # Make edited dict
            for key in EnzymeDict.keys():
                if key not in wantedkeys:
                    del(EnzymeDict[key])
            # Replace dict
            FoundEnzymes[Enzyme] = EnzymeDict
        return FoundEnzymes

    def parse_list(self, enzyme_list, partial=False):
        'Parses a list of enzyme names or search keys, and returns matching enzymes.'
        output = {}
        for key in enzyme_list:
            try:
                enzymedicts = self.search(key, partial)
                for enzymed in enzymedicts:
                    enzymename = enzymed['name'].lower()
                    output[enzymename] = enzymed
            except KeyError as e:
                self.report(e)
            except ValueError as e:
                self.report(e)
        return output

    def search(self, searchkey, partial=False):
        'Wraps _search; if searchkey ends in *, wildcard mode is enabled.'
        if searchkey[-1] == "*":
            return self._search(searchkey.strip().strip("*"), True)
        else:
            return self._search(searchkey.strip(), partial)

    def _search(self, searchkey, partial=False):
        '''Takes a search string which is matched to enzyme names.
        Returns list of matched enzymes.
        If wildcard is True, all matches are returned. Else, ValueError is raised if
        more than one value is found. If no match is found, KeyError is raised.'''
        hits = []
        if searchkey in self.enzymes:
            hits.append(searchkey)
        else:
            for enzyme_name in self.enzymes:
                if searchkey.lower() in enzyme_name.lower():
                    hits.append(enzyme_name)
        if hits:
            if partial:
                # Promiscuous behaviour: output all partial/matches.
                returned_enzymes = []
                for hit in hits:
                    returned_enzymes.append(self.enzymes[hit])
                return returned_enzymes
            elif len(hits)>1:
                raise ValueError("Several matches found for query "+searchkey+": "+str(hits)) 
            else:
                # Normal behaviour: output enzyme data.
                return [self.enzymes[hits[0]]]
        else:
            raise KeyError("No matches found for "+searchkey+".")

def main():
    import argparse
    Parser = argparse.ArgumentParser(
               description=(
                 "enzymedirectory - A tool for searching for restriction "
                 "enzyme sites by name and creating json-formatted lists "
                 "of IUPAC enzyme sites for parsing with other tools."),
               epilog=(
                 "enzymedirectory is part of the PySplicer toolkit by "
                 "Cathal Garvey, available here: "
                 "https://www.github.com/cathalgarvey/pysplicer")
               )

    Parser.add_argument("-a", "--append", default=None, type=str, metavar='File-to-Append', 
                        help=
                        "Append returned target sites to a specified existing list.")
    Parser.add_argument("-f", "--file", action="store_true", help=
                         ("Specify a file instead of an enzyme name. "
                          "File should consist of enzyme names on separate lines, "
                          "with optional comments on lines starting with a '#' symbol."))
    Parser.add_argument("input", default='', type=str, metavar = 'Enzyme/File', help=
                         ('Input to work with; an enzyme name or (if using -f) a filename.'))
    Parser.add_argument("-i", "--iupac-only", action='store_true', default=False, 
                        help=
                         ("Strip non-IUPAC symbols (i.e. ^()?) from output sequences."))
    Parser.add_argument("-s", "--strip-wildcards", action="store_true", default=False,
                        help=
                         ("Strip flanking wildcards from output, if present."))
    Parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Output verbose information.")
    Parser.add_argument("-l", "--logfile", default='', type=str, metavar = 'Logfile', help=
                        "Append full output (errors included) to file")
    Parser.add_argument("-p", "--partial-match", default=False, action='store_true', help=
                        ("Return all partial matches to the search string(s)."
                        "e.g. for all B.subtilis enzymes, could use command flags"
                        """'-p "bsu"'."""))
    Args = vars(Parser.parse_args())
    with open("enzymes.json") as EnzymeRepo:
        Enzymes = json.loads(EnzymeRepo.read())
    EnzymeDir = EnzymeDirectory(Enzymes, Args['verbose'], Args['logfile'])

    # Default mode: Get enzyme sites.
    if Args['file']:
        FoundEnzymes = EnzymeDir.sites_from_file(Args['input'],
                          iupac_only=Args['iupac_only'],
                          strip_wildcards=Args['strip_wildcards'],
                          partial_matches = Args['partial_match'] )
    else:
        FoundEnzymes = EnzymeDir.enzyme_sites([Args['input']],
                          iupac_only=Args['iupac_only'],
                          strip_wildcards=Args['strip_wildcards'],
                          partial_matches = Args['partial_match'] )

    if FoundEnzymes:
        # Got enzymes to use. Print, or append to existing file.
        if Args['append']:
            with open(Args['append']) as AppendFile:
                try:
                    old_list = json.loads(AppendFile.read())
                except ValueError:
                    print("Specified append file does not appear to be valid JSON list.",
                          file=sys.stderr)
                    sys.exit(1)
            old_list.extend(FoundEnzymes)
            new_list = EnzymeDir.uniquify(old_list)
            with open(Args['append'],'w') as AppendFile:
                AppendFile.write(json.dumps(new_list))
        else:
            print(json.dumps(FoundEnzymes))
    else:
        # No matches found
        print("No matches found. Try -v for more information.", file=sys.stderr)

if __name__ == "__main__":
    main()
