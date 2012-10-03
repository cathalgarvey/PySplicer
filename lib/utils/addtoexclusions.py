#!/usr/bin/env python3
# Quick script to append DNA/RNA sequences to enzyme exclusion/ignore lists.
import json
import argparse
import sys

Parser = argparse.ArgumentParser(description="Simple script to add string values to a JSON list file. It is included as part of the PySplicer toolkit to offer a trivial way to add arbitrary sites to exclusion lists before using these lists for genetic optimisation. For example, if you want to add an obscure enzyme recognition sequence not found in the enzyme database included with PySplicer, you can use this script to add it easily, rather than trying to add the site to the database.", epilog="By Cathal Garvey. Part of the PySplicer toolkit, released under the GNU General Public License. Code can be found here: https://www.github.com/cathalgarvey/pysplicer")
Parser.add_argument("file", type=str, help="JSON file to append to.")
Parser.add_argument("value", type=str, help="Value to append to the file.")
Parser.add_argument("-u", "--force-unique", default=False, action='store_true',
                    help="Edits list after appending value to remove all list duplicates.")
Args = Parser.parse_args()

def uniquify(seq):
    '''This is a fast, order-preserving function for removing duplicates 
    from a sequence/list.
    Found here: http://www.peterbe.com/plog/uniqifiers-benchmark'''
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]

def open_and_append():
    with open(Args.file) as AppendFile:
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

try:
    with open(Args.file) as AppendFile:
        try:
            json_list = json.loads(AppendFile.read())
            assert isinstance(json_list, list)
        except ValueError:
            print("Specified append file does not appear to be valid JSON file.",
                  file=sys.stderr)
            sys.exit(1)
        except AssertionError:
            print(("JSON contents are not a simple list."
                   " This script can't handle other data types, sorry."),
                   file=sys.stderr)
except IOError:
    with open(Args.file, 'w') as NewFile:
        NewFile.write(json.dumps([Args.value]) + "\r\n")
    sys.exit(0)

json_list.append(Args.value)
if Args.force_unique: json_list = uniquify(json_list)

try:
    with open(Args.file,'w') as AppendFile:
        AppendFile.write(json.dumps(json_list) + "\r\n")
except IOError:
    print("Cannot write to file; do you have write access?", file=sys.stderr)
    sys.exit(1)
