PySplicer
=========
A frequency/bias codon optimisation system with early NGG avoidance and optional
IUPAC match avoidance.

By Cathal Garvey, released under the GNU Affero General Public License.

Installation
------------
PySplicer is pure python and requires no compilation or building nonsense.

The easiest way to install pysplicer is using pip-3.2, if you have it installed.
easy-install may also work to install pysplicer without headaches.

The next best thing is to download the latest version from the github repository,
enter the directory through the terminal, and call "python3 setup.py install"
with administrator/sudo/root privileges.

Usage
-----
PySplicer consists of a Python Module and a python script for terminal usage.
For usage information on the script, after installing this package with pip-3.2
or setup.py, try typing "pysplicer --help" into the terminal.

What is PySplicer?
------------------
PySplicer is a Free Software implementation of a codon optimisation method where
codons are selected based on a frequency table, usually to match target host
frequencies but also (where possible) to match empirically determined
high-expression tables which usually give better expression results.

PySplicer does not (yet) implement all best practices for codon optimisation,
but compared to most available Free Software codon optimisation programs, it is
pretty good. While most programs still make use of the somewhat-discredited
"best pick" method, where the most common codon in "highly expressed" genes is
preferentially selected wherever possible, the frequency-matching approach appears
to deliver better results in general, and a frequency matching approach that
biases towards codons whose tRNAs remain charged under starvation conditions
appears to give even better results.

PySplicer can be directed to avoid DNA/RNA sequences in the final output, which
it attempts to accomplish firstly be generating large numbers of candidates and
walking through them to avoid such sites, and then by attempting to substitute
synonymous codons at each such site to remove it. Subsequences to avoid can be
given in full extended IUPAC notation, so that AWGS can refer to AAGG, AAGC, ATGG,
or ATGC. It *may* be relevant to note that all sequences are converted to DNA
internally prior to use for the sake of internal consistency,

PySplicer also maps simple hairpin structures and attempts to remove these.
However, this procress is extremely intensive (DNA/RNA structure prediction
is NP-Hard, and my implementation is pretty dumb), so expect PySplicer to be
very, very slow at optimising sequences due to this. However, this structure
avoidance/removal process can have a significant and positive impact on
successful gene expression, particularly when structures in the 5' UTR or the
first few codons of a CDS are eliminated.

Using Your Own Tables
---------------------
If a codon usage table is not available, a JSON-formatted file containing a
codon usage table can be specified. The format should be a JSON Object (akin to
a python dict with the same syntax) containing single-letter amino keys, which
point to dicts containing Codon keys, which each point to float values representing
relative frequencies. Keys are all uppercase DNA, not RNA. For example:
{"A": {"GCA": 0.15, "GCC": 0.18, "GCG": 0.51, "GCT": 0.15},
 "C": {"TGC": 0.65, "TGT": 0.35}...}

To simplify this process, a utility script is included, "cud-to-pysplicer", which
will accept a filename of a Codon Usage Database (CUD) frequency table and will
translate this table to a pysplicer-format JSON table. It accepts an optional
"-t" switch which specifies with codon table to use; it defaults to the "standard"
table. Calling cut-to-pysplicer with the "--help" switch will print usage
information and also a list of possible translation tables to specify.
