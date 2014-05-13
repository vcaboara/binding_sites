#!/usr/bin/python -u
"""General description/objective.
Developed using Python 2.7.3

Binding Sites
-------------

The title of the files are important as they correspond to the part of the gene
that we're interested in.

Examples:
dmel-4-exon-r5.57.fasta.gz            dmel-4-intron-r5.57.fasta.gz           
dmel-4-five_prime_UTR-r5.57.fasta.gz  dmel-4-three_prime_UTR-r5.57.fasta.gz

If you browse through any one of them, you'll see lines looking like this:
    >FBtr0089105 type=three_prime_untranslated_region; \
	    loc=4:complement(588922..588992); name=CG1970-RA; \
	    MD5=9890392c617f7db0c6b335e5fcb3bd1f; length=71; \
	    parent=FBgn0039909; release=r5.57; species=Dmel;
    TTGTTTTTAGAGAATGTGTGTCTTCTCTGTACGAAACTGCTTAAATATAT
    AAATATAAATACGTTCAAAAT

Given a set search string (i.e. TGTGTGTGTGTG) or a text file with a list of
    search string's:
    1) Search those 4 files for the string(s)
    2) When it finds a match output the following to a text file:
    [part of gene (e.g. exon/intron/etc.)] [pattern matched] \
	    [gene that it matches] [length of sequence matched]

The "gene that it matches" is the parent=FBgn0039909 value seen above.

The complete genome will be just like these examples except it will be much
larger and called dmel-all-exon-r5.57.fasta.gz etc.

The source for all these files is located at:
    ftp://ftp.flybase.net/releases/current/dmel_r5.57/fasta/
"""

import re
import os
import csv
import gzip
import argparse

header = ["Part of gene", "Pattern Matched", "Gene Matched",
    "Length of Sequence Matched"]

def GetFileContents(filePath):
    """Return the contents of the gzipped file.
    """
    with gzip.open(filePath) as f:
	content = f.read()
	
    return content

def WriteHeader(output):
    """Write header for the csv if necessary.
    """
    # Determine if file for results exists
    if not os.path.exists(output):
	# Assume the header does not exist
	with open(output, "wb") as csvFile:
	    writer = csv.writer(csvFile, delimiter="\t",
		    quotechar="|", quoting=csv.QUOTE_MINIMAL)
	    writer.writerow(header)
	    

def SearchForBindingSite(filePath, searchString):
    """Search the file for the specified string.
    """
    content = GetFileContents(filePath)
    results = []

    for data in [c for c in content.split(">") if c]:
	seq = ''.join(data.split("\n")[1:]) 
	if searchString not in seq:
	    continue
	
	header = dict([d.strip().split("=")
	    for d in data.split("\n")[0].split(";")[1:-1]])
	results.append([searchString, header.get("parent"),
	    header.get("length")])

    return results

def FindBindingSites(args):
    """Find binding sites and write appropriate information to the output file.
    """
    searchString = args["searchString"]
    if os.path.exists(searchString):
	with open(args["searchString"]) as f:
	    searchString = [s for s in f.read().split("\n") if s]
    else:
	searchString = searchString.split()
    
    output = args["output"]
    if output is None:
	output = "%(species)s-bs_finder_results.csv" % args
    
    WriteHeader(output)
    
    for gene in args["genes"]:
	args["gene"] = gene
	p = os.path.join(args["path"],
		"%(species)s-%(chromosome)s-%(gene)s-%(version)s.fasta.gz"
		% args)
	
	if "ftp://" in args["path"]:
	    # XXX FTP location
	    print("FTP location: %s" % p)
	elif os.path.exists(p):
	    for ss in searchString:
		print("Searching %s for %s..." % (p, ss))
		result = SearchForBindingSite(p, ss)
		
		with open(output, 'ab') as csvFile:
		    writer = csv.writer(csvFile, delimiter="\t", quotechar="|",
			    quoting=csv.QUOTE_MINIMAL)
		    for r in result:
			writer.writerow([gene] + r)
	else:
	    print("Could not locate %s" % p)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
	    description="Script to identify matching genetic binding sites. " +
	    "[species]-[chromosome #|'all']-[gene]-[version].fasta.gz: " +
	    "[gene] will be taken from the list of genes.")
    parser.add_argument("searchString", action="store",
	    help="Search string to use. Can be a plain text file with " +
	    "strings on each line.")
    parser.add_argument("-p", "--path", action="store", default=".",
	    help="Path to the files (default is '.').")
    parser.add_argument("-o", "--output", action="store", default=None,
	    help="Name for the outputfile.")
    parser.add_argument("-s", "--species", action="store", default="dmel",
	    help="Species to use examine.")
    parser.add_argument("-g", "--genes", nargs="+", type=str,
	    default=["exon", "five_prime_UTR", "intron", "three_prime_UTR"],
	    help="List of genes parts to check, i.e. exon intron ...")
    parser.add_argument("-c", "--chromosome", action="store", default="4",
	    help="Chromosome # or 'all'")
    parser.add_argument("-v", "--version", action="store", default="r5.57",
	    help="Release version, i.e. 'r5.57'.")

    args = vars(parser.parse_args())
    
    FindBindingSites(args)
