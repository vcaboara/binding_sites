#!/usr/bin/python -u
"""General description/objective.
Developed using Python 2.7.3

Binding Sites
-------------

The title of the files are important as they correspond to the part of the gene
that we're interested in.

Examples:
dmel-4-exon-r5.56.fasta.gz            dmel-4-intron-r5.56.fasta.gz           
dmel-4-five_prime_UTR-r5.56.fasta.gz  dmel-4-three_prime_UTR-r5.56.fasta.gz

If you browse through any one of them, you'll see lines looking like this:
    >FBtr0089105 type=three_prime_untranslated_region; \
	    loc=4:complement(588922..588992); name=CG1970-RA; \
	    MD5=9890392c617f7db0c6b335e5fcb3bd1f; length=71; \
	    parent=FBgn0039909; release=r5.56; species=Dmel;
    TTGTTTTTAGAGAATGTGTGTCTTCTCTGTACGAAACTGCTTAAATATAT
    AAATATAAATACGTTCAAAAT

Given a set regex (i.e. (?x)(TG){4,} ) or a text file with a list of regex's:
    1) Search those 4 files for the regex
    2) When it finds a match output the following to a text file:
    [part of gene (e.g. exon/intron/etc.)] [pattern matched] \
	    [gene that it matches] [length of sequence matched]

The "gene that it matches" can be snagged from parent=FBgn0039909 as seen above.

The complete genome will be just like these examples except it will be much
larger and called dmel-all-exon-r5.56.fasta.gz etc.

The source for all these files is located at:
    ftp://ftp.flybase.net/releases/current/dmel_r5.56/fasta/
"""

import re
import os
import csv
import gzip
import argparse

def GetFileContents(filePath):
    """Return the contents of the gzipped file.
    """
    with gzip.open(filePath) as f:
	content = f.read()
	
    return content
    
def GetRegEx(regex):
    """Return the regex object(s) to use.
    """
    if os.path.exists(regex):
	# RegEx provided is a file
	with open(regex) as f:
	    myREs = [re.compile(tempRE) for tempRE in f.read().split("\n")]
    else:
	myREs = [re.compile(regex)]
	
    return myREs

def SearchForBindingSite(filePath, regex):
    """Search for the file content for the regex.
    """
    fileContent = GetFileContents(filePath)
    res = regex.search(fileContent)
    match = []
    
    if res is None:
	print("%s not found in %s." % (res.pattern, filePath))
	return None
    
    for groupNo in xrange(len(res.groups())):
	# XXX Code needs to be altered for groupNo?
	data = res.string[:res.string.find(">", res.end())].split(">")[-1]
	dataDict = dict([d.strip().split("=") for d in data.split(";")[1:-1]])
	match.append([res.group(groupNo), dataDict['parent'],
	    dataDict['length']])
    
    return match

def FindBindingSites(args):
    """Find binding sites and write appropriate information to the output file.
    """
    myRegExs = GetRegEx(args['regex'])
    
    output = args['output']
    if output is None:
	output = "%(species)s-bs_finder_results.csv" % args
    
    for gene in args['genes']:
	args['gene'] = gene
	p = os.path.join(args['path'],
		"%(species)s-%(chromosome)s-%(gene)s-%(version)s.fasta.gz"
		% args)
	
	if os.path.exists(p):
	    for regex in myRegExs:
		result = SearchForBindingSite(p, regex)
		
		if result:
		    with open(output, 'ab') as csvfile:
			writer = csv.writer(csvfile, delimiter='\t',
				quotechar='|', quoting=csv.QUOTE_MINIMAL)
			for r in result:
			    writer.writerow([gene] + r)
	# XXX elif FTP:
	else:
	    print("Could not locate %s" % p)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
	    description="Script to identify matching genetic binding sites. " +
	    "[species]-[chromosome #|'all']-[gene]-[version].fasta.gz: " +
	    "[gene] will be taken from the list of genes.")
    parser.add_argument("path", action="store", help="Path to the files.")
    parser.add_argument("regex", action="store", help="RegEx to use. " +
	    "Can be a plain text file with regexs on each line.")
    parser.add_argument("-o", "--output", action="store", default=None,
	    help="Name for the outputfile.")
    parser.add_argument("-s", "--species", action="store", default="dmel",
	    help="Species to use examine.")
    parser.add_argument("-c", "--chromosome", action="store", default="4",
	    help="Chromosome # or 'all'")
    parser.add_argument("-g", "--genes", nargs="+", type=str,
	    default=["exon", "five_prime_UTR", "intron", "three_prime_UTR"],
	    help="List of genes parts to check, i.e. exon intron ...")
    parser.add_argument("-v", "--version", action="store", default="r5.56",
	    help="Release version, i.e. 'r5.56'.")

    args = vars(parser.parse_args())
    
    FindBindingSites(args)
