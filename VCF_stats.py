#!/usr/bin/env python3
# coding: utf-8

"""
Usage: VCF_stats.py -c <INT> [-h]

  [Options]
    -c, --callable <INT>                       Number of callable sites
    -h, --help                                 Show this message
"""

import os
import sys
import collections
import itertools
from docopt import docopt

args = docopt(__doc__)
callable_sites = int(args["--callable"])

dxy_dict = collections.defaultdict(float)
het_dict = collections.defaultdict(float)

for line in sys.stdin:
	line = line.strip()
	elements = line.split()

	if len(elements[3].split(",")) == 1: # only biallelic SNPs

		genotypes = {} # collect genotypes for this SNP
		for i in range(4, len(elements)):
			ID, GT = elements[i].split("=")
			genotypes[ID] = set()
			genotypes[ID].add(GT[0])
			genotypes[ID].add(GT[2])

		for combo in itertools.combinations(genotypes.keys(), 2): # do pairwise comparisons for dxy
			if len(genotypes[combo[0]].union(genotypes[combo[1]])) == 1: # no differences
				pass
			elif len(genotypes[combo[0]]) == 1 and len(genotypes[combo[1]]) == 1: # fixed diff
				dxy_dict[frozenset(combo)] += 1
			else: # non-fixed diff
				dxy_dict[frozenset(combo)] += 0.5

		for ID in genotypes: # check which are heterozygous
			if len(genotypes[ID]) == 2:
				het_dict[ID] += 1

print("###indA\tindB\tFst\tdxy")

for comparison in dxy_dict: # print pairwise stats
	dxy = dxy_dict[comparison]
	mean_het = 0
	taxa = [taxon for taxon in comparison]
	mean_het += het_dict[taxa[0]]
	mean_het += het_dict[taxa[1]]
	mean_het = mean_het / 2
	Fst = (dxy - mean_het) / (dxy + mean_het)
	print("{}\t{}\t{}\t{}".format(taxa[0], taxa[1], round(Fst, 5), dxy / callable_sites))

print("###ind\theterozygosity")

for ID in genotypes: # print heterozygosity
	print("{}\t{}".format(ID, het_dict[ID] / callable_sites))

sys.exit()

