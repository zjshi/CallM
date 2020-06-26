from __future__ import division

import sys, os
import vcf, argparse, operator

import numpy as np

from time import time

def parse_args():
	""" Return dictionary of command line arguments
	"""
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS)
	parser.add_argument('--dist', type=str, dest='dist_path', required=True,
		help="""Path to core-genome SNPs""")
	parser.add_argument('--tag-list', type=str, dest='tag_path', required=True,
		help="""Path to core-genome SNPs""")
	parser.add_argument('--out-dir', type=str, dest='out_dir', required=True,
		help="""Path to output directory""")
	parser.add_argument('--max-sites', dest='max_sites', default=float('inf'), type=int,
		help="""Number of sites to process from input (use all)""")
	parser.add_argument('--max-dist', dest='max_d', default=0.001, type=float,
		help="""Minimum r2 for identifying linked SNPs and loci (0.81)""")

	return vars(parser.parse_args())

def read_tags(tag_path):
	tag_map = dict()
	with open(tag_path, 'r') as fh:
		for line in fh:
			tag_genome = line.rstrip()
			tag_map[tag_genome] = 0
	
	return tag_map

def calc_tag_weights(tag_map, dist_path):
	sys.stderr.write("[clustering] start\n")

	with open(dist_path, 'r') as fh:
		for line in fh:
			items = line.rstrip().split('\t')
			genome1, genome2, d = items[0], items[1], float(items[2])

			if genome1 >= genome2:
				#sys.stderr.write("{} {}\n".format(genome1, genome2))
				continue
			# sys.stderr.write("{} {}\n".format(genome1, genome2))

			if genome1 in tag_map and genome2 in tag_map:
				tag_map[genome1] += d
				tag_map[genome2] += d

	sys.stderr.write("[clustering] done\n")

	return tag_map

def id_centroid(tag_map):
	centroid = None

	for tag in tag_map.keys():
		if centroid is None:
			centroid = tag
		else:
			if tag_map[tag] < tag_map[centroid]:
				centroid = tag

	return centroid

def main():
	args = parse_args()
	
	tag_path = args["tag_path"]
	dist_path = args["dist_path"]


	tag_map = read_tags(tag_path)
	tag_map = calc_tag_weights(tag_map, dist_path)

	centroid = id_centroid(tag_map)

	sys.stderr.write("\n{}Done!\n".format(centroid))


if __name__ == "__main__":
	main()
