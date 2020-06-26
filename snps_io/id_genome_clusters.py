from __future__ import division

import sys, os
import vcf, argparse, operator

import numpy as np

from time import time

class GenomeCluster:
	def __init__(self, max_d=0.001):
		self.max_d  = max_d

		self.genomes = dict()
		self.links = dict()

		self.tag_genome = None

	def size(self):
		return len(genomes)

	def add(self, genome1, genome2, d):
		if genome1 not in self.genomes:
			self.genomes[genome1] = 1
		else:
			self.genomes[genome1] = self.genomes[genome1] + 1

		if genome2 not in self.genomes:
			self.genomes[genome2] = 1
		else:
			self.genomes[genome1] = self.genomes[genome1] + 1

		link = "{}|{}".format(genome1, genome2)
		if genome1 > genome2:
			link = "{}|{}".format(genome2, genome1)

		#sys.stderr.write("{} {}\n".format(genome1, genome2))

		if link not in self.links:
			self.links[link] = d
		else:
			#assert True
			assert self.links[link] == d
			#sys.stderr.write("{} {}\n".format(genome1, genome2))
			

	def merge(self, cluster2, genome1, genome2, d):
		if self.contains(genome1) and cluster2.contains(genome2):
			self.add(genome1, genome2, d)
			for link in cluster2.links.keys():
				genomes = link.split("|")
				self.add(genomes[0], genomes[1], cluster2.links[link])
		else:
			sys.exit("{} is not in cluster1 or {} is not cluster2, nullify the basis for merging".format(genome1, genome2))

	def is_empty(self):
		return len(self.genomes.keys()) == 0

	def contains(self, genome):
		return genome in self.genomes

	def id_tag_genome(self):
		# no snps
		if len(self.genomes) == 0:
			sys.exit("\nError: no genomes on cluster: cannot id tag genome\n")
		# one snp
		elif len(self.genomes) == 1:
			self.tag_genome = self.genomes.keys()[0]
		else:
			tmp_min = 0
			tmp_genome = None
			for genome in self.genomes.keys():
				if self.genomes[genome] > tmp_min:
					tmp_min = self.genomes[genome]
					tmp_genome = genome
			self.tag_genome = tmp_genome

		return self.tag_genome

	def fmtout(self):
		sorted_tuples = sorted(self.genomes.items(), key=operator.itemgetter(1), reverse=True)
		sorted_genomes = [genome_tuple[0] for genome_tuple in sorted_tuples]

		return "* {} {}".format(self.tag_genome, " ".join(sorted_genomes))

	def fmtout_all(self):
		fmt_str = "{}\n".format(self.fmtout())

		for link in self.links.keys():
			fmt_str += "- {} {}\n".format(link, self.links[link])

		return fmt_str

def parse_args():
	""" Return dictionary of command line arguments
	"""
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS)
	parser.add_argument('--dist', type=str, dest='dist_path', required=True,
		help="""Path to core-genome SNPs""")
	parser.add_argument('--out-dir', type=str, dest='out_dir', required=True,
		help="""Path to output directory""")
	parser.add_argument('--max-sites', dest='max_sites', default=float('inf'), type=int,
		help="""Number of sites to process from input (use all)""")
	parser.add_argument('--max-dist', dest='max_d', default=0.001, type=float,
		help="""Minimum r2 for identifying linked SNPs and loci (0.81)""")

	return vars(parser.parse_args())

def search_genome_clusters(dist_path, args):
	sys.stderr.write("[clustering] start\n")

	genome_clusters = []
	genome_lookup = dict()

	with open(dist_path, 'r') as fh:
		for line in fh:
			items = line.rstrip().split('\t')
			genome1, genome2, d = items[0], items[1], items[2]

			if genome1 >= genome2:
				#sys.stderr.write("{} {}\n".format(genome1, genome2))
				continue
			# sys.stderr.write("{} {}\n".format(genome1, genome2))

			if genome1 not in genome_lookup and genome2 not in genome_lookup:
				new_cluster = GenomeCluster(float(args['max_d']))
				new_cluster.add(genome1, genome2, d)
				genome_clusters.append(new_cluster)
				genome_lookup[genome1] = len(genome_clusters) - 1
				genome_lookup[genome2] = len(genome_clusters) - 1
			elif genome1 in genome_lookup and genome2 not in genome_lookup:
				cluster_indx = genome_lookup[genome1]
				genome_lookup[genome2] = cluster_indx
				genome_clusters[cluster_indx].add(genome1, genome2, d)
			elif genome1 not in genome_lookup and genome2 in genome_lookup:
				cluster_indx = genome_lookup[genome2]
				genome_lookup[genome1] = cluster_indx
				genome_clusters[cluster_indx].add(genome1, genome2, d)
			else:
				if genome_lookup[genome1] == genome_lookup[genome2]:
					pass
				else:
					cluster_indx1 = genome_lookup[genome1]
					cluster_indx2 = genome_lookup[genome2]
					genome_clusters[cluster_indx1].merge(genome_clusters[cluster_indx2], genome1, genome2, d)

					for genome in genome_clusters[cluster_indx2].genomes:
						genome_lookup[genome] = cluster_indx1

					genome_clusters[cluster_indx2] = None

	sys.stderr.write("[clustering] done\n")
	sys.stderr.write("[clustering] {} genomes have been included in clusters\n".format(len(genome_lookup.keys())))

	good_clusters = verify_clusters(genome_clusters, genome_lookup)

	for gcluster in good_clusters:
		print gcluster.fmtout_all()

	return good_clusters

def verify_clusters(genome_clusters, genome_lookup):
	for genome in genome_lookup.keys():
		assert genome in genome_clusters[genome_lookup[genome]].genomes

	good_clusters = []
	for i, cluster in enumerate(genome_clusters):
		if cluster is not None:
			cluster.id_tag_genome()
			good_clusters.append(cluster)

			for genome in cluster.genomes:
				assert genome_lookup[genome] == i

	return good_clusters

def build_genome_blocks(args):
	dist_path = args["dist_path"]

	genome_clusters = search_genome_clusters(dist_path, args)

def main():
	args = parse_args()
	args['out_dir'] = args['out_dir'].rstrip('/')
	genome_blocks = build_genome_blocks(args)

	sys.stderr.write("\nDone!\n")


if __name__ == "__main__":
	main()
