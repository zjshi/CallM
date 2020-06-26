from __future__ import division

import gc
import numpy as np

"""
Class for encapsulating a single Multiple Sequence Alignment(MSA)
"""
class Alignment:
	def __init__(self):
		self.desc = '' #not used for now
		# chrom has identifier like 'cluster123', which mark a division on the core-genome, which is a conserved region found cross all reference genomes
		self.chrom = ''
		self.nseqs = 0
		self.ncols = 0 # number of sites
		self.seqs = [] # actual sequence for each sample
		self.sample_ids = []

		"""
		attributes below were described in function update()
		"""
		self.char_mat = []
		self.local_pos = []
		self.count_mat = []
		self.freq_mat = []

		self.ref_alleles = []
		self.alt_alleles = []
		self.third_alleles = []
		self.forth_alleles = []

		self.ref_prob_mat = []
		self.alt_prob_mat = []
		self.third_prob_mat = []
		self.forth_prob_mat = []

		self.sample_presence = []
		self.ref_freqs = []
		self.alt_freqs = []
		self.third_freqs = []
		self.forth_freqs = []

		self.prevalence = []
		self.aligned_pctg = []

	def update(self):
		assert len(self.seqs) > 1
		self.nseqs = len(self.seqs)
		self.ncols = len(self.seqs[0].seq)
		self.sample_ids = [seq.id for seq in self.seqs]

		"""
		generate char matrix from the aligned sequences

		example:
		sample1: ATCG
		sample2: ATGG
		sample3: ATGC

		the char matrix is the transpose of [[A, T, G, C],[A, T, G, G],[A, T, G, C]]
		"""
		self.char_mat = np.array([np.fromstring(seq.seq, dtype='c') for seq in self.seqs])

		#print self.char_mat.shape
		#print self.char_mat.nbytes

		"""
		count A, T, G, C, N and - for each site on the sequences
		local_pos stores all local positions of each site on this core-genome division (or alignment)
		"""
		As = (self.char_mat == 'A').sum(0)
		Ts = (self.char_mat == 'T').sum(0)
		Gs = (self.char_mat == 'G').sum(0)
		Cs = (self.char_mat == 'C').sum(0)
		Ns = (self.char_mat == 'N').sum(0)
		Gaps = (self.char_mat == '-').sum(0)

		self.local_pos = np.arange(len(seq.seq))

		self.count_mat = np.array([As, Ts, Gs, Cs, Ns, Gaps])

		"""
		char_template: complete set of chars for each site on the sequences, from which the ref allele and alt allele will be selected
		"""
		char_template = np.array([
			np.repeat('A', self.count_mat.shape[1]),
			np.repeat('T', self.count_mat.shape[1]),
			np.repeat('G', self.count_mat.shape[1]),
			np.repeat('C', self.count_mat.shape[1])
			# np.repeat('N', self.count_mat.shape[1])
			# np.repeat('-', self.count_mat.shape[1])
		])

		"""
		sorting the counts of different chars at each site,
		then select the indices of chars whose counts are in top 2,
		then using the indices of chars to select the chars,
		then the top 1 char (the char with highest count) at each site is considered as ref allele for now,
		then the char with the second highest count at each site is considered as alt allele for now,
		for those sites that have only one char, the counts of other chars are zeroes, so theorectically any of them can be selected by program, but in reality '-' will be selected, no exception was found.
		"""
		count_inds_mat = self.count_mat[0:4,:].argsort(axis=0)
		top2_inds = count_inds_mat[-4:,]
		top2_char_mat = np.choose(top2_inds, char_template)
		self.ref_alleles = top2_char_mat[3,:]
		self.alt_alleles = top2_char_mat[2,:]

		self.third_alleles = top2_char_mat[1,:]
		self.forth_alleles = top2_char_mat[0,:]

		"""
		frequency matrix has the same shape as the char matrix
		it is initialize to have None only,
		in the end, it is used to store the presence/absence of ref/alt allele for each site cross all samples with the following rules:
		- presence of ref allele: 1
		- absence of ref allele: 0
		- presence of N or -: None
		"""
		self.freq_mat = np.repeat(np.int8(-1), self.char_mat.shape[0]*self.char_mat.shape[1]).reshape(self.char_mat.shape)


		"""
		these two masks have the same shape as the frequence matrix
		"""
		# ref_mask = ((self.char_mat == self.ref_alleles) & (self.char_mat != '-') & (self.char_mat != 'N'))
		# alt_mask = ((self.char_mat == self.alt_alleles) & (self.char_mat != '-') & (self.char_mat != 'N'))

		ref_mask = (self.char_mat == self.ref_alleles)
		alt_mask = (self.char_mat == self.alt_alleles)
		third_mask = (self.char_mat == self.third_alleles)
		forth_mask = (self.char_mat == self.forth_alleles)
		"""
		such operation is possible because, I guess, numpy store matrix in a gigantic 1D array.
		"""
		self.freq_mat[ref_mask] = 0
		self.freq_mat[alt_mask] = 1
		self.freq_mat[third_mask] = 2
		self.freq_mat[forth_mask] = 3

		self.ref_prob_mat = np.repeat(np.int8(0), self.char_mat.shape[0]*self.char_mat.shape[1]).reshape(self.char_mat.shape)
		self.alt_prob_mat = np.repeat(np.int8(0), self.char_mat.shape[0]*self.char_mat.shape[1]).reshape(self.char_mat.shape)
		self.third_prob_mat = np.repeat(np.int8(0), self.char_mat.shape[0]*self.char_mat.shape[1]).reshape(self.char_mat.shape)
		self.forth_prob_mat = np.repeat(np.int8(0), self.char_mat.shape[0]*self.char_mat.shape[1]).reshape(self.char_mat.shape)



		self.ref_prob_mat[ref_mask] = 1
		self.alt_prob_mat[alt_mask] = 1
		self.third_prob_mat[third_mask] = 1
		self.forth_prob_mat[forth_mask] = 1

		ref_counts = ref_mask.sum(0)
		alt_counts = alt_mask.sum(0)
		third_counts = third_mask.sum(0)
		forth_counts = forth_mask.sum(0)
		"""
		the sample presence here only sum up the counts of A, T, G, C for each site, leave out the N and -,
		it facilitate the calculation of prevalence of the site on certain sequence alignment position
		ref and alt allele frequencies were calculated with denominator of sample_presence rather than the number of sample.
		naturally, there is another route to calculate them.
		"""
		self.sample_presence = self.count_mat[0:4,:].sum(0)
		self.prevalence = self.sample_presence/self.nseqs

		zero_mask = (self.sample_presence == 0)
		self.sample_presence[zero_mask] = 1
		self.ref_freqs = ref_counts/self.sample_presence
		self.alt_freqs = alt_counts/self.sample_presence
		self.third_freqs = third_counts/self.sample_presence
		self.forth_freqs = forth_counts/self.sample_presence

		self.ref_freqs[zero_mask] = 0
		self.alt_freqs[zero_mask] = 0
		self.third_freqs[zero_mask] = 0
		self.forth_freqs[zero_mask] = 0

		self.sample_presence[zero_mask] = 0

		unaligned_masks = (self.count_mat[4:6,:] != 0)
		self.aligned_pctg = 1 - (unaligned_masks.sum(0) / self.nseqs)
