import numpy as np
import sys
import random
import math
from operator import itemgetter
from collections import Counter, defaultdict

def isolate_motifs(mini_seqs, motif_length):
	motifs = [0 for subseq in mini_seqs]
	for i, (subseq, index) in enumerate(mini_seqs):
		motifs[i] = subseq[index: index + motif_length]
	return motifs

def PSSM(motifs, motif_length, background = 0.25):
	PSSM = [{'A':0.0, 'C':0.0, 'G':0.0, 'T':0.0} for i in range(motif_length)]
	pseudo = 1
	sigma  = 4
	for i in range(motif_length - 1):
		column = ""
		for motif in motifs:
			column += motif[i]
		freqs = Counter(column)
		for base, n in freqs.items():
			f = (n + pseudo) / (len(motifs) + sigma*pseudo)
			PSSM[i][base] = str(math.log2(f / background))
	return PSSM

def the_matrix(mini_seqs, motif_length, background):
	motifs = isolate_motifs(mini_seqs, motif_length)
	P = PSSM(motifs, motif_length, background)
	return P, motifs

def set_initial_positions(seqs, motif_length, seed):
	indices = [0 for seq in seqs]
	random.seed(seed)
	for i, seq in enumerate(seqs):
		indices[i] = random.randint(1, len(seq) - motif_length)
	return list(map(list, zip(seqs, indices)))

def FASTA_reader(f_ptr=None):
	f_ptr = f_ptr if f_ptr else sys.stdin.readlines()
	lines = [line.rstrip() for line in f_ptr]
	trash = [line for line in lines if line[0] == '>']
	seqs = ["" for i in range(len(trash))]
	j, current_sequence = 0, ""
	for line in lines[1:]:
		if line[0] == '>':
			seqs[j] = current_sequence
			current_sequence = ""
			j += 1
			continue
		current_sequence += line
	seqs[j] = current_sequence
	return seqs

def best_motif(star, PSSM, motif_length):
	best_motif = ["", 0, 0.0] # [motif, start, score]
	for i in range(0, len(star) - motif_length, motif_length):
		candidate_motif = star[i: i+motif_length]
		score = 0
		for pos, base in enumerate(candidate_motif):
			score += float(PSSM[pos][base])
		#print("{0:} had a score of {1:}".format(candidate_motif, score))
		best_motif = [candidate_motif, i, score] if score > best_motif[2] else best_motif
	return best_motif

def printMSA(motifs):
	for i, motif in enumerate(motifs):
		print("({0:})  {1:}".format(i, motif))
	return

def printPSSM(PSSM):
	final_PSSM = defaultdict(list)
	for column_dict in PSSM:
		for base, logodd in column_dict.items():
			final_PSSM[base].append(logodd)
	for base, row in final_PSSM.items():
		print("{}:  {}".format(base, [(math.floor(float(logodd) * 100)/100) for logodd in row]))
	return

def gibbs(iseqs, motif_length, background=0.25):
	old_motif = ["", 0, 0.0] # [motif, start, score]
	while True:
		print("-" * 50)
		star = random.randint(0, len(iseqs)-1)
		no_star = [iseqs[i] for i in range(len(seqs)) if i != star]
		PSSM, motifs = the_matrix(no_star, motif_length, background)
		new_motif = best_motif(iseqs[star][0], PSSM, motif_length)
		MSA = motifs.append(new_motif[0])
		printMSA(motifs)
		printPSSM(PSSM)
		if new_motif[0] == old_motif[0]:
			old_motif = new_motif
			break
		print("CURRENT BEST: ", new_motif[0])
		iseqs[star][1] = new_motif[1]
		old_motif = new_motif
	return old_motif

if __name__ == '__main__':
	seed = 75
	motif_length = 6
	background   = 0.25
	seqs = FASTA_reader()
	iseqs = set_initial_positions(seqs, motif_length, seed)
	m = gibbs(iseqs, motif_length, background)
	print("FINAL BEST: ", m[0])
