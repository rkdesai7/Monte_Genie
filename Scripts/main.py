#main script
from optimized_genie import sequence, snRNP, u1, u5
from fitness import calculate_fitness
import gzip
import sys
import argparse
import numpy as np
import csv
import pandas as pd

parser = argparse.ArgumentParser(description="Return all possible isoforms based on a randomization algorithm")
parser.add_argument("sequences", type=str, help="Path to fasta files of sequences you want to decode")
parser.add_argument("--u1num", type=int, default = 5, help="Maximum amount of u1's that can bind")
parser.add_argument("--u5num", type=int, default=5, help="Maximum amount of u5's that can bind")
parser.add_argument("--z", type=int, default=10, help="Iteration parameter")
parser.add_argument("--acceptor_pwm", type=str, default = "models/acc.pwm", help="Path to a .pwm file containing acceptor site probabilities")
parser.add_argument("--donor_pwm", type=str, default = "models/don.pwm", help="Path to a .pwm file containing donor site probabilities")

arg = parser.parse_args()
def readfasta(filename):
	
	name = None
	seqs = []
	
	fp = None
	if filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-': fp = sys.stdin
	else: fp = open(filename)
	
	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()
	

def display_exon(transcript):
	""" Convert to ideal display """
	indexes = np.array(transcript)[:,1]
	start = indexes[0]
	temp = indexes[0]
	end = None
	storage = []
	for index, value in enumerate(indexes):
		if start == value: continue
		if index == len(indexes)-1:
			end = indexes[index]
			storage.append((start, end))
		if int(value) == int(temp) + 1: 
			temp = value
			continue
		if (value != int(temp) + 1):
			end = indexes[index - 1]
			storage.append([start, end])
			start = value
			temp = value
	
	#Output
	for i in storage:
		print("exon", i[0], i[1])

def display_intron(transcript):
	introns = []
	indexes = np.array(transcript)[:,1]
	prev = int(indexes[0])
	for index, value in enumerate(indexes):
		value = int(value)
		if prev == value: continue
		if value == prev + 1:
			prev = value
		else:
			introns.append([prev+1, value-1])
			prev = value
	#for i in introns:
		#print(i[0], i[1])
		
	return introns
#Create Objects
def main(u1arg, u5arg, acc_arg, don_arg, seq_arg, z_arg):
	
	test_sequences = readfasta(seq_arg)
	curr_sequence = next(test_sequences)
	curr_sequence = curr_sequence[1]
	seq_name = curr_sequence[0]
	
	#Run simulation
	introns = []
	for j in range(1000):
		
		#Initialize
		u1s = []
		u5s = []
		for i in range(u1arg):
			u1s.append(u1(don_arg))
		for i in range(u5arg):
			u5s.append(u5(acc_arg))
		seq = sequence(curr_sequence, seq_name, z_arg)
		
		#Run
		i = 1
		final = 0
		while final < 50:
			seq.one_iteration(i)
			for u in u1s: u.one_iteration(seq)
			for u in u5s: u.one_iteration(seq)
			i += 1
			if seq.curr == len(curr_sequence):
				final += 1
		
		#display
		intron = display_intron(seq.transcript)
		introns.extend(intron)
	
	#Calculate fitness
	real_path = seq_arg.replace(".fa", ".gff3")
	fit = calculate_fitness(introns, real_path)
	print(fit)
	

if __name__ == "__main__":
	main(arg.u1num, arg.u5num, arg.acceptor_pwm, arg.donor_pwm, arg.sequences, arg.z) 
		
	

	
