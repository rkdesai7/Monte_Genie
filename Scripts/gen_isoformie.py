#main script
from genie2 import sequence, snRNP, u1, u5
from fitness import calculate_fitness, collect_pred_data, count2prob
import gzip
import sys
import argparse
import numpy as np
import csv
import os

parser = argparse.ArgumentParser(description="Return all possible isoforms based on a randomization algorithm")
parser.add_argument("sequences", type=str, help="Path to fasta files of sequences you want to decode")
parser.add_argument("--u1num", type=int, default = 5, help="Maximum amount of u1's that can bind")
parser.add_argument("--u5num", type=int, default=5, help="Maximum amount of u5's that can bind")
parser.add_argument("--z", type=int, default=10, help="Iteration parameter")
parser.add_argument("--acceptor_pwm", type=str, default = "models/acc.pwm", help="Path to a .pwm file containing acceptor site probabilities")
parser.add_argument("--donor_pwm", type=str, default = "models/don.pwm", help="Path to a .pwm file containing donor site probabilities")
parser.add_argument("--runs", type=int, default=1000, help="How many isoforms you want to generate (how many times you want to run the simulation")

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
	
test_sequences = readfasta(arg.sequences)
curr_sequence = next(test_sequences)
curr_sequence = curr_sequence[1]
seq_name = curr_sequence[0]

#Run simulation
introns = []
isoforms = []
count = 0
for j in range(arg.runs):
	
	#Initialize
	u1s = []
	u5s = []
	for i in range(arg.u1num):
		u1s.append(u1(arg.donor_pwm))
	for i in range(arg.u5num):
		u5s.append(u5(arg.acceptor_pwm))
	seq = sequence(curr_sequence, seq_name, arg.z)
	
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
	intron = seq.compile_intron()
	count += seq.detect_intrasplicing()
	isoforms.append(intron)
	introns.extend(intron)

#Calculate Probabilites
pred = collect_pred_data(introns)
probs = count2prob(pred)
total = len(introns)
#Calculate fitness
real_path = arg.sequences.replace(".fa", ".gff3")
fit = calculate_fitness(introns, real_path)
print(fit)

#write output to a file
out_name = os.path.splitext(os.path.basename(arg.sequences))[0]
with open(f"{out_name}.txt", 'w') as f:
	f.write(f"Number of Intraspliced Introns: {count}, {total}\n")
	f.write("Fitness:")
	f.write(str(fit))
	f.write("\nSplicing Events:\n")
	for i in seq.splicing_events:
		f.write(str(i))
		f.write("\n")
	f.write("Frequencies:\n")
	for key,value in probs.items():
		statement = str(key) + ": " + str(value)
		f.write(statement)
		f.write("\n")
		
	


		
	

	
