import os
import argparse
import subprocess

parser = argparse.ArgumentParser(description="Find distribution of predicted isoforms and calculate accuracy")
parser.add_argument("gene", type=str, help="Name of gene to run on, ie.'ch.1_0'")
parser.add_argument("--u1num", type=int, default = 5, help="Maximum amount of u1's that can bind")
parser.add_argument("--u5num", type=int, default=5, help="Maximum amount of u5's that can bind")
parser.add_argument("--z", type=int, default=10, help="Iteration parameter")
arg = parser.parse_args()
 
def collect_data(real_path = f"{arg.gene}.gff3"):
	"""Run 10k samples and get data from gff file, output file"""
	
	#with open("args.txt", "w") as file:
		#for i in range(1, 10001):
			#f.write(f"../smallgenes/{arg.gene}.fa {arg.u1num} {arg.u5num} {arg.z}\n")
	command = f"seq 10000 | parallel -j 5 --results results_dir python3 main.py ../smallgenes/{arg.gene}.fa --u1num {arg.u1num} --u5num {arg.u5num} --z {arg.z}")
    os.system(command)
	os.system("cat results_dir/*/stdout > outs.txt")
	os.getpid()
	#Predicted data
	pred_path = "outs.txt"
	line_counts = {}
	with open(pred_path, "r") as file:
		for line in file:
			line = line.strip().replace('"', '')
			start, end = map(int, line.split())
			line_counts[(start, end)] = line_counts.get((start, end), 0) + 1
	pred = [(start, end, count) for (start, end), count in line_counts.items()]
	
	#Real data
	real = []
	with open(f"../smallgenes/{real_path}", "r") as file:
		for line in file:
			if "splice" in line:
				vals = line.split()
				real.append((float(vals[3]), float(vals[4]), float(vals[5])))
	return real, pred
	
def count2prob(introns):
	"""Get probabilities from counts"""
	outs = {}
	total = sum([x[2] for x in introns])
	for beg, end, n in introns:
		outs[(beg, end)] = n/total
	
	return outs

def calc_distance(rp, pp):
	"""Calculate Manhattan Distance"""
	fitness = 0
	for coor in rp.keys() | pp.keys():
		if coor in rp and coor in pp:
			fitness += abs(rp[coor]-pp[coor])
		elif coor in pp:
			fitness += pp[coor]
		else:
			fitness += rp[coor]
	return fitness
				
	
real, pred = collect_data()
rp = count2prob(real)
pp = count2prob(pred)
fitness = calc_distance(rp, pp)
print("hi")
print(fitness)
