import csv
import os
import sys
import time

from pathlib import Path

folder = Path("/home/rkdesai7/Documents/Monte_Genie/smallgenes")
data = []
for file in folder.glob("*.fa"):
	with open(file, "r") as f:
		seq = "".join(line.strip() for line in f if not line.startswith(">"))
	start = time.perf_counter()
	os.system(f"python3 gen_isoformie.py {file} --runs 1")
	end = time.perf_counter()
	data.append((len(seq), end-start))

filename = 'runtimes.csv'
with open(filename, mode='w', newline='') as file:
	writer = csv.writer(file)
	writer.writerow(['Length', 'Runtime'])
	writer.writerows(data)
	
	
