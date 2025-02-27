#Calculates fitness

def collect_real_data(real_path):
	"""Get real data from .gff file"""
	real = []
	with open(f"../smallgenes/{real_path}", "r") as file:
		for line in file:
			if "splice" in line:
				vals = line.split()
				real.append((float(vals[3]), float(vals[4]), float(vals[5])))
	return real
	
def collect_pred_data(introns):
	"""Gets frequencies from introns start and stop list"""
	freq_dict = {}
	for intron in introns:
		pair = tuple(intron)
		if pair in freq_dict:
			freq_dict[pair] += 1
		else:
			freq_dict[pair] = 1
	result = [(beg, end, count) for (beg, end), count in freq_dict.items()]
	return result
	
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

def calculate_fitness(pred_introns, real_data_path):
	"""Calculate fitness of model given outputs"""
	real = collect_real_data(real_data_path)
	pred = collect_pred_data(pred_introns)
	rp = count2prob(real)
	pp = count2prob(pred)
	fitness = calc_distance(rp, pp)
	return fitness
