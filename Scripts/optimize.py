import argparse
import json
import os
import random
import main
import subprocess

def get_fitness(wormie):
	u5num = wormie['genotype']['--u5num']
	u1num = wormie['genotype']['--u1num']
	rez = wormie['genotype']['--rez']
	result=subprocess.run(['python3', 'fitness.py', f'{arg.gene}', '--u1num', f'{u1num}', '--u5num', f'{u5num}', '--z', f'{rez}'], capture_output=True, text=True)
	return result.stdout

def random_wormie():
	wormie = {
				'genotype':{
					'--u5num': random.randint(1, arg.num_snRNP),
					'--u1num': random.randint(1, arg.num_snRNP),
					'--rez': random.randint(1, arg.resolution),
				},
				'fitness': None
			}
			
	return wormie
	
def mate(p1, p2, mut):
	child = {
			'genotype': {},
			'fitness': None
			}
	weight = {'--u5num', '--u1num', '--rez'}
	for k in weight:
		if random.random() < 0.5: child['genotype'][k] = p1['genotype'][k]
		else:					  child['genotype'][k] = p2['genotype'][k]
		if random.random() < mut: 
			new = random_wormie()
			child['genotype'][k] = new['genotype'][k];
	return child
	
parser = argparse.ArgumentParser(description='Parameter Optimization Program')
parser.add_argument('gene', type=str, help="name of gene to finetune on ie. ch.1_0")
parser.add_argument('--pop', required=False, type=int, default=100, help='population size [%(default)i]')
parser.add_argument('--gen', required=False, type=int, default=100, help='generations [%(default)i]')
parser.add_argument('--die', required=False, type=float, default=0.5, help='fraction that die each gen [%(default).2f]')
parser.add_argument('--mut', required=False, type=float, default=0.1, help='mutation frequency [%(default).2f]')
parser.add_argument('--seed', required=False, type=int, help='random seed')
parser.add_argument('--name', required=False, type=str, default='', help='name the output')
parser.add_argument('--num_snRNP', required = False, type=int, default=20, help='range for snRNP values')
parser.add_argument('--resolution', required=False, type=int, default=40, help='range for sequence resolution (how often splicing and transcripition occurs)')
parser.add_argument('--verbose', action='store_true', help='show progress')
arg = parser.parse_args()

#Initialize
if arg.seed: random.seed(arg.seed)
pop = []
for i in range(arg.pop): pop.append(random_wormie())
for wormie in pop: wormie['fitness'] = get_fitness(wormie)

#Evolve population
half = int(len(pop)*arg.die)
for g in range(arg.gen):
	pop = sorted(pop, key=lambda item: item['fitness'])
	if arg.verbose: print(f'generation: {g}, fitness: {pop[0]["fitness"]}')
	
	#mate
	children = []
	for i in range(half, len(pop)):
		p1 = random.randint(0, half)
		p2 = random.randint(0, half)
		pop[i] = mate(pop[p1], pop[p2], arg.mut)
		children.append(pop[i])
		
	#fitness
	for child in children: child['fitness'] = get_fitness(child)

#Final report
pop = sorted(pop, key=lambda item: item['fitness'])
best = pop[0]
print(best["fitness"], end='\t')
for prop, val in best['genotype'].items():
	print(f'{val:.4f}', end='\t')
print(arg.name)
