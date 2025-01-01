#main script
from genie import sequence, snRNP, u1, u5

parser = argparse.ArgumentParser(description="Return all possible isoforms based on a randomization algorithm")
parser.add_argument("sequence", type=str, help="Path to fasta files of sequences you want to decode")
parser.add_argument("acceptor_pwm", type=str, help="Path to a .pwm file containing acceptor site probabilities")
parser.add_argument("donor_pwm", type=str, help="Path to a .pwm file containing donor site probabilities")
parser.add_argument("--u1num", type=int, default = 1, help="Maximum amount of u1's that can bind")
parser.add_argument("--u5num", type=int, default=1, help="Maximum amount of u5's that can bind")
parser.add_argument("--z", type=int, default=1, help="Iteraion parameter")
arg = parser.parse_args()

#Read FASTA File
#test_sequence = readfasta(arg.sequence)

#storage of muliptle snRNP's
#option1: list
u1s = []
for i in range(arg.u1num):
	a_u1 = u1(arg.donor_pwm)
	u1s.append(a_u1)
	
#option2: dict
u1s = {}
for i in range(arg.u1num):
	key = f"u1_{i}"
	u1s[key] = u1(arg.donor_pwm)
	
iters = 0
my_seq = sequence(test_sequence, "genie")
while True:
	my_seq.one_iteration(iters)
	for i in range(arg.u1num):
		curr_u1 = u1s[i]
		curr_u1.one_iteration(my_seq)
	for i in range(arg.u5num):
		curr_u5 = u5s[i]
		curr_u5.one_iteration(my_seq)

	
	
