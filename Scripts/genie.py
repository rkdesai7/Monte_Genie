import random
import math
import argparse
import sys
import gzip
import statistics as stats

from calculate_intron_prob import intron_probability

class sequence:
	def __init__(self, sequence, name, resolution=10):
		""" Sequence Object 
		- sequence: entire sequence that we are working with
		- name: name of the sequence
		- curr: current index of sequence that has been transcribed
		- transcript: the part of the sequence that has been transcribed and where splicing may occur
		- bindings: storage of what snRNPs are bound to parts of the transcript
		- resolution: how many iterations until a new base is transcribed
		- u1_bound_ids: collection of u1s currently bound to the transcript
		- u5_bound_ids: collection of u5s currently bound to the transcript
		- splicing_events: storage of splicing events that occur"""
		self.sequence = sequence
		self.name = name
		self.curr = 0
		self.transcript = []
		self.bindings = [] 
		self.resolution = resolution
		self.u1_bound_ids = []
		self.u5_bound_ids = []
		self.splicing_events = []
	def transcribe(self):
		""" Transcribe one base pair in the sequence """
		new_base = self.sequence[self.curr]
		self.transcript.append(new_base)
		self.curr += 1
		self.bindings.append(0)
		return self.curr
	def seq_available(self, size):
		""" Checks which regions of the sequence are available for binding based on the size"""
		available_starts = []
		if len(self.transcript) < size:
			return available_starts
		for i in range(len(self.bindings)):
			end = i + size
			full = self.bindings[i:end]
			if len(full) < size:
				continue
			if self.bindings[i] != 0:
				continue
			elif all(x == 0 for x in full):
				available_starts.append(i)
		return available_starts	
	def extract_region(self, start, end):
		""" Get specific section of the sequence based on the start and end index"""
		interest_region = self.transcript[start: end]
		return interest_region
	def splice(self):
		""" Selects which regions can be spliced at a given time """
		#For c. elegans: (30, 60)
		min_splice_num = 10
		max_splice_num = 20
		#check that sequence has at least one u1 and u5 bound to it
		if (len(self.u1_bound_ids) > 0) and (len(self.u5_bound_ids) > 0):
			#search for u1s
			u1_curr = None
			for u1_index, value in enumerate(self.bindings):
				if value == 0: continue
				if value[0] == u1_curr: continue
				elif isinstance(value, tuple) and (value[0][0:2] == "u1"):
					u1_curr = value[0]
					start_range = u1_index + 5 + min_splice_num
					end_range = u1_index + 5 + max_splice_num
					temp_binds = self.bindings[start_range:end_range]
					#see if there is a u5 within the right range of the u1
					u5_best = None
					u5_best_power = 0
					u5_best_start_index = None
					u5_curr = None
					for u5_index, value in enumerate(temp_binds):
						if value == 0: continue
						if value[0] == u5_curr: continue
						elif isinstance(value, tuple) and (value[0][0:2] == "u5"):
							u5_curr = value[0]
							#calculate 'likelihood' of splicing - prob of u5 binding for now, edit later
							power = self.power(value[1], #Redo using apc code
							if power > u5_best_power:
								u5_best_power = power
								u5_best = value[0]
								u5_best_start_index = u5_index + u1_index + 5 + min_splice_num
					if u5_best != None:
						self.cut(u1_index, u5_best_start_index, u1_curr, u5_best)
	def power(self, binding_prob, length):
		###Use apc code to calculate intron probabilities for this part
	def cut(self, u1_start, u5_start, u1_id, u5_id):
		""" Performs the splicing of the select region """
		cut_sequence = ''.join(self.transcript[u1_start:u5_start+6])
		self.splicing_events.append((cut_sequence, u1_start, u5_start, u1_id, u5_id, self.curr))
		del self.bindings[u1_start:u5_start+6]
		del self.transcript[u1_start:u5_start+6]
	def bind(self, start, end, prob, id, size):
		""" Performs binding """
		if size==5:
			self.u1_bound_ids.append(id)
		elif size==6:
			self.u5_bound_ids.append(id)
		self.bindings[start: end] = [(id, prob)] * (end - start)
	def unbind(self, start, size, id):
		""" Performs unbinding """
		end = start+size
		self.bindings[start: start+size] = [0]*(end - start)
		if size == 5:
			self.u1_bound_ids.remove(id)
		elif size == 6:
			self.u5_bound_ids.remove(id)
	def one_iteration(self, iter_num):
		""" Sequence behavior for each iteration """
		if iter_num//self.resolution == 0: #transcribe one base every t iterations
			self.transcribe()
			self.splice()

class snRNP:
	def __init__(self, pwm, size, number):
		""" The snRNP object superclass 
		- pwm_path: path to position weight matrix
		- pwm: the position weight matrix as a list
		- bind_start: the start index of where the snRNP is bound to, if bound
		- bindtime: how long the snRNP has left to be bound for, at 0 it is not bound
		- size: the size of the snRNP
		- type: u1 or u5?
		- id: unique id of snRNP: u{1 or 5}_{number}
		- prob: probability that it would bind to the specific region of the sequence"""
		self.pwm_path = pwm
		self.pwm = self.read_pwm()
		self.bind_start = None
		self.bindtime = 0
		self.size = size #5 for u1s, 6 for u5s
		self.type = self.which_type()
		self.id = self.type + "_" + str(number)
		self.prob = None
	def which_type(self):
		if self.size == 5:
			return "u1"
		elif self.size == 6:
			return "u5"
	def read_pwm(self):
		""" Takes in path to acceptor pwm file and converts into a list """
		pwm = []
		with open(self.pwm_path, 'r') as f:
			for line in f:
				if line[0] == "%":
					continue
				line = line.strip()
				probs = list(map(float, line.split()))
				pwm.append(probs)
		return pwm
	def bind(self, sequence, start_index = None): 
		""" Perform an instance of randomly binding the snRNP to the sequence """
		if start_index == None:
			start_index = random.randint(-2, sequence.curr-1)
		print(start_index)
		#extract region
		available_starts = sequence.seq_available(self.size)
		if (start_index in available_starts) and (start_index >= 0) and (sequence.curr >= self.size):
			#Calculate probability of binding to that region
			end_index = start_index + self.size
			bind_seq = sequence.extract_region(start_index, end_index)
			mappings = {"A": 0, "C": 1, "G": 2, "T": 3}
			x = 0
			self.prob = 1
			for i in bind_seq:
				index = mappings[i]
				curr_prob = self.pwm[x][index]
				x += 1
				self.prob = self.prob*curr_prob
			if self.prob > 0:
				sequence.bind(start_index, end_index, self.prob, self.id, self.size)
				self.bind_start = start_index
				self.bindtime = 10 #how many iterations it is bound for 
			
	def unbind(self, sequence):
		""" Unbind the snRNP from the sequence """
		sequence.unbind(self.bind_start, self.size, self.id)
		self.bind_start = None
		self.prob = None
		self.bindtime = 0
	def one_iteration(self, sequence):
		""" Perform snRNP behaviour for each iteration of the algorithm """
		if self.bindtime > 0: self.bindtime -= 1
		if self.bindtime == 0:
			unbind(sequence)
			bind(sequence)
	
class u5(snRNP):
	number = 0
	def __init__(self, pwm, size = 6):
		""" u5 object, passes u5 pwm and size of 6 to snRNP superclass """
		super().__init__(pwm, size, u5.number)
		u5.number += 1
		
class u1(snRNP):
	number = 0
	def __init__(self, pwm, size = 5):
		""" u1 object, passes u1 pwm and size of 5 to snRNP superclass """
		super().__init__(pwm, size, u1.number)
		u1.number += 1
		
