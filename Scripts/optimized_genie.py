import random
import math
import sys
import statistics as stats
import numpy as np
from collections import deque


class sequence:
    def __init__(self, sequence, name, resolution=5):
        self.sequence = sequence
        self.name = name
        self.curr = 0
        self.transcript = []
        self.bindings = []
        self.resolution = resolution
        self.u1_bound_ids = set()  
        self.u5_bound_ids = set()  
        self.splicing_events = []
        self.bound_u1s = {}  
        self.bound_u5s = {}  
        self.len_probabilities = self._load_length_probabilities()
        
    def _load_length_probabilities(self):
        """Load length probabilities once during initialization"""
        probabilities = "models/intron.len"
        data = []
        try:
            with open(probabilities, 'r') as f:
                lines = f.readlines()
            for line in lines[1:]:  # Skip header
                data.append(float(line.strip()))
        except FileNotFoundError:
            print("Could not find {probabilities}")
            sys.exit()
        return data
    
    def transcript_length(self):
        """Return length of transcript"""
        return len(self.transcript)
    
    def transcribe(self):
        """ Transcribe one base pair in the sequence """
        if self.curr < len(self.sequence):
            new_base = self.sequence[self.curr]
            self.transcript.append([new_base, self.curr])
            self.curr += 1
            self.bindings.append(0)
        return self.curr
    
    def seq_available(self, size):
        """ Find available binding regions """
        available_starts = []
        transcript_len = len(self.transcript)
        
        # Early return if transcript is too short
        if transcript_len < (size + 100):
            return available_starts
        
        # Use a sliding window approach starting from the specified index
        for i in range(100, transcript_len - size + 1):
            # Check if binding site is available (all positions are 0)
            if self.bindings[i] == 0:
                region_free = True
                for j in range(i, i + size):
                    if self.bindings[j] != 0:
                        region_free = False
                        break
                if region_free:
                    available_starts.append(i)
        
        return available_starts
    
    def extract_region(self, start, end):
        """ Get specific section of the sequence based on the start and end index"""
        return self.transcript[start:end]
    
    def splice(self):
        """ 
        Selects which regions can be spliced at a given time
        """
        min_splice_num = 40
        max_splice_num = 2000
        
        # Early return if no binding
        if not self.u1_bound_ids or not self.u5_bound_ids:
            return
        
        # Process each bound U1 snRNP based on their bind_start attribute
        for u1_id, u1_snrp in self.bound_u1s.items():
            if u1_snrp.bind_start is None:
                continue
                
            u1_index = u1_snrp.bind_start
            u1_prob = u1_snrp.prob
            
            start_range = u1_index + 5 + min_splice_num
            end_range = min(u1_index + 5 + max_splice_num, len(self.bindings))
            
            potential_u5s = []
            
            # Check each u5 snRNP
            for u5_id, u5_snrp in self.bound_u5s.items():
                if u5_snrp.bind_start is None:
                    continue
                    
                u5_index = u5_snrp.bind_start
                u5_prob = u5_snrp.prob
                
                # Check if it's in the appropriate range from this u1
                if start_range <= u5_index < end_range:
                    potential_u5s.append([u5_id, u5_prob, u1_index, u5_index])
            
            if potential_u5s:
                u5_best_start_index, u5_best = self.best_splice(potential_u5s)
                self.cut(u1_index, u5_best_start_index, u1_id, u5_best)
    
    def best_splice(self, potential_u5s):
        """Optimized selection of most likely u5 splice site"""
        splicing_weights = np.zeros(len(potential_u5s))
        
        for i, value in enumerate(potential_u5s):
            length = (value[3] + 6) - value[2]
            acc_prob = value[1]
            splicing_weight = self.power(length, acc_prob)
            splicing_weights[i] = splicing_weight
            
        # Normalize weights
        total = np.sum(splicing_weights)
        if total <= 0:
            # Fallback if all weights are zero or negative
            result_index = 0
        else:
            splicing_weights_norm = splicing_weights / total
            result_index = np.random.choice(len(potential_u5s), p=splicing_weights_norm)
            
        u5_best_start_index = potential_u5s[int(result_index)][3]
        u5_best = potential_u5s[int(result_index)][0]
        
        return u5_best_start_index, u5_best
    
    def power(self, length, acc_prob):
        """Returns the probability using cached length data"""
        if length <= 0 or length > len(self.len_probabilities):
            len_prob = -10  # Default value for out-of-range lengths
        else:
            len_prob = math.log(self.len_probabilities[length - 1])
            
        acc_prob = math.log(acc_prob)
        return len_prob + acc_prob
    
    def cut(self, u1_start, u5_start, u1_id, u5_id):
        """ Performs the splicing of the select region """
        cut_items = self.transcript[u1_start:u5_start+6]
        begin = cut_items[0][1]
        end = cut_items[-1][1]
        cut_sequence = ''.join(x[0] for x in cut_items)
        
        self.splicing_events.append((cut_sequence, begin, end, u1_id, u5_id, self.curr))
        del self.bindings[u1_start:u5_start+6]
        del self.transcript[u1_start:u5_start+6]
        
        # Update bound snRNPs that are affected by the splice
        for snrp_id, snrp in list(self.bound_u1s.items()):
            if snrp.bind_start is not None and u1_start <= snrp.bind_start <= u5_start + 6:
                # The snRNP binding site was deleted during splicing
                snrp.bind_start = None
                snrp.bindtime = 0
                snrp.prob = None
                self.u1_bound_ids.remove(snrp_id)
                
        for snrp_id, snrp in list(self.bound_u5s.items()):
            if snrp.bind_start is not None and u1_start <= snrp.bind_start <= u5_start + 6:
                # The snRNP binding site was deleted during splicing
                snrp.bind_start = None
                snrp.bindtime = 0
                snrp.prob = None
                self.u5_bound_ids.remove(snrp_id)
                
        # Update bind_start for snRNPs after the splice point
        for snrp_id, snrp in list(self.bound_u1s.items()):
            if snrp.bind_start is not None and snrp.bind_start > u5_start + 6:
                # Adjust the position to account for deleted region
                snrp.bind_start -= (u5_start + 6 - u1_start)
                
        for snrp_id, snrp in list(self.bound_u5s.items()):
            if snrp.bind_start is not None and snrp.bind_start > u5_start + 6:
                # Adjust the position to account for deleted region
                snrp.bind_start -= (u5_start + 6 - u1_start)
    
    def bind(self, start, end, prob, id, size, snrp_obj=None):
        """ 
        Performs binding and records the snRNP object 
        """
        if size == 5:
            self.u1_bound_ids.add(id)
            if snrp_obj:
                self.bound_u1s[id] = snrp_obj
        elif size == 6:
            self.u5_bound_ids.add(id)
            if snrp_obj:
                self.bound_u5s[id] = snrp_obj
                
        self.bindings[start:end] = [(id, prob)] * (end - start)
    
    def unbind(self, start, size, id):
        """ Performs unbinding """
        end = start + size
        self.bindings[start:start+size] = [0] * (end - start)
        
        if size == 5:
            self.u1_bound_ids.remove(id)
            if id in self.bound_u1s:
                del self.bound_u1s[id]
        elif size == 6:
            self.u5_bound_ids.remove(id)
            if id in self.bound_u5s:
                del self.bound_u5s[id]
    
    def one_iteration(self, iter_num):
        """ Sequence behavior for each iteration """
        self.transcribe()
        if iter_num % self.resolution == 0:
            self.splice()


class snRNP:
    def __init__(self, pwm, size, number):
        self.pwm_path = pwm
        self.pwm = self.read_pwm()
        self.bind_start = None
        self.bindtime = 0
        self.size = size
        self.type = "u1" if size == 5 else "u5"
        self.id = f"{self.type}_{number}"
        self.prob = None
        
        # Precompute nucleotide mappings
        self.mappings = {"A": 0, "C": 1, "G": 2, "T": 3}
    
    def read_pwm(self):
        """ Takes in path to acceptor pwm file and converts into a list """
        pwm = []
        try:
            with open(self.pwm_path, 'r') as f:
                for line in f:
                    if line[0] == "%":
                        continue
                    line = line.strip()
                    probs = list(map(float, line.split()))
                    pwm.append(probs)
            return pwm
        except FileNotFoundError:
            print(f"Warning: Could not find {self.pwm_path}")
            return [[0.25, 0.25, 0.25, 0.25]] * self.size  # Default fallback
    
    def bind(self, sequence, start_index=None): 
        """ Perform an instance of binding the snRNP to the sequence """
        available_starts = sequence.seq_available(self.size)
        if not available_starts:
            return
            
        if start_index is None or start_index not in available_starts:
            if not available_starts:
                return
            start_index = random.choice(available_starts)
        
        # Calculate binding probability
        end_index = start_index + self.size
        bind_seq = sequence.extract_region(start_index, end_index)
        
        # Optimize probability calculation
        self.prob = 1.0
        for x, base_info in enumerate(bind_seq):
            base = base_info[0]
            if base not in self.mappings:
                continue  # Skip invalid bases
            index = self.mappings[base]
            self.prob *= self.pwm[x][index]
        
        if self.prob > 0:
            # Pass self as the snRNP object to be stored
            sequence.bind(start_index, end_index, self.prob, self.id, self.size, self)
            self.bind_start = start_index
            self.bindtime = 50
    
    def unbind(self, sequence):
        """ Unbind the snRNP from the sequence """
        if self.bind_start is not None:
            sequence.unbind(self.bind_start, self.size, self.id)
            self.bind_start = None
            self.prob = None
            self.bindtime = 0
    
    def one_iteration(self, sequence):
        """ Perform snRNP behaviour for each iteration of the algorithm """
        if self.bindtime > 0:
            self.bindtime -= 1
        if self.bindtime == 0:
            self.unbind(sequence)
            self.bind(sequence)


class u5(snRNP):
    number = 0
    def __init__(self, pwm, size=6):
        super().__init__(pwm, size, u5.number)
        u5.number += 1
        

class u1(snRNP):
    number = 0
    def __init__(self, pwm, size=5):
        super().__init__(pwm, size, u1.number)
        u1.number += 1
