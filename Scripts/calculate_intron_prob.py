import matplotlib.pyplot as plt
import statistics as stats
import numpy as np

from scipy.stats import gaussian_kde, lognorm

file_name = "/home/rkdesai7/Documents/Monte_Genie/Data/1pct_elegans.gff3"
class intron_probability():
	def __init__(self, gff_path, intron_length):
		"""Returns the probability that a intron is of the length intron_length"""
		self.file_path = gff_path
		self.intron_length = intron_length
		self.intron_lengths = self.extract_lengths()
	def extract_lengths(self):
		""" Gather the lengths of all the introns from a .gff file """
		intron_lengths = []
		val_starts = ['X', 'V', 'IV', 'III', 'II', 'I']
		total_lines = 0
		with open(self.file_path, 'r') as file:
			for line in file:
				total_lines += 1
				data = line.split()
				if data[0] not in val_starts: continue
				if data[2] == "intron":
					intron_len = int(data[4]) - int(data[3])
					if intron_len > 200: continue
					intron_lengths.append(intron_len)
		return intron_lengths
	def calculate_prob(self, method):
		"""Returns probability based method specified. Current options: Gaussian KDE, Log Normal, Frequency Based"""
		if method == "Gaussian KDE": return self.gaussian()
		if method == "Log Normal": return self.log_normal()
		if method == "Frequency Based": return self.frequency_based()
	def gaussian(self):
		"""Calculate probability using the gaussian kde method from scipy"""
		kde = gaussian_kde(self.intron_lengths, bw_method ='scott')
		prob = kde.evaluate([self.intron_length])
		return prob
	def log_normal(self):
		"""Calculates probability using the log normal distribution"""
		mean, sd = self.calc_summary_stats()
		distribution = lognorm(sd, scale=np.exp(mean))
		prob = distribution.pdf(self.intron_length)
		return prob
	def frequency_based(self):
		freq = 0
		for i in self.intron_lengths:
			if i == self.intron_length: freq += 1
		freq = freq/len(self.intron_lengths)
		return freq
	def calc_summary_stats(self):
		"""Returns mean and standard deviation of intron length data"""
		mean = sum(self.intron_lengths)/len(self.intron_lengths)
		sd = stats.stdev(self.intron_lengths)
		return mean, sd

	
intron = intron_probability(file_name, 300)
intron1 = intron_probability(file_name, 60)
print(intron.calculate_prob("Gaussian KDE"), intron1.calculate_prob("Gaussian KDE"))

