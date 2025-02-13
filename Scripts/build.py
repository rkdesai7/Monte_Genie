data = "../smallgenes/ch.1_403.gff3"
real_data = []
with open(data, "r") as file:
	for line in file:
		if "splice" in line:
				name, prog, type, start, end, freq = map(float, line.split())
				real.append((start, end, freq))
			
print(real_data)

line_counts = {}
pred_path = "outs.txt"
with open(pred_path, "r") as file:
	for line in file:
	#strip whitespace and remove quotes
		line = line.strip().replace('"', '')
		start, end = map(int, line.split())
		line_counts[(start, end)] = line_counts.get((start, end), 0) + 1
pred = [(start, end, count) for (start, end), count in line_counts.items()]

			
print(pred)
