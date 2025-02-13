from genie import sequence, snRNP, u5, u1

#Inputs for test purposes:TTTCAG
test_sequence = "AGCTGTAAGTAAGTGAAGATATTTCAGTCTCGCGCTAGGGCGTAAGTCTAAAGCTGAGGTCAAAAG"
u1_pwm_path = "models/don.pwm"
u5_pwm_path = "models/acc.pwm"

###Look at the sequence properties first:
#First initialize the sequence object and let's look at it's name
my_seq = sequence(test_sequence, "genie")
print("SEQUENCE PROPERTIES")
print(f"My full sequence {my_seq.name}:\n", my_seq.sequence)
#Run 40 rounds of transcription, and check the object properties
for i in range(40):
	my_seq.transcribe()
print("Transcribed Sequence: \n", my_seq.transcript,"\nNumber of bases transcribed:\n", my_seq.curr,"\nBindings: \n", my_seq.bindings)
#See which regions of the sequence are available for binding
print("\nAvailable u1 binding starts:\n", my_seq.seq_available(5))
print("\nAvailable u5 binding starts:\n", my_seq.seq_available(6))

#Let's create a u1 and u5 snRNP
print("\n=======================================")
print("CREATE NEW U1 and U5 snRNP")
my_u1 = u1(u1_pwm_path)
print("U1 PWM:", my_u1.pwm)
my_u5 = u5(u5_pwm_path)
print("U5 PWM:", my_u5.pwm)
print("The respecitce id's are:", my_u1.id, my_u5.id)
#Let's look at the position weight matrix of the u1:
print("\nThe u1 position weight matrix is:\n", my_u1.pwm)

#Let's force the u1 to bind at the 4th position, and the u5 to bind to the 21nd positions:
my_u1.bind(my_seq, 4)
my_u5.bind(my_seq, 21)
print("\n=======================================")
print("BIND snRNPs")
print("\nOur sequence bindings:\n", my_seq.bindings)
print("\nOur u1 snRNP updated info:")
print("Binding Probability:\n", my_u1.prob)
print("Binding start:\n", my_u1.bind_start)
print("Binding Time:\n", my_u1.bindtime)
print("\nOur u5 snRNP updated info:")
print("Binding Probability:\n", my_u5.prob)
print("Binding start:\n", my_u5.bind_start)
print("Binding Time:\n", my_u5.bindtime)

#Let's see if any splicing can happen:
print("\n=======================================")
my_seq.splice()
print("RUN ONE INSTANCE OF SPLICING:")
print("\nNew post-splice sequence:\n", my_seq.transcript)
print("\nSplicing events, if any:\n", my_seq.splicing_events)




