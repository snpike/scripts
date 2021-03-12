import os.path
import sys

# RULES:

# From Jill:
# MOSFET syntax
# <name> <drain node> <gate node> <source node> <bulk/substrate/back node> <model name>

# All ne and ne5 should have the bulk/substrate node connected to agnd.   
# The pe should have transistor connections to d1p5V and the pe5 should have a5V, d5V etc.

# From Rick:
# There are four types of FETs used in the design named ne, ne5, pe and pe5.
# The ne is a 1.5V nmos fet; the pe is a 1.5V pmos fet; the ne5 is a 5V nmos
# fet and the pe5 is the 5V pmos fet. We need to make sure that the pe and
# ne fets are not used in a situation where any of the terminal voltages can exceed
# 1.5V. The fet gates are typically driven by other fets, so to make sure the
# fet gates cannot exceed 1.5V requires looking upstream to make sure the
# driving fets are not capable of driving above 1.5V. The source and drain connections
# also need to be checked, as well as the back gate connection. In all cases
# the ne and ne5 back gates need to be connected to a ground, usually agnd.
# The back gates of the pe fets should typically be connected to the 1.5V power
# rail, while the back gates of the pe5 fets should be connected to a 5V power rail.



# Get the path to the file and make sure it exists')
filepath = input('Where\'s the file to check? ').strip()
if '/' not in filepath:
	filepath = './' + filepath
while not os.path.exists(filepath):
	filepath = input('Sorry, I couldn\'t find that file. Try again or type "q" to quit: ').strip()
	if filepath=='q':
		sys.exit(0)


# The output will be placed in the same directory as the input
infile = open(filepath, 'r')
extension = -(filepath[::-1].index('.') + 1)
outfile = open(filepath[:extension] + '_errors.txt', 'w')

# Go through line by line and check each rule. If a rule is broken, fix it.
cell_name = ''
for line in infile:
	splitline = line.split()
	if splitline[1]=='Cell':
		cell_name = splitline[-1]
	if line[:1]=='xm':
		circuit = splitline[0]
		fet = splitline[5]
		back = splitline[4]
		source = splitline[3]
		gate = splitline[2]
		drain = splitline[1]
		if fet=='ne' or fet=='ne5':
			if back != 'agnd':
				outfile.write('Error at cell ' + cell_name + ' circuit ' + circuit + '\n')
		if fet=='pe':
			if back != 'agnd':
				outfile.write('Error at cell ' + cell_name + ' circuit ' + circuit + '\n')



# 	add text at bottom of file
# 	     .model PMOS PMOS
# 	     .model NMOS NMOS
# 	     .END
outfile.write('.model PMOS PMOS\n')
outfile.write('.model NMOS NMOS\n')
outfile.write('.END\n')

infile.close()
outfile.close()
print('Done! Please double check the new file for errors.')