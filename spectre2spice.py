import os.path
import sys

# RULES:

# rm lines starting with global
# rm lines starting with include
# rm all "(" and ")"
# convert strings \+ to +, \* to * , \- to -
# rm c= and r=
# rm text resistor and capacitor
# replace nmos with NMOS
# replace pmos with PMOS
# replace subckt with .SUBCKT
# replace ends with .ENDS
# replace instances starting with I with XI
# unwrap lines that wrap to next line
# add text at bottom of file
#      .model PMOS PMOS
#      .model NMOS NMOS
#      .END

badslash = ['\\+', '\\-', '\\*']
badequals = [' c=', ' r=']
badstrings = ['(', ')', ' resistor', ' capacitor']
capitalize = ['nmos', 'pmos']
capandperiod = ['ends', 'subckt']

# Get the path to the file and make sure it exists')
filepath = input('Where\'s the file to convert? ').strip()
if '/' not in filepath:
	filepath = './' + filepath
while not os.path.exists(filepath):
	filepath = input('Sorry, I couldn\'t find that file. Try again or type "q" to quit: ').strip()
	if filepath=='q':
		sys.exit(0)
if filepath[-7:] != 'spectre':
	response = input('This doesn\'t look like a spectre file. Continue? (y/n) ').strip()
	while response not in ['y','yes','n','no']:
		response = input('y/n? ')
	if response in ['n', 'no']:
		sys.exit(0)



# The output will be placed in the same directory as the input
infile = open(filepath, 'r')
extension = -filepath[::-1].index('.')
outfile = open(filepath[:extension] + 'spc', 'w')

# Go through line by line and check each rule. If a rule is broken, fix it.
lastbreak = False
for line in infile:

	# 	rm lines starting with global
	# 	rm lines starting with include
	if (line[:6]!='global') and (line[:7]!='include'):

		if line == 'simulator lang=spectre\n':
			line = 'simulator lang=hspice\n'

		# Get rid of tabs if we just unwrapped
		if lastbreak:
			line = line.strip() + '\n'

		# 	convert strings \+ to +, \* to * , \- to -
		for x in badslash:
			while x in line:
				for i in range(len(line)):
					if line[i:i+2] == x:
						line = line[:i] + line[i+1:]
						break

		# 	rm c= and r=
		for x in badequals:
			while x in line:
				for i in range(len(line)):
					if line[i:i+3] == x:
						line = line[:i+1] + line[i+3:]
						break

		# 	rm all "(" and ")"
		# 	rm text resistor and capacitor
		for x in badstrings:
			while x in line:
				for i in range(len(line)):
					if line[i:i+len(x)] == x:
						line = line[:i] + line[i+len(x):]
						break

		# 	replace nmos with NMOS
		# 	replace pmos with PMOS
		for x in capitalize:
			while x in line:
				for i in range(len(line)):
					if line[i:i+len(x)] == x:
						line = line[:i] + line[i:i+len(x)].upper() + line[i+len(x):]
						break

		# 	replace subckt with .SUBCKT
		# 	replace ends with .ENDS
		for x in capandperiod:
			while x in line:
				for i in range(len(line)):
					if line[i:i+len(x)] == x:
						line = line[:i] + '.' + line[i:i+len(x)].upper() + line[i+len(x):]
						break

		# 	replace instances starting with I with XI
		if 'I' in line:
			if line.strip()[0]=='I':
				line = line[:line.index('I')] + 'X' + line[line.index('I'):]

		# 	unwrap lines that wrap to next line
		if ('\\\n') in line:
			line = line[:-2]
			lastbreak = True
		else:
			lastbreak = False

		outfile.write(line)

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