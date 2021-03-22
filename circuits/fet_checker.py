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

class node:
	def __init__(self, owner=None, links=[], flavor=None):
		self.owner = owner
		self.links = links
		self.flavor = flavor

	def add_link(self, node):
		if node is not None:
			self.links.append(node)

	def get_owner(self):
		return self.owner

	def find_power(self, link):
		i = 0
		root_power = None
		while root_power is None:
			if self.links[i] != link:
				root = self.links[i].find_power(link=self)
		return root_power

class rail(node):
	def __init__(self, owner=None, links=[], flavor):
		node.__init__(self, owner, links, flavor)

	def find_power(self, link):
		return self.flavor

class ground(node):
	def __init__(self, owner=None, links=[], flavor):
		node.__init__(self, owner, links, flavor)

	def find_power(self, link):
		return None

# class d1p5v(rail):
# 	def __init__(self, owner=None, links=[]):
# 		rail.__init__(self, owner, links, flavor='d1p5v')

# class d5v(rail):
# 	def __init__(self, owner=None, links=[]):
# 		rail.__init__(self, owner, links, flavor='d5v')

# class ground(node):
# 	def __init__(self, owner=None, links=[], flavor):
# 		node.__init__(self, owner, links, flavor)

class mosfet:
	def __init__(self, subckt=None, name=None, model=None):
		self.subckt = subckt
		self.name = name
		self.model = model
		
		self.drain = None
		self.gate = None
		self.source = None
		self.bulk = None

	def set_drain(self, node):
		drain = node
		if self.source is not None:
			self.drain.add_link(self.source)
			self.source.add_link(self.drain)

	def set_gate(self, node):
		gate = node

	def set_source(self, node):
		source = node
		if self.drain is not None:
			self.drain.add_link(self.source)
			self.source.add_link(self.drain)

	def set_bulk(self, node):
		bulk = node

	def get_gate(self):
		return self.gate

	def get_nodes(self):
		return([self.drain, self.gate, self.source, self.bulk])

class subckt:
	def __init__(self, name=None, nodes = {}, mosfets = {}):
		self.name=name
		self.nodes=nodes
		self.mosfets=mosfets

	def add_node(self, node, node_name):
		self.nodes[node_name] = node

	def add_mosfet(self, mosfet):
		self.mosfets{mosfet.name} = mosfet
		return mosfet

	def get_nodes(self):
		return self.nodes

if __name__=='__main__':

	cell_dict = {}

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
			cell = subckt(cell_name)
		if splitline[0]=='.subckt':
			for x in splitline[2:]:
				if x in ['agnd', 'dgnd']:
					cell.add_node(ground(owner=cell, flavor=x), node_name=x)
				if x in ['d1p5v', 'd5v']:
					cell.add_node(rail(owner=cell, flavor=x), node_name=x)
				else:
					cell.add_node(node(owner=cell), node_name=x)
		if line[:1]=='xm':
			new_fet = cell.add_mosfet(mosfet, subckt=cell, name=splitline[0], model=splitline[5])
			for i in range(splitline[1:5]):
				if splitline[1:5][i] not in cell.get_nodes():
					cell.add_node(node(owner=cell), splitline[1:5][i])
				if i ==0:
					new_fet.set_drain(cell.get_nodes()[splitline[1:5][i]])
				if i ==1:
					new_fet.set_gate(cell.get_nodes()[splitline[1:5][i]])
				if i ==2:
					new_fet.set_source(cell.get_nodes()[splitline[1:5][i]])
				if i ==3:
					new_fet.set_bulk(cell.get_nodes()[splitline[1:5][i]])
		if line[:1]=='xi':
			if splitline[-1] in cell_dict:
				


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