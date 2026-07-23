import sys
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt
import copy
from pyvis.network import Network
import inspect
import types

#sys.path.append('..')
#from problems.problem_definition import *

"""
#A directed, weighted graph
class DWGraph:
def __init__(self):
		self.graph = {}

	def add_node(self, node):
		if node not in self.graph:
			self.graph[node] = {}

	def add_edge(self, u, v, weight):
		# Ensure both vertices exist in our graph
		self.add_node(u)
		self.add_node(v)
		
		# Add the edge from u to v with the designated weight
		self.graph[u][v] = weight

	def get_neighbors(self, node):
		return self.graph.get(node, {})

	def get_weight(self, u, v):
		return self.graph.get(u, {}).get(v, None)

	def display(self):
		for node in self.graph:
			connections = ", ".join([f"{neighbor}(w:{weight})" for neighbor, weight in self.graph[node].items()])
			print(f"{node} -> {connections}")
"""

class DependencyGraph:
	def __init__(self):
		self.dg = nx.DiGraph() #private
		#self.allowable_categories = ["EVENT","VAR"]
	
	def add_node(self, node):
		self.dg.add_node(node, category=_category, subsystem=_subsystem)

	def add_dir_edge(self, u, v):
		self.dg.add_edge(u,v)
		
	def add_nodes_inputs(self, node, input_list):
		for in_node in input_list:
			self.add_dir_edge(in_node, node)

	def get_inputs(self, node):
		print("Inputs to",node,":",list(self.dg.predecessors(node)))
		return self.dg.predecessors(node)
		
	def get_dependencies(self, node):
		print("Dependencies of",node,":",list(self.dg.successors(node)))
		return self.dg.successors(node)

	def printout(self):
		print(f"NODES")
		for node in self.dg.nodes(data=False):
			print(node) 
		print(f"EDGES")
		for edge in self.dg.edges:
			print(edge[0],"->", edge[1]) 
			
	def print_subsystems(self):
		#print(self.dg.nodes(data=True))
		each_subsystem = [node_data["subsystem"] for node,node_data in self.dg.nodes(data=True)]
		print(list(set(each_subsystem)))
		
	def display(self, layout_type):		
		#set node colors and shapes
		node_list = []
		for node, node_data in self.dg.nodes(data=True):
			node_list.append(node)
		
		#draw
		if layout_type == "circle":
			pos = nx.circular_layout(self.dg)
		if layout_type == "spring":
			pos = nx.spring_layout(self.dg)
		nx.draw_networkx_nodes(self.dg, pos,
			nodelist=node_list,
			node_color="lightblue"
			)
		nx.draw_networkx_edges(self.dg, pos, edge_color="gray", arrowsize=10)
		nx.draw_networkx_labels(self.dg, pos, font_size=5)
		plt.show()


class SystemDependencyGraph:
	##############################################
	###Define the graph structure
	##############################################
	def __init__(self, node_input_dict, timestep=1):
		self.__dg = nx.DiGraph() #private
		self.__dt = timestep
		self.__steps = -1
		self.__special_args = ['t','dt','prev_val']
		for node, input_list in node_input_dict.items():
			self.__add_nodes_inputs(node, input_list)
	
	def __add_node(self, node):
		#parse category
		if "[" in node and "]" in node and not ("(" in node and ")" in node):
			_category = "EVENT" 
		elif "(" in node and ")" in node:
			_category = "VAR" 
		else:
			_category = "VAR" 
			#print("error parsing new node",node)
			#sys.exit()
			
		if node in self.__special_args:
			print("Error: Node named",node,"not allowed.")
			sys.exit()
		
		#parse subsystem
		if node.find(":") != -1:
			_subsystem = node.split(":")[0]
		else:
			_subsystem = ""
			#print("error parsing new node",node)
			#sys.exit()
		
		self.__dg.add_node(node, category=_category, subsystem=_subsystem, fn=None, init=0, val_timeseries=[])

	def __add_dir_edge(self, u, v):
		self.__add_node(u)
		self.__add_node(v)
		self.__dg.add_edge(u,v)
		
	def __add_nodes_inputs(self, node, input_list):
		for in_node in input_list:
			self.__add_dir_edge(in_node, node)
			
	##############################################
	###Define the graph content
	##############################################	
	def is_graph_fully_specified(self, doPrint=True, throwError=False):
		fully_specified = True
		specs = self.get_specifications(doPrint=False)
		unspecified = []

		for node,spec in specs:
			if spec == None:
				fully_specified = False
				unspecified.append(node)
			else:
				continue
		
		if fully_specified:
			return True
		else:		
			if throwError or doPrint:
				print("Error: unspecified nodes")
				for node in unspecified:
					print(node,"is unspecified")
			if throwError:
				sys.exit()
			else:
				return False
	
	#TODO implement a way for constant values to be specified; this will make PP and GP easier to implement
	def specify_node(self, node, fn, arg_edge_mapping_dict=None, get_history=[]):	
		if not arg_edge_mapping_dict:
			arg_edge_mapping_dict = {}
		if node in self.__dg:
			self.__dg.nodes[node]["fn"] = fn
			if get_history:
				self.__dg.nodes[node]["get_history"] = get_history
			fn_args = inspect.signature(fn).parameters.keys() #ordered list of fn parameters
			
			defined_args = set(arg_edge_mapping_dict.keys())
			defined_in_nodes = sorted(arg_edge_mapping_dict.values())
			
			#sanity check: don't duplicate things in the spec, why would you do that?
			#if len(defined_args) != len(set(defined_args)):
			#	print("Error: the specified arguments",defined_args,"contain repeated arguments")
			#if len(defined_in_nodes) != len(set(defined_in_nodes)):
			#	print("Error: the specified in nodes",defined_in_nodes,"contain repeated node names")
			
			#sanity check: make sure defined arguments are valid
			if not defined_args.issubset(set(fn_args)):
				print("Error: the specified arguments",defined_args,"includes arguments that aren't among the arguments",set(fn_args)," for",fn,"on node",node)
				sys.exit()
				
			#sanity check: make sure all edges (and only the edges) connect to an argument
			in_nodes = sorted(self.get_inputs(node))
			if not defined_in_nodes == in_nodes:
				print("Error: the input nodes",defined_in_nodes,"specified for",fn,"on node",node,"disagree with the graph structure:",in_nodes)
				sys.exit()
				
			#for all nodes, allow for special arguments that don't refer to a node
			for special_arg in self.__special_args:
				if special_arg in fn_args:
					#supply the current epoch
					arg_edge_mapping_dict[special_arg] = special_arg
			
			for i,arg in enumerate(fn_args):
				#enumerating through fn_args because those are guaranteed to be in correct order
				#for each node, add date defining that argument i of its function corresponds to a particular incident node
				self.__dg.nodes[node]["arg"+str(i+1)+"_node"] = arg_edge_mapping_dict[arg]
		else:
			print("Error: node",node,"not recognized")
			sys.exit()
			
	def set_initial_values(self, node_val_dict):
		for node, init in node_val_dict.items():
			self.__dg.nodes[node]["init"] = init

	##############################################
	###Execute the graph
	##############################################	
	def init_simulation(self):
		self.is_graph_fully_specified(doPrint=False, throwError=True)
			
		#TODO catch other errors
		
		#reset time, set all node val_timeseries to initial values
		self.__steps = 0
		for node in self.__dg.nodes():
			initial_value = self.__dg.nodes[node]["init"]
			self.__dg.nodes[node]["val_timeseries"] = [initial_value]
		
	def simulation_step(self):
		if self.__steps == -1:
			print("Error: simulation has not been initialized")
			sys.exit()
			
		times = [self.__dt*t for t in range(self.__steps+1)]
			
		for node in self.__dg.nodes():
			#identify the function inputs
			fn = self.__dg.nodes[node]["fn"]
			fn_arg_names = inspect.signature(fn).parameters.keys()
			arguments = []
						
			for i,arg in enumerate(fn_arg_names):
				#Handle special arguments manually here
				if arg == 't':
					arguments.append(self.__dt*self.__steps) #supply the current epoch
				elif arg == 'dt':
					arguments.append(self.__dt) #supply the timestep
				elif arg == 'prev_val':
					arguments.append(self.__dg.nodes[node]["val_timeseries"][self.__steps]) #supply the value returned by this function at t-1
				else:
					incident_node_i = self.__dg.nodes[node]["arg"+str(i+1)+"_node"]
					if incident_node_i in self.__dg.nodes[node].get("get_history",[]): #return [] if nothing is defined for get_history
						#Handle cases when you want history from a node
						val_i = self.__dg.nodes[incident_node_i]["val_timeseries"][:self.__steps+1] #from 0:t
					else:
						#Handle default case
						val_i = self.__dg.nodes[incident_node_i]["val_timeseries"][self.__steps] #at time t
					arguments.append(val_i)
					
			#print(arguments) #debug
				
			#call the node's function with those inputs as arguments
			new_val = fn(*arguments)
			
			#append the new value to the timeseries
			self.__dg.nodes[node]["val_timeseries"].append(new_val)
		
		self.__steps += 1
		time = self.__dt * self.__steps
		
		#prepare for print:
		values = [0]*len(self.__dg.nodes())
		for i,(node,value) in enumerate(self.get_node_vals()):
			if self.__dg.nodes[node]["category"]=="EVENT":
				values[i] = [node, sum(self.__dg.nodes[node]["val_timeseries"])]
			else:
				values[i] = [node,value]
		return time, values

	##############################################
	###Access the graph
	##############################################
	def get_inputs(self, node, doPrint=False):
		if doPrint:
			print("Inputs to",node,":",list(self.__dg.predecessors(node)))
		return list(self.__dg.predecessors(node))
		
	def get_dependencies(self, node, doPrint=False):
		if doPrint:
			print("Dependencies of",node,":",list(self.__dg.successors(node)))
		return list(self.__dg.successors(node))
	
	def get_specifications(self, doPrint=True):
		specs = []
		for node in self.__dg.nodes():
			spec = self.get_specification(node, doPrint=doPrint)
			specs.append([node, spec])
			
		return specs
	
	def get_specification(self, node, doPrint=False):
		if self.__dg.nodes[node]["fn"] == None:
			if doPrint:
				print(node, "=== no function specified")
		else:
			if doPrint:
				print(node, "===", self.__dg.nodes[node]["fn"])
				fn_args = inspect.signature(self.__dg.nodes[node]["fn"]).parameters.keys()
				for i,arg in enumerate(fn_args):
					if arg not in self.__special_args:
						print('\tArgument',arg,'from node',self.__dg.nodes[node]["arg"+str(i+1)+"_node"])
		
		return self.__dg.nodes[node]["fn"]
		
	def get_node_vals(self, doPrint=False):
		data = []
		for node in self.__dg.nodes():
			val = self.__dg.nodes[node]["val_timeseries"][-1]
			if doPrint:
				print(node,val)
			data.append([node,val])
		return data
		
	def generate_specification_template(self):
		for node in self.__dg.nodes():
			in_nodes = self.get_inputs(node)
			print('DG.specify_node("'+str(node)+'", TODO', end="")
			if in_nodes:
				print(", {")
				for in_node in in_nodes:
					print('\tTODO : "'+str(in_node)+'",')
				print("})")
			else:
				print(")")
			print("")
	
	##############################################
	###Visualize the graph
	##############################################
	def printout(self):
		print(f"NODES")
		for node in self.__dg.nodes(data=False):
			print(node) 
		print(f"EDGES")
		for edge in self.__dg.edges:
			print(edge[0],"->", edge[1]) 
			
	def print_subsystems(self):
		#print(self.__dg.nodes(data=True))
		each_subsystem = [node_data["subsystem"] for node,node_data in self.__dg.nodes(data=True)]
		print(list(set(each_subsystem)))
		
	def display_simple(self, bipartite_nodes=None):		
		#set node colors and shapes
		node_color = []
		node_shape = []
		node_list = []
		for node, node_data in self.__dg.nodes(data=True):
			color = '#ADD8E6' if node_data["category"]=="VAR" else '#FF6F61'
			#shape = "o" if node_data["category"]=="VAR" else "s"
			node_list.append(node)
			node_color.append(color)
			#node_shape.append(shape)
		
		#draw
		is_planar, _ = nx.check_planarity(self.__dg)
		if bipartite_nodes != None:
			pos = nx.bipartite_layout(self.__dg, bipartite_nodes)
		elif is_planar:
			pos = nx.planar_layout(self.__dg)
		else:
			pos = nx.circular_layout(self.__dg)
		nx.draw_networkx_nodes(self.__dg, pos,
			nodelist=node_list,
			node_color=node_color,
			)
		nx.draw_networkx_edges(self.__dg, pos, edge_color="gray", arrowsize=20)
		nx.draw_networkx_labels(self.__dg, pos, font_size=10)
		plt.show()

	def display(self, layout_type):
		#first, make a copy of the network that adds hidden edges among all shared subsystems
		#inefficient, thats ok
		draw_dg = copy.deepcopy(self.__dg)
		hidden_edges = []
		for node,node_data in self.__dg.nodes(data=True):
			for vnode,vnode_data in self.__dg.nodes(data=True):
				#i.e. create a hidden suprious edge between two disconnected different nodes with the same subsystem
				if node_data["subsystem"] == vnode_data["subsystem"] and node != vnode and not self.__dg.has_edge(node,vnode) and not self.__dg.has_edge(vnode,node):
					draw_dg.add_edge(node,vnode)
					draw_dg.add_edge(vnode,node)
					hidden_edges.append((node,vnode))
		
		#set edge colors, hiding the hidden ones
		edge_color = []
		for u, v in draw_dg.edges():
			if (u, v) in hidden_edges or (v, u) in hidden_edges:
				edge_color.append('none')  # 'none' keyword makes it invisible
			else:
				edge_color.append("gray")  # Visible color
		
		#set node colors and shapes
		node_color = []
		node_shape = []
		node_list = []
		for node, node_data in draw_dg.nodes(data=True):
			color = '#ADD8E6' if node_data["category"]=="VAR" else '#FF6F61'
			shape = "o" if node_data["category"]=="VAR" else "s"
			node_list.append(node)
			node_color.append(color)
			node_shape.append(shape)
		
		#draw
		if layout_type == "spring":
			pos = nx.spring_layout(self.__dg)
		else:
			pos = nx.kamada_kawai_layout(draw_dg)
		nx.draw_networkx_nodes(draw_dg, pos,
			nodelist=node_list,
			node_color=node_color,
			)
		nx.draw_networkx_edges(draw_dg, pos, edge_color=edge_color, arrowsize=10)
		nx.draw_networkx_labels(draw_dg, pos, font_size=5)
		plt.show()
		
	def display_interactive(self):
		net = Network(height="800px", width="100%", bgcolor="#FFFFFF", font_color="black", directed=True)
		#net.set_options("""
		#const options = {
		#  "physics": {
		#	"forceAtlas2Based": {
		#	  "springLength": 100,
		#	  "springConstant": 0.375,
		#	  "damping": 0.18,
		#	  "avoidOverlap": 0.8
		#	},
		#	"minVelocity": 0.75,
		#	"solver": "forceAtlas2Based"
		#  }
		#}
		#""")
		net.options.edges.smooth = False
		
		draw_dg = copy.deepcopy(self.__dg)
		# Sanitize function attributes
		for node, data in draw_dg.nodes(data=True):
			for key, value in list(data.items()):
				# If the attribute is a function or lambda, convert it to a string
				if isinstance(value, (types.FunctionType, types.BuiltinFunctionType)):
					data[key] = str(value.__name__)
		
		for node,node_data in draw_dg.nodes(data=True):
			#node_data['title'] = f"Hover info for Node {node}"
			node_data['color'] = '#ADD8E6' if node_data["category"]=="VAR" else '#FF6F61'
		
		net.from_nx(draw_dg)
		for edge in net.get_edges():
			edge["color"] = '#808080'

		# 4. Turn on the interactive physics configuration UI panel
		net.show_buttons(filter_=['physics'])
		
		# 5. Save and open the file in your browser
		net.write_html("network.html", open_browser=True)
	
	
	
	
if __name__ == "__main__":
	#unit testing
	debug_init = {
		"B" :["A"],
		"C" : ["A","B"],
	}
	DG_debug = SystemDependencyGraph(debug_init, timestep=0.01)
	DG_debug.get_inputs("A", doPrint=True)
	DG_debug.get_inputs("B", doPrint=True)
	DG_debug.get_inputs("C", doPrint=True)
	
	#DG_debug.display_interactive()
	#DG_debug.get_specifications()
	#DG_debug.is_graph_fully_specified()
	
	def A_fn(t):
		return t
	
	def B_fn(a_in, t):
		return a_in*2
	
	def C_fn(b_in, a_in, t):
		return a_in + b_in + t
	
	DG_debug.specify_node("A", A_fn, {})
	DG_debug.specify_node("B", B_fn, {"a_in":"A"})
	DG_debug.specify_node("C", C_fn, {"a_in":"A","b_in":"B"})
	DG_debug.set_initial_values({
		"A":0,
		"B":0,
		"C":0
	})
	
	DG_debug.get_specifications()
	
	#run sim
	print("t=0")
	DG_debug.init_simulation()
	DG_debug.get_node_vals(doPrint=True)
	while True:
		time, _ = DG_debug.simulation_step()
		print("t="+str(time))
		DG_debug.get_node_vals(doPrint=True)
		input()