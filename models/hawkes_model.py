import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#sys.path.append('..')
#from problems.problem_definition import *

#A directed, weighted graph
class DWGraph:
def __init__(self):
        self.graph = {}

    def add_vertex(self, vertex):
        """Adds a new vertex to the graph if it doesn't exist."""
        if vertex not in self.graph:
            self.graph[vertex] = {}

    def add_edge(self, u, v, weight):
        """Adds a weighted edge between vertex u and vertex v."""
        # Ensure both vertices exist in our graph
        self.add_vertex(u)
        self.add_vertex(v)
        
        # Add the edge from u to v with the designated weight
        self.graph[u][v] = weight

    def get_neighbors(self, vertex):
        """Returns a dictionary of neighbors and their edge weights."""
        return self.graph.get(vertex, {})

    def get_weight(self, u, v):
        """Returns the weight of the edge between u and v, or None if it doesn't exist."""
        return self.graph.get(u, {}).get(v, None)

    def display(self):
        """Prints the adjacency list representation of the graph."""
        for vertex in self.graph:
            connections = ", ".join([f"{neighbor}(w:{weight})" for neighbor, weight in self.graph[vertex].items()])
            print(f"{vertex} -> {connections}")

class HGEM_constantkernel:
	def __init__(self, node_list, edge_list, delta_t):
		#check to make sure the graph is fully connected
		
		#check to make sure there are no redundant edges
		
		#check to make sure mu>=0, beta>=0, g
		
		#create an internal graph with named nodes and beta edges
		
		#create dict for mu,gamma accessible by node name
		self.mu
		self.gamma
		
		#create dict for time-series of N,lambda accessible by node name
		self.N
		self.lambda_rate
		
		#define delta_t is the shortest timescale for the problem
		self.delta_t
		
		#track total time elapsed as an internal state
		self.t = 0
	
	#return a new HGEM object that only includes these nodes, 
	#and the edges between these nodes
	def subset(self, node_list):
	
	#simulate a short amount of time
	#maintain state of N counts so that the simulation can be picked up again
	def simulate_step(self):
		#elapse time
		self.t += self.delta_t
	
		#return all lambda
		#return all N
	
	#call simulate_step several times
	def simulate(self, time_elapse):
	
	def get_rates(self):
		#return all lambda
		
	def get_counts(self):
		#return all N
	
	#set all N counts to 0, set time to 0
	def reset(self):
		self.t = 0
	
	
	
	
if __name__ == "__main__":
    #unit testing
	
	node_list = [ #name, mu, gamma
		["mde", 1, 0],
		["overtemp", 0, 0],
		["overvolt", 0, 0],
		["stall", 0, 0],
	]
		
	edge_list = [ #"source_node", "target_node", beta
		["mde","mde",0]
		["overtemp","mde",0]
		["overvolt","mde",0]
		["stall","mde",0]
	]
	
	hgem_test = HGEM_constantkernel(node_list,edge_list)