#!/usr/bin/env python3

# usage parentPath.py -g blunt.PackedGraph.pg -m mat.node.kmer.count -p pat.node.kmer.count

from bdsg.bdsg import PackedGraph
import argparse 
import numpy as np
import json

def parseArguments():
	# take I/O arguments from command line
	parser=argparse.ArgumentParser()
	parser.add_argument('-g','--input_pg')
	# parser.add_argument('-m','--mat_alns')
	# parser.add_argument('-p','--pat_alns')
	parser.add_argument('-j','--input_json')

	args=parser.parse_args()
	inputGraph=args.input_pg
	inputJson=args.input_json
	return inputGraph, inputJson

class snrlManager:
	'''class that loads and contains the snarl dictionaries '''
	def __init__ (self,inputJson):
		''' class variables:
			 dictionaries that hold the dictionaries of snarl starts and ends'''
		self.inputJson=inputJson
		self.sd={\
			'parentSnarls':{},\
			'snarls':{},\
			'biAllelic_snarls':{},\
			'directed_acyclic_net_graph':{},\
			'parentSnarlsEnds':{},\
			'snarlsEnds':{},\
			'biAllelic_snarlsEnds':{},\
			'directed_acyclic_net_graphEnds':{}}
		self.children=[]
		# self.parentSnarls={}
		# self.snarls={}						# complicated snarls
		# self.biAllelic_snarls={}				# 
		# self.directed_acyclic_net_graph={}   	# tips reported as alts ( make a size requirement?)
		# #ends
		# self.parentSnarlsEnds={}
		# self.snarlsEnds={}						# complicated snarls
		# self.biAllelic_snarlsEnds={}				# 
		# self.directed_acyclic_net_graphEnds={}
		
		
	# input method
	def snarlsInput(self):
		'''load snarls.json into their respective snarl type dictionaries'''
		snarl_json=[json.loads(line) for line in open(self.inputJson,'r')]
		for sn in snarl_json:
		# bi-allelic snarls are type: 1
			if sn.get('type',False):
				if sn.get('type')==1:
					self.sd['biAllelic_snarls'][int(sn.get('start').get('node_id'))] = {\
						'visited' : False, \
						'backward': sn.get('start').get('backward', False),\
						'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward', False)},\
						'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
						'child' : sn.get('parent',False)}
					self.sd['biAllelic_snarlsEnds'][int(sn.get('end').get('node_id'))] = {\
						'visited' : False, \
						'backward': sn.get('end').get('backward', False),\
						'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward',False)},\
						'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
						'backward' : sn.get('end').get('backward', False), \
						'child' : sn.get('parent',False)}
				else:
					print('Snarl Type not 1:',sn.get('type')) 
			elif sn.get('directed_acyclic_net_graph',False):
				self.sd['directed_acyclic_net_graph'][int(sn.get('start').get('node_id'))] = {\
					'visited' : False, \
					'backward': sn.get('start').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward', False)},\
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'child' : sn.get('parent',False)}
				self.sd['directed_acyclic_net_graphEnds'][int(sn.get('end').get('node_id'))] = {\
					'visited' : False, \
					'backward': sn.get('end').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')), 'backward': sn.get('start').get('backward',False)}, \
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'child' : sn.get('parent',False)}
			else:
				self.sd['snarls'][int(sn.get('start').get('node_id'))] = {\
					'visited' : False, \
					'backward': sn.get('start').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward', False)},\
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'child' : sn.get('parent',False)}
				self.sd['snarlsEnds'][int(sn.get('end').get('node_id'))] = {\
					'visited' : False, \
					'backward': sn.get('end').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward',False)},\
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'child' : sn.get('parent',False)}
			if sn.get('parent',False):
				# print('parent',sn.get('parent'),'Parent Start',sn.get('parent').get('start'),sn.get('parent').get('end'))
				if int(sn.get('parent').get('start').get('node_id')) in list(self.sd['parentSnarls'].keys()):
				# 	# print(parentSnarls[int(sn.get('parent').get('start').get('node_id'))]['child'])
					self.sd['parentSnarls'][int(sn.get('parent').get('start').get('node_id'))]['child'].append((int(sn.get('start').get('node_id')),int(sn.get('end').get('node_id')),sn.get('end').get('backward', False)))
				else:
					self.sd['parentSnarls'][int(sn.get('parent').get('start').get('node_id'))] = {\
						'visited' : False, \
						'start': sn.get('parent').get('start').get('node_id'), \
						'end': sn.get('parent').get('end').get('node_id'), \
						'backward' : sn.get('parent').get('end').get('backward', False), \
						'child' : [(int(sn.get('start').get('node_id')),int(sn.get('end').get('node_id')),sn.get('end').get('backward', False))]}
					self.sd['parentSnarlsEnds'][int(sn.get('parent').get('end').get('node_id'))] = {\
						'visited' : False, \
						'start': sn.get('parent').get('start').get('node_id'), \
						'end': sn.get('parent').get('end').get('node_id'), \
						'backward' : sn.get('parent').get('end').get('backward', False), \
						'child' : [(int(sn.get('start').get('node_id')),int(sn.get('end').get('node_id')),sn.get('end').get('backward', False))]}
					self.children.append(int(sn.get('start').get('node_id')))

		for p in list(self.sd['parentSnarls'].keys()):
			if p in self.children:
				print('PARENT IS ALSO CHILD',p,self.sd['parentSnarls'][int(p)])
				self.sd['parentSnarls'].pop(p)
				if int(p) in list(self.sd['snarls'].keys()):
					print('p in snarls')
				if p in list(self.sd['directed_acyclic_net_graph'].keys()):
					print('p in dag')
		print('check again',list(self.sd['parentSnarls'].keys()))

		if len(list(self.sd['snarls'].keys()))>0:
			print('snarls here',list(self.sd['snarls'].keys()))
		return self.sd

class graphPathFinder:
	'''class that loads and contains the snarl dictionaries '''
	def __init__ (self,sd,inputGraph):
		# snarl dictionary and inputGraph File
		self.sd=sd
		self.inputGraph=inputGraph
		# make a PackedGraph object
		self.graph=PackedGraph() 
		# hap 1 list of handles in the path
		self.hap1Path=[]
		self.hap1NodePath=[]
		self.hap1LeftNodePath=[]
		self.hap1RightNodePath = []
		# hap 2 list of handles in the path
		self.NodePath=[]
		# snarlWalk lists to explore complicated areas of the graph
		self.pSnarlWalk=[]
		self.snarlWalk=[]
		self.snarlEnd=-1
		self.snarlStart=-1
		# set max distance we will walk into a Parent snarl w/o turning around for the other way
		self.maxParentSnarlDepth = 20
		# no parent path variables  TODO make this the same as hap1 I think 
		self.path =[]
		self.pathLeft=[]
		self.pathRight=[]
		self.plen=0
		# next parent node
		self.nextParentStartNode = -1
		self.currentParentStartNode = -1
		self.parentSnarlDepth = 0
		self.returnAndFlipNode = False
		self.nextSnarlStartNode=0
		self.parentSnarlStartKeyList=[]
		# list to get through parent snarl
		self.neighbors=[]
		self.tipWalk=[]
		self.startSnarlOver=False
		self.snarlStartNodeDegrees=-1
		self.snarlAttempts=0
		self.noParentSnarlEnd=0
		self.dagOtherEnd=0
		self.dagEnd=0
		self.dagDict={}
		self.dagSide=''
		#stats
		self.biAllelicLengths=[]
		self.tipLengths=[]

	def deserializeGraph(self):
		# deserialize the packed graph
		self.graph.deserialize(inputGraph)
		# Q? make a path in the graph instead of the list of graph_ids? 
		# A: no probably too slow? also I'll be altering the origional gfa file 
		# path = graph.create_path_handle("hap1path")

	# def chainParents(self):
		# order the parent snarls in the forward (left -> right) orientation for forward traversal



def writeOutPathCSV(hap1NodePath):
	# open the path.csv and write out the path with colors for bandage viewing 
	outpathCSV=open(args.input_pg[:-3]+".path.csv",'w')
	outpathCSV.write("Node,Path,Color\n")
	# parent1='hg03'
	parent2='hg04'
	# parent1Color = "light sea green"
	parent2Color = 'tomato'

	for n in hap1NodePath:
		outpathCSV.write( str(n) +",M,"+parent2Color+"\n")

def revComp(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

# TODO make it a list of handles that are orientated 
def writeContig(hap,contigCount,primary):
	# output all the sequence from the path to a fasta file ( 80 characters)
			# Don't forget to check orientations with e.g. graph.get_is_reverse(here)
	charMax = 70
	charCount = 0
	priLength = 0
	if primary:
		# for now I'm opening 'w'. if I loop over all connected components then append 'a'
		file = open(args.input_pg[:-3]+'.pri.fa','w')
	else:
		file = open(args.input_pg[:-3]+'.alt.fa','w')

	file.write('>'+args.input_pg[:-3]+'_'+str(contigCount)+'\n')
	# change to handles instead of graph ids 
	for node in hap:
		seq = graph.get_sequence( graph.get_handle( int(node) ) )
		priLength+=len(seq)
		# if the node is backwards use the rev compliment 
		if graph.get_is_reverse( graph.get_handle(int(node)) ):
			seq = revComp(seq)

		writeList=[]
		writeString=""
		for base in seq:
			# file.write(base)
			writeString=writeString+base
			charCount+=1
			if charCount >= charMax:
				writeString=writeString+'\n'
				# file.write(writeString)
				writeList.append(writeString)
				writeString=""
				charCount=0
		if len(writeString)>0:
			writeString=writeString+'\n'
			writeList.append(writeString)
			writeString=""
		for s in writeList:
			file.write(s)
		writeList.clear()

	file.write('\n')
	file.close()
	return priLength

def getStats(plen,pri_nodes,alt_nodes,dbls_pri,dbls_alt):
	### open the stats.csv and add the header
	global tipLengths
	global biAllelicLengths
	
	outSTAT=open(args.input_pg[:-3]+".stats.csv",'w')
	outSTAT.write("chunkNum,chunk_Len,NumNodes,parent_snarls,num_snarls,covered_snarls,num_dang,num_bi_allelics,hetSeq_Len,bi_allelic_avgLen,bi_allelic_medianLen,total_seq_inTip,PrimaryPath_len,Primary_Nodes,Alt_Nodes,Pri_dups,Alt_dups,NotIncluded_Nodes,NotInc_Nodes/TotalNodes,bps_not_inPri,missingParentSnarls,covered/totParentSnarls,plen/totalLen\n")
	outSTAT.write(args.input_pg[6:][:-3]+",")
	outSTAT.write(str(graph.get_total_length())+",")
	outSTAT.write(str(graph.get_node_count())+",")
	outSTAT.write(str(len(parentSnarls))+",")
	outSTAT.write(str(len(snarls))+","+str(len(parentSnarlStartKeyList))+",")
	outSTAT.write(str(len(directed_acyclic_net_graph))+",")
	outSTAT.write(str(len(biAllelic_snarls))+",")
	if len(biAllelicLengths) > 0:
		outSTAT.write(str(int(np.sum(biAllelicLengths)))+",")
		outSTAT.write(str(int(np.mean(biAllelicLengths)))+",")
		outSTAT.write(str(int(np.median(biAllelicLengths)))+",")
	else:
		outSTAT.write('0,0,0,')
	if len(tipLengths) > 0:
		outSTAT.write(str(int(np.sum(tipLengths)))+",")
	else:
		outSTAT.write('0,')
	outSTAT.write(str(plen)+",")
	outSTAT.write(str(pri_nodes)+",")
	outSTAT.write(str(alt_nodes)+",")
	outSTAT.write(str(dbls_pri)+",")
	outSTAT.write(str(dbls_alt)+",")
	outSTAT.write(str( (graph.get_node_count()-(pri_nodes+alt_nodes) ) )+",")
	outSTAT.write(str( ((pri_nodes+alt_nodes)/graph.get_node_count()) )+",")
	outSTAT.write(str( (graph.get_total_length()-(plen) ) )+",")
	outSTAT.write(str(len(snarls)-len(parentSnarlStartKeyList))+",")
	if len(snarls)>0:
		outSTAT.write(str(len(parentSnarlStartKeyList)/len(snarls))+",")
	else:
		outSTAT.write('0,')
	outSTAT.write(str( (plen)/(graph.get_total_length() ) )+"\n")

	outSTAT.close()
# function to visit each graph handle and determine it's neighbors parent of origin
# visit a handle in the graph


lis=[]
def makelis(y):
	# make a list of all reachable nodes
	global lis
	lis.append(int(graph.get_id( y)))
	return True

def checkAllReachableNodes(other):
	# follow edges in both directions adding nodes to the lis
	lis.clear()	
	graph.follow_edges(other,False,makelis) 
	graph.follow_edges(graph.flip(other),False,makelis)

def checkEndReachable(thisNode,endNode):
	# make a list of all reachable nodes return end node if it is reachable
	lis.clear()
	checkAllReachableNodes(thisNode)
	print(graph.get_id(thisNode),'made lis:', lis, endNode, type(endNode), endNode in lis )
	# if the end of the snarl is in that list - go there 
	if endNode in lis:
		lis.clear()
		return endNode
	else:
		lis.clear()
		return -1

def checkDirection(nodeID,theList):
	## TODO change it to a count in the list?
	# if pSnarlWalk.count(graph.get_id(other)) <3:
	# if we haven't visited this node before add it to the walk
	if nodeID not in theList:
		# print(theList,'\t',nodeID)
		theList.append(nodeID)
	else:
		# if we have gone this way get out and turn around
		print(nodeID,'in list wrong way')
		return True 
	return False

def runFunc(node,dic,fun):
	# if the node is backwards flip it, with an included double check for going in the wrong direction
	print('doing runFunc',graph.get_id(node))
	if dic[graph.get_id(node)]['backward']:
		print('flip node!', graph.get_id(node), str(fun),'degree', graph.get_degree(graph.flip(node),False))
		# Funciton returns False if we find the next parent and all is good
		if graph.follow_edges(graph.flip(node),False,fun):
			# if for some reason in that direction we never found the end, or had already gone that way
			# it returns True and we will flop and go the other way looking for the next parent
			print('rtIntoSnarl rt True, un flipping node')
			graph.follow_edges(node,False,fun)
	else:
		# if it's not backwards then just go to the 'right'
		print('going to the rt, no flip')
		if graph.follow_edges(node,False,fun):
			print('rt no good flipping')
			graph.follow_edges(graph.flip(node),False,fun)
	return False

def runSideFun(ode,theDict,side,fun):
	# determine if snarl is backwards, run it and have an alternative flip if that direction isn't in the snarl
	if theDict[graph.get_id(ode)][side]['backward']:
		print('flippin from side fun',theDict[graph.get_id(ode)])
		checkAllReachableNodes(ode)
		print('all reachable',lis)
		lis.clear()
		if graph.follow_edges(graph.flip(ode),False,fun):
			print('un flipping the side node in fun now')
			graph.follow_edges(ode,False,fun)
			
	else:
		print('running right from side fun',theDict[graph.get_id(ode)])
		checkAllReachableNodes(ode)
		print('all reachable',lis)
		lis.clear()
		if graph.follow_edges(ode,False,fun):
			print('\n that didnt work: flipping the side fun node', ode)
			graph.follow_edges(graph.flip(ode),False,fun)
	return False

def navigateSnarl(snStart, dictForSnarl,dictForSnarlEnd,fun,endFun,thisPath,side):
	# decide how to proceed at snarl edge: 'start' or 'end'
	global snarlEnd
	global snarlStart
	snrlEnd = snarlEnd
	snStart =snarlStart

	def runSnrlFun(snNode,dictSn,side,fun):
		# determine if snarl is backwards, run it and have an alternative flip if that direction isn't in the snarl
		if dictSn[snNode][side]['backward']:
			print('running right')
			if graph.follow_edges(graph.flip(graph.get_handle(snNode)),False,fun):
				print('flipping the snarl node')
				graph.follow_edges(graph.get_handle(snNode),False,fun)
				
		else:
			print('running flip')
			if graph.follow_edges(graph.get_handle(snNode),False,fun):
				print('unflipping the snarl node')
				graph.follow_edges(graph.flip(graph.get_handle(snNode)),False,fun)
		return False

	print('side',side,'snarlEnd',snarlEnd,type(snarlEnd))
	print('path',thisPath, snarlEnd not in thisPath)
	# if we haven't already passed the end of the snarl go through it
	if snrlEnd not in thisPath:
		print('so what is up', fun)
		print('rt into snarl? end:', snrlEnd)
		runSnrlFun(snStart,dictForSnarl,side,fun)
	else: 
		# if we have already passed through the snarl set the current node to the snarl end and keep going
		snrlEnd = snStart # graph.get_id(snStart)

	#go from snrlEnd
	if snrlEnd in list(dictForSnarlEnd.keys()):
		runSnrlFun(snrlEnd,dictForSnarlEnd,side,endFun)
	elif snrlEnd in list(dictForSnarl.keys()):
		runSnrlFun(snrlEnd,dictForSnarl,side,endFun)

def visit_right_intoParentsnarl_neighbor(other):
	# go from the left 'Start' node to the 'end' rt node
	global currentParentStartNode
	global lis
	global pSnarlWalk
	global returnAndFlipNode
	
	def checkSnarlDepth():
		# increment the counter: max snarl depth is 20
		global maxParentSnarlDepth
		global pSnarlWalk

		if len(pSnarlWalk) > maxParentSnarlDepth:
			print('snarl is too deep (>20)', len(pSnarlWalk))
			# pSnarlWalk.clear()
			return True
		else:
			return False

	print('\n',graph.get_id(other),'current:',currentParentStartNode,'parent end',parentSnarls[ currentParentStartNode ]['end'])
	print('psnarlWalk',pSnarlWalk)

	# stop if we are more than 20 nodes in 
	if checkSnarlDepth():
		print('do we return True here?')
		returnAndFlipNode = True
		return True
	# if we have visited this node clear the walk and get out
	if checkDirection(graph.get_id(other),pSnarlWalk):
		# pSnarlWalk.clear()
		returnAndFlipNode =  True
		return True

	# if the node is already in the haploid path, clear the walk and get out
	if graph.get_id(other) in hap1NodePath:
		print('node in path',graph.get_id(other))#, hap1NodePath)
		pSnarlWalk.clear()
		# returnAndFlipNode =  True
		return True
	# if we haven't reached the end of the snarl keep looking
	elif str(graph.get_id(other)) != parentSnarls[ currentParentStartNode ]['end']:
		# make a list of all reachable nodes
		checkAllReachableNodes(other)
		# graph.follow_edges(other,False,makelis) 
		# graph.follow_edges(graph.flip(other),False,makelis)
		print('made lis:', lis,  parentSnarls[ currentParentStartNode ]['end'], int(parentSnarls[ currentParentStartNode ]['end']) in lis )
		# if the end of the snarl is in that list - go there 
		if int(parentSnarls[ currentParentStartNode ]['end']) in lis:
			lis.clear()
			print('choose node go into Parent Snrl')
			print('end walk on',parentSnarls[currentParentStartNode]['end'],pSnarlWalk,'\n')
			# 
			visit_right_intoParentsnarl_neighbor(graph.get_handle(int(parentSnarls[ currentParentStartNode ]['end'])))
			# graph.follow_edges(graph.get_handle(int(parentSnarls[ currentParentStartNode ]['end'])),False,visit_right_intoParentsnarl_neighbor)
		else:
			lis.clear()
			graph.follow_edges(other,False,visit_right_intoParentsnarl_neighbor) 
	# if we reach the end of the snarl stop the walk. 
	elif str(graph.get_id(other)) == parentSnarls[ currentParentStartNode ]['end']:
		print('end walk on',parentSnarls[list(parentSnarls.keys())[0]]['end'],pSnarlWalk,'\n')

	if not (returnAndFlipNode):
		# store the snarl walk in the hapPath
		for n in pSnarlWalk: 
			if n in hap1NodePath:
				print('wrong way?',n)
				pSnarlWalk.clear()
				# returnAndFlipNode =  True
				return True
			else:
				print(n,'append')
				hap1NodePath.append(n)
		pSnarlWalk.clear()
	
	# return True allows for both sides of bubble to be traversed
	#, returnAndFlipNode)
	if len(pSnarlWalk) > 0:
		if graph.get_id(other) == pSnarlWalk[0]:
			if returnAndFlipNode:
				print(graph.get_id(other),'rt True')
				returnAndFlipNode=False
				pSnarlWalk.clear()
				return True
	else:
		print(graph.get_id(other), 'or are we here and end up returning False')
		return False

def visit_handle_checkParent(here):
	global snarls
	global snarlsEnds

	# finds all the paths through the snarl - but 'returns True' so it checks all paths
	if checkDirection(graph.get_id(here),hap1NodePath+hap1LeftNodePath):
		return True

	print('handle',graph.get_id(here), 'backward',parentSnarls[graph.get_id(here)]['backward'])
	# if graph.get_id(here) not in hap1NodePath:
	hap1NodePath.append(graph.get_id(here))
	
	global nextParentStartNode
	def visit_outOfsnarlsEnd_neighbor(other):
		global snarlEnd
		global snarlStart
		global dagOtherEnd
		global dagEnd
		global dagDict
		global dagSide
		global hap1NodePath
		# going now from start of snarl to the left
		print('here \t',graph.get_id(other), (graph.get_id(other) in hap1NodePath+hap1RightNodePath),'lis',lis)


		## 5:01 do something with the lis of reachable nodes to keep from going back into the snarl or whatever


		# print('check tip')
		if ( graph.get_degree(other,True) > 0):
			if (graph.get_degree(other,False) > 0):
				print('\t not a tip:',graph.get_id(other))
			else:
				print('go back and try another edge?', print('it is a tip:',graph.get_id(other)))
				return True
		 # something broken here
		if checkDirection(graph.get_id(other),hap1NodePath+hap1RightNodePath):
			print('wrong direction?', graph.get_id(other))
			return True

		def check_tipDegree_unvisitedNode(node,dEnd,dStart):
			print('ends',dStart,dEnd, 'thisNode',graph.get_id(node))
			if ( graph.get_degree(node,True) > 0):
				if (graph.get_degree(node,False) > 0):
					print('\t not a tip:',graph.get_id(node))
					if (graph.get_id(node)==dEnd):# | (graph.get_id(node)==dStart):
						if graph.get_id(node) not in hap1NodePath+hap1RightNodePath:
							tipWalk.append(graph.get_id(node))
							print('\tappend tipWalk',tipWalk)
							for t in tipWalk:
								hap1NodePath.append(t)
							tipWalk.clear()
							graph.follow_edges(node,False,visit_outOfsnarlsEnd_neighbor)
							return True

			return False



		def avoidTip(node):
			global dagOtherEnd
			global dagEnd
			global dagDict
			global dagSide
			global tipWalk
			if dagSide=='end':
				dEND=dagOtherEnd             # an id int
				dENDint = dagOtherEnd        # the end int id
				dStart=dagEnd                # the handle
			else:
				dENDint=graph.get_id(dagEnd) # the end int id
				dEND=dagEnd 		         # the handle
				dStart=dagOtherEnd           # an id int
			# dStart,dEND = find_D_startAndEnd(node)

			# make a list of all reachable nodes return end node if it is reachable
			foundEnd = checkEndReachable(node,dENDint) 
			print(graph.get_id(node),'found End',foundEnd, 'dEND',dENDint)
			if foundEnd > 0: 
				foundEnd=-1
				tipWalk.append(graph.get_id(node))
				print()
				if dENDint not in hap1NodePath+hap1RightNodePath:
					print('append tipWalk to path and go out of snarl')
					tipWalk.append(dEND)
					for t in tipWalk:
						hap1NodePath.append(t)
					tipWalk.clear()
					# error here?
					print('side',dagSide,'dEND',dENDint,'dStart',dStart,'dagDict',dagDict.keys(),'\n',hap1NodePath)
					if dagSide=='end':
						runSideFun(dStart,dagDict,dagSide,visit_outOfsnarlsEnd_neighbor)
					else:
						runSideFun(dEND,dagDict,dagSide,visit_outOfsnarlsEnd_neighbor)
					
					############################
					# 	if dagDict[dStart][side]['backward']:
					# 		print('flip')
					# 		graph.follow_edges(graph.get_handle(dEND),False,visit_outOfsnarlsEnd_neighbor)
	
					# if dagDict[dEND][side]['backward']:
					# 	print('flip')
					# 	graph.follow_edges(graph.get_handle(dEND),False,visit_outOfsnarlsEnd_neighbor)

			else:
				foundEnd=-1
				print(graph.get_id(node),'dEND',dEND, 'dStart', dStart)
				if check_tipDegree_unvisitedNode(node,dEND,dStart):
					print('tip Degree unvisitedNode rt False')
					return False
				else:
					if ( graph.get_degree(node,True) > 0):
						if (graph.get_degree(node,False) > 0):
							print('dag but not a simple tip?',graph.get_id(node), graph.get_id(node) not in hap1NodePath)
							print(hap1NodePath)
							##  NOT WORKING 
							if graph.get_id(node) not in hap1NodePath: #+hap1RightNodePath:
								print('append Node', graph.get_id(node))
								tipWalk.append(graph.get_id(node))
								graph.follow_edges(node,False,avoidTip)
							else:
								print('keep looking',tipWalk)
								tipWalk[:-1]#.clear() 
								return True
					else:
						# keep looking for a node that isn't a tip
						return True
			# else:
			# 	return True

		def visit_right_intosnarl_neighbor(snarlNode):
			global snarls
			global snarlsEnds
			global startSnarlOver
			global snarlStartNodeDegrees
			global snarlAttempts
			global snarlWalk
			global snarlEnd 
			global snarlStart

			print('checking direction of:',graph.get_id(snarlNode))
			if checkDirection(graph.get_id(snarlNode),hap1NodePath):
				# pSnarlWalk.clear()
				print('are we in here?')
				startSnarlOver =  True
				return True
			# Make this usable for both snarls and snarlsEnd dictionaries
			d={}
			side=''
			print('snarlNode',graph.get_id(snarlNode),'snarlEnd',snarlEnd)
			if snarlStart in list(snarls.keys()):
				d=snarls
				side='end'
			elif snarlStart in list(snarlsEnds.keys()):
				d=snarlsEnds
				side='start'
			
			# print('side',side,'dict',d.keys())
			# going into the parent snarl and exploring all paths to the other side
			# store the node name and parental kmer count if there is one
			print('\t',graph.get_id(snarlNode),'end',d[snarlStart][side][side+'_id'])
			
			# if we reach the end of the snarl stop the walk. 
			if graph.get_id(snarlNode) != d[snarlStart][side][side+'_id']:
				# limit the number of times  a node is visited in a snarl : 2 max rn
				if snarlWalk.count(graph.get_id(snarlNode)) <3:
					# append this node to the walk and check all edges for the end of the snarl
					snarlWalk.append(graph.get_id(snarlNode))
					#### change this to foundEnd = checkEndReachable(snarlNode,int(d[snarlStart][side][side+'_id']))
					# if foundEnd>0:
						#foundEnd=-1
						#visit_right_intosnarl_neighbor(graph.get_handle(int(d[snarlStart][side][side+'_id'])))
					#else:
						#foundEnd=-1
						#graph.follow_edges(snarlNode,False,visit_right_intosnarl_neighbor)
					graph.follow_edges(snarlNode,False,makelis)
					graph.follow_edges(graph.flip(snarlNode),False,makelis) # 
					print('lis',lis)
					# if the end of the snarl is in that list - go there 
					if int(d[snarlStart][side][side+'_id']) in lis:
						print('lis',lis)
						lis.clear()
						visit_right_intosnarl_neighbor(graph.get_handle(int(d[snarlStart][side][side+'_id'])))
					else:
						lis.clear()
						graph.follow_edges(snarlNode,False,visit_right_intosnarl_neighbor)
					####
				else:
					# changed to FALSE 6/8
					print('rt true?',graph.get_id(snarlNode),snarlWalk)
					startSnarlOver=False
					snarlWalk.clear()
					return False
			else:
				print('end walk',d[snarlStart][side],'\n')
				snarlWalk.append(graph.get_id(snarlNode))
				for s in snarlWalk:
					hap1NodePath.append(s)
				# hap1NodePath = hap1NodePath+snarlWalk
				snarlWalk.clear()
				print('hap1Nodes\n',hap1NodePath)

				return False
			# TODO: this doesn't work 
			print('getting here? startOver?',startSnarlOver,'attmps',snarlAttempts,'degrees',snarlStartNodeDegrees)
			# if startSnarlOver:
			# 	if snarlAttempts <= snarlStartNodeDegrees:
			# 		snarlAttempts+=1
			# 		return True
			# 	else:
			# 		snarlAttempts=0
			# 		snarlStartNodeDegrees=-1
			# 		return False
			# else:
			return False
			
			# return True allows for all sides of bubble to be traversed
			# and can get stuck - never finding the end?
			# return False

		# stop conditions:
		## if we are hitting the start of another parent snarl stop
		if graph.get_id(other) in list(parentSnarls.keys()):
			
			global nextParentStartNode 
			print('parent snarl stop?',graph.get_id(other),'next parentStart', nextParentStartNode)
			nextParentStartNode= graph.get_id(other)
			print('in def next parent node?',nextParentStartNode,'\n')
			return False
		# if we run out of graph stop

		# Go if:
		## if we are going into a bi-allelic bubble
			#1 consider both sides and compare parental kmers?
		### TODO save lengths of both sides of the bubbles - for N50/ len distribution
		elif graph.get_id(other) in list(biAllelic_snarls.keys()):
			# print('\t',graph.get_id(other))
			hap1NodePath.append(graph.get_id(other))
			graph.follow_edges(other,False,visit_outOfsnarlsEnd_neighbor)

		## if we are entering a DAG proceed until we find the end of the DAG 
		# in the case it is a tip checking both sides should work
		elif graph.get_id(other) in list(directed_acyclic_net_graph.keys()):
			print('\t DAG',graph.get_id(other) , directed_acyclic_net_graph[graph.get_id(other)])
			dagDict = directed_acyclic_net_graph
			dagSide='end'
			dagEnd=other
			dagOtherEnd=directed_acyclic_net_graph[graph.get_id(other)]['end']['end_id']
			if dagOtherEnd not in hap1NodePath:
				hap1NodePath.append(graph.get_id(other))
				graph.follow_edges(other,False,avoidTip)
				# avoidTip(other)
			else:
				# hap1NodePath.append(graph.get_id(other))
				graph.follow_edges(other,False,visit_outOfsnarlsEnd_neighbor)

		elif graph.get_id(other) in list(directed_acyclic_net_graphEnds.keys()):
			print('\t DAG end', graph.get_id(other), directed_acyclic_net_graphEnds[graph.get_id(other)])
			dagDict = directed_acyclic_net_graphEnds
			dagSide='start'
			dagEnd=other
			dagOtherEnd=directed_acyclic_net_graphEnds[graph.get_id(other)]['start']['start_id']
			print('dagEnd',graph.get_id(other),'dagOtherEnd',dagOtherEnd, 'side',dagSide)
			if dagOtherEnd not in hap1NodePath:
				hap1NodePath.append(graph.get_id(other))
				graph.follow_edges(other,False,avoidTip)
				# avoidTip(other)
			else:
				graph.follow_edges(other,False,visit_outOfsnarlsEnd_neighbor)
		
		## Entering a snarl
		elif graph.get_id(other) in list(snarls.keys()):
			hap1NodePath.append(graph.get_id(other))

			print('\t SNARL',graph.get_id(other), 'end', snarls[graph.get_id(other)]['end']['end_id'], snarls[graph.get_id(other)]['start']['backward'], graph.get_degree(other,False) ) 
			
			snarlEnd = snarls[graph.get_id(other)]['end']['end_id']
			snarlStart = graph.get_id(other)

			navigateSnarl(other,  snarls, snarlsEnds, visit_right_intosnarl_neighbor, visit_outOfsnarlsEnd_neighbor, hap1NodePath+hap1RightNodePath,'end')

		elif graph.get_id(other) in list(snarlsEnds.keys()):
			print('append',graph.get_id(other))
			hap1NodePath.append(graph.get_id(other))
			
			print('\t SNARL',graph.get_id(other), 'end', snarlsEnds[graph.get_id(other)]['start']['start_id'], snarlsEnds[graph.get_id(other)]['start']['backward'], graph.get_degree(other,False) ) 
			snarlEnd = snarlsEnds[graph.get_id(other)]['start']['start_id']
			snarlStart = graph.get_id(other)
			
			navigateSnarl(other, snarlsEnds, snarls, visit_right_intosnarl_neighbor, visit_outOfsnarlsEnd_neighbor, hap1NodePath+hap1RightNodePath,'start')

		# if we are just in none of those situations
			# print('no snarl?',graph.get_id(other), hap1NodePath)
		elif graph.get_id(other) not in hap1NodePath:
			print('just append it and go out of snarls end again')
			hap1NodePath.append(graph.get_id(other))
			graph.follow_edges(other,False,visit_outOfsnarlsEnd_neighbor)
		else:
			nextParentStartNode=-1
		return False

	### Second Pass
	# if this is second pass, go out of snarl ( left )
	# we know it is the second pass because the right hand pass has been stored in RightNodePath list
	if len(hap1RightNodePath)>0:
		# do this for chunk 86 -- or make a flag about the direction of first snarl when we started
		if not (runFunc(here,parentSnarls,visit_outOfsnarlsEnd_neighbor)):
			runFunc(graph.flip(here),parentSnarls,visit_outOfsnarlsEnd_neighbor)
		# This is the end, get out of the function
		return False

	### visit_right_intoParentsnarl_neighbor
	else: # this is the first pass to the right into the snarl
		print('into Parent Snarl neighbor')
		print('is this the end?',parentSnarls[graph.get_id(here)]['end'], currentParentStartNode)
		psnarlEndFound = checkEndReachable(here,int(parentSnarls[graph.get_id(here)]['end']))
		if psnarlEndFound > 0:
			visit_right_intoParentsnarl_neighbor(graph.get_handle(psnarlEndFound))
		else:
			# find all neighbors to the right - into the parent snarl
			runFunc(here,parentSnarls,visit_right_intoParentsnarl_neighbor)

	# print the path so far L->R through the first parent snarl:
	print(graph.get_id(here),'hap1 path L->R thru',graph.get_id(here),'-',parentSnarls[graph.get_id(here)]['end'],' Parent Snarl',hap1NodePath,'\n\t continue from ',parentSnarls[graph.get_id(here)]['end'])

	### visit_outOfsnarlsEnd_neighbor
	# find path to the right -> after the snarl
	# if it is not a tip/end of the graph
	psnarlEnd = parentSnarls[graph.get_id(here)]['end']
	if ( graph.get_degree(graph.get_handle(int(psnarlEnd)),True) > 0):
		if (graph.get_degree(graph.get_handle(int(psnarlEnd)),False) > 0):
			# after going through the snarl go out the right side of the Snarl 
			# until another parent snarl is found, or we run out of graph
			if int(psnarlEnd) not in hap1NodePath:
				print('append snarlEnd',psnarlEnd)
				hap1NodePath.append(int(psnarlEnd))
			runFunc(graph.get_handle(int(psnarlEnd)),parentSnarlsEnds,visit_outOfsnarlsEnd_neighbor)

	return nextParentStartNode 

def visit_handle_checkSnarls(here, snDict):
	global path
	global pathRight
	global pathLeft
	# global dagOtherEnd
	# global dagEnd
	# finds all the paths through the snarl - but 'returns True' so it checks all paths
	if checkDirection(graph.get_id(here),path+pathLeft):
		return True

	# here 
	print('handle',graph.get_id(here))#, 'backward',snDict[graph.get_id(here)], graph.get_is_reverse(here))
	
	global nextSnarlStartNode
	def visit_outOfsnarlsEnd_neighbor2(other):
		global dagOtherEnd
		global dagEnd
		# going now from start of snarl to the left
		print('here \t',graph.get_id(other), (graph.get_id(other) in path+pathRight))

		 # something broken here
		if checkDirection(graph.get_id(other),path+pathRight):
			print('wrong direction?', graph.get_id(other))
			return False


		def check_tipDegree_unvisitedNode(node,dEnd,dStart):
			print('ends',dStart,dEnd, 'thisNode',graph.get_id(node))
			if ( graph.get_degree(node,True) > 0):
				if (graph.get_degree(node,False) > 0):
					print('\t not a tip:',graph.get_id(node))
					if (graph.get_id(node)==dEnd) | (graph.get_id(node)==dStart):
						if graph.get_id(node) not in path+pathRight:
							tipWalk.append(graph.get_id(node))
							print('\tappend tipWalk',tipWalk)
							for t in tipWalk:
								path.append(t)
							tipWalk.clear()
							graph.follow_edges(node,False,visit_outOfsnarlsEnd_neighbor2)
							return True

			return False

		def avoidTip2(node):
			global dagOtherEnd
			global dagEnd
			global tipWalk
			dEND=dagOtherEnd
			dStart=dagEnd
			print(graph.get_id(node),'dEND',dEND, 'dStart', dStart, 'other',graph.get_id(other))
			
			if check_tipDegree_unvisitedNode(node,dEND,dStart):
				return False
			else:
				if ( graph.get_degree(node,True) > 0):
					if (graph.get_degree(node,False) > 0):
						print('dag but not a simple tip?')
						if graph.get_id(node) not in path+pathRight:
							tipWalk.append(graph.get_id(node))
							# look for the end of the DAG
							foundEnd = checkEndReachable(node,dEND)
							if foundEnd > 0: # if we found the end go there
								if avoidTip2(graph.get_handle(foundEnd)):
									tipWalk.append(foundEnd)
									print('it returned True, but we found the end tip',tipWalk)
									for t in tipWalk:
										path.append(t)
									graph.follow_edges(graph.get_handle(foundEnd),False,visit_outOfsnarlsEnd_neighbor2)
									return False
							else:            #otherwise just follow the next available node
								print('going back into avoidTip2')
								graph.follow_edges(node,False,avoidTip2)
						else:
							# keep looking
							tipWalk.clear() 
							return True
					else:
						return True
				else:
					return True
			return False

		def visit_right_intosnarl_neighbor2(snarlNode):
			# going into the parent snarl and exploring all paths to the other side
			# store the node name and parental kmer count if there is one
			print('\t',graph.get_id(snarlNode),'end',snarls[graph.get_id(other)]['end']['end_id'])
			global startSnarlOver
			global snarlStartNodeDegrees
			global snarlAttempts
			global snarlWalk
			
			# if we reach the end of the snarl stop the walk. 
			if graph.get_id(snarlNode) != snarls[graph.get_id(other)]['end']['end_id']:
				# print('whats happening?',graph.get_id(snarlNode),snarlWalk)
				# if graph.get_id(snarlNode) not in snarlWalk:
				if snarlWalk.count(graph.get_id(snarlNode)) <3:
					snarlWalk.append(graph.get_id(snarlNode))

					graph.follow_edges(snarlNode,False,makelis) # 
					print('lis',lis)
					# if the end of the snarl is in that list - go there 
					if int(snarls[graph.get_id(other)]['end']['end_id']) in lis:
						print('lis',lis)
						lis.clear()
						visit_right_intosnarl_neighbor2(graph.get_handle(int(snarls[graph.get_id(other)]['end']['end_id'])))
					else:
						lis.clear()
						graph.follow_edges(snarlNode,False,visit_right_intosnarl_neighbor2)
				else:

					print('rt true?',graph.get_id(snarlNode),snarlWalk)
					startSnarlOver=True
					snarlWalk.clear()
					# changed this 6/8 from True
					return False
			else:
				print('end walk',snarls[graph.get_id(other)]['end'],'\n')
				snarlWalk.append(graph.get_id(snarlNode))
				for s in snarlWalk:
					path.append(s)
				# hap1NodePath = hap1NodePath+snarlWalk
				snarlWalk.clear()
				print('path\n',path)

				return False
			# TODO: this doesn't work 
			print('getting here? startOver?',startSnarlOver,'attmps',snarlAttempts,'degrees',snarlStartNodeDegrees)
			if startSnarlOver:
				if snarlAttempts <= snarlStartNodeDegrees:
					snarlAttempts+=1
					return True
				else:
					snarlAttempts=0
					snarlStartNodeDegrees=-1
					return False
			else:
				return False
			
			# return True allows for all sides of bubble to be traversed
			# and can get stuck - never finding the end?
			# return False
		
		def navigateSnarl(currentOther, dictForSnarl):
			# decide how to proceed at snarl edge: 'start' or 'end'
			global path
			global snarlStartNodeDegrees
			snarlEnd = dictForSnarl[graph.get_id(other)]['end']['end_id']
			# if we haven't already passed the end of the snarl go through it
			if snarlEnd not in path:
				snarlStartNodeDegrees=graph.get_degree(other,False)
				print('rt into snarl?',graph.follow_edges(other,False,visit_right_intosnarl_neighbor2))
			else: 
				# if we have already passed through the snarl set the current node to the snarl end and keep going
				snarlEnd = graph.get_id(other)
			# check the orientation of the end node of the snarl and go from there out of snarl
			if dictForSnarl[graph.get_id(other)]['end']['backward']:
				graph.follow_edges(graph.flip(graph.get_handle(int(snarlEnd))),False,visit_outOfsnarlsEnd_neighbor2)
			else:
				graph.follow_edges(graph.get_handle(int(snarlEnd)),False,visit_outOfsnarlsEnd_neighbor2)

		# stop conditions:
		#  - if we are hitting the start of another parent snarl stop
		#  - if we run out of graph stop
		### Entering a snarl
		if graph.get_id(other) in list(snarls.keys()):
			path.append(graph.get_id(other))
			print('\t SNARL',graph.get_id(other), 'end', snarls[graph.get_id(other)]['end']['end_id'], snarls[graph.get_id(other)]['start']['backward'], graph.get_degree(other,False) ) 
			navigateSnarl(other, snarls)
		# sometimes we can hit the 'end' of a snarl first
		elif graph.get_id(other) in list(snarlsEnds.keys()):
			path.append(graph.get_id(other))
			print('\t SNARL',graph.get_id(other), 'end', snarlsEnds[graph.get_id(other)]['start']['start_id'], snarlsEnds[graph.get_id(other)]['end']['backward'], graph.get_degree(other,False) ) 
			snarlEnd = snarlsEnds[graph.get_id(other)]['start']['start_id']
			navigateSnarl(other, snarlsEnds)

		## if we are going into a bi-allelic bubble
			#TODO: 1 consider both sides and compare parental kmers?
		elif graph.get_id(other) in list(biAllelic_snarls.keys()):
			path.append(graph.get_id(other))
			graph.follow_edges(other,False,visit_outOfsnarlsEnd_neighbor2)

		## if we are entering a DAG proceed until we find the end of the DAG avoing the tip
		elif graph.get_id(other) in list(directed_acyclic_net_graph.keys()):
			dagEnd=graph.get_id(other)
			dagOtherEnd=directed_acyclic_net_graph[graph.get_id(other)]['end']['end_id']
			print('\t DAG start',graph.get_id(other), dagOtherEnd )
			if dagOtherEnd not in path+pathRight:
				path.append(graph.get_id(other))
				graph.follow_edges(other,False,avoidTip2)
			else:
				graph.follow_edges(other,False,visit_outOfsnarlsEnd_neighbor2)

		elif graph.get_id(other) in list(directed_acyclic_net_graphEnds.keys()):
			dagEnd=graph.get_id(other)
			dagOtherEnd=directed_acyclic_net_graphEnds[graph.get_id(other)]['start']['start_id']
			print('\t DAG end',graph.get_id(other), dagOtherEnd )
			if dagOtherEnd not in path+pathRight:
				path.append(graph.get_id(other))
				graph.follow_edges(other,False,avoidTip2)
			else:
				graph.follow_edges(other,False,visit_outOfsnarlsEnd_neighbor2)
		


		# if we are just in none of those situations
			# print('no snarl?',graph.get_id(other), hap1NodePath)
		elif graph.get_id(other) not in path:
			path.append(graph.get_id(other))
			graph.follow_edges(other,False,visit_outOfsnarlsEnd_neighbor2)
		else:
			nextSnarlStartNode=-1
		return False

	### Second Pass
	# if this is second pass, go out of snarl ( left )
	# we know it is the second pass because the right hand pass has been stored in RightNodePath list
	if len(pathRight)>0:
		# do this for chunk 86 -- or make a flag about the direction of first snarl when we started
		if not (runFunc(here,snDict,visit_outOfsnarlsEnd_neighbor2)):
			runFunc(graph.flip(here),snDict,visit_outOfsnarlsEnd_neighbor2)

		# This is the end, get out of the function
		return False

	### visit_outOfsnarlsEnd_neighbor2 with start
	else: # this is the first pass to the right into the snarl
		# don't use runFun here because we don't want to follow any edges at first here
		if snDict[graph.get_id(here)]['backward']:
			if visit_outOfsnarlsEnd_neighbor2(graph.flip(here)):
				print('unflip',graph.get_id(here))
				visit_outOfsnarlsEnd_neighbor2(here)

		else:
			if visit_outOfsnarlsEnd_neighbor2(here):
				print('unflip',graph.get_id(here))
				visit_outOfsnarlsEnd_neighbor2(graph.flip(here))
		# find all neighbors to the right - into the parent snarl
		# runFunc(here,snDict,visit_outOfsnarlsEnd_neighbor2)

	# print the path so far L->R through the first parent snarl:
	print(graph.get_id(here),'hap1 path L->R thru',graph.get_id(here),'-',snDict[graph.get_id(here)]['end']['end_id'],' section ',path,'\n\t continue from ',snDict[graph.get_id(here)]['end']['end_id'])# 'ends dict:',snDictEnds[graph.get_id(here)])

	### visit_outOfsnarlsEnd_neighbor2
	# find path to the right -> after the snarl
	# if it is not a tip/end of the graph
	thisSnarlEnd = snDict[graph.get_id(here)]['end']['end_id']
	thisSnarlEndBackwards = snDict[graph.get_id(here)]['end']['backward']
	if ( graph.get_degree(graph.get_handle(int(thisSnarlEnd)),True) > 0):
		if (graph.get_degree(graph.get_handle(int(thisSnarlEnd)),False) > 0):
			# after going through the snarl go out the right side of the Snarl 
			# until another parent snarl is found, or we run out of graph
			print('pSnarlEnd',thisSnarlEnd)
			if thisSnarlEnd not in path:
				print('append snarlEnd')
				path.append(int(thisSnarlEnd))
			if thisSnarlEndBackwards: 
				if visit_outOfsnarlsEnd_neighbor2(graph.flip(graph.get_handle(thisSnarlEnd))):
					visit_outOfsnarlsEnd_neighbor2(graph.get_handle(thisSnarlEnd))
			else:
				if visit_outOfsnarlsEnd_neighbor2(graph.get_handle(thisSnarlEnd)):
					visit_outOfsnarlsEnd_neighbor2(graph.flip(graph.get_handle(thisSnarlEnd)))

			# runFunc(graph.get_handle(int(psnarlEnd)),snDictEnds,visit_outOfsnarlsEnd_neighbor2)
			# print(snDictEnds)

	# now check those neighbors for parent asignment
	# neighborParents = labelNeighborParents(neighbors)

	# write node to bandage csv for coloring nodes when viewing gfa
	# writeOutColorLabeledBandageCSV(neighborParents)

	return nextSnarlStartNode

def noParentPath(node):
	# there are no parents no snarls
	global path
	global lis
	foundNewNode=False
	print('in noParentPath', graph.get_id(node))
	if graph.get_id(node) not in path:
		print('append',graph.get_id(node))
		path.append(graph.get_id(node))
	# 	graph.follow_edges(node,False,noParentPath)
	# else:
		lis.clear()
		graph.follow_edges(node,False,makelis) 
		graph.follow_edges(graph.flip(node),False,makelis)
		print(graph.get_id(node),'lis',lis,path)
		for n in lis:
			if n not in path:
				print('choose',n)
				lis.clear()
				path.append(n)
				foundNewNode=True
				graph.follow_edges(graph.get_handle(n),False,noParentPath)
				break
		lis.clear()
	if not foundNewNode:
		return False
	return True

# here is main()
# getParentKmers()

def traverseThroughParentSnarls():

	global currentParentStartNode
	global nextParentStartNode
	global hap1Path
	global hap1NodePath
	global hap1RightNodePath
	global hap1LeftNodePath
	global parentSnarlStartKeyList

	parentSnarlStartKeyList=list(parentSnarls.keys())
	if len(parentSnarlStartKeyList)>0:
		currentParentStartNode = parentSnarlStartKeyList[0]
		parentSnarls[currentParentStartNode]['visited']=True
	firstSnarltoFlip=currentParentStartNode

	if currentParentStartNode != -1:
		# if there is a parent snarl, visit it's start node and traverse to the right
		visit_handle_checkParent( graph.get_handle( currentParentStartNode, False ) )
		print('removing currentParentStartNode',currentParentStartNode,'from parent list', parentSnarlStartKeyList, '\n nextParent',nextParentStartNode)
		parentSnarlStartKeyList.remove(currentParentStartNode)

	# if the next start node was found from that traversal, visit it's start node and travese to the right
	while nextParentStartNode != -1:
		# start again at next parent start node
		print('nextParentStartNode',nextParentStartNode, '\n  list of parent starts:',parentSnarlStartKeyList)
		if nextParentStartNode in parentSnarlStartKeyList:
			print('removing nextParentStartNode', nextParentStartNode)
			parentSnarlStartKeyList.remove( nextParentStartNode )
		else:
			break
		# store next parent node as current node
		currentParentStartNode=nextParentStartNode
		# 
		if visit_handle_checkParent( graph.get_handle( currentParentStartNode ) ):
			print('flip handle and go again!', currentParentStartNode) 
			visit_handle_checkParent( graph.flip(graph.get_handle( currentParentStartNode ) ) )

	if len(parentSnarlStartKeyList) > 0:
		print('went throuh next parent and parentSnarlStartKeyList >1\n',parentSnarlStartKeyList)

	# store right path and go the other way
	hap1RightNodePath= [m for m in hap1NodePath]
	hap1NodePath.clear()
	print('hap1NodePath',hap1NodePath,'\nhap1RightNodePath',hap1RightNodePath)
	if len(list(parentSnarls.keys())) > 0:
		print('\n\ngoing left from first parent', firstSnarltoFlip)


		if visit_handle_checkParent( graph.flip(graph.get_handle( firstSnarltoFlip ) )):
			print('un-flip handle and go again?')
			visit_handle_checkParent( graph.get_handle( firstSnarltoFlip ))

	# organize the lists 
	hap1LeftNodePath=[o for o in hap1NodePath[::-1]]
	hap1NodePath.clear()
	# remove 1st node from left path and reverse it before adding L and R together
	hap1Path=hap1LeftNodePath[1:][::-1]+hap1RightNodePath
	# check for doubles
	dbls_hap1 = set([x for x in hap1Path if hap1Path.count(x) > 1])
	print('dbls',dbls_hap1)
	return dbls_hap1

def noParents():
	global whichDict
	global pathStartNode
	global pathRight
	global pathLeft
	global path 
	global plen
	# MAYBE: if there are no biallelic bubbles - just take the whole graph
	# 		 if there are only biallelic bubbles subtract hap1/dad from each other

	goLeft=True	
	### after all that if there is still no path
	if len(list(snarls.keys())) > 0:
		print('\n\t running snarl Key')
		whichDict=snarls
		if visit_handle_checkSnarls( graph.get_handle( list(snarls.keys())[0] ), snarls):
			print('\t flip snarlkey',list(snarls.keys())[0])
			visit_handle_checkSnarls( graph.flip(graph.get_handle( list(snarls.keys())[0] ) ), snarls)
	elif len(list(directed_acyclic_net_graph)) > 0:
		print('\n\t running dag Key')
		whichDict=directed_acyclic_net_graph
		if visit_handle_checkSnarls( graph.get_handle( list(directed_acyclic_net_graph.keys())[0] ), directed_acyclic_net_graph ):
			visit_handle_checkSnarls( graph.flip(graph.get_handle( list(directed_acyclic_net_graph.keys())[0] ) ), directed_acyclic_net_graph)
		
	elif len(list(biAllelic_snarls)) > 0:
		print('\n\t running bi-allelic Key')
		whichDict=biAllelic_snarls
		if visit_handle_checkSnarls( graph.get_handle( list(biAllelic_snarls.keys())[0] ), biAllelic_snarls ):
			visit_handle_checkSnarls( graph.flip(graph.get_handle( list(biAllelic_snarls.keys())[0] ) ), biAllelic_snarls)

	else:
		goLeft=False
		print('runnin no parent path')
		graph.for_each_handle(noParentPath)

	# store path and go the other way
	pathStartNode = path[0]
	pathRight = [p for p in path]
	print('path Rt',pathRight)
	path.clear()
	if goLeft:
		print('\n\ngoing left from first parent', pathStartNode)
		if visit_handle_checkSnarls( graph.flip(graph.get_handle(pathStartNode)), whichDict ):
			print('okay, unflipping', pathStartNode )
			visit_handle_checkSnarls( graph.get_handle(pathStartNode), whichDict )

		pathLeft = [p for p in path]
		path.clear()
	fullPath = pathLeft[1:][::-1]+pathRight
	pathDbls = set([x for x in fullPath if fullPath.count(x)>1])
	print('pathDbls',pathDbls,len(fullPath))
	writeOutPathCSV(fullPath)
	# write primary.fa and calc length of nodes on path
	plen = writeContig(fullPath,args.input_pg[6:][:-3],True)
	getStats(plen,len(fullPath),0,len(pathDbls),0)


def main():
	''' creates and runs the object methods to create kmers '''
	# read in the data with a read in file object
	# read = readInFile()
	# k,dna = read.read()
	# # create kmer object and create kmers 
	# mk = makeKmers(k, dna)
	# mk.kmers()

	''' create snarl and graphPath objects
		run object methods to deseralize graph and traverse it 
		output: 2 haplotype paths ''' 
	inputGraph, inputJson = parseArguments()
	# create an instance of the snarl manager obj and input the snarl file into dictionaries
	sm = snrlManager(inputJson)
	allSnarlDict = sm.snarlsInput()

	# pass the snarl dict and input graph to the graphPathManager object
	gpf = graphPathFinder(allSnarlDict,inputGraph)



	# dbls_hap1 = traverseThroughParentSnarls()

	# # write either the no parent path or the hap1Path through the snarls to CSV
	# whichDict={}

	# if len(hap1Path)==0:
	# 	noParents()
	# else:
	# 	# print('hap1Path',hap1Path)
	# 	writeOutPathCSV(hap1Path)
	# 	# write primary.fa and calc length of nodes on path
	# 	plen = writeContig(hap1Path,args.input_pg[6:][:-3],True)
	# 	# print out stats to csv
	# 	getStats(plen,len(hap1Path),0,len(dbls_hap1),0)

	# print('Path bps / total:',(plen)/(graph.get_total_length() ))

	# # outCSV.close()
	# outpathCSV.close()

	# starts from one known node
	# visit_handle(graph.get_handle(50373))#195))


if __name__ == "__main__":

    main()