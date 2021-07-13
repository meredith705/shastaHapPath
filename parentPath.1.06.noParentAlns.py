#!/usr/bin/env python3

# usage parentPath.py -g blunt.PackedGraph.pg -j blunt.PackedGraph.snarls.json -m mat.node.kmer.count -p pat.node.kmer.count

from bdsg.bdsg import PackedGraph
import argparse 
import numpy as np
import json

SNARL_TYPE_ULTRABUBBLE:1

def parseArguments():
	''' take I/O arguments from command line '''
	parser=argparse.ArgumentParser()
	parser.add_argument('-g','--input_pg')
	# parser.add_argument('-m','--mat_alns')
	# parser.add_argument('-p','--pat_alns')
	parser.add_argument('-j','--input_snarl_filename')

	args=parser.parse_args()
	inputGraph=args.input_pg
	inputSnarlFilename=args.input_snarl_filename
	return inputGraph, inputSnarlFilename

class SnarlManager:
	'''class that loads and contains the snarl dictionaries '''
	def __init__ (self,inputSnarlFilename):
		''' class variables:
			 dictionaries that hold the dictionaries of snarl starts and ends'''
		self.inputSnarlFilename=inputSnarlFilename
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

	# TODO this is a problem because parents can start on the same node.
	# TODO snarls can also start and end on the same node.. Work on this chunk_12
	def snarlsInput(self):
		'''load snarls.json into their respective snarl type dictionaries'''
		snarl_json=[json.loads(line) for line in open(self.inputSnarlFilename,'r')]
		for sn in snarl_json:
		# bi-allelic snarls are type: 1
			if sn.get('parent',False):
				# make sure that the start node has a smaller node id than the end
				snarlStart = int(sn.get('parent').get('start').get('node_id'))
				snarlStartBkwd = sn.get('parent').get('start').get('backward',False)
				snarlEnd = int(sn.get('parent').get('end').get('node_id'))
				snarlEndBkwd = sn.get('parent').get('end').get('backward', False)

				### I'm not orientating the top level snarl start and end nodes, just leaving them as is
				# if snarlStart > snarlEnd:
				# 	# if the start is larger than the end, swich them and toggle the backwards flags
				# 	tempSnStart = snarlStart
				# 	snarlStart = snarlEnd
				# 	snarlEnd = tempSnStart
				# 	print('switched: snarlStart',snarlStart,'snarlEnd',snarlEnd)
				# 	# toggle the backwards flags - I don't think so _ I think snarls just gets this wrong and flips the start/end labels
				# 	if snarlStartBkwd:
				# 		snarlStartBkwd = False
				# 	else:
				# 		snarlStartBkwd = True
				# 	if snarlEndBkwd:
				# 		snarlEndBkwd = False
				# 	else:
				# 		snarlEndBkwd = True

				# TODO this is a problem because parents can start on the same node. check if end nodes are the same?
				# set it up so that smaller numerical value snarl boundry nodes are the start and larger id nodes are the end
				if snarlStart in list(self.sd['parentSnarls'].keys()):
				# 	# print(parentSnarls[int(sn.get('parent').get('start').get('node_id'))]['child'])
					self.sd['parentSnarls'][snarlStart]['child'].append((sn.get('start'),sn.get('end')))
					self.children.append(int(sn.get('start').get('node_id')))
					self.children.append(int(sn.get('end').get('node_id')))
				else:
					
					self.sd['parentSnarls'][snarlStart] = {'unvisited' : True,\
						'start': { 'start_id' : snarlStart, 'backward': snarlStartBkwd},\
						'end': {'end_id' : snarlEnd,'backward' : snarlEndBkwd},'child' : [(sn.get('start'),sn.get('end'))]}
					self.sd['parentSnarlsEnds'][snarlEnd] = {'unvisited' : True, \
						'start': {'start_id' : snarlStart,'backward': snarlStartBkwd}, \
						'end': {'end_id' : snarlEnd,'backward' :snarlEndBkwd}, 'child' : [(sn.get('start'),sn.get('end'))]}
					self.children.append(int(sn.get('start').get('node_id')))
					self.children.append(int(sn.get('end').get('node_id')))
			
			if sn.get('type',False):
				if sn.get('type')==	1:#SNARL_TYPE_ULTRABUBBLE:
					self.sd['biAllelic_snarls'][int(sn.get('start').get('node_id'))] = {\
						'unvisited' : True, \
						'backward': sn.get('start').get('backward', False),\
						'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward', False)},\
						'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
						'child' : sn.get('parent',False)}
					self.sd['biAllelic_snarlsEnds'][int(sn.get('end').get('node_id'))] = {\
						'unvisited' : True,\
						'backward': sn.get('end').get('backward', False),\
						'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward',False)},\
						'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
						'backward' : sn.get('end').get('backward', False), \
						'child' : sn.get('parent',False)}
				else:
					print('Snarl Type not 1:',sn.get('type')) 
			elif sn.get('directed_acyclic_net_graph',False):
				self.sd['directed_acyclic_net_graph'][int(sn.get('start').get('node_id'))] = {\
					'unvisited' : True, \
					'backward': sn.get('start').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward', False)},\
					'start_self_reachable' : sn.get('start_self_reachable',False),\
					'start_end_reachable' : sn.get('start_end_reachable',False),\
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'end_self_reachable' : sn.get('end_self_reachable',False),\
					'child' : sn.get('parent',False)}
				self.sd['directed_acyclic_net_graphEnds'][int(sn.get('end').get('node_id'))] = {\
					'unvisited' : True, \
					'backward': sn.get('end').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')), 'backward': sn.get('start').get('backward',False)}, \
					'start_self_reachable' : sn.get('start_self_reachable',False),\
					'start_end_reachable' : sn.get('start_end_reachable',False),\
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'end_self_reachable' : sn.get('end_self_reachable',False),\
					'child' : sn.get('parent',False)}
			else:
				self.sd['snarls'][int(sn.get('start').get('node_id'))] = {\
					'unvisited' : True, \
					'backward': sn.get('start').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward', False)},\
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'start_end_reachable' : sn.get('start_end_reachable',False),\
					'child' : sn.get('parent',False)}
				self.sd['snarlsEnds'][int(sn.get('end').get('node_id'))] = {\
		
					'unvisited' : True, \
					'backward': sn.get('end').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward',False)},\
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'start_end_reachable' : sn.get('start_end_reachable',False),\
					'child' : sn.get('parent',False)}

		# remove top level snarls from the list when they are also children ( making them lower level snarls)
		for p in list(self.sd['parentSnarls'].keys()):		
			if p in self.children:
				print('PARENT IS ALSO CHILD',p,self.sd['parentSnarls'][int(p)])
				self.sd['parentSnarlsEnds'].pop(int(self.sd['parentSnarls'][int(p)]['end']['end_id']))
				self.sd['parentSnarls'].pop(p)
		# remove top level parent snarls from the lower level snarls when the start and end id's match		
		for p in list(self.sd['parentSnarls'].keys()):
			# if p in keys:
			# pop from snarl dicts? I think so
			for s in ['snarls','snarlsEnds']:
				if p in list(self.sd[s].keys()):
					print('do the snarls match start and end? ',self.sd[s][int(p)])
					print(self.sd['parentSnarls'][int(p)])
					if self.sd[s][int(p)]['end']['end_id'] == self.sd['parentSnarls'][int(p)]['end']['end_id']:
						print('pop',self.sd[s].pop(int(p),s))
						print('pop',self.sd[s].pop(int(self.sd['parentSnarls'][int(p)]['end']['end_id']),s))


		print('check again',list(self.sd['parentSnarls'].keys()),'ends', self.sd['parentSnarlsEnds'].keys()  )

		if len(list( self.sd['snarls'].keys() )) >0:
			print('snarls here',list(self.sd['snarls'].keys()))
		else:
			print('snarls len:',len(list( self.sd['snarls'].keys() )))
		return self.sd

class graphPathFinder:
	'''class that loads and contains the snarl dictionaries '''
	def __init__ (self,sd,inputGraph):
		# snarl dictionary and inputGraph File
		self.sd=sd
		self.inputGraph=inputGraph
		# make a PackedGraph object
		self.graph=PackedGraph() 
		self.topLevelChain=[]
		self.tlEndNode=-1
		self.currentNode=0
		# hap 1 list of handles in the path
		self.hap1Path=[]
		self.hap1PathIds=[]
		self.hap1FinalWalk=[]
		# boolean flag for the final backward traversal
		self.goingBackward = False
		# walks inbetween top level snrls
		self.walkFromEnd_toStart=[]
		# snarls
		self.whichSnDic=''
		self.exit=0
		self.chosenSnarlNodes=0

		# graph stats
		self.hap1_pathLenInBps=[]
		self.hap1_doubledHandles=[]
		# sanarl stats
		self.homNodeLens=[]			# list of hom nodes inbetween snarl boundries
		self.hetNodeAvgLens=[]      # avg len of pairs of biallelic nodes
		self.bluntifiedNodeLens=[]	# snarl boundry nodes, artifact of bluntification
		# self.biAllelicLengths=[]
		self.dagWalkLens=[]			# DAG walk lengths (>3 node DAGs)
		self.dagTipLens=[]			# lens of DAG tips (3 node DAGs)
		# reachable node list
		self.lis=[]
		

	def writeOutPathCSV(self,p,singleNode):
		# open the path.csv and write out the path with colors for bandage viewing 
		# TODO: Output the other haploid path
		outpathCSV=open(self.inputGraph[:-3]+".path.csv",'w')
		outpathCSV.write("Node,Path,Color\n")
		# parent1='hg03'
		parent2='hg04'
		# parent1Color = "light sea green"
		parent2Color = 'tomato'
		if singleNode:
			for n in p:
				outpathCSV.write( str(self.graph.get_id(n)) +",M,"+parent2Color+"\n")
		else:
			for walk in p:
				for n in walk:
					outpathCSV.write( str(self.graph.get_id(n)) +",M,"+parent2Color+"\n")
		outpathCSV.close()

	#  TODO find faster fasta output method
	def writeContig(self, hap1Path):
		# output all the sequence from the path to a fasta file ( 80 characters)
				# Don't forget to check orientations with e.g. graph.get_is_reverse(here)
		charMax = 70
		charCount = 0
		priLength = 0
		# if primary:
			# for now I'm opening 'w'. if I loop over all connected components then append 'a'
		file = open(self.inputGraph[:-3]+'.pri.fa','w')
		# else:
		# 	file = open(self.inputGraph[:-3]+'.alt.fa','w')

		file.write('>'+self.inputGraph[:-3]+'\n')
		# change to handles instead of graph ids 
		for node in hap1Path:
			seq = self.graph.get_sequence( node )
			priLength+=len(seq)


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

	def writeOneContig(self, hap1Path):
		# output all the sequence from the path to a fasta file ( 80 characters)
				# Don't forget to check orientations with e.g. graph.get_is_reverse(here)
		charMax = 70
		charCount = 0
		priLength = 0
		# if primary:
			# for now I'm opening 'w'. if I loop over all connected components then append 'a'
		file = open(self.inputGraph[:-3]+'.pri.fa','w')
		# else:
		# 	file = open(self.inputGraph[:-3]+'.alt.fa','w')

		file.write('>'+self.inputGraph[:-3]+'\n')
		# change to handles instead of graph ids 

		seq = self.graph.get_sequence( hap1Path )
		priLength+=len(seq)


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

	def getStats(self):
		''' open the stats.csv, add the header, and all stats values'''		
		outSTAT=open(self.inputGraph[:-3]+".stats.csv",'w')
		outSTAT.write("chunkNum,chunk_Len,NumNodes,num_parent_snarls,num_snarls,num_dang,num_bi_allelics,Sum_of_Avg_hetSeq_Len,hetSeqAvgLen_avgLen,hetSeqAvgLen_medianLen,total_seq_inTip,hap1Path_len,num_hap1Path_Nodes,Alt_Nodes,duplicate_nodes,NotIncluded_Nodes,NotInc_Nodes/TotalNodes,bps_not_inPri,hap1len/totalLen\n")
		#chunkNum
		outSTAT.write(self.inputGraph[6:][:-3]+",")
		#chunk_Len
		outSTAT.write(str(self.graph.get_total_length())+",")
		#NumNodes
		outSTAT.write(str(self.graph.get_node_count())+",")
		#num_parent_snarls
		outSTAT.write(str(len(self.sd['parentSnarls']))+",")
		#num_snarls, num_dang_snarls, num_biallelic
		outSTAT.write(str(len(self.sd['snarls']))+",")#+str(len(parentSnarlStartKeyList))+",")
		outSTAT.write(str(len(self.sd['directed_acyclic_net_graph']))+",")
		outSTAT.write(str(len(self.sd['biAllelic_snarls']))+",")
		
		# hetSeq_Len,bi_allelic_avgLen,bi_allelic_medianLen
		if len(self.hetNodeAvgLens) > 0:
			outSTAT.write(str(int(np.sum(self.hetNodeAvgLens)))+",")
			outSTAT.write(str(int(np.mean(self.hetNodeAvgLens)))+",")
			outSTAT.write(str(int(np.median(self.hetNodeAvgLens)))+",")
		else:
			outSTAT.write('0,0,0,')

		# total_seq_inTip
		if len(self.dagTipLens) > 0:
			outSTAT.write(str(int(np.sum(self.dagTipLens)))+",")
		else:
			outSTAT.write('0,')
		# hap1Path_len
		hap1_length = np.sum(np.array(self.hap1_pathLenInBps))
		print('hap1_length',hap1_length)
		outSTAT.write(str(hap1_length)+",")
		# num_hap1Path_Nodes
		outSTAT.write(str(len(self.hap1FinalWalk))+",")
		outSTAT.write("0,")# alt nodes
		#duplicate nodes ~ kinda different now with handles 
		outSTAT.write(str(len(self.hap1_doubledHandles))+",")
		# NotIncluded_Nodes
		notIncludedNodes = (self.graph.get_node_count()-(len(self.hap1FinalWalk)))
		outSTAT.write(str( notIncludedNodes  )+",")
		# NotInc_Nodes/TotalNodes
		outSTAT.write(str( (notIncludedNodes/self.graph.get_node_count()) )+",")
		# bps_not_inPri
		outSTAT.write(str( (self.graph.get_total_length()-(hap1_length) ) )+",")
		# missingParentSnarls
		# outSTAT.write(str(len(snarls)-len(parentSnarlStartKeyList))+",")
		# covered/totParentSnarls
		# if len(snarls)>0:
		# 	outSTAT.write(str(len(parentSnarlStartKeyList)/len(snarls))+",")
		# else:
		# 	outSTAT.write('0,')
		# path_len/totalLen (bps/bps)
		outSTAT.write(str( hap1_length/self.graph.get_total_length() ) +"\n")

		outSTAT.close()

	def calcHetAvgLen(self,node):
		''' get length of both sides of a biallelic bubble and return mean length of the bubbles sides '''
		# fill the lis with the forward reachable nodes from this boundry node
		self.fwdReachableNodes(node,False)
		ls = np.array([self.graph.get_length(l) for l in self.lis])
		print('biallelic mean',np.mean(ls))
		# append the mean length to the list of all biallelic bubble lengths
		self.hetNodeAvgLens.append(np.mean(ls))
	
	def calcTipLen(self,node):
		''' check degree and find len of tip '''
		def itIsATip(tNode):
			''' if either side of this node has an edge degree of 0, it is a tip
				returns True if the node has 0 degrees on either side '''
			# note the graph often ends in a tip
			if self.graph.get_degree(tNode,False)==0:
				print(self.graph.get_id(tNode),'is a tip')
				return True
			elif self.graph.get_degree(tNode,True)==0:
				print(self.graph.get_id(tNode),'is a tip')
				return True
			else:
				return False
		# from a list of all fwd reachable nodes store the length of any tips
		self.fwdReachableNodes(node,False)
		ls = np.array([self.graph.get_length(l) for l in self.lis if itIsATip(l)])
		for l in ls:
			print('tip Len',l)
			self.dagTipLens.append(l)

	def deserializeGraph(self):
		''' deserialize the packed graph into memory and begin chaining parents together'''
		self.graph.deserialize(self.inputGraph)
		self.chooseWalkingMethod()

	def chooseWalkingMethod(self):
		''' depending on the shape of the graph and the number of parent snarls , choose a walking method '''
		print('node count',self.graph.get_node_count())
		if self.graph.get_node_count()==1:
			self.graph.for_each_handle(lambda y: self.hap1Path.append(y))
			print('THe one node hap1Path:\n')
			for h in self.hap1Path:
				print(self.graph.get_id(h))
				if self.hap1FinalWalk.count(h) == 0:
					self.hap1FinalWalk.append(h)
					self.hap1_pathLenInBps.append(self.graph.get_length(h))
			 # output the path to a csv with appropriate colors for bandage
			self.writeOutPathCSV(self.hap1Path,True)
			self.getStats()
			print(np.sum(np.array(self.hap1_pathLenInBps)), (np.sum(np.array(self.hap1_pathLenInBps))/self.graph.get_total_length()))
			self.writeContig(self.hap1FinalWalk)
		elif len(self.sd['parentSnarls'])==0:
			self.noParentPath()
		else:
			# chain parents together
			self.chainParents()
	
	def fwdReachableNodes(self,other,bothEnds):
		''' follow edges in fwd direction adding nodes to the lis, 
			bothEnds is a boolean to look at both sides of the handle, False = check both sides'''
		def add2List(y):
			self.lis.append(y)
			return True
		self.lis.clear()	
		self.graph.follow_edges(other,False,add2List) 

		if (bothEnds==False):
			self.graph.follow_edges(other,True,add2List) 

	def goOppDirFromLeftmostTopLevelSnStart(self):
		'''go 'left' from most 'left' start snarl - the one that hasn't been visited yet'''
		unVisitedTopSns=[x for x in list(self.sd['parentSnarls'].keys()) if self.sd['parentSnarls'][x]['unvisited']]
		print('going opposite from first start, unVisitedTopSns',unVisitedTopSns, self.graph.get_id(self.currentNode))
		# set the start node to visited and flip the directional flag
		#TODO change this incase ther are no parent snarls 
		if len(unVisitedTopSns)>0:
			self.sd['parentSnarls'][unVisitedTopSns[0]]['unvisited']=False
			self.goingBackward=True
		else:
			print('\nWe lost the origional parent start node!\n')

		#flip the only unvisited start node and go from there till end
		if len(unVisitedTopSns)>0:
			# flip the boolean handle orientation to go the out of the snarl in the backward direction
			snarlHandleDirection = self.sd['parentSnarls'][unVisitedTopSns[0]]['start']['backward']
			handleDir = True
			if snarlHandleDirection:
				handleDir = False
			# set the graph handle to the opposite of whatever the snarls file told us origionally
			self.currentNode=self.graph.get_handle(unVisitedTopSns[0],handleDir)
			self.fwdReachableNodes(self.currentNode,True)
			nodesInWalk = [ l for l in self.lis if self.graph.get_id(l) in self.hap1PathIds]
			print('lis',[self.graph.get_id(l) for l in self.lis],'self.hap1PathIds',self.hap1PathIds)
			print('nodesInWalk',[self.graph.get_id(n) for n in nodesInWalk]) 
			# if there are reachable nodes that are in the current walk - then we shouldn't have flipped the start node
			if len(nodesInWalk)>0:
				self.currentNode=self.graph.get_handle(unVisitedTopSns[0],snarlHandleDirection)
				self.fwdReachableNodes(self.currentNode,True)
				print('lis',[self.graph.get_id(l) for l in self.lis],'self.hap1PathIds',self.hap1PathIds)
			# get left side of graph
			self.findNextTLs_Start()

	def nextNode(self,nextN):
		'''moves pointer to first reachable node in the graph'''
		print('nextNode?',self.graph.get_id(nextN))
		self.currentNode=nextN
		return False # rt true didn't stop the flow at degree>1

	def chooseNextNode(self,nextN):
		'''moves pointer to each reachable node until the chosen self.exit node is reached
			self.graph.follow_edges(self.currentNode,False,self.chooseNextNode)
		'''
		print('choose nextNode?',self.graph.get_id(nextN), int(self.graph.get_id(nextN)) == self.exit)
		if self.graph.get_id(nextN) == self.exit:
			self.currentNode=nextN
			# self.exit=0
			return False
		else: 
			# flipedHandleNextN=self.get_handle()
			return True

	def noParentPath(self):
		''' there are no top level snarls and we only have either lower level 
			snarls, dag snarls, or biallelic snarls 
			 '''
		startOn = 0
		def getTheHandle(node):
			''' get a handle from the graph'''
			self.currentNode = node
			return False
		print('no tlparent Path',)
		# look for snarls to anchor ourselves to the graph as a place to start
		for s in ['snarlsEnds','directed_acyclic_net_graphEnds','biAllelic_snarlsEnds']:
			if len(self.sd[s])>0:
				for key in list(self.sd[s].keys()):
					print(s,key,self.sd[s][key]['start'])
					if self.graph.get_degree(self.graph.get_handle(self.sd[s][key]['end']['end_id'],self.sd[s][key]['end']['backward']),False) > 0:
						node_id = self.sd[s][key]['end']['end_id']
						orientation = self.sd[s][key]['end']['backward']
						self.currentNode=self.graph.get_handle(node_id,orientation)
						print('startNode',self.graph.get_id(self.currentNode) )
						break
					if self.graph.get_degree(self.graph.get_handle(self.sd[s][key]['start']['start_id'],self.sd[s][key]['start']['backward']),False) > 0:
						node_id = self.sd[s][key]['start']['start_id']
						orientation = self.sd[s][key]['start']['backward']
						self.currentNode=self.graph.get_handle(node_id,orientation)
						print('startNode',self.graph.get_id(self.currentNode) )
						break
		if self.currentNode != 0:
			print('picked a node to start at',self.graph.get_id(self.currentNode) )
			self.fwdReachableNodes(self.currentNode,True)
			startOn = self.currentNode
			self.findNextTLs_Start()
			
			# flip the boolean handle orientation to go the out of the snarl in the backward direction
			print('now go he opposite way')
			handleDir = True
			if orientation:
				handleDir = False
			self.currentNode=self.graph.get_handle(node_id,handleDir)
			self.fwdReachableNodes(self.currentNode,True)
			nodesInWalk = [ l for l in self.lis if self.graph.get_id(l) in self.hap1PathIds]

			print('nodesInWalk',[self.graph.get_id(n) for n in nodesInWalk]) 
			# if we already reached nodes that are in the current walk - then we shouldn't have flipped the start node
			if len(nodesInWalk)>0:
				self.currentNode=self.graph.get_handle(node_id,orientation)
			# get left side of graph
			self.findNextTLs_Start()

			for walk in self.hap1Path:
				for h in walk:
					print(self.graph.get_id(h))
					if self.hap1FinalWalk.count(h) == 0:
						self.hap1FinalWalk.append(h)
						self.hap1_pathLenInBps.append(self.graph.get_length(h))
				 # output the path to a csv with appropriate colors for bandage
				self.writeOutPathCSV(self.hap1FinalWalk,True)
				self.getStats()
			self.writeContig(self.hap1FinalWalk)

		else:		
			# pick a random handle for the current node and go from there

			self.graph.for_each_handle(getTheHandle)
			print('\t\npick random handle', self.graph.get_id(self.currentNode))
			self.fwdReachableNodes(self.currentNode,True)
			if len(self.lis) > 0:
				startOn = self.currentNode
				self.findNextTLs_Start()
			else:
				self.graph.for_each_handle(getTheHandle)
				print('\tpick another random handle', self.graph.get_id(self.currentNode))
				startOn = self.currentNode
				self.findNextTLs_Start()

				# error here for chunk 107
			for w in self.hap1Path:
				# this is a n error because its graph.get_id([handle])
				for h in w:
					print(self.graph.get_id(h))
					if self.hap1FinalWalk.count(h) == 0:
						self.hap1FinalWalk.append(h)
						self.hap1_pathLenInBps.append(self.graph.get_length(h))
			print('double nodes here?:', [self.graph.get_id(l) for l in self.hap1FinalWalk])
			### TODO: check that the length of path isn't longer than total chunk length
			# like chunk 283
			if self.graph.get_node_count() == 2:
				self.hap1FinalWalk = self.hap1FinalWalk[0:2]
				self.hap1_pathLenInBps=[]
				for h in self.hap1FinalWalk:
					self.hap1_pathLenInBps.append(self.graph.get_length(h))
			print('double nodes here?:', [self.graph.get_id(l) for l in self.hap1FinalWalk])
			# output the path to a csv with appropriate colors for bandage
			self.writeOutPathCSV(self.hap1FinalWalk,True)
			self.getStats()
			
			self.writeContig(self.hap1FinalWalk)

		return False

	def checkIsNodeSnBoundry(self,node):
		''' return true if node is a key in any of the snarl dictionaries (start/end)'''

		def addNode_advancePointer():
			''' add the handle to the walk, move pointer usint nextNode '''
			self.walkFromEnd_toStart.append(self.currentNode)
			self.hap1PathIds.append(self.graph.get_id(self.currentNode))
			self.graph.follow_edges(self.currentNode,False,self.nextNode)

		#bi-allelic bubble Is Navigated right here
		# don't do this - fix it to be like the rest of the sn types
			#TODO: 1 consider both sides and compare parental kmers?
			# 	   2 make these biAllelic boundry bits functions ^
		if self.graph.get_id(node) in list(self.sd['biAllelic_snarls'].keys()):
			print('\t bi allelic start',self.graph.get_id(node),int(self.sd['biAllelic_snarls'][int(self.graph.get_id(node))]['end']['end_id']) )
			# if self.sd['biAllelic_snarls'][int(self.graph.get_id(node))]['end']['end_id'] not in :
			# only for the start node? or check fwd degree count to make sure?
			self.calcHetAvgLen(node) 
			self.exit = self.sd['biAllelic_snarls'][int(self.graph.get_id(node))]['end']['end_id']
			# walk forward until we land on the opposide sn boundry node
			# while int(self.graph.get_id(self.currentNode)) != int(self.sd['biAllelic_snarls'][int(self.graph.get_id(node))]['end']['end_id']):
			# 	# TODO use parent info to assign hap
			# 	addNode_advancePointer()
			# move to that boundry node
			# addNode_advancePointer()
			return 'biAllelic_snarls'
		elif self.graph.get_id(node) in list(self.sd['biAllelic_snarlsEnds'].keys()):
			self.exit = self.sd['biAllelic_snarlsEnds'][int(self.graph.get_id(node))]['start']['start_id']
			print('\t bi allelic end',self.graph.get_id(node), 'start',self.exit) 
			# Since we move backward once check direction and move accordingly
			if self.goingBackward:
				print('going from end to start of biAlic')
				self.calcHetAvgLen(node)
				# while int(self.graph.get_id(self.currentNode)) != self.sd['biAllelic_snarlsEnds'][int(self.graph.get_id(node))]['start']['start_id']:
				# 	# this is where I use parent info to assign hap
				# 	addNode_advancePointer()
				# # addNode_advancePointer()
			else:
				#we could hit end of bi allelic even though we are going forward
				if self.sd['biAllelic_snarlsEnds'][int(self.graph.get_id(node))]['start']['start_id'] in self.hap1PathIds:
					print('going out of biAlic')
					addNode_advancePointer()
					return False
				# else:
				# 	while int(self.graph.get_id(self.currentNode)) != self.sd['biAllelic_snarlsEnds'][int(self.graph.get_id(node))]['start']['start_id']:
				# 		#TODO this is where I use parent info to assign hap
				# 		## THIS IS MESSING UP CUNK 12 I THINK WHEN SNARL END IS ALSO BI ALLELIC END
				# 		addNode_advancePointer()

				# 	addNode_advancePointer()
			return 'biAllelic_snarlEnds'

		#DAG set the end of the DAG and return to avoid the tip
		if self.graph.get_id(node) in list(self.sd['directed_acyclic_net_graph'].keys()):
			if self.sd['directed_acyclic_net_graph'][self.graph.get_id(node)]['unvisited']:
				self.exit = self.sd['directed_acyclic_net_graph'][self.graph.get_id(node)]['end']['end_id']
				print('\t DAG start',self.graph.get_id(node),'DAG end:',self.exit )
				return 'directed_acyclic_net_graph'
		elif self.graph.get_id(node) in list(self.sd['directed_acyclic_net_graphEnds'].keys()):
			if self.sd['directed_acyclic_net_graphEnds'][self.graph.get_id(node)]['unvisited']:
				self.exit = self.sd['directed_acyclic_net_graphEnds'][self.graph.get_id(node)]['start']['start_id']
				print('\t DAG end',self.graph.get_id(node),'DAG exit:',self.exit)
				return 'directed_acyclic_net_graphEnds'
		#snarls 
		if self.graph.get_id(node) in list(self.sd['snarls'].keys()): 
			if self.sd['snarls'][self.graph.get_id(node)]['unvisited']:
				print('\t sn start',self.graph.get_id(node) )
				self.exit = self.sd['snarls'][self.graph.get_id(node)]['end']['end_id']
				return 'snarls'
		# sometimes we can hit the 'end' of a snarl first
		elif self.graph.get_id(node) in list(self.sd['snarlsEnds'].keys()):
			if self.sd['snarlsEnds'][self.graph.get_id(node)]['unvisited']:
				print('\t sn end',self.graph.get_id(node) )
				self.exit = self.sd['snarlsEnds'][self.graph.get_id(node)]['start']['start_id']
				return 'snarlsEnds'

		return False

	def biallelicWalk(self, llSnDict):
		''' walks across the bi allelic bubble
			TODO add parent kmer aligns and all
		'''
		biallelicWalk = []
		# boundry node
		biallelicWalk.append(self.currentNode)
		# bi-allelic node
		self.graph.follow_edges(self.currentNode,False,self.nextNode)
		biallelicWalk.append(self.currentNode)
		# other boundry node, unless it is a in/del bi allelic like in chunk_25
		if self.graph.get_id(self.currentNode) != self.exit:
			self.graph.follow_edges(self.currentNode,False,self.chooseNextNode)
			biallelicWalk.append(self.currentNode)

		if self.graph.get_id(self.currentNode) != self.exit:
			print('\n\tsomething wrong with biAllelic?', [self.graph.get_id(b) for b in biallelicWalk])
		else:
			for d in biallelicWalk:
				print( self.graph.get_id(d) )
				self.walkFromEnd_toStart.append(d)
				self.hap1PathIds.append(self.graph.get_id(d))
					# 		#TODO this is where I use parent info to assign hap
					# 		## THIS IS MESSING UP CUNK 12 I THINK WHEN SNARL END IS ALSO BI ALLELIC END
					# 		addNode_advancePointer()
		return

	def avoidTip(self, llSnDict):
		''' move through a DAG avoiding walking into the tip '''
		dagWalk = []	
		dagWalkLens = []
		self.calcTipLen(self.currentNode)
		self.sd[llSnDict][self.graph.get_id(self.currentNode)]['unvisited']=False
		# send start self reachable to see if we will need to flip the orientation to look for start
		self.fwdReachableNodes(self.currentNode,self.sd[llSnDict][self.graph.get_id(self.currentNode)]['start_self_reachable'] )
		# do this to avoid DAG tips when exit node is preset
		endInList = [ l for l in self.lis if self.exit == self.graph.get_id(l)] 
		print('endinList',[self.graph.get_id(l) for l in endInList], self.exit)
		# if it is a tip the list will have just the exit node in it
		if len(endInList) ==1:
			dagWalk.append(self.currentNode)
			# move to one of the nodes in lis 
			self.graph.follow_edges(self.currentNode,False,self.chooseNextNode)
			print('\tavoiding the tip, going to:',self.graph.get_id(self.currentNode))
			
		print('\tBoundry node:',self.graph.get_id(self.currentNode), 'endInList',self.exit,len(endInList),'lis',[self.graph.get_id(i) for i in self.lis])
		
		while len(endInList) < 1:
			#TODO calc tip len or whatever we have
			# if self.lis > 1 : check degree of both and get tip len
			print('\t dagNode:',self.graph.get_id(self.currentNode),'dagWalk',[self.graph.get_id(i) for i in dagWalk])
			# stop after 2 visits to any node
			if dagWalk.count(self.currentNode) >2:#< 3:
				print('node in dagWalk too many times',[self.graph.get_id(i) for i in dagWalk])
				prevId = 0
				for i in dagWalk:
					if i == prevId:
						dagWalk.remove(i)
					prevId = i
				if i == prevId:
					dagWalk.remove(i)
				print('removed suquential repeats',[self.graph.get_id(i) for i in dagWalk])
				break
			# if this is the second time we have passed this node in a dag, it is an inversion
			elif dagWalk.count(self.currentNode) == 2:
				print('second visit, flip node orientation looking for end')
				dagWalk.append(self.currentNode)
					# move to one of the nodes in lis 
				self.graph.follow_edges(self.currentNode,False,self.nextNode)
				self.fwdReachableNodes(self.currentNode,True)
				endInList = [ l for l in self.lis if self.exit == self.graph.get_id(l)]	
			elif len(self.lis)==0:
				print('are we out of graph ?')
				break			
			else:

				dagWalk.append(self.currentNode)
				dagWalkLens.append(self.graph.get_length(self.currentNode))
				print('add ',self.graph.get_id(self.currentNode),'to dagWalk')
				# move to one of the nodes in lis 
				self.graph.follow_edges(self.currentNode,False,self.nextNode)
				self.fwdReachableNodes(self.currentNode,True)
				if len(self.lis) ==0:
					if len(dagWalk)<2:
						print('flip the node and follow edges')
						self.graph.follow_edges(self.graph.flip(self.currentNode),False, self.nextNode)
						self.fwdReachableNodes(self.currentNode,True)
				endInList = [ l for l in self.lis if self.exit == self.graph.get_id(l)]
				print(self.graph.get_id(self.currentNode),'endInList',self.exit,len(endInList),'lis',[self.graph.get_id(i) for i in self.lis])
				if len(endInList) == 1:
					dagWalk.append(self.currentNode)


		# If we haven't reached it, go to end and append the end of the dag		
		if self.exit not in self.hap1PathIds:
			print('go to the end. exit',self.exit)
			self.graph.follow_edges(self.currentNode,False,self.chooseNextNode)
			print('are we at the end of the dag?',self.graph.get_id(self.currentNode))
			if self.graph.get_id(self.currentNode) == self.exit:
				dagWalk.append(self.currentNode)
			# move currentNode pointer to the next reacheable node
			# self.graph.follow_edges(self.currentNode,False,self.nextNode)
		else:
			print('just keep going')
			self.graph.follow_edges(self.currentNode,False,self.nextNode)
		
		for d in dagWalk:
			print( self.graph.get_id(d) )
			self.walkFromEnd_toStart.append(d)
			self.hap1PathIds.append(self.graph.get_id(d))
		dagWalkLens.append(np.sum(np.array(dagWalkLens)))

		print('final dag Walk',[self.graph.get_id(i) for i in dagWalk])
		return True

	def snarlWalk(self, llSnDict):
		''' move through a DAG avoiding walking into the tip '''
		def moveToUnvisitedNode(node):
			print('choosing nextNode',[self.graph.get_id(i) for i in self.chosenSnarlNodes], 'node',self.graph.get_id(node), )
			chosenSnarlIds = [self.graph.get_id(i) for i in self.chosenSnarlNodes]
			if node in self.chosenSnarlNodes:
				self.currentNode=node
				print('moved to:',self.graph.get_id(self.currentNode))
				# self.exit=0
				return False
			else: 
				# flipedHandleNextN=self.get_handle()
				return True

		snarlWalk = []	
		snarlWalkIds = []
		snarlWalkLens = []
		self.sd[llSnDict][self.graph.get_id(self.currentNode)]['unvisited']=False
		self.fwdReachableNodes(self.currentNode,self.sd[llSnDict][self.graph.get_id(self.currentNode)])
		# do this to avoid DAG tips when exit node is preset
		endInList = [ l for l in self.lis if self.exit == self.graph.get_id(l)] 
		# if it is a tip the list will have just the exit node in it
		if len(endInList)==1:
			snarlWalk.append(self.currentNode)
			snarlWalkIds.append(self.graph.get_id(self.currentNode))
			# move to the exit node using chooseNextNode
			self.graph.follow_edges(self.currentNode,False,self.chooseNextNode)
			print('\tcan reach the end, going to:',self.graph.get_id(self.currentNode))
			
		print('\tBoundry node:',self.graph.get_id(self.currentNode), 'endInList',self.exit,len(endInList),'lis',[self.graph.get_id(i) for i in self.lis])

		print('exit',self.exit)
		while len(endInList)<1:
			# TODO limit the # of times a node is in the list
			if snarlWalk.count(self.currentNode) > 2:
				print('node in snarlWalk too many times',[self.graph.get_id(i) for i in snarlWalk])
				prevId = 0
				for i in snarlWalk:
					if i == prevId:
						snarlWalk.remove(i)
					prevId = i
				if i == prevId:
					snarlWalk.remove(i)
				print('removed suquential repeats',[self.graph.get_id(i) for i in snarlWalk])
				break

			elif snarlWalkIds.count(self.graph.get_id(self.currentNode)) > 2:
				print('node in snarlWalk too many times',[self.graph.get_id(i) for i in snarlWalk])

				break

			snarlWalk.append(self.currentNode)
			snarlWalkIds.append(self.graph.get_id(self.currentNode))
			print('snarl Walk',snarlWalkIds)
			self.fwdReachableNodes(self.currentNode,False)
			self.chosenSnarlNodes = [l for l in self.lis if l not in snarlWalk]
			if len(self.chosenSnarlNodes) > 0:
				# self.chosenSnarlNode=self.graph.get_id(unusedNodes[-1])
				# print('chosen node', self.chosenSnarlNode)
				prvNode = self.currentNode
				self.graph.follow_edges(self.currentNode,False,moveToUnvisitedNode)

				if self.currentNode==prvNode:
					self.graph.follow_edges(self.currentNode,False,self.nextNode)
				# self.fwdReachableNodes(self.currentNode,False)
				# snarlWalk.append(self.currentNode)
				# snarlWalkIds.append(self.graph.get_id(self.currentNode))
			else:
				self.graph.follow_edges(self.currentNode,False,self.nextNode)
				# self.fwdReachableNodes(self.currentNode,False)
				# snarlWalk.append(self.currentNode)
				# snarlWalkIds.append(self.graph.get_id(self.currentNode))
			self.fwdReachableNodes(self.currentNode,False)
			endInList = [ l for l in self.lis if self.exit == self.graph.get_id(l)]
			print(self.graph.get_id(self.currentNode),'endInList',self.exit,len(endInList),'lis',[self.graph.get_id(i) for i in self.lis])

		# If we haven't reached it, go to end and append the end of the dag		
		if self.exit not in [self.hap1PathIds+snarlWalkIds]:
			print('go to the end. exit',self.exit)
			if self.currentNode not in snarlWalk and self.graph.get_id(self.currentNode) != self.exit:
				snarlWalk.append(self.currentNode)
				snarlWalkIds.append(self.graph.get_id(self.currentNode))
			self.graph.follow_edges(self.currentNode,False,self.chooseNextNode)
			print('are we at the end of the snarl?',self.exit,self.graph.get_id(self.currentNode), self.graph.get_id(self.currentNode)==self.exit)
			print('exit',self.exit)
			if self.graph.get_id(self.currentNode) == self.exit:
				snarlWalk.append(self.currentNode)
				snarlWalkIds.append(self.graph.get_id(self.currentNode))
			# move currentNode pointer to the next reacheable node
			# self.graph.follow_edges(self.currentNode,False,self.nextNode)
		else:
			print('just keep going')
			self.graph.follow_edges(self.currentNode,False,self.nextNode)

		for s in snarlWalk:
			print( self.graph.get_id(s) )
			self.walkFromEnd_toStart.append(s)
			self.hap1PathIds.append(self.graph.get_id(s))

	def navigateSnarl(self, llSnDict):
		''' call snarl walk functions
			TODO maybe i don't need to check if it is the 'ends' dict now that I use self.exit
		 '''
		# DAGs
		if llSnDict[0] == 'd':
			if llSnDict[-4:] == 'Ends':	
				print('navigate DAG end',self.goingBackward)
				self.avoidTip(llSnDict)
				print('is this the right dag?',list(self.sd[llSnDict[:-4]].keys()) )
				self.sd[llSnDict[:-4]][self.exit]['unvisited']=False
			else:
				# run avoid tip
				self.avoidTip(llSnDict)
				print('is this the right dag?',self.sd[str(llSnDict+'Ends')][self.exit])
				self.sd[str(llSnDict+'Ends')][self.exit]['unvisited']=False
		# lower level snarls
		elif llSnDict[0] == 's':
			if llSnDict[-4:] == 'Ends':	
				if self.exit in self.hap1PathIds:
					print('follow edges')
					self.graph.follow_edges(self.currentNode,False,self.nextNode)
				else:
					print('going to walk snarl')
					self.snarlWalk(llSnDict)
					# self.sd[llSnDict[:-4]][self.exit]['unvisited']=False
			else:
				print(llSnDict)
				self.snarlWalk(llSnDict)
				self.sd[str(llSnDict+'Ends')][self.exit]['unvisited']=False
		# biallelic bubbles 
		elif llSnDict[0] == 'b':
			self.biallelicWalk(llSnDict)
		return

	def findNextTLs_Start(self):
		''' go from the end of a top level snarl to either: 
				1) the next top level snarl start
				2) the end of the graph
			passing through all other snarls along the way '''
		def walkInbetweenTopLevelSnarls():
			'''Add handles to the path list from the End of one snarl until 
				the start of another top level snarl is reached or there is no more graph'''
			unVisitedTopSnrls=[x for x in list(self.sd['parentSnarls'].keys()) if self.sd['parentSnarls'][x]['unvisited']]
			unVisitedTopSnrlsEnds=[x for x in list(self.sd['parentSnarlsEnds'].keys()) if self.sd['parentSnarlsEnds'][x]['unvisited']]
			print('unVisitedTopSns',unVisitedTopSnrls)
			# while currentNode isn't a top level snarl key
			while True:
				print('walk inbetween now Node:',self.graph.get_id(self.currentNode), 'degree sum',self.graph.get_degree(self.currentNode,False) + self.graph.get_degree(self.currentNode,True))
				# if we hit a start of a next top level snarl stop 
				if self.graph.get_id(self.currentNode) in unVisitedTopSnrls:
					# append the start node of the next top level snarl for varifacation later
					self.walkFromEnd_toStart.append(self.currentNode)
					self.hap1PathIds.append(self.graph.get_id(self.currentNode))
					if len(unVisitedTopSnrls)>1:
						self.sd['parentSnarls'][int(self.graph.get_id(self.currentNode))]['unvisited']=False
					print('found next ParentStart',self.hap1PathIds)
					break
				if self.graph.get_id(self.currentNode) in unVisitedTopSnrlsEnds:
					# append the start node of the next top level snarl for varifacation later
					self.walkFromEnd_toStart.append(self.currentNode)
					self.hap1PathIds.append(self.graph.get_id(self.currentNode))
					self.sd['parentSnarlsEnds'][int(self.graph.get_id(self.currentNode))]['unvisited']=False
					print('found next ParentEnd',self.graph.get_id(self.currentNode))
					break

				# if this is a snarl boundry
				#check snarl dicts
				snBndryDict = self.checkIsNodeSnBoundry(self.currentNode)
				if snBndryDict:
					# these boundries will be the small bluntified nodes - I don't think I want to skew the distribution like that. # self.homNodeLens.append(self.graph.get_length(self.currentNode)) maybe add it to it's own len distribution
					self.bluntifiedNodeLens.append(self.graph.get_length(self.currentNode))
					print('\twhich dict?',snBndryDict)
					if self.exit not in self.hap1PathIds:
						print('\tnavigate snarl')
						self.navigateSnarl(snBndryDict)

					else:
						print('\talready went through this',self.graph.get_id(self.currentNode))
						# go through hom graph portion record length and append to walk 
						self.homNodeLens.append(self.graph.get_length(self.currentNode))
						self.walkFromEnd_toStart.append(self.currentNode)
						self.hap1PathIds.append(self.graph.get_id(self.currentNode))
						# move currentNode pointer to the next reacheable node
						self.graph.follow_edges(self.currentNode,False,self.nextNode)
						# fill list with reachable nodes
						self.fwdReachableNodes(self.currentNode,False)						
					# break
				
				# if this is the end of the graph stop
				if (self.graph.get_degree(self.currentNode,False) + self.graph.get_degree(self.currentNode,True) ) <2:
					print('at the end', self.graph.get_id(self.currentNode))
					if self.currentNode not in self.walkFromEnd_toStart:
						self.homNodeLens.append(self.graph.get_length(self.currentNode))
						self.walkFromEnd_toStart.append(self.currentNode)
						self.hap1PathIds.append(self.graph.get_id(self.currentNode))
					break
				
				
				else:
					if self.hap1PathIds.count(self.graph.get_id(self.currentNode)) < 2:
						# go through hom graph portion record length and append to walk 
						self.homNodeLens.append(self.graph.get_length(self.currentNode))
						self.walkFromEnd_toStart.append(self.currentNode)
						self.hap1PathIds.append(self.graph.get_id(self.currentNode))
						# move currentNode pointer to the next reacheable node
						self.graph.follow_edges(self.currentNode,False,self.nextNode)
						# fill list with reachable nodes
						self.fwdReachableNodes(self.currentNode,False)
						# this python handle is broken: the method looks for a handle as the 2nd argument instead of a boolean. Slacked Adam: 6/14
						# print(self.graph.traverse_edge_handle((self.currentNode,self.lis[0]),False))
					else:
						break
			
			# append the walk from tl snarl end to the next tl snarl start, or EOG: end of graph
			self.hap1Path.append(self.walkFromEnd_toStart)
			self.hap1PathIds.append(self.graph.get_id(self.currentNode))
			print('\nnow',self.graph.get_id(self.currentNode),'lis:',[self.graph.get_id(l) for l in self.lis], self.hap1PathIds)
			print('this is the still walking hap1Path:\nlen hap1 path',len(self.hap1Path))
			for hap in self.hap1Path:
				print('walk')
				for h in hap:
					print(self.graph.get_id(h))
			self.writeOutPathCSV(self.hap1Path,False)
		
		print('\nfindNextTLs_Start here \t',self.graph.get_id(self.currentNode))
		self.fwdReachableNodes(self.currentNode,True)
		print( 'can go to',[ self.graph.get_id(f) for f in self.lis] )	
		self.walkFromEnd_toStart=[]
		self.walkFromEnd_toStart.append(self.currentNode)
		self.hap1PathIds.append(self.graph.get_id(self.currentNode))
		# when going the opposite direction from start - check that we step onto a unwallked node 
		nodesNOTInWalk = [ l for l in self.lis if self.graph.get_id(l) not in self.hap1PathIds]
		print('nodesNOTInWalk',[self.graph.get_id(n) for n in nodesNOTInWalk])
		if len(self.lis) == nodesNOTInWalk:
			self.graph.follow_edges(self.currentNode,False,self.nextNode)
		elif len(nodesNOTInWalk) > 0:
				self.exit = self.graph.get_id(nodesNOTInWalk[0])
				self.graph.follow_edges(self.currentNode,False,self.chooseNextNode)
		else:
			print('so are there no untouched nodes here?', [self.graph.get_id(l) for l in self.lis])

		# if len(self.lis) < 1:
		# 	print('no more graph')
		# 	self.walkFromEnd_toStart=[]
		# 	self.walkFromEnd_toStart.append(self.currentNode)
		# 	self.hap1Path.append(self.walkFromEnd_toStart)
		# 	self.hap1PathIds.append(self.graph.get_id(self.currentNode))
		# 	return True
		
			
		walkInbetweenTopLevelSnarls()
		return True

	def chainParents(self):
		''' create a chain of top level parent snarls (TLs) in the forward (left -> right) orientation by walking regions inbetween the top level snarls.
		From a list of walks: [ [TLs_1End->TLs_2St],[TLs_2End->EOG_1],[TLs_1St->EOG_2] ] 
		Order the walks ( reversing the [TLs_1St->EOG_2] walk ), 
		TODO: figure out the snarl walk, and output the contig.
		'''
		if len(list(self.sd['parentSnarlsEnds'].keys())) > 0: 
			# go until next parent start ( even if it is the current parent's own start)
			for eNode in list(self.sd['parentSnarlsEnds'].keys()):
				# we might end up at a end node instead of a start... then we neeed to figure out how to start - maybe check current node dict?
				# store endNode and set it as visited
				self.tlEndNode = int(eNode)#list(self.sd['parentSnarlsEnds'].keys())[0])
				tlstartNode = self.sd['parentSnarlsEnds'][self.tlEndNode]
				print('\n start chainng at: end',self.tlEndNode,'start',self.sd['parentSnarlsEnds'][self.tlEndNode]['start']['start_id']) 
				
				# set the top level boundary node to visited
				self.sd['parentSnarlsEnds'][self.tlEndNode]['unvisited']=False
				# set graph pointer to the End of the snarl and orient it according to snarl file
				# self.graph.get_handle(handle,'backward'=False) => handle in 'forward' orientation
				self.currentNode=self.graph.get_handle(self.tlEndNode,self.sd['parentSnarlsEnds'][self.tlEndNode]['end']['backward'])
				
				# walk to next TLs
				self.findNextTLs_Start()
		# after all the snarls are traveled then go find the one Start node that is unvisited and go backward ('left')
		self.goOppDirFromLeftmostTopLevelSnStart()
		
		#TODO order the walks so that it is eog -> -> eog in hap1FinalWalk

		print('THIS IS THE hap1Path:\n')
		for hap in self.hap1Path:
			print('walk')
			for h in hap:
				print(self.graph.get_id(h))
				if self.hap1FinalWalk.count(h) == 0:
					self.hap1FinalWalk.append(h)
					self.hap1_pathLenInBps.append(self.graph.get_length(h))
				else:
					self.hap1_doubledHandles.append(h)

		print('\n hap1_pathLenInBps:',self.hap1_pathLenInBps)

		

		# output the path to a csv with appropriate colors for bandage
		self.writeOutPathCSV(self.hap1Path,False)
		self.getStats()
		self.writeContig(self.hap1FinalWalk)
		print('\ntipLens',self.dagTipLens,'\nhomLens',self.homNodeLens,'\nhetLens',self.hetNodeAvgLens,'\nbluntedNodes',self.bluntifiedNodeLens,'\nwalkLens',self.dagWalkLens,'\ntipLens',self.dagTipLens	)
		print(np.sum(np.array(self.hap1_pathLenInBps)), (np.sum(np.array(self.hap1_pathLenInBps))/self.graph.get_total_length()))


def main():
	''' create snarl and graphPath objects
		run object methods to deseralize graph and traverse it 
		output: 2 haplotype paths ''' 
	inputGraph, inputSnarlFilename = parseArguments()
	# create an instance of the snarl manager obj and input the snarl file into dictionaries
	# TODO run this inside of graphPathFinder
	sm = SnarlManager(inputSnarlFilename)
	allSnarlDict = sm.snarlsInput()

	# pass the snarl dict and input graph to the graphPathManager object
	gpf = graphPathFinder(allSnarlDict,inputGraph)
	gpf.deserializeGraph()


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

	# starts from one known node
	# visit_handle(graph.get_handle(50373))#195))


if __name__ == "__main__":

    main()

