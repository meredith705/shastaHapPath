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

	def snarlsInput(self):
		'''load snarls.json into their respective snarl type dictionaries'''
		snarl_json=[json.loads(line) for line in open(self.inputJson,'r')]
		keys=[]
		for sn in snarl_json:
		# bi-allelic snarls are type: 1
			if sn.get('parent',False):
				# print('parent',sn.get('parent'),'Parent Start',sn.get('parent').get('start'),sn.get('parent').get('end'))
				if int(sn.get('parent').get('start').get('node_id')) in list(self.sd['parentSnarls'].keys()):
				# 	# print(parentSnarls[int(sn.get('parent').get('start').get('node_id'))]['child'])
					self.sd['parentSnarls'][int(sn.get('parent').get('start').get('node_id'))]['child'].append((int(sn.get('start').get('node_id')),int(sn.get('end').get('node_id')),sn.get('end').get('backward', False)))
				else:
					self.sd['parentSnarls'][int(sn.get('parent').get('start').get('node_id'))] = {\
						'unvisited' : True, \
						'start': { 'start_id' : int(sn.get('parent').get('start').get('node_id')),\
										'backward': sn.get('parent').get('start').get('backward',False)}, \
						'end': {'end_id' : sn.get('parent').get('end').get('node_id'), \
										'backward' : sn.get('parent').get('end').get('backward', False)}, \
						'child' : [(int(sn.get('start').get('node_id')),int(sn.get('end').get('node_id')),sn.get('end').get('backward', False))]}
					self.sd['parentSnarlsEnds'][int(sn.get('parent').get('end').get('node_id'))] = {\
						'unvisited' : True, \
						'start': {'start_id' : int(sn.get('parent').get('start').get('node_id')),\
										'backward': sn.get('parent').get('start').get('backward',False)}, \
						'end': {'end_id' : sn.get('parent').get('end').get('node_id'), \
										'backward' : sn.get('parent').get('end').get('backward', False)}, \
						'child' : [(int(sn.get('start').get('node_id')),int(sn.get('end').get('node_id')),sn.get('end').get('backward', False))]}
					self.children.append(int(sn.get('start').get('node_id')))
			elif sn.get('type',False):
				if sn.get('type')==1:
					keys.append(int(sn.get('start').get('node_id')))
					keys.append(int(sn.get('end').get('node_id')))
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
				keys.append(int(sn.get('start').get('node_id')))
				keys.append(int(sn.get('end').get('node_id')))
				self.sd['directed_acyclic_net_graph'][int(sn.get('start').get('node_id'))] = {\
					'unvisited' : True, \
					'backward': sn.get('start').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward', False)},\
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'child' : sn.get('parent',False)}
				self.sd['directed_acyclic_net_graphEnds'][int(sn.get('end').get('node_id'))] = {\
					'unvisited' : True, \
					'backward': sn.get('end').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')), 'backward': sn.get('start').get('backward',False)}, \
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'child' : sn.get('parent',False)}
			else:
				keys.append(int(sn.get('start').get('node_id')))
				keys.append(int(sn.get('end').get('node_id')))
				self.sd['snarls'][int(sn.get('start').get('node_id'))] = {\
					'unvisited' : True, \
					'backward': sn.get('start').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward', False)},\
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'child' : sn.get('parent',False)}
				self.sd['snarlsEnds'][int(sn.get('end').get('node_id'))] = {\
					'unvisited' : True, \
					# 'backward': sn.get('end').get('backward', False),\
					'start' : { 'start_id' : int(sn.get('start').get('node_id')),'backward': sn.get('start').get('backward',False)},\
					'end' : {'end_id' : int(sn.get('end').get('node_id')),'backward' : sn.get('end').get('backward', False)}, \
					'child' : sn.get('parent',False)}


		for p in list(self.sd['parentSnarls'].keys()):
			if p in self.children:
				print('PARENT IS ALSO CHILD',p,self.sd['parentSnarls'][int(p)])
				self.sd['parentSnarlsEnds'].pop(int(self.sd['parentSnarls'][int(p)]['end']['end_id']))
				self.sd['parentSnarls'].pop(p)
			if int(p) in keys:
				print('p in keys')
				for s in ['snarls','snarlsEnds','directed_acyclic_net_graph','directed_acyclic_net_graphEnds','biAllelic_snarls','biAllelic_snarlsEnds']:
					print('pop',self.sd[s].pop(int(p),s))
					print('pop',self.sd[s].pop(int(self.sd['parentSnarls'][p]['end']['end_id']),s))		

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
		# boolean flag for the final backward traversal
		self.goingBackward = False
		# walks inbetween top level snrls
		self.walkFromEnd_toStart=[]
		# snarls
		self.whichSnDic=''
		self.exit=0

		# stats
		self.homNodeLens=[]			# list of hom nodes inbetween snarl boundries
		self.hetNodeAvgLens=[]      # avg len of pairs of biallelic nodes
		self.bluntifiedNodeLens=[]	# snarl boundry nodes, artifact of bluntification
		# self.biAllelicLengths=[]
		self.dagWalkLens=[]			# DAG walk lengths (>3 node DAGs)
		self.dagTipLens=[]			# lens of DAG tips (3 node DAGs)
		# reachable node list
		self.lis=[]

		#########################
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
		self.dagEnd=0
		self.dagDict={}
		self.dagSide=''
		

	def writeOutPathCSV(self,p):
		# open the path.csv and write out the path with colors for bandage viewing 
		outpathCSV=open(self.inputGraph[:-3]+".path.csv",'w')
		outpathCSV.write("Node,Path,Color\n")
		# parent1='hg03'
		parent2='hg04'
		# parent1Color = "light sea green"
		parent2Color = 'tomato'
		for walk in p:
			for n in walk:
				outpathCSV.write( str(self.graph.get_id(n)) +",M,"+parent2Color+"\n")
		outpathCSV.close()

	# TODO make it a list of handles that are orientated 
	def writeContig(self,hap,contigCount,primary):
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

	def getStats(self,plen,pri_nodes,alt_nodes,dbls_pri,dbls_alt):
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

	def calcHetAvgLen(self,node):
		''' get length of both sides of a biallelic bubble and rt avg'''
		self.fwdReachableNodes(node,False)
		ls = np.array([self.graph.get_length(l) for l in self.lis])
		print('biallelic mean',np.mean(ls))
		self.hetNodeAvgLens.append(np.mean(ls))
	
	def calcTipLen(self,node):
		''' check degree and find len of tip '''
		def itIsATip(tNode):
			if self.graph.get_degree(tNode,False)==0:
				print(self.graph.get_id(tNode),'is a tip')
				return True
			elif self.graph.get_degree(tNode,True)==0:
				print(self.graph.get_id(tNode),'is a tip')
				return True
			else:
				return False
		self.fwdReachableNodes(node,False)
		ls = np.array([self.graph.get_length(l) for l in self.lis if itIsATip(l)])
		for l in ls:
			print('tip Len',l)
			self.dagTipLens.append(l)

	def deserializeGraph(self):
		# deserialize the packed graph
		self.graph.deserialize(self.inputGraph)
		# Q? make a path in the graph instead of the list of graph_ids? 
		# A: no probably too slow? also I'll be altering the origional gfa file 
		# path = graph.create_path_handle("hap1path")
		self.chainParents()
	
	def fwdReachableNodes(self,other,bothEnds):
		''' follow edges in fwd direction adding nodes to the lis'''
		def add2List(y):
			self.lis.append(y)
			return True
		self.lis.clear()	
		self.graph.follow_edges(other,False,add2List) 
		if bothEnds:
			self.graph.follow_edges(other,True,add2List) 

	def goOppDirFromLeftmostTopLevelSnStart(self):
		'''go 'left' from most 'left' start snarl - the one that hasn't been visited yet'''
		unVisitedTopSns=[x for x in list(self.sd['parentSnarls'].keys()) if self.sd['parentSnarls'][x]['unvisited']]
		print('unVisitedTopSns',unVisitedTopSns)
		# set the start node to visited and flip the directional flag
		#TODO change this incase ther are no parent snarls 
		if len(unVisitedTopSns)>0:
			self.sd['parentSnarls'][unVisitedTopSns[0]]['unvisited']=False
			self.goingBackward=True

		#flip the only unvisited start node and go from there till end
		if len(unVisitedTopSns)>0:
			handleDir = True
			if self.sd['parentSnarls'][unVisitedTopSns[0]]['start']['backward']:
				handleDir = False
			# set the graph handle to the opposite of whatever the snarls file told us origionally
			self.currentNode=self.graph.get_handle(unVisitedTopSns[0],handleDir)
			# get left side of graph
			self.findNextTLs_Start()

	def nextNode(self,nextN):
		'''moves pointer to first reachable node in the graph'''
		print('nextNode?',self.graph.get_id(nextN))
		self.currentNode=nextN
		return False # rt true didn't stop the flow at degree>1

	def chooseNextNode(self,nextN):
		'''moves pointer to each reachable node until the chosen self.exit node is reached'''
		print('nextNode?', int(self.graph.get_id(nextN)) == self.exit)
		if self.graph.get_id(nextN) == self.exit:
			self.currentNode=nextN
			self.exit=0
			return False
		else: 
			return True

	def checkIsNodeSnBoundry(self,node):
		''' return true if node is a key in any of the snarl dictionaries (start/end)'''

		def addNode_advancePointer():
			''' add the handle to the walk, move pointer usint nextNode '''
			self.walkFromEnd_toStart.append(self.currentNode)
			self.hap1PathIds.append(self.graph.get_id(self.currentNode))
			self.graph.follow_edges(self.currentNode,False,self.nextNode)

		#snarls
		if self.graph.get_id(node) in list(self.sd['snarls'].keys()): 
			print('\t sn start',self.graph.get_id(node) )
			self.exit = self.sd['snarls'][self.graph.get_id(node)]['end']['end_id']
			return 'snarls'
		# sometimes we can hit the 'end' of a snarl first
		elif self.graph.get_id(node) in list(self.sd['snarlsEnds'].keys()):
			print('\t sn end',self.graph.get_id(node) )
			self.exit = self.sd['snarls'][self.graph.get_id(node)]['start']['start_id']
			return 'snarlsEnds'

		#bi-allelic bubble Is Navigated right here
			#TODO: 1 consider both sides and compare parental kmers?
			# 	   2 make these biAllelic boundry bits functions ^
		elif self.graph.get_id(node) in list(self.sd['biAllelic_snarls'].keys()):
			print('\t bi allelic start',self.graph.get_id(node),int(self.sd['biAllelic_snarls'][int(self.graph.get_id(node))]['end']['end_id']) )
			# if self.sd['biAllelic_snarls'][int(self.graph.get_id(node))]['end']['end_id'] not in :
			# only for the start node? or check fwd degree count to make sure?
			self.calcHetAvgLen(node) 
			while int(self.graph.get_id(self.currentNode)) != int(self.sd['biAllelic_snarls'][int(self.graph.get_id(node))]['end']['end_id']):
				# TODO use parent info to assign hap
				addNode_advancePointer()
			addNode_advancePointer()
			return 'biAllelic_snarls'
		elif self.graph.get_id(node) in list(self.sd['biAllelic_snarlsEnds'].keys()):
			print('\t bi allelic end',self.graph.get_id(node) )
			# Since we move backward once check direction and move accordingly
			if self.goingBackward:
				print('going from end to start of biAlic')
				self.calcHetAvgLen(node)
				while int(self.graph.get_id(self.currentNode)) != self.sd['biAllelic_snarlsEnds'][int(self.graph.get_id(node))]['start']['start_id']:
					# this is where I use parent info to assign hap
					addNode_advancePointer()
				addNode_advancePointer()
			else:
				#we could hit end of bi allelic even though we are going forward
				if self.sd['biAllelic_snarlsEnds'][int(self.graph.get_id(node))]['start']['start_id'] in self.hap1PathIds:
					print('going out of biAlic')
					addNode_advancePointer()
					return False
				else:
					while int(self.graph.get_id(self.currentNode)) != self.sd['biAllelic_snarlsEnds'][int(self.graph.get_id(node))]['start']['start_id']:
						#TODO this is where I use parent info to assign hap
						addNode_advancePointer()
					addNode_advancePointer()
			return 'biAllelic_snarlEnds'

		#DAG proceed until we find the end of the DAG avoing the tip
		# TODO make DAG toy
			''' Current Errors chunk_4.pg
			nextNode? False
	         DAG start 36949
			ok now which dict? directed_acyclic_net_graph
			endInList [<bdsg.handlegraph.handle_t object at 0x7fc817ade330>]
			go to the end. dagWalk: [<bdsg.handlegraph.handle_t object at 0x7fc817ade2b0>]
			nextNode? False
			         DAG start 36949
			ok now which dict? directed_acyclic_net_graph
			endInList [<bdsg.handlegraph.handle_t object at 0x7fc817ade370>]
			go to the end. dagWalk: [<bdsg.handlegraph.handle_t object at 0x7fc817ade2b0>]
			nextNode? False
			         DAG start 36949
			'''
		elif self.graph.get_id(node) in list(self.sd['directed_acyclic_net_graph'].keys()):
			print('\t DAG start',self.graph.get_id(node) )
			self.exit = self.sd['directed_acyclic_net_graph'][self.graph.get_id(node)]['end']['end_id']
			return 'directed_acyclic_net_graph'
		elif self.graph.get_id(node) in list(self.sd['directed_acyclic_net_graphEnds'].keys()):
			print('\t DAG end',self.graph.get_id(node))
			self.exit = self.sd['directed_acyclic_net_graphEnds'][self.graph.get_id(node)]['start']['start_id']
			return 'directed_acyclic_net_graphEnds'
		else:
			return False

	def avoidTip(self, llSnDict):
		''' move through a DAG avoiding walking into the tip '''
		# dagEnd = self.sd[llSnDict][int(self.graph.get_id(self.currentNode))]['end']
		# datStart = self.sd[llSnDict][int(self.graph.get_id(self.currentNode))]['start']
		dagWalk = []	
		dagWalkLens = []
		self.calcTipLen(self.currentNode)
		self.fwdReachableNodes(self.currentNode,True)
		# do this to avoid DAG tips when exit node is preset
		endInList = [ l for l in self.lis if self.exit == self.graph.get_id(l)] 
		print('endInList',self.exit,len(endInList))
		while len(endInList) < 1:
			#TODO calc tip len or whatever we have
			# if self.lis > 1 : check degree of both and get tip len
			dagWalkLens.append(self.graph.get_length(self.currentNode))
			print('\t dagNode:',self.graph.get_id(self.currentNode),'dagWalk',dagWalk)
			if dagWalk.count(self.currentNode) < 3:
				dagWalk.append(self.currentNode)
			else:
				print('node in dagWalk too many times',dagWalk)
				break
			# move to one of the nodes in lis 
			self.graph.follow_edges(self.currentNode,False,self.nextNode)
			self.fwdReachableNodes(self.currentNode,True)
			endInList = [ l for l in self.lis if self.exit == self.graph.get_id(l)]

		# append the node that connects to the end of dag
		dagWalk.append(self.currentNode)
		# go to end and append the end of the dag
		
		if self.exit not in self.hap1PathIds:
			print('go to the end. exit',self.exit)
			self.graph.follow_edges(self.currentNode,False,self.chooseNextNode)
			dagWalk.append(self.currentNode)
		else:
			print('just keep going')
			self.graph.follow_edges(self.currentNode,False,self.nextNode)
		
		for d in dagWalk:
			self.walkFromEnd_toStart.append(d)
			self.hap1PathIds.append(self.graph.get_id(d))
		dagWalkLens.append(np.sum(np.array(dagWalkLens)))



		return True

	def navigateSnarl(self, llSnDict):
		''' call snarl walk functions '''
		# DAGs
		if llSnDict[0] == 'd':
			if llSnDict[-4:] == 'Ends':	
				print('navigate DAG end',self.goingBackward)
				self.avoidTip(llSnDict)
			else:
				# run avoid tip

				self.avoidTip(llSnDict)
		# TODO: lower level snarls
		elif llSnDict[0] == 's':
			print(llSnDict)

	def findNextTLs_Start(self):
		''' go from the end of a top level snarl to either: 
				1) the next top level snarl start
				2) the end of the graph
			passing through all other snarls along the way '''
		def walkInbetweenTopLevelSnarls():
			'''Add handles to the path list from the End of one snarl until 
				the start of another top level snarl is reached or there is no more graph'''
			self.walkFromEnd_toStart=[]
			unVisitedTopSnrls=[x for x in list(self.sd['parentSnarls'].keys()) if self.sd['parentSnarls'][x]['unvisited']]
			print('unVisitedTopSns',unVisitedTopSnrls)
			# while currentNode isn't a top level snarl key
			while True:
				# if we hit a start of a next top level snarl stop 
				if self.graph.get_id(self.currentNode) in unVisitedTopSnrls:
					# append the start node of the next top level snarl for varifacation later
					self.walkFromEnd_toStart.append(self.currentNode)
					self.hap1PathIds.append(self.graph.get_id(self.currentNode))
					self.sd['parentSnarls'][int(self.graph.get_id(self.currentNode))]['unvisited']=False
					print('found next ParentStart')
					break
				
				#check snarl dicts
				snBndryDict = self.checkIsNodeSnBoundry(self.currentNode)
				if snBndryDict:
					# these boundries will be the small bluntified nodes - I'm not going to skew the distribution like that. # self.homNodeLens.append(self.graph.get_length(self.currentNode)) maybe add it to it's own len distribution
					self.bluntifiedNodeLens.append(self.graph.get_length(self.currentNode))
					print('ok now which dict?',snBndryDict)
					self.navigateSnarl(snBndryDict)
					# break
				# if this is the end of the graph stop
				elif self.graph.get_degree(self.currentNode,False)<1:
					self.homNodeLens.append(self.graph.get_length(self.currentNode))
					self.walkFromEnd_toStart.append(self.currentNode)
					self.hap1PathIds.append(self.graph.get_id(self.currentNode))
					break
				else:
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
			
			print('\nnow',self.graph.get_id(self.currentNode),'lis:',self.lis)
			# append the walk from tl snarl end to the next tl snarl start, or EOG: end of graph
			self.hap1Path.append(self.walkFromEnd_toStart)
			self.hap1PathIds.append(self.graph.get_id(self.currentNode))
			print('THIS IS THE hap1Path:\n')
			for hap in self.hap1Path:
				print('walk')
				for h in hap:
					print(self.graph.get_id(h))
			self.writeOutPathCSV(self.hap1Path)
		
		print('\nhere \t',self.graph.get_id(self.currentNode))
		self.fwdReachableNodes(self.currentNode,False)
		print('can go to',[ self.graph.get_id(f) for f in self.lis] )
			
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
				# store endNode and set it as visited
				self.tlEndNode = int(eNode)
				print('end',self.tlEndNode,'start',self.sd['parentSnarlsEnds'][self.tlEndNode]['start']['start_id']) 
				
				# set the top level boundary node to visited
				self.sd['parentSnarlsEnds'][self.tlEndNode]['unvisited']=False
				# set graph pointer to the End of the snarl and orient it according to snarl file
				# self.graph.get_handle(handle,'backward'=False) => handle in 'forward' orientation
				self.currentNode=self.graph.get_handle(self.tlEndNode,self.sd['parentSnarlsEnds'][self.tlEndNode]['end']['backward'])
				# walk to next TLs
				self.findNextTLs_Start()
		# after all the snarls are traveled then go find the one Start node that is unvisited and go backward ('left')
		self.goOppDirFromLeftmostTopLevelSnStart()
		# output the path to a csv with appropriate colors for bandage
		print('THIS IS THE hap1Path:\n')
		for hap in self.hap1Path:
			print('walk')
			for h in hap:
				print(self.graph.get_id(h))
		self.writeOutPathCSV(self.hap1Path)
		print('\ntipLens',self.dagTipLens,'\nhomLens',self.homNodeLens,'\nhetLens',self.hetNodeAvgLens,'\nbluntedNodes',self.bluntifiedNodeLens,'\nwalkLens',self.dagWalkLens,'\ntipLens',self.dagTipLens	)



def main():
	''' create snarl and graphPath objects
		run object methods to deseralize graph and traverse it 
		output: 2 haplotype paths ''' 
	inputGraph, inputJson = parseArguments()
	# create an instance of the snarl manager obj and input the snarl file into dictionaries
	sm = snrlManager(inputJson)
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