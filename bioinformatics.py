import sys
import collections
import os
import operator

def generateReadpair(genome,kmerLength,distance):
	pairs = []
	for i in range(len(genome)-2*kmerLength-distance+1):
		secondKmerLoc = i+kmerLength+distance
		pairs.append(genome[i:i+kmerLength]+'|'+genome[secondKmerLoc:secondKmerLoc+kmerLength])
	return pairs

# genome = input()
# kmerLength = int(input())
# distance = int(input())
# pairs = generateReadpair(genome,kmerLength,distance)
# pairs.sort()
# print (pairs)

class Node:
	def __init__(self,value,outDegree,inDegree):
		self.value = value
		self.inDegree = inDegree
		self.outDegree = outDegree
	
	def __str__(self):
		return (self.value + ' '+ str(self.inDegree) + ' ' + str(self.outDegree))
	
class deBruijnGraph:
	def __init__(self):
		self.edges = {}
		self.size = 0
		self.nodes = {} #key = nodeValue, val = Node object
	
	def createOrUpdateNode(self,value,outDegree,inDegree):
		if value not in self.nodes:
			node = Node(value,outDegree,inDegree)
			self.nodes[value]  = node
		else:
			self.updateNodes(value,outDegree,inDegree)
	
	def updateNodes(self,value,outDegree,inDegree):
		self.nodes[value].outDegree = self.nodes[value].outDegree + outDegree
		self.nodes[value].inDegree = self.nodes[value].inDegree + inDegree
		if self.nodes[value].outDegree==0 and self.nodes[value].inDegree==0 :
			del self.nodes[value]
		
	def createEdge(self,fromNode,toNode):
		if fromNode not in self.edges:
			self.edges[fromNode] = []
		self.edges[fromNode].append(toNode)			
		#increment count of edges
		self.size = self.size + 1
		#modify nodes
		self.createOrUpdateNode(fromNode,1,0)
		self.createOrUpdateNode(toNode,0,1)
		
	def removeEdge(self,fromNode,toNode):
		self.edges[fromNode].remove(toNode)
		if len(self.edges[fromNode])==0:
			del self.edges[fromNode]
		#decrement edges count
		self.size = self.size - 1
		#update nodes
		self.updateNodes(fromNode,-1,0)
		self.updateNodes(toNode,0,-1)
		
	def getStartAndEndNode(self):	
		startNode = None
		endNode = None
		for node in self.nodes.values():
			effectiveDegree = node.outDegree-node.inDegree
			if effectiveDegree > 0:
				startNode = node
			elif effectiveDegree < 0:
				endNode = node
		return startNode.value,endNode.value
	
	def oneInOneOutNode(self,value):
		if self.nodes[value].inDegree == 1 and self.nodes[value].outDegree == 1:
			return True
		return False
	
	def findOneInOneOutCycle(self,value):
		currentNode = self.edges[value][0]
		path = [value]
		while currentNode != value:
			path.append(currentNode)
			print(path)
			currentNode = self.edges[currentNode][0]
		return path
	
	def findCycle(self,startNode,unvisitedEdges):
		path = [startNode]
		prevNode = startNode
		currentNode = ''
		while currentNode != startNode:
			currentNode = unvisitedEdges[prevNode][0]
			path.append(currentNode)
			unvisitedEdges[prevNode].remove(currentNode)
			if not unvisitedEdges[prevNode]:
				del unvisitedEdges[prevNode]
			prevNode = currentNode
		return path
		
	
	def findEulerianPath(self):
		pathStart,pathEnd = self.getStartAndEndNode()
		self.createEdge(pathEnd,pathStart)
		unvisitedEdges = self.edges
		path = collections.deque([])
		start = pathStart
		while unvisitedEdges:
			cycle = self.findCycle(start,unvisitedEdges)
			#rotate existing path to start from start
			if path:
				amountToRotate = path.index(start)
				path.rotate(-amountToRotate)
			#remove last element from cycle to convert to path
			cycle.pop()
			path.extend(cycle)
			#finding node in cycle with unvisitedEdges
			for key in cycle:
				if key in unvisitedEdges:
					start = key
					break
		# rotating path to start from startPath
		amountToRotate = path.index(pathStart)
		path.rotate(-amountToRotate)
		return path
				
### Gene		path.append(key)rate text from kmers
### Logic: connect end and start node. Keep finding cycles till all edges covered.
def createGenomeFromKmers():
	graph = deBruijnGraph()
	for kmer in sys.stdin:
		graph.createEdge(kmer.strip()[:-1], kmer.strip()[1:])
	path = graph.findEulerianPath()
	genome = path[0]
	for i in range(1,len(path)):
		genome = genome + path[i][-1]
	return genome
	
# print (createGenomeFromKmers())
	
def populateGraph(graph,pairs):
	for pair in pairs:
		kmers = pair.split('|')
		firstKmer = kmers[0]
		secondKmer = kmers[1]
		fromNode = firstKmer[:-1]+'|'+secondKmer[:-1] # prefixes
		toNode = firstKmer[1:]+'|'+secondKmer[1:] # suffixes
		graph.createEdge(fromNode,toNode)		
	
def isNextKmer(genome,kmerLength,kmer,pos):
	mismatchFound = False
	for i in range(kmerLength):
		if genome[pos + i ] != '' and genome[pos+i] !=kmer[i]:
			mismatchFound = True
	return not mismatchFound
		
# place 2 kmers d dist apart, and keep adding pair in next node if prefix of new kmer matches suffix of existing kmer 
def createGenomeFromGraph(graph,kmerLength,distance):
	startNode,endNode = graph.getStartAndEndNode()
	genome = ['']*(graph.size+2*kmerLength+distance) #deduced from example
	#adding start node
	genome[0:kmerLength] = startNode.split('|')[0]
	genome[kmerLength+distance:2*kmerLength+distance] = startNode.split('|')[1]
	
	iterator = 1
	currentNode = startNode
	while graph.edges:
		for toNode in graph.edges[currentNode]:
			kmers = toNode.split('|')
			firstKmer = kmers[0]
			secondKmer = kmers[1]
			if (isNextKmer(genome,kmerLength,firstKmer,iterator) and 
			isNextKmer(genome,kmerLength, secondKmer,iterator+kmerLength+distance)):
				genome[iterator:iterator+kmerLength] = firstKmer
				secondKmerLoc = iterator+kmerLength+distance
				genome[secondKmerLoc:secondKmerLoc+kmerLength] = secondKmer
				#prepare for next iteration
				graph.removeEdge(currentNode,toNode)
				currentNode = toNode
				iterator = iterator + 1
				break
	return genome
	
def createGenomeFromReadPairs(pairs,kmerLength,distance):
	graph = deBruijnGraph()
	populateGraph(graph,pairs)
	genome = createGenomeFromGraph(graph,kmerLength-1,distance+1) #Prefixes and suffixes are (k-1,d+1)-mers
	return genome
		
			
# line  = input().split(' ')		
# kmerLength = int(line[0])
# distance = int(line[1])
# pairs = []
# for pair in sys.stdin:
	# pairs.append(pair.strip())
# print(''.join(createGenomeFromReadPairs(pairs,kmerLength,distance)))


### Generate contigs from a collection of imperfect coverage reads (IMPortant)
def generateContigs():
	graph = deBruijnGraph()
	contigs = []
	# need list for finding standalone 1in1out nodes for finding cycles
	oneInOneOutList = {} # key is node and valye is boolean indicating if the node has been visited or not
	
	for read in sys.stdin:
		graph.createEdge(read.strip()[:-1],read.strip()[1:])
		
	for node in graph.edges:
		if not graph.oneInOneOutNode(node) and graph.nodes[node].outDegree > 0 : #starting Node
			for toNode in graph.edges[node]: # each tonode is a separate Path => separate contig
				currentContig = node + toNode[-1]
				currentNode = toNode
				while graph.oneInOneOutNode(currentNode):
					oneInOneOutList[currentNode] = 1
					currentNode = graph.edges[currentNode][0]
					currentContig = currentContig + currentNode[-1]
				contigs.append(currentContig)
		elif node not in oneInOneOutList and graph.oneInOneOutNode(node):
			oneInOneOutList[node] = 0
	#find disconnected 1-in-1-out cycles
	for node in oneInOneOutList:
		if oneInOneOutList[node]==0:
			path = graph.findOneInOneOutCycle(node)
			currentContig = []
			for pathNode in path:
				currentContig.append(pathNode)
				oneInOneOutList[pathNode] = 1
			contigs.append(currentContig)
	return contigs

# contigs  = generateContigs()
# contigs.sort()
# print (contigs)

def getRNACodonMap():
	#Read RNA codon translation file
	rnaCodonMap = {}
	f= open('RNA_codon_table_1.txt','r')
	lines = f.readlines()
	for line in lines:
		values = line.strip().split(' ')
		if len(values) > 1:
			rnaCodonMap[values[0]] = values[1]
		else:
			rnaCodonMap[values[0]] = ''
	return rnaCodonMap
	
### Convert RNA to Amino acid sequence
def convertRNAToAmino(rnaString):
	#Read RNA codon translation file
	rnaCodonMap = getRNACodonMap()
	i = 0
	output = ''
	while i < len(rnaString):
		codon = rnaString[i:i+3]
		output = output + rnaCodonMap[codon]
		i = i+3
	return output

# rnaString = input()
# print(convertRNAToAmino(rnaString))

### Find DNA substrings that encode a given amino sequence
def getReverseComplement(dnaString):
	reverComplement = ''
	for i in range(len(dnaString)):
		char = dnaString[i]
		replaceChar = ''
		if char == 'A':
			replaceChar = 'T'
		elif char == 'T':
			replaceChar ='A'
		elif char == 'G':
			replaceChar ='C'
		elif char == 'C':
			replaceChar ='G'
		reverComplement= reverComplement + replaceChar
	return reverComplement[::-1] # returns reverse of string Eg. reverse complement of GGCCAT is AUGGCC
	

def getSubstring(pos,rnaString,aminoString,rnaCodonMap):
	i = pos
	subString = ''
	aminoiterator = 0
	while aminoiterator < len(aminoString) and i< len(rnaString)-2:
		codon = rnaString[i:i+3]
		if rnaCodonMap[codon] == aminoString[aminoiterator]:
			subString = subString + codon
			i = i + 3
			aminoiterator = aminoiterator + 1
		else:
			break
	if aminoiterator == len(aminoString):
		return subString
	else:
		return ''

def getSubstringList(dnaString,aminoString,rnaCodonMap):
	subStringList = []
	#Convert DNA to RNA
	rnaString = dnaString.replace("T","U")
	i =0
	while i< len(rnaString)-3*len(aminoString) + 1:
		codon = rnaString[i:i+3]
		if rnaCodonMap[codon] == aminoString[0]:
			subString = getSubstring(i,rnaString,aminoString,rnaCodonMap)
			if subString:
				subStringList.append(subString)
		i = i+1
	# Convert RNA strings to DNA strings
	dnaSubStringList= []
	for string in subStringList:
		dnaSubStringList.append(string.replace("U","T"))
	return dnaSubStringList
		
def findDNASeqForAminoSeq(dnaString,aminoString):
	#Read RNA codon translation file
	rnaCodonMap = getRNACodonMap()
	dnaSubStringList = getSubstringList(dnaString,aminoString,rnaCodonMap)
	#Get reverse complement of dnaString
	reverseComplement = getReverseComplement(dnaString)
	reverseSubstringList = getSubstringList(reverseComplement,aminoString,rnaCodonMap)
	for string in reverseSubstringList:
		dnaSubStringList.append(getReverseComplement(string))
	return dnaSubStringList

# dnaString = input()
# aminoString = input()
# print(findDNASeqForAminoSeq(dnaString,aminoString))


	

	
# peptide = input()
# spectrum = generatetheoreticalSpectrum(peptide)
# print (*spectrum) # prints list space separated

### Finding how many peptide sequence exist with given mass (IMPortant) - DP
def countPeptideWithGivenMass(parentMass):
	peptideCounts = [0]*(parentMass+1)
	peptideCounts[0] = 1
	atomicMasses = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]
	for possibleMass in range(atomicMasses[0],parentMass+1):
		for aminoMass in atomicMasses:
			if possibleMass >= aminoMass:
				peptideCounts[possibleMass] = peptideCounts[possibleMass] + peptideCounts[possibleMass - aminoMass]
				if possibleMass == 243 and peptideCounts[possibleMass - aminoMass]:
					print(possibleMass-aminoMass,aminoMass)
	return peptideCounts[parentMass]
#
# parentMass = int(input())
# print(countPeptideWithGivenMass(parentMass))


### Calculating number of subpeptides in a linear peptide of length n
# from below, fastest solution is N(N+1)/2
# def getSubpeptideCount(length):
	# peptideCount = 1 # 1 empty string
	# for i in range(length):
		# for j in range(i,length):
			# peptideCount = peptideCount + 1
	# return peptideCount

# length = int(input())


### Sample Input: 0 113 128 186 241 299 314 427
### Sample output: 186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186
class cycloPeptideSeq:
	aminoMasses = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]
	
	def getAtomicMasses(self):
		atomicMasses = {}
		# for v in self.aminoMasses:
		# 	atomicMasses[v] = []
		f = open('integer_mass_table.txt','r')
		lines = f.readlines()
		for line in lines:
			values = line.split(' ')
			atomicMasses[values[0]] = int(values[1])
		return atomicMasses

	def __init__(self,experimentalSpectrum):
		self.spectrum = experimentalSpectrum
		self.parentMass = max(experimentalSpectrum)
		self.validAminos = [] # stores only those aminoAcids that are in spectrum. Used when expanding peptidelist
		for amino in cycloPeptideSeq.aminoMasses:
			if amino in self.spectrum:
				self.validAminos.append(amino)

				
	### Generating theoretical spectrum - collection of masses of all subpeptides
	def generatetheoreticalSpectrum(self,peptide,isCyclic):
		# Generate prefix masses
		prefixMass = [0]
		for i in range(len(peptide)):
			prefixMass.append(peptide[i]+prefixMass[i])
		#Generating theoretical spectrum
		peptideMass = prefixMass[-1]
		theoreticalSpectrum = [0]
		for i in range(len(prefixMass)-1):
			for j in range(i+1,len(prefixMass)):
				subPeptideMass = prefixMass[j]-prefixMass[i]
				theoreticalSpectrum.append(subPeptideMass)
				if i> 0 and j<len(prefixMass)-1 and isCyclic:  # pnly need to subtract middle combinations(without 1st n last elements) to get circular combinations
					theoreticalSpectrum.append(peptideMass - subPeptideMass) # for cyclic cases eg. LEQN - EQ = mass of(NL)
		theoreticalSpectrum.sort()
		return theoreticalSpectrum
	
	# Append a single amino acid to end of existing peptides
	def expand(self,peptides):
		newPeptides = []
		for peptideMassPair in peptides:
			for amino in self.validAminos: # only append those amino acids that are in spectrum
				existingAmino = peptideMassPair[0]
				existingTotal = peptideMassPair[1]
				newAmino = str(amino)
				if existingAmino: #for first step when existing amino is empty
					newAmino = existingAmino+'-'+str(amino)
				newPeptides.append((newAmino,existingTotal+amino))
		return newPeptides
		
	def getSequence(self):
		validPeptides = [('',0)]
		while validPeptides:
			expandedPeptides = self.expand(validPeptides)
			validPeptides = list(expandedPeptides) # Cant remove from validpeptides, if for on validpeptides
			for peptideMassPair in expandedPeptides:
				totalMass = peptideMassPair[1]
				peptide = peptideMassPair[0]
				# if found a peptide which can be fit for spectrum
				if totalMass == self.parentMass:
					peptideMassesList = list(map(int,peptide.split('-')))
					theoreticalSpectrum = self.generatetheoreticalSpectrum(peptideMassesList,1) 
					if theoreticalSpectrum == self.spectrum:
						print (peptide)
					validPeptides.remove(peptideMassPair)
				else: 
					#Check of theoretical spectrum of linear peptide is consistent with spectrum
					peptideMassesList = list(map(int,peptide.split('-')))
					theoreticalSpectrum = self.generatetheoreticalSpectrum(peptideMassesList,0) 
					inconsistent = 0
					for mass in theoreticalSpectrum:
						if mass not in self.spectrum:
							inconsistent = 1
					if inconsistent:
						validPeptides.remove(peptideMassPair)	

	def getSeqScore(self,peptide):
		aminoMasses = self.getAtomicMasses()
		peptideMassesList = []
		for i in range(len(peptide)):
			peptideMassesList.append(aminoMasses[peptide[i]])
		theoreticalSpectrum = self.generatetheoreticalSpectrum(peptideMassesList,1)
		score  = 0
		for mass in self.spectrum:
			if mass in theoreticalSpectrum:
				theoreticalSpectrum.remove(mass)
				score = score + 1
		return score
		
				
# line = input()
# values = line.split(' ')
# experimentalSpectrum = []
# for value in values:
	# experimentalSpectrum.append(int(value.strip()))
# instance1 = cycloPeptideSeq(experimentalSpectrum)	
# instance1.getSequence()		


### Score cyclopeptides against experimental spectrum
peptide = input().strip()
line = input()
experimentalSpectrum = []
values = line.split(' ')
for value in values:
	experimentalSpectrum.append(int(value.strip()))
instance1 = cycloPeptideSeq(experimentalSpectrum)
print (instance1.getSeqScore(peptide))
