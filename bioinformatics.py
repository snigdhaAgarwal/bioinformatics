import sys
import collections
from operator import itemgetter
import decimal
from collections import namedtuple

def generateReadpair(genome, kmerLength, distance):
    pairs = []
    for i in range(len(genome)-2*kmerLength-distance+1):
        secondKmerLoc = i+kmerLength+distance
        pairs.append(genome[i:i+kmerLength]+'|' +
                     genome[secondKmerLoc:secondKmerLoc+kmerLength])
    return pairs

# genome = input()
# kmerLength = int(input())
# distance = int(input())
# pairs = generateReadpair(genome,kmerLength,distance)
# pairs.sort()
# print (pairs)


class Node:
    def __init__(self, value, outDegree, inDegree):
        self.value = value
        self.inDegree = inDegree
        self.outDegree = outDegree

    def __str__(self):
        return (self.value + ' ' + str(self.inDegree) + ' ' + str(self.outDegree))

class Edge:
    def __init__(self,fromNode,toNode,weight):
        self.fromNode = fromNode
        self.toNode = toNode
        self.weight = weight

class normalGraph:
    def __init__(self):
        self.edges = {}
        self.size = 0
        self.nodes = {}  # key = nodeValue, val = Node object

    def createOrUpdateNode(self, value, outDegree, inDegree):
        if value not in self.nodes:
            node = Node(value, outDegree, inDegree)
            self.nodes[value] = node
        else:
            self.updateNodes(value, outDegree, inDegree)

    def updateNodes(self, value, outDegree, inDegree):
        self.nodes[value].outDegree = self.nodes[value].outDegree + outDegree
        self.nodes[value].inDegree = self.nodes[value].inDegree + inDegree
        # if self.nodes[value].outDegree == 0 and self.nodes[value].inDegree == 0:
        #     del self.nodes[value] #TODO: FIX later fro debruijn graph

    def createEdge(self, fromNode, toNode):
        if fromNode not in self.edges:
            self.edges[fromNode] = []
        self.edges[fromNode].append(toNode) #Edge(fromNode,toNode,0))
        # increment count of edges
        self.size = self.size + 1
        # modify nodes
        self.createOrUpdateNode(fromNode, 1, 0)
        self.createOrUpdateNode(toNode, 0, 1)


    def removeEdge(self, fromNode, toNode):
        self.edges[fromNode].remove(toNode)
        if len(self.edges[fromNode]) == 0:
            del self.edges[fromNode]
        # decrement edges count
        self.size = self.size - 1
        # update nodes
        self.updateNodes(fromNode, -1, 0)
        self.updateNodes(toNode, 0, -1)

    def getTopologicalSortsRecursive(self, visitedNodes, topologies, resultSoFar):
        allVisited = True
        for index in range(len(visitedNodes)):
            if visitedNodes[index] == 0:
                allVisited = False

        if allVisited:
            topologies.append(resultSoFar)
            return

        for index in range(len(visitedNodes)):
            node = index + 1

            ''' Pick a node that is unvisited and has inDegree = 0'''
            if visitedNodes[index] == 0 and node in self.nodes and self.nodes[node].inDegree == 0:

                ''' Remove all outgoing edges from node and  mark currentNode as visited'''
                removedEdges = [] # storing removed edges for adding back later
                while node in self.edges:
                    toNode = self.edges[node][0]
                    self.removeEdge(node,toNode)
                    removedEdges.append((node,toNode))
                visitedNodes[index] = 1

                resultSoFar = resultSoFar + str(node)

                self.getTopologicalSortsRecursive(visitedNodes,topologies,resultSoFar)

                ''' Backtrack steps: Add all edges back, mark node as unvisited, remove node from result'''
                for pair in removedEdges:
                    self.createEdge(pair[0], pair[1])
                visitedNodes[index] = 0
                resultSoFar = resultSoFar[:-1]



class deBruijnGraph(normalGraph):
    def __init__(self):
        normalGraph.__init__(self)

    def getStartAndEndNode(self):
        startNode = None
        endNode = None
        for node in self.nodes.values():
            effectiveDegree = node.outDegree-node.inDegree
            if effectiveDegree > 0:
                startNode = node
            elif effectiveDegree < 0:
                endNode = node
        return startNode.value, endNode.value

    def oneInOneOutNode(self, value):
        if self.nodes[value].inDegree == 1 and self.nodes[value].outDegree == 1:
            return True
        return False

    def findOneInOneOutCycle(self, value):
        currentNode = self.edges[value][0]
        path = [value]
        while currentNode != value:
            path.append(currentNode)
            print(path)
            currentNode = self.edges[currentNode][0]
        return path

    def findCycle(self, startNode, unvisitedEdges):
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
        pathStart, pathEnd = self.getStartAndEndNode()
        self.createEdge(pathEnd, pathStart)
        unvisitedEdges = self.edges
        path = collections.deque([])
        start = pathStart
        while unvisitedEdges:
            cycle = self.findCycle(start, unvisitedEdges)
            # rotate existing path to start from start
            if path:
                amountToRotate = path.index(start)
                path.rotate(-amountToRotate)
            # remove last element from cycle to convert to path
            cycle.pop()
            path.extend(cycle)
            # finding node in cycle with unvisitedEdges
            for key in cycle:
                if key in unvisitedEdges:
                    start = key
                    break
        # rotating path to start from startPath
        amountToRotate = path.index(pathStart)
        path.rotate(-amountToRotate)
        return path

# Gene		path.append(key)rate text from kmers
# Logic: connect end and start node. Keep finding cycles till all edges covered.


def createGenomeFromKmers():
    graph = deBruijnGraph()
    for kmer in sys.stdin:
        graph.createEdge(kmer.strip()[:-1], kmer.strip()[1:])
    path = graph.findEulerianPath()
    genome = path[0]
    for i in range(1, len(path)):
        genome = genome + path[i][-1]
    return genome

# print (createGenomeFromKmers())


def populateGraph(graph, pairs):
    for pair in pairs:
        kmers = pair.split('|')
        firstKmer = kmers[0]
        secondKmer = kmers[1]
        fromNode = firstKmer[:-1]+'|'+secondKmer[:-1]  # prefixes
        toNode = firstKmer[1:]+'|'+secondKmer[1:]  # suffixes
        graph.createEdge(fromNode, toNode)


def isNextKmer(genome, kmerLength, kmer, pos):
    mismatchFound = False
    for i in range(kmerLength):
        if genome[pos + i] != '' and genome[pos+i] != kmer[i]:
            mismatchFound = True
    return not mismatchFound

# place 2 kmers d dist apart, and keep adding pair in next node if prefix of new kmer matches suffix of existing kmer


def createGenomeFromGraph(graph, kmerLength, distance):
    startNode, endNode = graph.getStartAndEndNode()
    genome = ['']*(graph.size+2*kmerLength+distance)  # deduced from example
    # adding start node
    genome[0:kmerLength] = startNode.split('|')[0]
    genome[kmerLength+distance:2*kmerLength+distance] = startNode.split('|')[1]

    iterator = 1
    currentNode = startNode
    while graph.edges:
        for toNode in graph.edges[currentNode]:
            kmers = toNode.split('|')
            firstKmer = kmers[0]
            secondKmer = kmers[1]
            if (isNextKmer(genome, kmerLength, firstKmer, iterator) and
                    isNextKmer(genome, kmerLength, secondKmer, iterator+kmerLength+distance)):
                genome[iterator:iterator+kmerLength] = firstKmer
                secondKmerLoc = iterator+kmerLength+distance
                genome[secondKmerLoc:secondKmerLoc+kmerLength] = secondKmer
                # prepare for next iteration
                graph.removeEdge(currentNode, toNode)
                currentNode = toNode
                iterator = iterator + 1
                break
    return genome


def createGenomeFromReadPairs(pairs, kmerLength, distance):
    graph = deBruijnGraph()
    populateGraph(graph, pairs)
    # Prefixes and suffixes are (k-1,d+1)-mers
    genome = createGenomeFromGraph(graph, kmerLength-1, distance+1)
    return genome


# line  = input().split(' ')
# kmerLength = int(line[0])
# distance = int(line[1])
# pairs = []
# for pair in sys.stdin:
    # pairs.append(pair.strip())
# print(''.join(createGenomeFromReadPairs(pairs,kmerLength,distance)))


# Generate contigs from a collection of imperfect coverage reads (IMPortant)
def generateContigs():
    graph = deBruijnGraph()
    contigs = []
    # need list for finding standalone 1in1out nodes for finding cycles
    # key is node and valye is boolean indicating if the node has been visited or not
    oneInOneOutList = {}

    for read in sys.stdin:
        graph.createEdge(read.strip()[:-1], read.strip()[1:])

    for node in graph.edges:
        # starting Node
        if not graph.oneInOneOutNode(node) and graph.nodes[node].outDegree > 0:
            # each tonode is a separate Path => separate contig
            for toNode in graph.edges[node]:
                currentContig = node + toNode[-1]
                currentNode = toNode
                while graph.oneInOneOutNode(currentNode):
                    oneInOneOutList[currentNode] = 1
                    currentNode = graph.edges[currentNode][0]
                    currentContig = currentContig + currentNode[-1]
                contigs.append(currentContig)
        elif node not in oneInOneOutList and graph.oneInOneOutNode(node):
            oneInOneOutList[node] = 0
    # find disconnected 1-in-1-out cycles
    for node in oneInOneOutList:
        if oneInOneOutList[node] == 0:
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
    # Read RNA codon translation file
    rnaCodonMap = {}
    f = open('RNA_codon_table_1.txt', 'r')
    lines = f.readlines()
    for line in lines:
        values = line.strip().split(' ')
        if len(values) > 1:
            rnaCodonMap[values[0]] = values[1]
        else:
            rnaCodonMap[values[0]] = ''
    return rnaCodonMap

# Convert RNA to Amino acid sequence


def convertRNAToAmino(rnaString):
    # Read RNA codon translation file
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

# Find DNA substrings that encode a given amino sequence


def getReverseComplement(dnaString):
    reverComplement = ''
    for i in range(len(dnaString)):
        char = dnaString[i]
        replaceChar = ''
        if char == 'A':
            replaceChar = 'T'
        elif char == 'T':
            replaceChar = 'A'
        elif char == 'G':
            replaceChar = 'C'
        elif char == 'C':
            replaceChar = 'G'
        reverComplement = reverComplement + replaceChar
    # returns reverse of string Eg. reverse complement of GGCCAT is AUGGCC
    return reverComplement[::-1]


def getSubstring(pos, rnaString, aminoString, rnaCodonMap):
    i = pos
    subString = ''
    aminoiterator = 0
    while aminoiterator < len(aminoString) and i < len(rnaString)-2:
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


def getSubstringList(dnaString, aminoString, rnaCodonMap):
    subStringList = []
    # Convert DNA to RNA
    rnaString = dnaString.replace("T", "U")
    i = 0
    while i < len(rnaString)-3*len(aminoString) + 1:
        codon = rnaString[i:i+3]
        if rnaCodonMap[codon] == aminoString[0]:
            subString = getSubstring(i, rnaString, aminoString, rnaCodonMap)
            if subString:
                subStringList.append(subString)
        i = i+1
    # Convert RNA strings to DNA strings
    dnaSubStringList = []
    for string in subStringList:
        dnaSubStringList.append(string.replace("U", "T"))
    return dnaSubStringList


def findDNASeqForAminoSeq(dnaString, aminoString):
    # Read RNA codon translation file
    rnaCodonMap = getRNACodonMap()
    dnaSubStringList = getSubstringList(dnaString, aminoString, rnaCodonMap)
    # Get reverse complement of dnaString
    reverseComplement = getReverseComplement(dnaString)
    reverseSubstringList = getSubstringList(
        reverseComplement, aminoString, rnaCodonMap)
    for string in reverseSubstringList:
        dnaSubStringList.append(getReverseComplement(string))
    return dnaSubStringList

# dnaString = input()
# aminoString = input()
# print(findDNASeqForAminoSeq(dnaString,aminoString))


# peptide = input()
# spectrum = generatetheoreticalSpectrum(peptide)
# print (*spectrum) # prints list space separated

# Finding how many peptide sequence exist with given mass (IMPortant) - DP
def countPeptideWithGivenMass(parentMass):
    peptideCounts = [0]*(parentMass+1)
    peptideCounts[0] = 1
    atomicMasses = [57, 71, 87, 97, 99, 101, 103, 113,
                    114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
    for possibleMass in range(atomicMasses[0], parentMass+1):
        for aminoMass in atomicMasses:
            if possibleMass >= aminoMass:
                peptideCounts[possibleMass] = peptideCounts[possibleMass] + \
                    peptideCounts[possibleMass - aminoMass]
                if possibleMass == 243 and peptideCounts[possibleMass - aminoMass]:
                    print(possibleMass-aminoMass, aminoMass)
    return peptideCounts[parentMass]
#
# parentMass = int(input())
# print(countPeptideWithGivenMass(parentMass))


# Calculating number of subpeptides in a linear peptide of length n
# from below, fastest solution is N(N+1)/2
# def getSubpeptideCount(length):
    # peptideCount = 1 # 1 empty string
    # for i in range(length):
    # for j in range(i,length):
    # peptideCount = peptideCount + 1
    # return peptideCount

# length = int(input())

class Peptide:

    def __init__(self, aminoList=[]):
        self.aminoList = aminoList
        self.mass = sum(aminoList)

    def __repr__(self):
        return str(self.aminoList)

# Sample Input: 0 113 128 186 241 299 314 427
# Sample output: 186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186


class cycloPeptideSeq:
    with open('resources/allAminoMasses') as f:
        aminoMasses = [int(word) for word in f.readline().split(',')]

    def getAtomicMasses(self):
        atomicMasses = {}
        # for v in self.aminoMasses:
        # 	atomicMasses[v] = []
        f = open('resources/integer_mass_table.txt', 'r')
        lines = f.readlines()
        for line in lines:
            values = line.split(' ')
            atomicMasses[values[0]] = int(values[1])
        return atomicMasses

    def __init__(self, experimentalSpectrum):
        self.spectrum = experimentalSpectrum
        self.parentMass = 1284#max(experimentalSpectrum)

    # Generating theoretical spectrum - collection of masses of all subpeptides
    def generatetheoreticalSpectrum(self, peptide, isCyclic):
        # Generate prefix masses
        prefixMass = [0]
        for i in range(len(peptide)):
            prefixMass.append(peptide[i]+prefixMass[i])
        # Generating theoretical spectrum
        peptideMass = prefixMass[-1]
        theoreticalSpectrum = [0]
        for i in range(len(prefixMass)-1):
            for j in range(i+1, len(prefixMass)):
                subPeptideMass = prefixMass[j] - prefixMass[i]
                theoreticalSpectrum.append(subPeptideMass)
                # only need to subtract middle combinations(without 1st n last elements) to get circular combinations
                if i > 0 and j < len(prefixMass)-1 and isCyclic:
                    # for cyclic cases eg. LEQN - EQ = mass of(NL)
                    theoreticalSpectrum.append(peptideMass - subPeptideMass)
        theoreticalSpectrum.sort()
        return theoreticalSpectrum

    # Append a single amino acid to end of existing peptides
    def expand(self, peptides):
        newPeptides = []
        for peptideMassPair in peptides:
            for amino in self.aminoMasses:
                existingAmino = peptideMassPair.aminoList
                newAmino = existingAmino + [amino]
                newPeptides.append(Peptide(aminoList=newAmino))
        return newPeptides

    def getSequence(self):
        validPeptides = [Peptide()]
        while validPeptides:
            expandedPeptides = self.expand(validPeptides)
            # Cant remove from validpeptides, if for loop on validpeptides
            validPeptides = list(expandedPeptides)
            for peptideMassPair in expandedPeptides:
                totalMass = peptideMassPair.mass
                peptideMassesList = peptideMassPair.aminoList
                # if found a peptide which can be fit for spectrum
                if totalMass == self.parentMass:
                    # peptideMassesList = list(map(int,peptide.split('-')))
                    theoreticalSpectrum = self.generatetheoreticalSpectrum(
                        peptideMassesList, True)
                    if theoreticalSpectrum == self.spectrum:
                        print(peptideMassesList, end='-')
                    validPeptides.remove(peptideMassPair)
                else:
                    # Check of theoretical spectrum of linear peptide is consistent with spectrum
                    # peptideMassesList = list(map(int,peptide.split('-')))
                    theoreticalSpectrum = self.generatetheoreticalSpectrum(
                        peptideMassesList, False)
                    inconsistent = 0
                    for mass in theoreticalSpectrum:
                        if mass not in self.spectrum:
                            inconsistent = 1
                    if inconsistent:
                        validPeptides.remove(peptideMassPair)

    def getLeaderboardSeq(self, leaderboardSize):
        '''Only topmost peptide with highest score is returned.
        At each step , leaderboard is expanded and trimmed to `leaderboardSize` keeping topmost `leaderboardSize` scores.'''

        leaderboard = [Peptide()]  # peptide,mass
        # will contain final peptide with highest score
        leaderPeptidePair = leaderboard[0]
        # contains all highest scored peptides
        leaderPeptides = [leaderPeptidePair]
        while len(leaderboard) > 0:
            expandedPeptides = self.expand(leaderboard)
            leaderboard = list(expandedPeptides)
            for peptideMassPair in expandedPeptides:
                mass = peptideMassPair.mass
                peptideMassesList = peptideMassPair.aminoList
                # if mass == self.parentMass:
                #     if self.getSeqScore(peptideMassesList, True) > self.getSeqScore(
                #             leaderPeptidePair.aminoList, True):
                #         leaderPeptidePair = peptideMassPair
                #         # If new leaderPeptide discovered empty leaderPeptidesList
                #         leaderPeptides.clear()
                #     # add all peptides with same score as leaderPeptide
                #     if self.getSeqScore(peptideMassesList, True) == self.getSeqScore(
                #             leaderPeptidePair.aminoList, True):
                #         leaderPeptides.append(peptideMassPair)
                # elif mass > self.parentMass:
                #     leaderboard.remove(peptideMassPair)

                ''' If exact peptide mass is not known'''
                if len(peptideMassesList)> 10:
                    leaderboard.remove(peptideMassPair)
                else:
                    if self.getSeqScore(peptideMassesList, True) > self.getSeqScore(
                            leaderPeptidePair.aminoList, True):
                        leaderPeptidePair = peptideMassPair
                        # If new leaderPeptide discovered empty leaderPeptidesList
                        leaderPeptides.clear()
                    if self.getSeqScore(peptideMassesList, True) == self.getSeqScore(
                            leaderPeptidePair.aminoList, True):
                        leaderPeptides.append(peptideMassPair)


            # BOunding step. Only keeping top N scored peptides
            leaderboard = self.trimLeaderboard(leaderboard, leaderboardSize)
        return leaderPeptides

    def trimLeaderboard(self, leaderboard, leaderboardSize):
        '''Bounding step where only top N scored peptides are retained. In the last position all peptides
        with that score are included.

        Logic is:

        ```
        Trim(Leaderboard, Spectrum, N, AminoAcid, AminoAcidMass)
        for j ← 1 to |Leaderboard|
            Peptide ← j-th peptide in Leaderboard
            LinearScores(j) ← LinearScore(Peptide, Spectrum)
        sort Leaderboard according to the decreasing order of scores in LinearScores
        sort LinearScores in decreasing order
        for j ← N + 1 to |Leaderboard|
            if LinearScores(j) < LinearScores(N)
                remove all peptides starting from the j-th peptide from Leaderboard
                return Leaderboard
        return Leaderboard
        ```
        '''

        peptideScoreList = []
        for peptideMassPair in leaderboard:
            peptideScoreList.append(
                (peptideMassPair, self.getSeqScore(peptideMassPair.aminoList, False)))
        sortedPeptides = list(
            sorted(peptideScoreList, key=itemgetter(1), reverse=True))

        lastIndex = -1
        for j in range(leaderboardSize, len(leaderboard)):
            if sortedPeptides[j][1] < sortedPeptides[leaderboardSize-1][1]:
                lastIndex = j
                break
        return [pmp for pmp, score in sortedPeptides[:lastIndex]]

    # converting Alphabetical peptide to masses list
    def getPeptideScore(self, peptide):
        aminoMasses = self.getAtomicMasses()
        peptideMassesList = []
        for i in range(len(peptide)):
            peptideMassesList.append(aminoMasses[peptide[i]])
        return self.getSeqScore(peptideMassesList, True)

    def getSeqScore(self, peptideMassesList, isCyclic=True):
        theoreticalSpectrum = self.generatetheoreticalSpectrum(
            peptideMassesList, isCyclic)
        score = 0
        for mass in self.spectrum:
            if mass in theoreticalSpectrum:
                theoreticalSpectrum.remove(mass)
                score = score + 1
        return score

    def getSpectrumConvolution(self):
        convolution = []
        for i in range(len(self.spectrum)):
            for j in range(i+1,len(self.spectrum)):
                difference = abs(self.spectrum[j] - self.spectrum[i])
                if difference!=0:
                    convolution.append(difference)
        return convolution

    def trimConvolution(self, convolutionSize, convolution):
        ''' Remove all values except between 57 and 200 and store the number of occurences of each'''
        filteredConvolution = {}
        for value in convolution:
            if value >= 57 and value <=200:
                if value in filteredConvolution:
                    filteredConvolution[value] = filteredConvolution[value] + 1
                else:
                    filteredConvolution[value] = 1

        ''' Sorting based on occurence. Using occurence as score keeping top N, including all values tied at last position'''
        sortedAminos = sorted(filteredConvolution.items(), key=lambda pair: (pair[1],pair[0]), reverse=True)
        lastIndex = -1
        for i in range(convolutionSize, len(sortedAminos)):
            if sortedAminos[i][1] < sortedAminos[convolutionSize-1][1]:
                lastIndex = i
                break
        aminoList = [amino for amino,count in sortedAminos[:lastIndex]]

        return aminoList

    def convolutionPeptideSequencing(self,convolutionSize, leaderBoardSize):
        ''' Generate convolution, trim it to required length,use it to perform leaderBoard seq '''
        convolution = self.getSpectrumConvolution()
        self.aminoMasses = self.trimConvolution(convolutionSize,convolution)
        self.aminoMasses = [99,128,113,147,97,186,137,163]#[99,128,113,147,97,114,163]#
        leaderPeptides = self.getLeaderboardSeq(leaderBoardSize)
        return leaderPeptides


def problem1():
    line = input()
    values = line.split(' ')
    experimentalSpectrum = []
    for value in values:
        experimentalSpectrum.append(int(value.strip()))
    instance1 = cycloPeptideSeq(experimentalSpectrum)
    instance1.getSequence()


def problem2():
    # Score cyclopeptides against experimental spectrum
    peptide = input().strip()
    line = input()
    experimentalSpectrum = []
    values = line.split(' ')
    for value in values:
        experimentalSpectrum.append(int(value.strip()))
    print("There are {0} values".format(len(values)))
    instance1 = cycloPeptideSeq(experimentalSpectrum)
    print(instance1.getPeptideScore(peptide))


def problem3():
    '''LeaderBoard Peptide Sequencing'''

    leaderboardSize = int(input())
    line = input()
    experimentalSpectrum = []
    values = line.split(' ')
    for value in values:
        experimentalSpectrum.append(int(value.strip()))
    instance1 = cycloPeptideSeq(experimentalSpectrum)
    leaderPeptides = instance1.getLeaderboardSeq(leaderboardSize)
    for pair in leaderPeptides:
        peptide = '-'.join([str(amino) for amino in pair.aminoList])
        print(peptide, end=' ')

def generateConvolution():
    ''' Generate the convolution of a spectrum= positive differences between all values in spectrum'''
    line = input()
    experimentalSpectrum = []
    values = line.split(' ')
    for value in values:
        experimentalSpectrum.append(int(value.strip()))
    instance1 = cycloPeptideSeq(experimentalSpectrum)
    convolution = instance1.getSpectrumConvolution()
    print(' '.join(str(x) for x in instance1.trimConvolution(20, convolution)))

def leaderBoardSeqWithConvolution():
    ''' Convolution is used to determine the possible amino acids making up the peptide
    to avoid other possible peptide solutions '''
    maxCharge = 1
    massOfCharge = decimal.Decimal(1.007)
    massDiscrepancy = decimal.Decimal(0.3)
    convolutionSize = int(input())
    leaderBoardSize = int(input())
    line = input()
    experimentalSpectrum = []
    values = line.split(' ')
    for value in values:
        realValue = maxCharge * decimal.Decimal(value.strip())-maxCharge*massOfCharge
        for i in range(int(realValue-massDiscrepancy), int(realValue+massDiscrepancy)+1):
            experimentalSpectrum.append(i)
    instance1 = cycloPeptideSeq(experimentalSpectrum)
    leaderPeptides = instance1.convolutionPeptideSequencing(convolutionSize,leaderBoardSize)
    for pair in leaderPeptides:
        peptide = ' '.join([str(int(amino)) for amino in pair.aminoList])
        print(peptide, end=',')

def numberOfWaysFromStartToEnd():
    rows = int(input())
    cols = int(input())
    dp = [[0 for x in range(cols)] for y in range(rows)]
    for i in range(rows):
        for j in range(cols):
            if i==0 or j==0:
                dp[i][j] = 1
            else:
                dp[i][j] = dp[i-1][j] + dp[i][j-1]
    print(dp[rows-1][cols-1])


''' Calculate minimum number of coins required to make a sum'''
def coinChangeProblem():
    coins = [int(x) for x in input().split(',')]
    totalValue = int(input())
    coinDP = [1000 for x in range(totalValue+1)] # contains min number of coins req to make all sums from 0 to totalValue
    for value in range(1,totalValue+1):
        if value in coins:
            coinDP[value] = 1
        else:
            for coin in coins:
                if coin < value:
                    coinDP[value] = min(coinDP[value], 1 + coinDP[value-coin])
    print(coinDP)


''' Get longest path length from 0,0 to row,col square '''
def manhattanTouristProblem():
    rows,cols = [int(x) for x in input().split(' ')]
    downWeights = [(0 for x in range(cols+1)) for y in range(rows)]
    line = input()
    rowNo = 0
    while '-' not in line:
        downWeights[rowNo] = [int(x) for x in line.split(' ')]
        rowNo = rowNo + 1
        line = input()
    rowNo = 0
    diagonalWeights = [(0 for x in range(cols)) for y in range(rows)]
    line = input()
    while '-' not in line:
        diagonalWeights[rowNo] = [int(x) for x in line.split(' ')]
        rowNo = rowNo + 1
        line = input()
    rowNo = 0
    rightWeights = [(0 for x in range(cols)) for y in range(rows+1)]
    for line in sys.stdin:
        rightWeights[rowNo] = [int(x) for x in line.split(' ')]
        rowNo = rowNo + 1
    maxWtMatrix = [[0 for x in range(cols+1)] for y in range(rows+1)]

    for i in range(1,rows+1):
        maxWtMatrix[i][0] = downWeights[i-1][0] + maxWtMatrix[i-1][0]
    for j in range(1,cols+1):
        maxWtMatrix[0][j] = rightWeights[0][j-1] + maxWtMatrix[0][j-1]
    for i in range(1,rows+1):
        for j in range(1,cols+1):
            maxWtMatrix[i][j] = max(maxWtMatrix[i-1][j]+downWeights[i-1][j], maxWtMatrix[i][j-1] +rightWeights[i][j-1],
                                    maxWtMatrix[i-1][j-1]+diagonalWeights[i-1][j-1])
    print(maxWtMatrix[rows][cols])


''' Count all topological orderings in a DAG'''
def count_topological_sorts():
    graph = normalGraph()
    graph.createEdge(1,2)
    graph.createEdge(2,3)
    graph.createEdge(2,5)
    graph.createEdge(2,8)
    graph.createEdge(3,4)
    graph.createEdge(3,7)
    graph.createEdge(5,6)

    visited = [0 for x in range(len(graph.nodes))]
    topologies = []
    graph.getTopologicalSortsRecursive(visited, topologies, '')
    print(len(topologies))

def LCS():
    string1 = input()
    string2 = input()
    LCSMatrix = [[(0, 'e') for x in range(len(string2) + 1)] for y in range(len(string1) + 1)]
    '''Leaving 0th row and col empty'''
    for i in range(1, len(string1) + 1):
        for j in range(1, len(string2) + 1):
            if string1[i - 1] == string2[j - 1]:
                LCSMatrix[i][j] = (LCSMatrix[i - 1][j - 1][0] + 1, '\\')
            else:
                maxValue = max(LCSMatrix[i - 1][j][0], LCSMatrix[i][j - 1][0])
                if maxValue == LCSMatrix[i - 1][j][0]:
                    LCSMatrix[i][j] = (LCSMatrix[i - 1][j][0], '^')
                elif maxValue == LCSMatrix[i][j - 1][0]:
                    LCSMatrix[i][j] = (LCSMatrix[i][j - 1][0], '<')

    '''
    To get LCS String start from bottom right corner and go according to direction symbols till reach 0,0
    reason to start from end: If start from 0,0 many dead ends along the way. But always straight path from end square
    if the direction is stored while creating grid
    '''

    i = len(string1)
    j = len(string2)
    LCSString = ''
    while i != 0 and j != 0:
        currentNode = LCSMatrix[i][j][1]
        if currentNode == '\\':
            LCSString = string1[i - 1] + LCSString
            i = i - 1
            j = j - 1
        elif currentNode == '^':
            i = i - 1
        else:
            j = j - 1

    print(LCSString)


''' 
Includes Global and local alignment solutions
'''
class LCSWithScoring:
    penalty_scores = []
    alphabetIndices = {}
    string1 = ''
    string2 = ''
    scoreTuple = namedtuple('scoreTuple','score direction')

    def __init__(self, penaltyResource):
        with open(penaltyResource) as f:
            alphabets = f.readline().split()
            for i in range(len(alphabets)):
                alphabetIndices[alphabets[i]] = i
            for line in f.readlines():
                values = [int(x) for x in line.split()]
                penalty_scores.append(values)

    def getPenaltyScores(self,character1,character2):
        index1 = self.alphabetIndices[character1]
        index2 = self.alphabetIndices[character2]
        return self.penalty_scores[index1][index2]

    ''' Indel and mismatch penalty along with different values for matching scenario. 
        Print the highest score and the related string alignment'''
    def GlobalLCSWithScoring(self):
        self.string1 = input()
        self.string2 = input()
        indelPenalty = -5

        LCSMatrix = [[self.scoreTuple(0, 'e')
                      for x in range(len(self.string2) + 1)]
                     for y in range(len(self.string1) + 1)]
        '''Filling first row'''
        for j in range(1,len(self.string2) + 1):
            LCSMatrix[0][j] = self.scoreTuple(LCSMatrix[0][j-1].score + indelPenalty,'--')
        ''' Filling first column'''
        for i in range(1,len(self.string1)+1):
            LCSMatrix[i][0] = self.scoreTuple(LCSMatrix[i-1][0].score + indelPenalty,'|')
        '''Filling rest of the matrix'''
        for i in range(1,len(self.string1) + 1):
            for j in range(1,len(self.string2) + 1):
                character1 = self.string1[i-1]
                character2 = self.string2[j-1]
                maxValue = -1000
                direction = 'e'

                diagonal = LCSMatrix[i-1][j-1].score + self.getPenaltyScores(character1,character2)
                if maxValue < diagonal:
                    maxValue = diagonal
                    direction = '\\'
                up = LCSMatrix[i-1][j].score + indelPenalty
                if maxValue < up:
                    maxValue = up
                    direction = '|'
                side = LCSMatrix[i][j-1].score + indelPenalty
                if maxValue < side :
                    maxValue = side
                    direction = '--'
                LCSMatrix[i][j] = self.scoreTuple(maxValue,direction)

        print(LCSMatrix[i][j].score)  #Printing score
        self.printStringAlignment(LCSMatrix)

    def printStringAlignment(self,LCSMatrix):
        shiftedString1 = ''
        shiftedString2 = ''
        i = len(self.string1)
        j = len(self.string2)
        while i != 0 or j != 0:
            currentDirection = LCSMatrix[i][j].direction
            if currentDirection == '\\':
                shiftedString1 = self.string1[i - 1] + shiftedString1
                shiftedString2 = self.string2[j-1] + shiftedString2
                i = i - 1
                j = j - 1
            elif currentDirection == '|':
                shiftedString1 = self.string1[i-1] + shiftedString1
                shiftedString2 = '-' + shiftedString2
                i = i - 1
            else:
                shiftedString1 = '-' + shiftedString1
                shiftedString2 = self.string2[j-1] + shiftedString2
                j = j - 1
        print(shiftedString1)
        print(shiftedString2)

    def localLCS(self):


''' 
    Starts at endNode and goes back recursively.
    for all edges leading to endNode:
        LongestPath(endNode) = max(EdgeWeight + LongestPath(fromNode))
    Intermediate results are stored in longestPath in form of :
        Node: longestPathTuple(LengthOfpath, Path)
'''
longestPathTuple = namedtuple('longestPathTuple','length path')
def longestPathDAGRecursive(startNode,endNode,adjMatrix, longestPath):
    if startNode == endNode or endNode not in adjMatrix:
        return
    max = -1
    for edgeTuple in adjMatrix[endNode]:
        fromNode = edgeTuple[0]
        edgeWeight = edgeTuple[1]
        if edgeWeight != -1:
            if fromNode not in longestPath:
                longestPathDAGRecursive(startNode,fromNode,adjMatrix,longestPath)

            if fromNode in longestPath:
                fromNodeTuple = longestPath[fromNode]
                if max < fromNodeTuple.length + edgeWeight:
                    max = fromNodeTuple.length + edgeWeight
                    endNodeTuple = longestPathTuple(max, fromNodeTuple.path + [endNode])
                    longestPath[endNode] = endNodeTuple


def longestPathDAG():
    startNode = int(input())
    endNode = int(input())
    adjMatrix = {} # toNode indexed matrix
    for line in sys.stdin:
        fromNode = int(line.split('->')[0])
        toNode = int(line.split('->')[1].split(':')[0])
        weight = int(line.split('->')[1].split(':')[1])
        if toNode in adjMatrix:
            adjMatrix[toNode].append((fromNode,weight))
        else:
            adjMatrix[toNode] = [(fromNode, weight)]

    longestPathDict = {startNode:longestPathTuple(0,[startNode])}
    longestPathDAGRecursive(startNode,endNode,adjMatrix, longestPathDict)

    print(longestPathDict[endNode].length)
    for node in longestPathDict[endNode].path:
        print(node, end='->')

def main():
    problem_name_to_fn_map = {
        'problem1': problem1,
        'problem2': problem2,
        'problem3': problem3
    }
    problem = LCSWithScoring('resources/BLOSUM62.txt').GlobalLCSWithScoring
    if len(sys.argv) > 1:
        problem = problem_name_to_fn_map[sys.argv[1]]
    problem()


if __name__ == '__main__':
    main()
