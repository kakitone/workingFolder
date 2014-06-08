import os
from itertools import groupby
from operator import itemgetter

###################################################### Helper Functions 
def extractMumData(folderName, fileName):
    #"Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
    f = open(folderName+fileName,'r')
    dataList = []
    
    for i in range(6):
        tmp = f.readline()

    while len(tmp) > 0:
        info = tmp.split('|')
        filterArr =  info[1].split()
        rdGpArr = info[-1].split('\t')
        firstArr = info[0].split()
        
        matchLenArr = info[2].split()
    
        matchLen1 = int(matchLenArr[0])
        matchLen2 = int(matchLenArr[1])    
        percentMatch = float(info[3])
        
        
        helperStart, helperEnd =  int( firstArr[0]), int( firstArr[1])
        readStart, readEnd =  int(filterArr[0]) , int(filterArr[1])
        
        helperName = rdGpArr[0].rstrip().lstrip()
        readName = rdGpArr[1].rstrip().lstrip()
        
        dataList.append([helperStart, helperEnd , readStart, readEnd,matchLen1,matchLen2,percentMatch,helperName,readName ])
    
        tmp = f.readline().rstrip()
                
    f.close()
    
    return dataList


def obtainLength(folderName, fileName):
    f = open(folderName+fileName, 'r')
    tmp = f.readline().rstrip()
    lenDic = {}
    tmplen = 0 
    tmpName = ""
    
    while len(tmp) > 0:
        
        if tmp[0] == '>':
            if tmplen != 0:
                lenDic[tmpName] =tmplen
                tmplen = 0
            tmpName = tmp[1:]
        else:
            tmplen += len(tmp)
        tmp  = f.readline().rstrip()
        
    lenDic[tmpName] =tmplen
    

    f.close()
    
    return lenDic
    
    
def putListToFileO(folderName, sourceFileName, targetFileName, myList):
    f = open(folderName + targetFileName+".txt", 'w')
    for eachitem in myList:
        f.write(eachitem+'\n')
        
    f.close()
    
    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "+folderName+targetFileName+".txt "+folderName+sourceFileName+" > "+folderName+targetFileName+".fasta"
    os.system(command)


def transformCoor(dataList):
    #"Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
    newList = []
    
    for eachitem in dataList:
        if eachitem[2] < eachitem[3]:
            newList.append(eachitem)
        else:
            tmpitem = eachitem
            tmp = tmpitem[2]
            tmpitem[2] = tmpitem[3]
            tmpitem[3] = tmp
            newList.append(tmpitem)
    
    return newList


def useMummerAlign(mummerLink, folderName, outputName, referenceName, queryName):
    command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName + outputName + " " + folderName+ referenceName +" "+ folderName+ queryName
    os.system(command)
    
    command  = mummerLink +"show-coords -r "+folderName+outputName+".delta > "+folderName+outputName+"Out"
    os.system(command)




def writeToFile_Double1(folderName, fileName1, fileName2, option = "contig"):

    f2 = open(folderName + fileName2, 'w')
    fOriginal = open(folderName + fileName1, 'r')
    
    readSet = []
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpRead) >0:
                readSet.append(tmpRead)
                tmpRead = ""
        else:
            tmpRead = tmpRead+ tmp
            
        tmp = fOriginal.readline().rstrip()
    readSet.append(tmpRead)
    
    print "len(readSet)", len(readSet)
    
    fOriginal.close()
    
    if option == "contig":
        header = ">Contig"
    else:
        header = ">Read"
    for eachcontig, dum in zip(readSet, range(len(readSet))):
        f2.write(header+ str(dum)+"_p\n")
        f2.write(eachcontig+'\n') 
        f2.write(header+ str(dum)+"_d\n")
        f2.write(reverseComplement(eachcontig)+'\n')
        
    f2.close()
    
def reverseComplement(myStr):
    myNewStr = myStr[::-1]
    myNewStr2 = ""
    for i in range(len(myNewStr)):
        if myNewStr[i] == 'A' or myNewStr[i] == 'a':
            myNewStr2 += 'T'
            
        elif  myNewStr[i] == 'T' or myNewStr[i] == 't':
            myNewStr2 +='A'
            
        elif  myNewStr[i] == 'C' or myNewStr[i] == 'c':
            myNewStr2 += 'G'
            
        elif  myNewStr[i] == 'G' or myNewStr[i] == 'g':
            myNewStr2 += 'C'
            
    return myNewStr2

def readConnectList(folderName, fileName):
    f = open(folderName+ fileName , 'r')
    tmp = f.readline().rstrip()
    
    dataList = []
    while len(tmp) > 0:
        index, connector, overlap = tmp.split(',')
        dataList.append([int(index), int(connector), int(overlap)])
        tmp = f.readline().rstrip()
        
    connectorList = [ [-1, -1] for i in range(len(dataList))]    
    for eachitem in dataList:
        connectorList[eachitem[0]] = [eachitem[1], eachitem[2]]
        
    f.close()

    return connectorList


def nameInEdgeList(name, myList):
    
    haveInserted = False
    for eachitem in myList:
        if eachitem[0] == name:
            haveInserted = True
    return haveInserted
                

def removeItem(myList,myname ):
    newList = []
    for eachitem in myList:
        if eachitem[0] != myname:
            newList.append(eachitem)
    #print len(newList)- len(myList)
    return newList




def quastEvaluate(folderName, quastLink, originalName, improvedNameList, referenceName ):
    
    # ./quast.py ~/git/myFinisher/finishingTool/S_cerivisea/contigs.fasta  ~/git/myFinisher/finishingTool/S_cerivisea/improved.fasta -R ~/git/myFinisher/finishingTool/S_cerivisea/reference.fasta
    header = quastLink + "quast.py"+ " "
    originalContigPath = folderName + originalName +" "
    improvedContigPath = "" 
    for eachname in improvedNameList:
        improvedContigPath = improvedContigPath +folderName+  eachname + " "
        
    
    referencePath = "-R "+folderName +referenceName +" "
    
    command = header + originalContigPath +improvedContigPath +referencePath 
    
    os.system(command)
    
    
    command = "cp "+"quast_results/latest/report.txt " + folderName + "assemblyAssessment.txt"
    os.system(command)

def compareGraphUnitTest(G, G2):
    assert(len(G.graphNodesList) == len(G2.graphNodesList))
    for index in range(len(G.graphNodesList)):
        assert(G.graphNodesList[index].nodeIndex == G2.graphNodesList[index].nodeIndex)
        assert(G.graphNodesList[index].nodeIndexList == G2.graphNodesList[index].nodeIndexList)
        assert(G.graphNodesList[index].overlapList == G2.graphNodesList[index].overlapList)
        assert( G.graphNodesList[index].listOfPrevNodes == G2.graphNodesList[index].listOfPrevNodes)
        assert( G.graphNodesList[index].listOfNextNodes == G2.graphNodesList[index].listOfNextNodes)
        assert(G.graphNodesList[index].visited == G2.graphNodesList[index].visited)


def loadContigsFromFile(folderName, fileName):
    f = open(folderName+fileName, 'r')
    tmp = f.readline().rstrip()
    dataDic = {}
    tmpSeq = ""
    tmpName = ""
    
    while len(tmp) > 0:
        
        if tmp[0] == '>':
            if len(tmpSeq) != 0:
                dataDic[tmpName] =tmpSeq
                tmpSeq = ""
            tmpName = tmp[1:]
        else:
            tmpSeq += tmp
        tmp  = f.readline().rstrip()
        
    dataDic[tmpName] =tmpSeq
    

    f.close()
    return dataDic


class seqGraphNode(object):
    def __init__(self, nodeIndex):
        self.nodeIndex = nodeIndex
        self.nodeIndexList = [nodeIndex]
        self.overlapList = []
        self.listOfPrevNodes = []
        self.listOfNextNodes = []
        self.visited = 0


class seqGraph(object):
    def __init__(self, numberOfNodes):
        self.graphNodesList = [seqGraphNode(i) for i in range(numberOfNodes)]
    
    def loadFromFile(self, folderName, fileName):
        
        f= open(folderName + fileName, 'r')
        numberOfNodes = 0 
        tmp = f.readline().rstrip()
        while (len(tmp) > 0):
            tmp = f.readline().rstrip()
            numberOfNodes += 1
        f.close()
        
        self.graphNodesList = [seqGraphNode(i) for i in range(numberOfNodes)]
        
        f= open(folderName + fileName, 'r')
        tmp = f.readline().rstrip()
        runningIndex = 0
        while (len(tmp) > 0):
            dataList = tmp.split(';')
            

            for i in range(6):
                if len(dataList[i]) >0:
                    myList = dataList[i].split(',')
                else:
                    myList = []

                    
                if i == 0:
                    self.graphNodesList[runningIndex].nodeIndex = int(dataList[i])
                elif i == 1:
                    self.graphNodesList[runningIndex].nodeIndexList = []
                    for eachitem in myList:
                        self.graphNodesList[runningIndex].nodeIndexList.append(int(eachitem))
                        
                elif i == 2:
                    self.graphNodesList[runningIndex].overlapList= []
                    for eachitem in myList:
                        self.graphNodesList[runningIndex].overlapList.append(int(eachitem))
                elif i == 3 or i == 4:
                    for eachitem in myList:
                        mydata = eachitem.split('-')
                        if i ==3 :
                            self.graphNodesList[runningIndex].listOfPrevNodes.append([int(mydata[0]), int(mydata[1])])
                        elif i == 4:
                            self.graphNodesList[runningIndex].listOfNextNodes.append([int(mydata[0]), int(mydata[1])])
                elif i == 5:
                    self.graphNodesList[runningIndex].visited = int(dataList[i])
            
                
            tmp = f.readline().rstrip()
            runningIndex  = runningIndex +1 
        f.close()
    
    def insertEdge(self, i, j, wt):
        if j != -1 and i!= -1:
            haveInserted = nameInEdgeList(j, self.graphNodesList[i].listOfNextNodes)

            if not haveInserted:      
                self.graphNodesList[i].listOfNextNodes.append([j,wt])
                
                
            haveInserted = nameInEdgeList(i, self.graphNodesList[j].listOfPrevNodes) 
            
            if not haveInserted:      
                self.graphNodesList[j].listOfPrevNodes.append([i,wt])
                


    def findStartList(self):
        self.myStartList = []
        for eachnode in self.graphNodesList:
            if len(eachnode.listOfPrevNodes) == 0:
                self.myStartList.append(eachnode)
    def reportEdge(self):
        for eachnode in self.graphNodesList:
            myName = eachnode.nodeIndex
            for eachnextnode in eachnode.listOfNextNodes:
                nextName = eachnextnode[0]
                wt = eachnextnode[1]
                if wt > 1000:
                    print str(myName)+"->" + str(nextName) +"{weight:" + str(int(wt/1000))+"}"
                else:
                    print str(myName)+"->" + str(nextName) +"{weight:" + str(0.5)+"}"

    def condense(self):
        for eachnode in self.graphNodesList:
            myname = eachnode.nodeIndex
            
            for eachnextnode in eachnode.listOfNextNodes:
                nextname = eachnextnode[0]
                wt = eachnextnode[1]
                
                if len(eachnode.listOfNextNodes) == 1 and len(self.graphNodesList[nextname].listOfPrevNodes) == 1:

                    originalPrev = self.graphNodesList[myname].listOfPrevNodes
                    
                    self.graphNodesList[nextname].nodeIndexList  = self.graphNodesList[myname].nodeIndexList +self.graphNodesList[nextname].nodeIndexList
                    self.graphNodesList[nextname].overlapList =  self.graphNodesList[myname].overlapList +[wt] +self.graphNodesList[nextname].overlapList
                    
                    haveInserted = nameInEdgeList(myname, self.graphNodesList[nextname].listOfPrevNodes)
                    if haveInserted :
                        self.graphNodesList[nextname].listOfPrevNodes = removeItem(self.graphNodesList[nextname].listOfPrevNodes, myname)
                    
                    for eachp in originalPrev:  
                        if not eachp in    self.graphNodesList[nextname].listOfPrevNodes:
                            self.graphNodesList[nextname].listOfPrevNodes.append(eachp)
                         
                    self.graphNodesList[myname].nodeIndexList  = [] 
                    self.graphNodesList[myname].overlapList = []
                    self.graphNodesList[myname].listOfNextNodes= []
                    self.graphNodesList[myname].listOfPrevNodes= []
                
                    for eachoriginalprev in originalPrev:
                        prevname = eachoriginalprev[0]
                        #if [myname,eachoriginalprev[1]] in self.graphNodesList[prevname].listOfNextNodes:
                        haveInserted = nameInEdgeList(myname, self.graphNodesList[prevname].listOfNextNodes)
                        if haveInserted:
                            self.graphNodesList[prevname].listOfNextNodes = removeItem(self.graphNodesList[prevname].listOfNextNodes,myname )
                            #self.graphNodesList[prevname].listOfNextNodes.remove([myname,eachoriginalprev[1]])
                        
                        self.graphNodesList[prevname].listOfNextNodes.append([nextname, eachoriginalprev[1]])
                    
                        
    def reportDummyUsefulNode(self):
        countUseful = 0
        countUseless = 0 
        for eachnode in self.graphNodesList:
            #print eachnode.listOfNextNodes
            if len(eachnode.nodeIndexList) > 0:
                countUseful += 1
            else:
                countUseless += 1
                
        print "useful, useless" ,countUseful, countUseless
      
    
    def checkSelfLoops(self):
        for eachnode in self.graphNodesList:
            myIndexList = eachnode.nodeIndexList
            for eachnext in eachnode.listOfNextNodes:
                if eachnext[0] in myIndexList:
                    print  "selfLoop"
        
    
    
    
    def saveToFile(self, folderName, fileName):
        ### Format : nodeIndex,nodeIndexList, overlapList, listOfPrevNodes, listOfNextNodes, visited
        # 3 ; [5, 3] ; [4757]; [[450, 3887], [123,5678]]  ; []; 0
        # 3 ; 5, 3 ; 4757; 450-3887,123-5678  ; ; 0
        
        f = open(folderName + fileName , 'w')
        for eachnode in self.graphNodesList:
            mystr = ""
            mystr =  mystr + str(eachnode.nodeIndex) + ';'
            
            for eachitem, index  in  zip(eachnode.nodeIndexList, range(len(eachnode.nodeIndexList))):
                if index != len(eachnode.nodeIndexList) - 1:
                    mystr= mystr + str(eachitem) + ','
                else: 
                    mystr = mystr + str(eachitem)
                    
            mystr = mystr + ';'
            
            for eachitem, index in  zip(eachnode.overlapList, range((len(eachnode.overlapList)))):
                if index != len(eachnode.overlapList)-1:
                    mystr = mystr + str(eachitem)+ ','
                else:
                    mystr = mystr + str(eachitem)
                
            mystr = mystr + ';'
            
            for eachitem, index in  zip(eachnode.listOfPrevNodes, range(len(eachnode.listOfPrevNodes))):
                if index != len(eachnode.listOfPrevNodes) - 1:
                    mystr = mystr  +   str(eachitem[0]) + '-'  + str(eachitem[1]) + ','
                else:
                    
                    mystr = mystr  +   str(eachitem[0]) + '-'  + str(eachitem[1]) 
            mystr = mystr + ';'
            
            for eachitem, index in  zip(eachnode.listOfNextNodes, range(len(eachnode.listOfNextNodes))):
                if index != len(eachnode.listOfNextNodes) - 1:
                    mystr = mystr  +   str(eachitem[0]) + '-'  + str(eachitem[1]) + ','
                else:
                    
                    mystr = mystr  +   str(eachitem[0]) + '-'  + str(eachitem[1]) 
            mystr = mystr + ';'
                

            mystr = mystr + str(eachnode.visited) + '\n'
            f.write(mystr)
            
        f.close()
          
    '''
    currentNode.nodeIndexList.extend(targetnextnode.nodeIndexList)
    targetnextnode.nodeIndexList = currentNode.nodeIndexList
    currentNode.nodeIndexList = deque()

    for eachprevnode in currentNode.listOfPrevNodes:
        eachprevnode.listOfNextNodes.append(targetnextnode)
        targetnextnode.listOfPrevNodes.append(eachprevnode)
        eachprevnode.listOfNextNodes.remove(currentNode)

    targetnextnode.listOfPrevNodes.remove(currentNode)                    
    currentNode.listOfPrevNodes = []
    currentNode.listOfNextNodes = []
    
    
    '''
       
###################################################### Key functions

### 0) Preprocess by removing embedded contigs (I: contigs.fasta ; O : noEmbed.fasta)


def removeEmbedded(folderName , mummerLink):
    print "removeEmbedded"
    thres = 3
    
    if False:
        useMummerAlign(mummerLink, folderName, "self", "contigs.fasta", "contigs.fasta")
    
    dataList = extractMumData(folderName, "selfOut")
    
    dataList= transformCoor(dataList)
    
    lenDic = obtainLength(folderName, 'contigs.fasta')
    
    removeList = []
    for eachitem in dataList:
        match1, match2, name1, name2 = eachitem[4],eachitem[5], eachitem[7], eachitem[8]
        
        if name1 != name2:
            l1, l2 = lenDic[name1], lenDic[name2]
            
            if abs(l1- match1) < thres and abs(l2-match2) > thres:
                removeList.append(name1)
            elif abs(l1- match1) > thres and abs(l2-match2) < thres:
                removeList.append(name2)
            elif abs(l1- match1) < thres and abs(l2- match2) < thres:
                print "Both shortembedd", eachitem
                
    
    
    nameList = []
    for eachitem in lenDic:
        nameList.append(eachitem)

    print len(nameList)
    
    for eachitem in removeList:
        if eachitem in nameList:
            nameList.remove(eachitem)
    print len(nameList)
    
    putListToFileO(folderName, "contigs.fasta", "noEmbed", nameList)


### 1) Fetch best successor and best predecessor for each contigs (I: noEmbed.fasta ;  O:  left_connect , right_connect )
def fetchSuccessor(folderName , mummerLink ): 
    
    print "fetchSuccessor"
    left_connect, right_connect = [],[] 
        
    print "Direct greedy"
    
    thres = 7
    minLen = 400
    #thres = 10
    #minLen = 200
    
    
    
    writeToFile_Double1(folderName, "contigs.fasta", "contigs_Double.fasta", "contig")
    
    fmyFile = open(folderName+ "contigs_Double.fasta", 'r')
    fSmaller = open(folderName+ "smaller_contigs_Double.fasta", 'w')

    tmp = fmyFile.readline().rstrip()
    maxSize = 50000

    myName = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            fSmaller.write(tmp+'\n')
            myName = tmp[1:]
        else:
            component = tmp[0:min(len(tmp), maxSize )] 
            countComp = len(component)
            fSmaller.write(component)
            
            component =tmp[max(0, len(tmp)-maxSize):len(tmp)]
            fSmaller.write(component)
            countComp = countComp + len(component)
            

            print "DebugName", myName, countComp
            fSmaller.write('\n')

        tmp = fmyFile.readline().rstrip()

    fSmaller.close()
    fmyFile.close()
    
    if False:
        useMummerAlign(mummerLink, folderName, "greedy", "smaller_contigs_Double.fasta", "smaller_contigs_Double.fasta")
        
        
    lengthDic = obtainLength(folderName, "smaller_contigs_Double.fasta") 
    
    dataSetRaw = extractMumData(folderName, "greedyOut")
    
    ### Format [ helperStart, helperEnd , readStart, readEnd,matchLen1,matchLen2,percentMatch,helperName,readName]
        
    dataSet = []
    
    for eachitem in dataSetRaw: 
        helperStart, helperEnd , readStart, readEnd, matchLen1, matchLen2, percentMatch, helperName, readName = eachitem 
        
        if helperName != readName and max(matchLen1, matchLen2) > minLen and readStart < readEnd  and min(helperStart,readStart) < thres and min(lengthDic[helperName]- helperEnd,  lengthDic[readName] - readEnd) < thres:
            conditionForMatch = True
        else:
            conditionForMatch = False

        if conditionForMatch :
            if helperStart < thres:
                
                dataSet.append((max(matchLen1, matchLen2), readName, helperName))
    
    dataSet.sort(reverse=True)
    
    numberOfContig = len(lengthDic)
    
    # [next_item, overlap_length]
    
    leftConnect = [[-1,-1] for i in range(numberOfContig)]
    rightConnect = [[-1,-1] for i in range(numberOfContig)]
    
    dataSet.sort(reverse=True, key = itemgetter(1))
    
    for key, items in groupby(dataSet, itemgetter(1)):
        #if key == "Contig217_d":
        #    print "dddd"
        maxVal = -1
        myName = key
        connectorName = "" 
        for eachsubitem in items:
            if eachsubitem[0] > maxVal:
                maxVal = eachsubitem[0]
                connectorName = eachsubitem[2]
        

        prefix = myName.split('_')
        suffix = connectorName.split('_')
        lengthOfOverlap = maxVal
        
        if prefix[1] == 'p':
            prefixContig = int(prefix[0][6:])*2 
        else:
            prefixContig = int(prefix[0][6:])*2 +1
        
        if suffix[1] == 'p':
            suffixContig = int(suffix[0][6:])*2 
        else:
            suffixContig = int(suffix[0][6:])*2 +1
            
        assert(rightConnect[prefixContig][0] == -1)
        rightConnect[prefixContig][0] = suffixContig
        rightConnect[prefixContig][1] = lengthOfOverlap
        

    dataSet.sort(reverse=True, key = itemgetter(2))
    
    for key, items in groupby(dataSet, itemgetter(2)):

        maxVal = -1
        myName = key
        connectorName = "" 
        for eachsubitem in items:
            if eachsubitem[0] > maxVal:
                maxVal = eachsubitem[0]
                connectorName = eachsubitem[1]
        

        prefix = connectorName.split('_')
        suffix = myName.split('_')
        lengthOfOverlap = maxVal
        
        if prefix[1] == 'p':
            prefixContig = int(prefix[0][6:])*2 
        else:
            prefixContig = int(prefix[0][6:])*2 +1
        
        if suffix[1] == 'p':
            suffixContig = int(suffix[0][6:])*2 
        else:
            suffixContig = int(suffix[0][6:])*2 +1
            
        assert(leftConnect[suffixContig][0] == -1)
        leftConnect[suffixContig][0] = prefixContig 
        leftConnect[suffixContig][1] = lengthOfOverlap
    
    
    ### Write to file: 
    f = open(folderName + 'rightConnect.txt', 'w')
    for eachitem,dummyIndex in zip(rightConnect, range(len(rightConnect))):
        f.write(str(dummyIndex) + ',' +str(eachitem[0])+','+str(eachitem[1])+'\n')
        
    f.close()
    
    f = open(folderName + 'leftConnect.txt', 'w')
    for eachitem,dummyIndex in zip(leftConnect, range(len(leftConnect))):
        f.write(str(dummyIndex) + ',' +str(eachitem[0])+','+str(eachitem[1])+'\n')
        
    f.close()

### 2) Form seqGraph (I: left_connect, right_connect ; O: startList, graphNodes )
def formSeqGraph(folderName , mummerLink ):
    print "formSeqGraph" 
    startList, graphNodes  = [], []
    
    rightConnect = readConnectList(folderName, "rightConnect.txt")
    leftConnect = readConnectList(folderName, "leftConnect.txt")
    
    numberOfNodes = len(rightConnect)
    
    G = seqGraph(numberOfNodes)
        
    for eachitem,i  in zip(rightConnect, range(len(rightConnect))):
        index = i
        connector, weight = eachitem
        G.insertEdge(index, connector,weight)
    
    for eachitem, i  in zip(leftConnect, range(len(leftConnect))):
        index = i
        connector, weight = eachitem
        G.insertEdge(connector, index, weight)
    

    
    G.condense()
    G.saveToFile(folderName, "condensedGraph.txt")
    G.checkSelfLoops()
    
    G2 = seqGraph(0)
    G2.loadFromFile(folderName, "condensedGraph.txt")
    
    compareGraphUnitTest(G, G2)
    #G.reportDummyUsefulNode()
    #G.reportEdge()
        

### 3) X-phased seqGraph (I: startList, graphNodes; O: startList, graphNodes )
def xPhased(folderName , mummerLink ):
    print "xPhased"



### 4) EC reduction (I: startList, graphNodes ; O: startList, graphNodes )
def ECReduction(folderName , mummerLink ):
    print "ECReduction" 
    

### 5) Read the contigs out (I: startList, graphNodes, ; O:improved.fasta, openZone.txt)
def readContigOut(folderName, mummerLink):
    
    print "readContigOut"
    
    G = seqGraph(0)
    G.loadFromFile(folderName, "condensedGraph.txt")
    
    myContigsDic = loadContigsFromFile(folderName, "contigs_Double.fasta")
    
    contigUsed = [False for i in range(len(G.graphNodesList)/2)]
     
    seqToPrint = []
    
    for eachnode in G.graphNodesList:
        if len(eachnode.nodeIndexList) > 0:
            tmpSeq = ""
            ### debug consistency of t/f
            ckList = []
            for dummy in eachnode.nodeIndexList:
                indexToAdd = dummy
                readNum = indexToAdd /2
                ckList.append(contigUsed[readNum])
                
            if len(ckList) >0 and not all(ckList) and any(ckList):
                print eachnode.nodeIndex, ckList
            
            ### end debug 
            
            for i in range(len(eachnode.nodeIndexList)):
                
                indexToAdd = eachnode.nodeIndexList[i]
                readNum = indexToAdd /2
                orientation = indexToAdd %2 

                if contigUsed[readNum] == False:
                    if i !=len(eachnode.nodeIndexList) -1:

                        overlapLen = eachnode.overlapList[i]
                        if orientation == 0:
                            tmpSeq = tmpSeq + myContigsDic['Contig'+str(readNum)+'_'+'p'][0:-overlapLen]
                        else:
                            tmpSeq = tmpSeq + myContigsDic['Contig'+str(readNum)+'_'+'d'][0:-overlapLen]
                    else:
                        if orientation == 0:
                            tmpSeq = tmpSeq + myContigsDic['Contig'+str(readNum)+'_'+'p']
                        else:
                            tmpSeq = tmpSeq + myContigsDic['Contig'+str(readNum)+'_'+'d']
                            
                    
                    contigUsed[readNum] = True
            if len(tmpSeq) > 0:
                seqToPrint.append(tmpSeq)
    
    
    fImproved = open(folderName + 'improved.fasta', 'w')
    for eachcontig, dummyIndex in zip(seqToPrint, range(len(seqToPrint))):
        fImproved.write(">Segkk"+str(dummyIndex)+'\n')
        fImproved.write(eachcontig+'\n')
         
    fImproved.close()

    


### 6) Fill gap(I: improved.fasta, openZone,txt ; O: improved2.fasta )
def fillGap(folderName , mummerLink):
    print "fillGap"


### 7) Compare with reference (I: improved.fasta, improved2.fasta, reference.fasta ; O : assembly assessment report )
def compareWithReference(folderName , mummerLink):
    print "compareWithReference"
    quastEvaluate(folderName, "quast-2.3/", originalName = "contigs.fasta", improvedNameList= ["noEmbed.fasta", "improved.fasta"] , referenceName= "reference.fasta" )
    
    
     

###################################################### Starting point
def mainFlow(folderName , mummerLink ):
    print "Go Bears! ! !" 
    #removeEmbedded(folderName , mummerLink)
    #fetchSuccessor(folderName , mummerLink )
    #formSeqGraph(folderName , mummerLink )
    
    xPhased(folderName , mummerLink )
    ECReduction(folderName , mummerLink )
    
    readContigOut(folderName, mummerLink)
    
    fillGap(folderName , mummerLink)
    
    compareWithReference(folderName , mummerLink)
    print "<3 Do cool things that matter <3"
    
    
    
folderName = "S_cerivisea/"
mummerLink = "MUMmer3.23/"

mainFlow(folderName,mummerLink)