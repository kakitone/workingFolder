import numpy as np
from scipy import weave
from itertools import groupby
from operator import itemgetter

def updateVote(newAlignedSeq1, newAlignedSeq2,tempVoteTable, istart, jstart, iend, jend ):

    temp = 1
        
    Lseq = min(len(newAlignedSeq1), len(newAlignedSeq2))
    cumCount = 0 
    insCumCount = 0
    
    for j in range(Lseq): 
        if newAlignedSeq2[j] == 0 :
            # format : [number of insert, base , count ] 
            subbase = newAlignedSeq1[j] -1 
            currentIns = [insCumCount, subbase, 1 ]
            existInIns = False
            for eachins, dummyindex in zip(tempVoteTable.insCount[cumCount + jstart-1], range(len(tempVoteTable.insCount[cumCount + jstart-1]))):
                #print "eachins" ,eachins
                if eachins[0:2] == currentIns[0:2]:
                    tempVoteTable.insCount[cumCount + jstart-1][dummyindex][2] = eachins[2]+1 
                    existInIns = True
                    #print "debug",tempVoteTable.insCount[cumCount + jstart-1][dummyindex][2] 
                    
            if not existInIns :
                tempVoteTable.insCount[cumCount -1  +jstart].append(currentIns)
            
            insCumCount = insCumCount + 1 
        elif newAlignedSeq1[j] == 0:
            tempVoteTable.delCount[cumCount+ jstart] += 1 
            
        elif newAlignedSeq2[j] == newAlignedSeq1[j]:
            tempVoteTable.confirmCount[cumCount + jstart] += 1 
            
        elif newAlignedSeq2[j]  != newAlignedSeq1[j]  :
            subbase = newAlignedSeq1[j] -1 
            tempVoteTable.subCount[cumCount+jstart][subbase] += 1 
            
            
        if newAlignedSeq2[j] != 0:
            if (j< len(newAlignedSeq2) -1 and  newAlignedSeq2[j+1] != 0) or j==len(newAlignedSeq2) -1:
                tempVoteTable.confirmNoIns[cumCount + jstart] += 1 
                
            cumCount = cumCount + 1 
            insCumCount = 0 
                
class voteTable(object):
    def __init__(self, index, eachlongread):
        self.index = index    
        self.longread = eachlongread

        
        
    def initVoteTable(self):
        L = len(self.longread)
        
        self.delCount = [0 for i in range(L)]
        self.confirmCount =[0 for i in range(L)]
        # 1, 2,3 ,4 - > 0, 1, 2, 3 
        self.subCount = [[0 , 0 ,0 ,0 ] for i in range(L)]
        
        self.insCount = [ [] for i in range(L) ]
        self.confirmNoIns= [0 for i in range(L)]
        
    def polished(self):
        #if len(self.longread) == 714:
        #    for i in range(len(self.longread)):
        #        print i, self.longread[i], self.delCount[i], self.subCount[i], self.confirmCount[i], self.insCount[i], self.confirmNoIns[i]
        
        print "Before", len(self.longread)
        returnString = ""
        insTotal, delTotal = 0,0
        for i in range(len(self.longread)):
            if (self.confirmCount[i]) >= self.delCount[i]:  
                if self.longread[i] == 1:
                    returnString += 'A'
                elif self.longread[i] == 2:
                    returnString += 'C'
                elif self.longread[i] == 3:
                    returnString += 'G'
                elif self.longread[i] == 4:
                    returnString += 'T'    
            else:
                delTotal += 1
                    
            maxScore = -1
            base = -1
            self.insCount[i].sort()
            
            for key, items in groupby(self.insCount[i], itemgetter(0)):
                maxScore = -1
                base = -1
                if key == 0:
                    for eachitem in items:
                        if eachitem[2] > maxScore:
                            maxScore = eachitem[2]
                            base = eachitem[1]
            
            if len(self.insCount[i]) > 0 and (self.confirmNoIns[i]) < maxScore:
                print "self.confirmNoIns[i], maxScore, self.insCount[i]", self.confirmNoIns[i], maxScore, self.insCount[i]
                if base == 0:
                    returnString += 'A'
                elif base == 1:
                    returnString += 'C'
                elif base == 2:
                    returnString += 'G'
                elif base == 3:
                    returnString += 'T'

                insTotal += 1
        print "After: length, insTotal, delTotal",  len(returnString), insTotal, delTotal
        return returnString
        
 


def SWAlignment(seq1 , seq2):
    score  = 0 

    wts = -10
    wti = -1 
    wtd = -1
    wtm=  1
    
    #print "wts, wti, wtd, wtm",wts, wti, wtd, wtm
    
    m = len(seq1) + 1
    n = len(seq2) +1
    
    #H = np.zeros(m*n, dtype = np.float64).reshape(m,n)
    #B = np.zeros(m*n, dtype = np.float64).reshape(m,n)
    H = np.zeros([m,n], dtype = np.float64)
    B = np.zeros([m,n], dtype = np.float64)

    # Assign weights 
    for i in range(m):
        H[i][0] = 0
        B[i][0] = 4
    for j in range(n):
        H[0][j]  =0
        B[0][j] = 4 
        
        
        
    seq1NP = np.zeros(m-1, dtype = np.float64)
    seq2NP = np.zeros(n-1, dtype = np.float64)
    
    for i in range(m-1):
        seq1NP[i] = seq1[i]
    for j in range(n-1) :
        seq2NP[j] = seq2[j]
        
 #   print "seq1NP, seq2NP"
 #   print seq1NP, seq2NP
 #   print "seq1, seq2"
 #   print seq1, seq2
    
    code =\
        """
        int i; 
        int j ;
        double w;

        for (i =1 ;i <m ; i++){
            for (j=1; j<n; j++){
                if (seq1NP[i-1] == seq2NP[j-1]){
                    w = wtm ;
                }
                else{
                    w= wts ;
                }
                
                    H2(i,j) = 0 ;
                    
                    if (H2(i,j) < H2(i-1,j-1) + w) {
                        H2(i,j) = H2(i-1,j-1) + w ; 
                    }

                    if (H2(i,j) < H2(i-1,j)+wtd ) {
                        H2(i,j) = H2(i-1,j)+wtd ;
                    }

                    if  (H2(i,j) < H2(i,j-1) + wti){
                        H2(i,j) = H2(i,j-1) + wti;
                    }

                    if (H2(i-1,j-1) + w == H2(i,j)){
                        B2(i,j) = 1 ; 
                    }
                    else if (H2(i-1,j)+wtd == H2(i,j)){
                        B2(i,j) = 2;
                    }
                    
                    else if  (H2(i,j-1) + wti == H2(i,j)){
                        B2(i,j) = 3;
                    }
                    else if (0 == H2(i,j)) {
                        B2(i,j) = 4 ;
                    }
            }
        }

        """
    
    weave.inline(code, ['H','B','m', 'n', 'wtm','wts','wti','wtd', 'seq1NP', 'seq2NP'])
    
    
#    for i in range(1, m):
#        for j in range(1, n):
#            if seq1[i-1] == seq2[j-1]:
#                w = wtm 
#            else:
#                w= wts
            
#            H[i][j] = max(H[i-1][j-1] + w, H[i-1][j]+wtd, H[i][j-1] + wti, 0)
#            if H[i-1][j-1] + w == H[i][j] :
#                B[i][j] = 1 
#            elif H[i-1][j]+wtd == H[i][j] :
#                B[i][j] = 2
#            elif  H[i][j-1] + wti == H[i][j] :
#                B[i][j] = 3
#            elif 0 == H[i][j]:
#                B[i][j] = 4 

    # Backtrack 
    alignedSeq1 = []
    alignedSeq2 = []
    
    bestindex = np.argmax(H)

 
    
    endi = bestindex / n
    endj = bestindex%n
    
    #print endi , endj
    score = H[endi][endj]
    
   # print "score", score 
   # print "B", B[1]
   # print "H", H
    
    tempi, tempj = endi , endj
    while (B[tempi][tempj] != 4):
        if B[tempi][tempj] == 1:
            alignedSeq1.append(seq1[tempi -1 ])
            alignedSeq2.append(seq2[tempj-1 ]) 
            tempi = tempi -1 
            tempj = tempj -1           

        elif B[tempi][tempj] == 2 :
            alignedSeq1.append(seq1[tempi-1 ])
            alignedSeq2.append(0) 
            tempi = tempi - 1

        elif B[tempi][tempj] == 3 :
            alignedSeq1.append(0)
            alignedSeq2.append(seq2[tempj-1 ])
            tempj = tempj -1

    
    starti, startj = tempi , tempj 
    
    returnalignedSeq1 = reverseString(alignedSeq1)
    returnalignedSeq2 = reverseString(alignedSeq2)
    
    #printSeq(returnalignedSeq1)
    #printSeq(returnalignedSeq2)
    return score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj


def reverseString(str):
    newstr = []
    for index in range(len(str)):
        newstr.append(str[len(str)-1-index])
    
    return newstr

def transformToDesiredForm(alignedSeq1, alignedSeq2):
    newAlignedSeq1, newAlignedSeq2 = [] ,[]

    Lseq = len(alignedSeq1)
    
    for i in range(Lseq):
        newAlignedSeq1.append(alignedSeq1[i])
        newAlignedSeq2.append(alignedSeq2[i])
    
    i =0 
    while (i <= len(newAlignedSeq1) -2 ):
        if newAlignedSeq1[i] == 0 and newAlignedSeq2[i] == 0 :
            newAlignedSeq1.pop(i)
            newAlignedSeq2.pop(i)
                    
        elif newAlignedSeq1[i] == newAlignedSeq2[i+1] and  newAlignedSeq1[i] == newAlignedSeq2[i+1] and  newAlignedSeq2[i] == 0:
            newAlignedSeq2[i] = newAlignedSeq1[i] 
            newAlignedSeq2[i+1] = 0
            i = i +1 
        elif newAlignedSeq2[i] == newAlignedSeq1[i+1] and  newAlignedSeq2[i] == newAlignedSeq1[i+1] and  newAlignedSeq1[i] == 0:
            newAlignedSeq1[i] = newAlignedSeq2[i] 
            newAlignedSeq1[i+1] = 0
            i = i +1 

        elif  newAlignedSeq1[i] == 0 and newAlignedSeq2[i+1] == 0:
            newAlignedSeq1[i] = newAlignedSeq1[i+1]
            newAlignedSeq1[i+1] =0 
            i = i +1 
        elif  newAlignedSeq2[i] == 0 and newAlignedSeq1[i+1] == 0:
            newAlignedSeq2[i] = newAlignedSeq2[i+1] 
            newAlignedSeq2[i+1] = 0   
            i = i +1          
        else: 
            i = i +1 

    return newAlignedSeq1, newAlignedSeq2



def polishing(eachlongread, helperList):
    
    tempVoteTable = voteTable(0, eachlongread) 
    tempVoteTable.initVoteTable()            
    
    for eachshortreaad in helperList:
        score, alignedSeq1, alignedSeq2, istart, jstart, iend, jend = SWAlignment(eachshortreaad, eachlongread )
        newAlignedSeq1, newAlignedSeq2 = transformToDesiredForm(alignedSeq1, alignedSeq2)
        updateVote(newAlignedSeq1, newAlignedSeq2,tempVoteTable,istart, jstart, iend, jend )
            
    cleanedLongRead = tempVoteTable.polished()
    
    return cleanedLongRead




