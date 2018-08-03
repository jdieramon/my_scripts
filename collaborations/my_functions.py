import collections
from random import choice


def getTwoSubst(seq,X,a, b, Y, c, d) : 
    '''It gets the set of sequences from a seq with 2 non-standard bases. 
    a/b  are the mixed bases for the 1st symbol(X)
    c/d  are the mixed bases for the 2nd symbol(Y)'''
    
    def setOneSubst(seq) : 
        motifs=[]
        optionX=[a,b]
        for itemb in optionX : 
            motifs.append(seq.replace(X,itemb))
        return motifs
        
    def setSecondSubst(seq) : 
        motifs=[]
        optionY=[c,d]
        for element in setOneSubst(seq) : 
            for itemy in optionY : 
                motifs.append(element.replace(Y, itemy))
        return motifs    
    
    setOneSubst(seq)
    return setSecondSubst(seq)

smre = 'ACCWAM'
getTwoSubst(smre, 'W', 'A', 'T', 'M', 'A', 'C')


def readGenome(filename):
    '''read and store the nucleotide seqs from a file'''
    genome = []
    seqs=''
    with open(filename, 'r') as f:
        # ignore header line with genome information
        for line in f:
            if line[0] == '>':
                genome.append(seqs) 
                seqs = ''
            if not line[0] == '>' : 
                    seqs += line.rstrip()
        genome.append(seqs) # add last sequence
        del genome[0]       # del first empty sequence
            
    return genome
    
def nameGenome(fileName) : 
    '''read and store the sequence name from a file'''
    seqName = []
    with open(fileName, 'r') as f:
        for line in f:
                if line[0] == '>':
                    seqName.append(line[1::].rstrip()) 
    return seqName
    
def naive(p,t): #pattern, text
    occurrences = []
    for i in range(len(t) - len(p) + 1): #every position in t where t can start w/o pass the end of t
        match = True
        for j in range(len(p)):        # loop over characters
            if not t[i+j] == p[j]:     # compare characters 
                match = False
                break #no sense comparing the rest of p because we find a mistmach
        if match:
            occurrences.append(i)
    return occurrences
    
    
    
    
def search(CRE, dataset) :  
    """Return list containing sequence and ocurrence of a CRE within a dataset"""
    's, sequence names containing CRE'
    'occurrences, position in the sequence'
    'genome, dataset'
    
    s = []
    occurrences = []
    for seqs in range(len(dataset)) :
        for motif in CRE : 
            if len(naive(motif, dataset[seqs])) >0 : # si encuentra el motivo :
                s.append(seqs)
                occurrences.append(naive(motif, dataset[seqs]))
                
    return s,occurrences
    
    

def countBases(sequence): 
    '''count total number of each base in a given sequence
    sequence must be a list of sequences'''
    count = collections.Counter()
    for seq in sequence:
        count.update(seq)
    return(count)
    
    
    

def weightedchoice(items): 
    return choice("".join(x * y for x, y in items))
    	

def getPromoter(length):
       '''Genera promotor con bases ponderadas'''
       DNA = ""
       for count in range(length):
            DNA+=weightedchoice([('C',15),('G',13),('A',35),('T',37)])  #AQUI tengo que introducir la freq. de cada base
       return DNA
    

        
def atLeastOneMotif(CRE, bp, nset):
    '''
    Return number of promoters in a given data set containing at least 1 CRE
    
    CRE, cis-regulatory element
    bp, promoter length (1500)
    nset, number of promoters 
    '''
    numSuccess = 0 #number of seqs with at least 1 motive
    
    for p in range(nset) : 
        promotor = getPromoter(bp)
        numMotives = 0
        
        for m in CRE : 
            if m in promotor : 
                numMotives +=1
                
                if numMotives >0:
                    numSuccess +=1
                    break
        
                        
    fracSuccess = numSuccess/float(nset)
    
    return fracSuccess

def runSimforOneMotif(CRE,bp,nset,numTrials) : 
    '''
    Run a MC simulation (numTrials). 
    Return a vector list with the number of expected promotores 
    in a given data set (nset) containing at least 1 CRE in each
    trial (numTrials)
    
    CRE, cis-regulatory element
    bp, promoter length (1500)
    nset, number promoters in the data set
    numTrials, number of simulations
    '''
    
    vals= []
    total = 0
    for i in range(1,numTrials+1) :
        ## Si quiero ir siguiendo las simulaciones
        #if i%200==0 : 
          #  print 'Simulacion',i+1
        total += atLeastOneMotif(CRE, bp, nset)*nset
        vals.append(round(total/i,3))
        
    return vals
    
    
    
    
    
def estimatedMotif(CRE, bp, numTrials,nset):
    ''' numTrials, numero promotores que se van a generar'''

    numSuccess = 0
    
    for i in range(numTrials) : 
        promotor = getPromoter(bp)

        for m in CRE : 
             numSuccess += len(naive(m,promotor)) #numero de motivos
               
    fracSuccess = numSuccess/float(numTrials)
        
    #print numSuccess, 'motivos'
    #print fracSuccess, 'motivos/promotor'
    #print 'At least once in',1/fracSuccess,'promoters'
    #print 'sobre un set de ',nset,'seqs, deberia haber encontrado',fracSuccess*nset,'motivos' 
    return fracSuccess*nset
    
    
    
def getTarget(CRE, bp, goal):
    
    ''' Get the number of promoter data set needed (goal) to find at least 1 CRE
    per promoter.
    
    bp = promoter length
    goal = number of promoters with at least 1 motif ''' 
    
    numSuccess = 0
    numTries = 0
    count = 0
    while numSuccess < goal : 
        promotor = getPromoter(bp)
        numTries +=1
        
        for m in CRE :
            if m in promotor : 
                count += 1
                if count >= goal : 
                    numSuccess= numTries
    
    return numSuccess
    
    
    
    
def runSim(CRE, bp, goal, numTrials):
    '''bp: bases del promotor (1500)
    goal (float) :numero de veces que se busca un motivo concreto
    numTrial : numero de simulaciones (2000)'''
    #val = [] #por si quiero llevar la cuenta del valor acumulado
    total = 0
    
    for i in range(1,numTrials+1):
        if i%200==0 : 
            print 'Sim.',i,'...',round(total/float(i),2) #veo parcial so far
        total += getTarget(CRE, bp, goal)
    #    val.append(total/float(i))
        
    
    #print val
    print 'Average number of tries =', total/float(numTrials)
    print 'P-value=',1/(total/float(numTrials))
    print''