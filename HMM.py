from math import *
import numpy as np
INF = float("inf")

def log_sum_exp(numlist):
    # using log-trick to compute log(sum(exp(x) for x in numlist))
    # mitigate the problem of underflow
    try:
        minx = min([x for x in numlist if x != -INF])
    except:
        return -INF    
    s = sum(exp(x-minx) for x in numlist)
    return minx + log(s) if s > 0 else -INF

class HMM:
    def __init__(self):
        self.alphabet = []
        # emission probability of the I states, one dictionary for each I in the model
        self.eI = [] 
        # emission probability of the M states, one dictionary for each M in the model
        # the first M state is called 'B'; it never emits anything so the associated dictionary is always empty
        self.eM = [{}] 
        # transition probability, one dictionary for each set of states (D,M,I) in the model
        self.t = [] 
        # number of matching states in the model, excluding B and E
        self.nstate = 0
    
    def load(self,hmmfile):
    # only load the first model in the given hmmfile if there are two or more models
        with open(hmmfile,'r') as fin:
            for line in fin:
                stream = line.strip().split()
                if stream[0] == "LENG":
                    self.nstate = int(stream[1])
                if stream[0] == "HMM": 
                    # read alphabet
                    self.alphabet = stream[1:]
                    # read transition order
                    stream = fin.readline().strip().split()
                    trans_order = [(y[0]+y[3]).upper() for y in stream]
                    # read the next line, if it is the COMPO line then ignore it and read one more line
                    stream = fin.readline().strip().split()
                    if stream[0] == "COMPO":
                        stream = fin.readline().strip().split()
                    # now the stream should be at the I0 state; read the emission of the I0 state 
                    e = {}
                    for (x,y) in zip(self.alphabet,stream):
                        e[x] = -float(y)
                    self.eI.append(e)    
                    # now the stream should be at the B state; read in the transition probability
                    stream = fin.readline().strip().split()
                    tB = {'MM':-INF,'MD':-INF,'MI':-INF,'IM':-INF,'II':-INF,'ID':-INF,'DM':-INF,'DI':-INF,'DD':-INF}
                    for x,y in zip(trans_order,stream):
                        tB[x] = -INF if y == '*' else -float(y)
                    self.t.append(tB)
                    break    
            
            for i in range(1,self.nstate+1):
                # read each set of three lines at a time
                stream = fin.readline().strip().split() # this one is the emission of the M state
                if float(stream[0]) != i:
                    print("Warning: incosistent state indexing in hmm file; expecting state "+ str(i) + "; getting state " + stream[0])
                e = {}
                for x,y in zip(self.alphabet,stream[1:]):
                    e[x] = -INF if y == "*" else -float(y) 
                self.eM.append(e)    
                # the next line is the emission of the I state
                stream = fin.readline().strip().split() 
                e = {}
                for x,y in zip(self.alphabet,stream):
                    e[x] = -INF if y == "*" else -float(y) 
                self.eI.append(e)                    
                # the next line contains the transition probs
                stream = fin.readline().strip().split() # this one is the transition prob
                tB = {'MM':-INF,'MD':-INF,'MI':-INF,'IM':-INF,'II':-INF,'ID':-INF,'DM':-INF,'DI':-INF,'DD':-INF}
                for x,y in zip(trans_order,stream):
                    tB[x] = -INF if y == '*' else -float(y)
                self.t.append(tB)
    
    def compute_llh(self,query):
        # compute likelihood of an aligned query to the HMM
        # return -inf if the query is not properly aligned
        j = 0
        prev_state = 'M'
        llh = 0

        for c in query:
            if c == '-': # gap -> deletion
                curr_state = 'D'
            elif c >= 'a' and c <= 'z': # lowercase -> insertion
                curr_state = 'I'
            elif c >= 'A' and c <= 'Z': # uppercase -> match
                curr_state = 'M'
            else: # encounter invalid symbol
                return -INF
            
            trans = prev_state + curr_state

            # penalize the transition
            if trans in self.t[j]:
                llh += self.t[j][trans]
            else: # encounter invalid transition
                return -INF 
            
            # transit: update the state index
            j += (curr_state != 'I') # move to the next state unless this is a move towards an 'I'
            if j > self.nstate: # reach end of the HMM chain but not end of the sequence
                return -INF

            # penalize the emission
            if curr_state == 'M':
                llh += self.eM[j][c]
            elif curr_state == 'I': 
                llh += self.eI[j][c.upper()]
            
            # update state
            prev_state = curr_state

        if j != self.nstate: # does not reach the E state at the end
            return float("-inf")
        
        # the last move: towards the 'E' state
        trans = prev_state + 'M'
        llh += self.t[j][trans]    
        return llh      
    
    def Viterbi(self,query):
        # implement the Viterbi algorithm to find the most probable path of a query sequence and its likelihood
        # return the log-likelihood and the path. If the likelihood is 0, return "$" as the path 
        
        n = len(query)
        m = self.nstate
        


        VM = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        VI = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        VD = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        arrows = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        
        Vscore = -INF
        aln = ""
        
        
        
        #import pdb; pdb.set_trace()
        # YOUR CODE HERE





        #for i in range(m):
        # first state first observations probability of 1
        VM[0][0] = 0


        # if statements
        #n and i have to do with queries
        # j and m have to do with states        
        #second one longer

        # each state I want to know where I came from
        arrowM = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        arrowI = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        arrowD = [[float("-inf") for j in range(m+2)] for i in range(n+1)]

        for i in range(1,n+1):
            # 0 , 1 
            #import pdb; pdb.set_trace()
            for j in range(m+1):
                # 0, 1, 2
                #import pdb; pdb.set_trace()
                #print('beep')
                if i>= 1 and j>=1:
                    option = [VM[i-1][j-1] + self.t[j-1]['MM'], VI[i-1][j-1] + self.t[j-1]['IM'], VD[i-1][j-1]  + self.t[j-1]['DM']]
                    VM[i][j] = self.eM[j][query[i-1]] + max(option)
                    arrowM[i][j] = max(enumerate(option), key=lambda x: x[1])[0]
                
                if i>=1:
                    option = [VM[i-1][j] + self.t[j]['MI'], VI[i-1][j] + self.t[j]['II'], VD[i-1][j] + self.t[j]['DI']]
                    VI[i][j] = self.eI[j][query[i-1]] + max(option)
                    arrowI[i][j] = max(enumerate(option), key=lambda x: x[1])[0]
                if j>=1:
                    option = [VM[i][j-1] + self.t[j-1]['MD'], VI[i][j-1] + self.t[j-1]['ID'], VD[i][j-1] + self.t[j-1]['DD']]
                    VD[i][j] = max( option)
                    arrowD[i][j] = max(enumerate(option), key=lambda x: x[1])[0]
                    #option = [max(VM[i][j]),max(VI[i][j]),max(VD[i][j])]
                # Make some arrow array here

                #import pdb; pdb.set_trace()
                #import pdb; pdb.set_trace()
        # some backtracking
                

            #option = [blosum + S[i-1][j-1][0], S[i-1][j][0] + gap*len(L2)*len(L1), S[i][j-1][0] + gap*len(L2)*len(L1)]  
                #max_index = max(enumerate(option), key=lambda x: x[1])[0]
        # reverse through the states
        #import pdb; pdb.set_trace()``
        
        # at VM
        # where did it come from? 
        #arrowM[i][j] (1,2,3)
        option = [VM[n][m] + self.t[m]['MM'],VI[n][m]+ self.t[m]['IM'],VD[n][m] + self.t[m]['DM']]
        VM[n][m+1] = max(option)
        arrow  = max(enumerate(option), key=lambda x: x[1])[0]
        arrow1 = arrow
        #import pdb; pdb.set_trace()
        aln = ""

        i = n
        j = m
        if arrow == 0:
            aln= aln + query[i-1]
            arrow = arrowM[i][j]
            
        elif arrow == 1:
            aln= aln + query[i-1].lower()
            arrow = arrowI[i][j]
        elif arrow == 2:
            arrow = arrowD[i][j]
            aln= aln + "-"
        

        #import pdb; pdb.set_trace()
        while not arrow== -INF:
            #import pdb; pdb.set_trace()
            if arrow == 0:
                j = j-1
                i = i-1
                if i>=0:
                    try:
                        if not arrow1 == 1: 
                            aln = aln + query[i-1]
                        elif not i==0:
                            aln = aln + query[i-1]
                    except:
                        pass
                #import pdb; pdb.set_trace()
                arrow = arrowM[i][j]
            elif arrow == 1:
                i = i-1
                if i>=0:
                    try:
                        aln= aln + query[i-1].lower()
                    except:
                        pass
                arrow = arrowI[i][j]
                #import pdb; pdb.set_trace()
            elif arrow ==2:
                j = j-1
                aln = aln + "-"
                arrow = arrowD[i][j]
        
        # take off final match witch represents begining 



        #import pdb; pdb.set_trace()
        #max_index = max(enumerate(option), key=lambda x: x[1])[0]

       # i = n-1
        #for letter in reversed(range(1,n+1)):
         #   option = [max(VM[j]),max(VI[j]),max(VD[j])]
          #  max_index = max(enumerate(option), key=lambda x: x[1])[0]
           #
            #if max_index == 0:
             #   aln = aln + query[i]
              #  i = i-1
            #if max_index == 1:
             #   aln = aln + query[i].lower()    
             #   i = i-1
           # if max_index == 2:
            #    aln = aln + "-"
                


        #VM[n][m+1] = max(VD[n][m] + self.t[m]['DM'],VM[n][m] + self.t[m]['MM'],VI[n][m]+ self.t[m]['IM'])
        #import pdb; pdb.set_trace()


        Vscore = VM[n][m +1]
        if Vscore == -INF:
            aln= '$'
        elif arrow1 ==0:
            aln = aln[::-1]
            aln = aln[1:]
        else:
            aln = aln[::-1]
        import pdb; pdb.set_trace()
        return Vscore,aln        


    def Forward(self,query):
        # implement the Forward algorithm to compute the marginal probability of a query sequence
        n = len(query)
        m = self.nstate
        
        FM = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        FI = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        FD = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        
        FM[0][0] = 0
        for i in range(1,n+1):
            # 0 , 1 
            #import pdb; pdb.set_trace()
            for j in range(m+1):
                # 0, 1, 2
                #import pdb; pdb.set_trace()
                #print('beep')
                if i>= 1 and j>=1:
                    #import pdb; pdb.set_trace()
                    try:
                        FM[i][j] = self.eM[j][query[i-1]] + log(exp(FM[i-1][j-1] + self.t[j-1]['MM']) + exp(FI[i-1][j-1] + self.t[j-1]['IM']) +  exp(FD[i-1][j-1] + self.t[j-1]['DM']))
                    except:
                        FM[i][j] = -INF
                if i>=1:
                    try:
                        FI[i][j] = self.eI[j][query[i-1]] + log(exp(FM[i-1][j] +self.t[j]['MI']) + exp(FI[i-1][j] + self.t[j]['II']) + exp(FD[i-1][j] + self.t[j]['DI']))
                    except:
                        FI[i][j] = -INF
                
                if j>=1:
                    try:
                        FD[i][j] = log(exp(FM[i][j-1] + self.t[j-1]['MD']) + exp(FI[i][j-1] + self.t[j-1]['ID']) + exp(FD[i][j-1] + self.t[j-1]['DD']))
                    except:
                        FD[i][j] = -INF
                #import pdb; pdb.set_trace()
                #import pdb; pdb.set_trace()
        # some backtracking

        #import pdb; pdb.set_trace()
        # YOUR CODE HERE!
        try:
            FM[n][m+1] = log(exp(FD[n][m] +self.t[m]['DM']) +exp(FM[n][m] +self.t[m]['MM']) + exp(FI[n][m] + self.t[m]['IM']))
        except:
            FM[n][m+1] = -INF
        Fscore = FM[n][m+1] 
        #import pdb; pdb.set_trace()       
        return Fscore
    
