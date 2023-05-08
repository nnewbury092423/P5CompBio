from math import *

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
        
        Vscore = -INF
        aln = "$"
       
        # YOUR CODE HERE
        
        Vscore = VM[n][m+1]
        return Vscore,aln        


    def Forward(self,query):
        # implement the Forward algorithm to compute the marginal probability of a query sequence
        n = len(query)
        m = self.nstate
        
        FM = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        FI = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        FD = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        
        # YOUR CODE HERE!

        Fscore = FM[n][m+1]        
        return Fscore
    
