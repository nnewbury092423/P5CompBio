#! /usr/bin/env python

from HMM import *
import argparse
import os

def read_FASTA(stream):
    '''
    This function reads a FASTA file from a given stream and returns a dictionary mapping identifiers to sequences
    '''
    seqs = {}; name = None; seq = ''
    for line in stream:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                assert len(seq) != 0, "Malformed FASTA"
                seqs[name] = seq
            name = l[1:]
            assert name not in seqs, "Duplicate sequence ID: %s" % name
            seq = ''
        else:
            seq += l
    assert name is not None and len(seq) != 0, "Malformed FASTA"
    seqs[name] = seq
    return seqs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument('-m', '--model', required=True, type=str, help="HMM model")
    #parser.add_argument('-q', '--query', required=True, type=str, help="Query sequences")
    #parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Write Viterbi and Forward scores to this file")

    #args = parser.parse_args()

    
    directory = os.fsencode('Data')
    queries = read_FASTA(open('queries.fas','r'))
    finalaln = ''
    with open('query_class','w') as fout:
        for q in queries:
            print(q)
            seq = queries[q]
            bestVscore= -INF
            bestFscore = -INF
            Vmodel = ''
            Fmodel = ''
            for file in os.listdir(directory):
                filename = os.fsdecode(file)
                print(filename)
                #import pdb; pdb.set_trace()
                #idk why i needed to do this, otherwise it wasnt changing models
                myHMM = 0
                myHMM = HMM()
                myHMM.load('Data/'+ filename)
                myHMM = HMM()
                myHMM.load('Data/'+ filename)    

                Vscore,aln = myHMM.Viterbi(seq)
                Vscore = round(Vscore,5)
                Fscore = round(myHMM.Forward(seq),5)
                print(Vscore)
                if Vscore >= bestVscore:
                    bestVscore = Vscore
                    Vmodel = filename
                    finalaln = aln
                if Fscore >= bestFscore:
                    bestFscore = Fscore
                    Fmodel = filename
            #import pdb; pdb.set_trace()
            fout.write('query name: '+ q + ' Assigned model: ' + Fmodel + " Best score(V,F): " + str(bestVscore) + ", " + str(bestFscore)+ " Aligned sequence: " +  finalaln + "\n")
