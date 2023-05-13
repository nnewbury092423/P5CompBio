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
    #parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument('-m', '--model', required=True, type=str, help="HMM model")
    #parser.add_argument('-q', '--query', required=True, type=str, help="Query sequences")
    #parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Write Viterbi and Forward scores to this file")

    #args = parser.parse_args()

    
    #import pdb; pdb.set_trace()
    directory = os.fsencode('Data')
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        print(filename)
        #import pdb; pdb.set_trace()
        myHMM = 0
        myHMM = HMM()
        myHMM.load('Data/'+filename)
        myHMM = HMM()
        myHMM.load('Data/'+filename)
    #myHMM.load(args.model)
        queries = read_FASTA(open('queries.fas','r'))
        with open('testout2','w') as fout:
            for q in queries:
            #import pdb; pdb.set_trace()
                print(q)
                seq = queries[q]
                Vscore,aln = myHMM.Viterbi(seq)
                Vscore = round(Vscore,5)
                Fscore = round(myHMM.Forward(seq),5)
                #import pdb; pdb.set_trace()
                print(Vscore)
                fout.write(q + " " + str(Vscore) + " " + str(Fscore) + " " + aln + "\n")
