'''            -- Jose Santiago Sanchez Fragoso --    2019

    proAmps   - Antiicrobial propetide cleavage site predictor and Antimicrobial peptides classificator
    
    DEPPENDENCIES:
                    - R < 3.x
                    - R package "Peptides"
                    - SciKitlearn
                    - rpy2
                    - Pandas
                    - Keras
    EXAMPLES:
                    import proAmps

                    #Cleavage site predictor
                    cs = proAmps.CSpredictor("test.fasta","/home/fragoso/machineLearning/sem4/cleaveSite/good_models/N_two_dense_concat_tr-val.h5")
                    cs.print_fasta()

                    #AMP predictor given mature sequences
                    ap = proAmps.AMPpredictor("test.fasta","/home/fragoso/machineLearning/sem4/deepLearning/LAST_acc.h5")
                    ap.predict_amps()
                    ap.print_fasta()

					#AMP predictor given AMP pro-peptides
					ap = proAmps.AMPpredictor("test.fasta","/home/fragoso/machineLearning/sem4/deepLearning/LAST_acc.h5")
					ap.pred_from_propeps("/home/fragoso/machineLearning/sem4/cleaveSite/good_models/N_two_dense_concat_tr-val.h5")
                    ap.print_fasta()
'''

import json, re, getopt, sys
import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
from keras import models


'''    CLASSES    '''
class AMPpredictor:
    def __init__(self,fasta,model):
        self.df = {}
        self.fasta = fasta
        self.model = model
    def pred_from_propeps(self,modelN):
        CS = CSpredictor(self.fasta,modelN)
        self.df = CS.df
        self.df["y_pred"], self.df["annot"] = pred_AMP(self.df["mature"],self.model)
        return self.df["y_pred"]
    def predict_amps(self):
        self.df["ids"],self.df["seqs"] = fasta_to_df(self.fasta)
        self.df["y_pred"], self.df["annot"] = pred_AMP(self.df["seqs"],self.model)
        return self.df["y_pred"]
    def print_fasta(self):
        try:
            if self.df["mature"]:
                for i in range(len(self.df["ids"])):
                    print(">" + self.df["ids"][i] + " | " + self.df['mature_annot'][i] 
                            + " | " + self.df["annot"][i]
                            + "\n" + self.df['mature'][i] )
        except:
            for i in range(len(self.df["ids"])):
                print(">" + self.df["ids"][i] + " | " + self.df["annot"][i]
                        + "\n" + self.df['seqs'][i] )

class CSpredictor:
    def __init__(self,fasta,modelN):
        self.df = {}
        self.df["ids"],self.df["seqs"] = fasta_to_df(fasta)
    #def predict_mature_seqs(self,modelN):
        self.clfN = models.load_model(modelN)
        self.df['mature'],self.df['mature_annot'] = pred_mature_from_list(self.df["seqs"],self.clfN)
    #    return self.df['mature'], self.df['mature_annot']
    def print_fasta(self):
        for i in range(len(self.df["ids"])):
            print(">" + self.df["ids"][i] + " | " + self.df['mature_annot'][i] 
                    + "\n" + self.df['mature'][i] )


'''    FUNCTIONS    '''
def encode_seqs(seqs):
    seqs = [x if len(x) < 102 else "" for x in seqs]
    characters = 'ACEDGFIHKMLNQPSRTWVY'
    token_index = dict(zip(characters,range(1, len(characters) + 1)))
    max_length = 102
    results = np.zeros((len(seqs), max_length, max(token_index.values()) + 1))
    for i, seqs in enumerate(seqs):
        for j, character in enumerate(seqs):
            index = token_index.get(character)
            results[i, j, index] = 1.
    return results

def pred_AMP(seqs,model):        
        pep = importr('Peptides')
        ohe_seqs = encode_seqs(seqs)
        chem_array = pd.DataFrame(calc(seqs,pep)).values
        chem_array -= chem_mean
        chem_array /= chem_std
        model = models.load_model(model)
        probs = model.predict([chem_array,ohe_seqs])
        annot = [" no-AMP " if x < 0.5 else "predicted AMP with a probabiity of: " + str(np.round(x,3)) for x in probs]
        y_pred = [0 if x < 0.5 else 1 for x in probs]
        return y_pred, annot

## take a sequence and make an np array with OHEing of the octamer, and the noralized position at the end of the vector
def seqs_to_ohe_octamers(seq):
    octas = []
    for i in range(len(seq)-7):
        octas.append(seq[i:i+8])
    characters = 'ACEDGFIHKMLNQPSRTWVY'
    token_index = dict(zip(characters,range(1, len(characters) + 1)))
    max_length = 8
    results = np.zeros((len(octas), max_length, max(token_index.values()) + 1))
    for i, octas in enumerate(octas):
        for j, character in enumerate(octas):
            index = token_index.get(character)
            results[i, j, index] = 1.
    results = results.reshape(results.shape[0],8*21)
    octas = []
    for i in range(len(results)):
        norm_p = i - 41.29971325271163
        norm_p /= 45.263206049562065 
        octas.append(np.append(results[i],norm_p))
    octas = np.array(octas)
    return octas

## makes prediction, returns the position of CS and the probaility 
def ohe_prediction(octas,model):
    y_pred = model.predict(octas)                             
    return np.argmax(y_pred),np.max(y_pred)  

## predict the mature sequence and return annotation    
def pred_mature(seq,model):
    pos,prob = ohe_prediction(seqs_to_ohe_octamers(seq),model)
    if prob >= 0.5:
        return seq[pos+4:],"N-terminus CS predicted at position: " + str(pos+4) + ", with a probability of: " + str(np.round(prob,3)) 
    else:
        return ""," No N-terminus CS predicted "

def pred_mature_from_list(seqs,model):
    mature, mature_annot = [],[]
    for seq in seqs:
        m,a = pred_mature(seq,model)
        mature.append(m)
        mature_annot.append(a)
    return mature, mature_annot

## calulate the chemechical propierties of a list of sequences. uses R Peptides package and rpy2 library
def calc(seqs,pep):
    df = {}
    df['mw'] = []
    df['tiny'] = []
    df['small'] = []
    df['aro'] = []
    df['nonP'] = []
    df['chrgd'] = []
    df['basic'] = []
    df['acidic'] = []
    df['charge'] = []
    df['pI'] = []
    df['aindex'] = []
    df['instaindex'] = []
    df['boman'] = []
    df['hmoment100'] = []
    df['hmoment160'] = []
    for i in seqs:
        df['mw'].append(pep.mw(i)[0])
        df['tiny'].append(pep.aaComp(i)[0][-9])
        df['small'].append(pep.aaComp(i)[0][-8])
        df['aro'].append(pep.aaComp(i)[0][-6])
        df['nonP'].append(pep.aaComp(i)[0][-5])
        df['chrgd'].append(pep.aaComp(i)[0][-3])
        df['basic'].append(pep.aaComp(i)[0][-2])
        df['acidic'].append(pep.aaComp(i)[0][-1])
        df['charge'].append(pep.charge(i,pH=0,pKscale = "EMBOSS")[0])
        df['pI'].append(pep.pI(i,pKscale = "EMBOSS")[0])
        df['aindex'].append(pep.aIndex(i)[0])
        df['instaindex'].append(pep.instaIndex(i)[0])
        df['boman'].append(pep.boman(i)[0])
        df['hmoment100'].append(pep.hmoment(i, angle = 100, window = 11)[0])
        df['hmoment160'].append(pep.hmoment(i, angle = 160, window = 11)[0])
    return df

chem_mean = np.array([1.03169026e+01, 7.88236788e+01, 1.00331206e+01, 1.47364910e+01,
       1.50721927e+00, 1.05481452e+01, 2.50533831e+01, 5.64771677e-01,
       5.32652144e-01, 3.99154949e+01, 7.15002626e+03, 5.59987213e+01,
       7.67410228e+00, 5.18110130e+01, 3.21122640e+01])

chem_std = np.array([6.32922159e+00, 3.14863376e+01, 5.29169301e+00, 6.87687902e+00,
       1.06299785e+00, 5.67288978e+00, 9.06730592e+00, 1.60777106e-01,
       1.51180417e-01, 2.13839052e+01, 2.97206228e+03, 9.74206925e+00,
       2.49596851e+00, 1.11045765e+01, 1.04242978e+01])

## take a fata and extract ids and sequences
def fasta_to_df(fasta):
    df = {}
    df['seqs'] = []
    df['ids'] = []
    ident = None
    seq = ""
    with open(fasta,"r") as myf:
        lines = myf.read().splitlines()
        while(not lines[-1].strip()):
            lines.pop()
        last_line = lines[-1]
        for line in range(len(lines)):
	    if len(lines[line]) < 1:
	        continue
            i = lines[line].split()[0]
            if re.search("^>",i):
                if ident and not len(seq) == 0:
                    if len(seq) < 8 :
                        None
                    else:
                        df['ids'].append(ident)
                        df['seqs'].append(seq)

                ident = i.split()[0].split(">")[1]
                seq = ""
            if not re.search("^>",i) and not re.search("[^ACEDGFIHKMLNQPSRTWVY]",i):

                seq = seq + i.strip()
            if line == len(lines) -1:
                df['ids'].append(ident)
                df['seqs'].append(seq)
    if len(df['seqs'][-1]) == 0:
        df['seqs'].pop()
        df['ids'].pop()
    return df['ids'],df['seqs']

