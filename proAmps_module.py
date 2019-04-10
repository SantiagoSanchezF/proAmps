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
                    ap = proAmps.AMPpredictor("test.fasta","proAmps.h5","proAmps_rf.pkl")
                    ap.predict_amps()
                    ap.print_fasta()

                    #AMP predictor given AMP pro-peptides
                    ap = proAmps.AMPpredictor("test.fasta","proAmps.h5","proAmps_rf.pkl")
                    ap.pred_from_propeps("/home/fragoso/machineLearning/sem4/cleaveSite/good_models/N_two_dense_concat_tr-val.h5")
                    ap.print_fasta()
'''

import json, re, getopt, sys
import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
from keras import models
from sklearn.externals import joblib


'''    CLASSES    '''
class AMPpredictor:
    def __init__(self,fasta,multi_rf,model,model_rf):
        self.df = {}
        self.fasta = fasta
        self.multi_rf = multi_rf
        self.model = model
        self.model_rf = model_rf
    def pred_from_propeps(self,modelN_dl,modelN_rf,modelN_svm,modelC_dl,modelC_rf,modelC_svm):
        CS = CSpredictor(self.fasta,modelN_dl,modelN_rf,modelN_svm,modelC_dl,modelC_rf,modelC_svm)
        self.df = CS.df
        self.df["y_pred"], self.df["annot"] = pred_AMP(self.df["mature"],self.multi_rf,self.model,self.model_rf)
        return self.df["y_pred"]
    def predict_amps(self):
        self.df["ids"],self.df["seqs"] = fasta_to_df(self.fasta)
        self.df["y_pred"], self.df["annot"] = pred_AMP(self.df["seqs"],self.multi_rf,self.model,self.model_rf)
        return self.df["y_pred"]
    def print_fasta(self):
        try:
            if self.df["mature"]:
                for i in range(len(self.df["ids"])):
                    print(">" + self.df["ids"][i] + " | " + self.df['mature_annot'][i] 
                            + "|" + self.df["annot"][i]
                            + "\n" + self.df['mature'][i] )
        except:
            for i in range(len(self.df["ids"])):
                print(">" + self.df["ids"][i] + " | " + self.df["annot"][i]
                        + "\n" + self.df['seqs'][i] )

class CSpredictor:
    def __init__(self,fasta,modelN_dl,modelN_rf,modelN_svm,modelC_dl,modelC_rf,modelC_svm):
        self.df = {}
        self.df["ids"],self.df["seqs"] = fasta_to_df(fasta)
        self.clfN = models.load_model(modelN_dl)
        self.modelN_rf = joblib.load(modelN_rf) 
        self.modelN_svm = joblib.load(modelN_svm) 
        self.clfC = models.load_model(modelC_dl)
        self.modelC_rf = joblib.load(modelC_rf) 
        self.modelC_svm = joblib.load(modelC_svm) 
        self.df['mature'],self.df['mature_annot'] = pred_mature_from_list(self.df["seqs"],self.clfN,self.modelN_rf,self.modelN_svm,self.clfC,self.modelC_rf,self.modelC_svm)
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

def pred_AMP(seqs,multi_rf,model,model_rf):        
        pep = importr('Peptides')
        ohe_seqs = encode_seqs(seqs)
        chem_array = pd.DataFrame(calc(seqs,pep)).values
        chem_array -= chem_mean
        chem_array /= chem_std
        model = models.load_model(model)
        model_multi_rf = joblib.load(multi_rf)
        model_rf = joblib.load(model_rf)
        pred_multi_rf = model_multi_rf.predict(chem_array )
        probs_dl = model.predict([chem_array,ohe_seqs])
        probs_dl = probs_dl.reshape(probs_dl.shape[0])
        probs_rf = np.array([x for x in zip(*model_rf.predict_proba(chem_array ))][1])
        probs = 0.61*probs_dl + 0.39*probs_rf 
        annot = ["non-AMP with a probabiity of: " + str(np.round(probs[x],3)) if probs[x] < 0.505 else "predicted AMP with a probabiity of: " + str(np.round(probs[x],3)) + " belonging to cluster:" + str(int(pred_multi_rf[x])) for x in range(len(probs))]
        y_pred = [0 if x < 0.5 else 1 for x in probs]
        return y_pred, annot

## take a sequence and make an np array with OHEing of the octamer, and the noralized position at the end of the vector
def seqs_to_ohe_octamers(seq,mean,std):
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
        norm_p = i - mean
        norm_p /= std
        octas.append(np.append(results[i],norm_p))
    octas = np.array(octas)
    return octas

## makes prediction, returns the position of CS and the probaility 
def ohe_prediction(octas,model_dl,model_rf,model_svm,wei): 
    y_pred_dl = model_dl.predict(octas) if model_dl else []
    y_pred_dl = y_pred_dl.reshape(y_pred_dl.shape[0]) if model_dl else []
    y_pred_rf = np.array([x for x in zip(*model_rf.predict_proba(octas))][1]) if model_rf else []
    y_pred_svm = np.array([x for x in zip(*model_svm.predict_proba(octas))][1]) if model_svm else  []
    y_pred = wei[0]*y_pred_dl + wei[1]*y_pred_rf + wei[2]*y_pred_svm if wei else None
    if wei:
        return np.argmax(y_pred),np.max(y_pred)
    else:
        return np.argmax([y_pred_dl,y_pred_rf,y_pred_svm][np.argmax([len(y_pred_dl),len(y_pred_rf),len(y_pred_svm)])]),np.max([y_pred_dl,y_pred_rf,y_pred_svm][np.argmax([len(y_pred_dl),len(y_pred_rf),len(y_pred_svm)])])
 

## predict the mature sequence and return annotation    
def pred_mature(seq,modelN,modelN_rf,modelN_svm,modelC,modelC_rf,modelC_svm):
    ohe_N = seqs_to_ohe_octamers(seq,43.98830865159782,82.4440045905732)
    ohe_C = seqs_to_ohe_octamers(seq,35.84697508896797, 56.168662322490874)
    #pos_N,prob_N  = ohe_prediction(ohe_N,modelN,modelN_rf,modelN_svm,(0.5,0.1,0.4))
    pos_N,prob_N  = ohe_prediction(ohe_N,0,0,modelN_svm,0)
    pos_C,prob_C  = ohe_prediction(ohe_C,modelC,modelC_rf,modelC_svm,(0.02, 0.89, 0.09))
    
   # pos_N,prob_N,pos_C,prob_C = ohe_prediction(seqs_to_ohe_octamers(seq),modelN,modelC)
    if prob_N >= 0.6 and len(seq[pos_N+4:]) > 4:
        if prob_C >= 0.6 and len(seq[pos_N+4:pos_C+4]) > 4:
            return seq[pos_N+4:pos_C+4],"N-terminus CS predicted at position: " + str(pos_N+4) + ", with a probability of: " + str(np.round(prob_N,3)) + "|" + "C-terminus CS predicted at position: " + str(pos_C+4) + ", with a probability of: " + str(np.round(prob_C,3))
        else:    
            return seq[pos_N+4:],"N-terminus CS predicted at position: " + str(pos_N+4) + ", with a probability of: " + str(np.round(prob_N,3)) + "| No C-terminus CS predicted"
    elif prob_C >= 0.6 and len(seq[:pos_C+4]) > 4:
        return seq[:pos_C+4],"No N-terminus CS predicted | C-terminus CS predicted at position: " + str(pos_C+4) + ", with a probability of: " + str(np.round(prob_C,3))
    else:
        return seq,"No N-terminus CS predicted | No C-terminus CS predicted"

def pred_mature_from_list(seqs,modelN,modelN_rf,modelN_svm,modelC,modelC_rf,modelC_svm):
    mature, mature_annot = [],[]
    for seq in seqs:
        m,a = pred_mature(seq,modelN,modelN_rf,modelN_svm,modelC,modelC_rf,modelC_svm)
        mature.append(m)
        mature_annot.append(a)
    return mature, mature_annot

## calulate the chemechical propierties of a list of sequences. uses R Peptides package and rpy2 library
def calc(seqs,pep):
    seqs = [x if len(x) < 102 else "WWWW" for x in seqs]
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

chem_mean = np.array([8.42277739e+00, 8.30251263e+01, 1.07009223e+01, 1.47836489e+01,
       1.35775798e+00, 9.21663941e+00, 2.32064232e+01, 5.21321870e-01,
       4.94085645e-01, 3.94308015e+01, 5.90153175e+03, 5.71550613e+01,
       8.09864228e+00, 5.13415030e+01, 3.18229247e+01])

chem_std = np.array([6.41429619e+00, 3.80489830e+01, 5.84315338e+00, 8.87892116e+00,
       1.37619010e+00, 6.96499326e+00, 1.06032731e+01, 1.75051327e-01,
       1.62239968e-01, 2.41712756e+01, 3.26236859e+03, 1.17307177e+01,
       2.73216063e+00, 1.27165049e+01, 1.19224287e+01])

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

