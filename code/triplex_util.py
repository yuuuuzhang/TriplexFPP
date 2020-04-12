import numpy as np
from Bio import SeqIO
import math
import pandas as pd
import csv
from keras.models import load_model
import pickle


def countnum(seq,nuacid):
    return len([1 for i in range(len(seq)) if seq.startswith(nuacid,i)])

def construct_kmer():
	ntarr = ("A","C","G","T")

	kmerArray = []


	for n in range(4):
		kmerArray.append(ntarr[n])

	for n in range(4):
		str1 = ntarr[n]
		for m in range(4):
			str2 = str1 + ntarr[m]
			kmerArray.append(str2)
#############################################
	for n in range(4):
		str1 = ntarr[n]
		for m in range(4):
			str2 = str1 + ntarr[m]
			for x in range(4):
				str3 = str2 + ntarr[x]
				kmerArray.append(str3)
#############################################
#change this part for 3mer or 4mer
	for n in range(4):
		str1 = ntarr[n]
		for m in range(4):
			str2 = str1 + ntarr[m]
			for x in range(4):
				str3 = str2 + ntarr[x]
				for y in range(4):
					str4 = str3 + ntarr[y]
					kmerArray.append(str4)
############################################
	for n in range(4):
		str1 = ntarr[n]
		for m in range(4):
			str2 = str1 + ntarr[m]
			for x in range(4):
				str3 = str2 + ntarr[x]
				for y in range(4):
					str4 = str3 + ntarr[y]
					for z in range(4):
						str5 = str4 + ntarr[z]
						kmerArray.append(str5)
####################### 6-mer ##############
	for n in range(4):
		str1 = ntarr[n]
		for m in range(4):
			str2 = str1 + ntarr[m]
			for x in range(4):
				str3 = str2 + ntarr[x]
				for y in range(4):
					str4 = str3 + ntarr[y]
					for z in range(4):
						str5 = str4 + ntarr[z]
						for t in range(4):
							str6 = str5 + ntarr[t]
							kmerArray.append(str6)
    
	return kmerArray

def kmer_encode(seq,kmerarray):
    result = np.zeros((len(seq),len(kmerarray)))
    for i in range(len(seq)):
        for j in range(len(kmerarray)):
            result[i,j] = seq[i].count(kmerarray[j])/len(seq[i])
    return result

def mer_sin(seq,nc_m,c_m,kmerarray,x):   
    
    l = len(seq)-x+1
    log_r = np.zeros((l))
    for i in range(l):
        tempseq = seq[i:i+x]
        idx = kmerarray.index(tempseq)
        Fc = c_m[int(idx)]
        Fnc = nc_m[int(idx)]
        if Fc==0 and Fnc==0:
            log_r[i]=0
        elif Fc==0 and Fnc!=0:
            log_r[i]=-1
        elif Fnc==0 and Fc!=0:
            log_r[i]=1
        else:
            log_r[i] = math.log(Fc/Fnc)
    miu = sum(log_r)/l
    
    return miu
   
def mer_score(seq,nc_m,c_m,kmerarray,x):
    miu = np.zeros((len(seq)))
    for i in range(len(seq)):
        miu[i] = mer_sin(seq[i],nc_m,c_m,kmerarray,x)
        
    miu0 = np.expand_dims(miu, axis=1)
    return miu0


def triplex_lncRNA(inputfile,datapath,outputname):
    
    records = list(SeqIO.parse(inputfile, "fasta")) 

    seq = []
    name=[]
    for j in range(len(records)):
        seq.append(records[j].seq)
        name.append(records[j].id)
        

    kmerArray = construct_kmer() 
    kmer1 = kmer_encode(seq,kmerArray[0:4])
    kmer2 = kmer_encode(seq,kmerArray[4:20])
    kmer3 = kmer_encode(seq,kmerArray[20:84])
    
    pos = pd.read_csv(datapath+'embed/mer_rnapos_mean.csv',header=None,delimiter = ',')
    pos = np.array(pos)
    neg = pd.read_csv(datapath+'embed/mer_rnaneg_mean.csv',header=None,delimiter = ',')
    neg = np.array(neg)
    
    merscore1 = mer_score(seq,pos[0:4],neg[0:4],kmerArray[0:4],1)
    merscore2 = mer_score(seq,pos[4:20],neg[4:20],kmerArray[4:20],2)
    merscore3 = mer_score(seq,pos[20:84],neg[20:84],kmerArray[20:84],3)
    merscore4 = mer_score(seq,pos[84:340],neg[84:340],kmerArray[84:340],4)
    merscore5 = mer_score(seq,pos[340:1364],neg[340:1364],kmerArray[340:1364],5)
    merscore6 = mer_score(seq,pos[1364:5460],neg[1364:5460],kmerArray[1364:5460],6)
    
    fea = np.concatenate((merscore1,merscore2,merscore3,merscore4,merscore5,merscore6,kmer1,kmer2,kmer3),axis=1)
    testdata = np.expand_dims(fea, axis=2) 
    bs = 256
    model = load_model(datapath+'embed/triplexlncRNA.h5')
    probs = model.predict(testdata)
    
    
    label=[]
    for i in range(len(probs)):
        if (probs[i,0]>0.49):
            label.append('Triplex')         
        else:
            label.append('Nontriplex')
    
    with open(datapath+'output/'+outputname, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(zip(name,list(probs[:,0]),label))
    
 
    return 

def triplex_DNA(inputfile,datapath,outputname):
    
    records = list(SeqIO.parse(inputfile, "fasta")) 

    seq = []
    name=[]
    for j in range(len(records)):
        seq.append(records[j].seq)
        name.append(records[j].id)
        

    kmerArray = construct_kmer() 
    kmer1 = kmer_encode(seq,kmerArray[0:4])
    kmer2 = kmer_encode(seq,kmerArray[4:20])
    kmer3 = kmer_encode(seq,kmerArray[20:84])
    
    pos = pd.read_csv(datapath+'embed/mer_dnapos_mean.csv',header=None,delimiter = ',')
    pos = np.array(pos)
    neg = pd.read_csv(datapath+'embed/mer_dnaneg_mean.csv',header=None,delimiter = ',')
    neg = np.array(neg)
    
    merscore1 = mer_score(seq,pos[0:4],neg[0:4],kmerArray[0:4],1)
    merscore2 = mer_score(seq,pos[4:20],neg[4:20],kmerArray[4:20],2)
    merscore3 = mer_score(seq,pos[20:84],neg[20:84],kmerArray[20:84],3)
    merscore4 = mer_score(seq,pos[84:340],neg[84:340],kmerArray[84:340],4)
    merscore5 = mer_score(seq,pos[340:1364],neg[340:1364],kmerArray[340:1364],5)
    merscore6 = mer_score(seq,pos[1364:5460],neg[1364:5460],kmerArray[1364:5460],6)
    
    fea = np.concatenate((merscore1,merscore2,merscore3,merscore4,merscore5,merscore6,kmer1,kmer2,kmer3),axis=1)
    testdata = np.expand_dims(fea, axis=2) 
    bs = 256
    model = load_model(datapath+'embed/triplexDNA.h5')
    probs = model.predict(testdata)
    
    
    label=[]
    for i in range(len(probs)):
        if (probs[i,0]>0.489):
            label.append('Triplex')         
        else:
            label.append('Nontriplex')
    
    with open(datapath+'output/'+outputname, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(zip(name,list(probs[:,0]),label))
    return

def cis_trans(inputfile,datapath,outputname):
    
    records = list(SeqIO.parse(inputfile, "fasta")) 

    seq = []
    name=[]
    for j in range(len(records)):
        seq.append(records[j].seq)
        name.append(records[j].id)
        

    kmerArray = construct_kmer() 
    
    pos = pd.read_csv(datapath+'embed/mer_cis_mean.csv',header=None,delimiter = ',')
    pos = np.array(pos)
    neg = pd.read_csv(datapath+'embed/mer_trans_mean.csv',header=None,delimiter = ',')
    neg = np.array(neg)
    
    merscore1 = mer_score(seq,pos[0:4],neg[0:4],kmerArray[0:4],1)
    merscore2 = mer_score(seq,pos[4:20],neg[4:20],kmerArray[4:20],2)
    merscore3 = mer_score(seq,pos[20:84],neg[20:84],kmerArray[20:84],3)
    merscore4 = mer_score(seq,pos[84:340],neg[84:340],kmerArray[84:340],4)
    merscore5 = mer_score(seq,pos[340:1364],neg[340:1364],kmerArray[340:1364],5)
    merscore6 = mer_score(seq,pos[1364:5460],neg[1364:5460],kmerArray[1364:5460],6)
    
    testdata = np.concatenate((merscore1,merscore2,merscore3,merscore4,merscore5,merscore6),axis=1)
    
    svc = pickle.load(open(datapath+'embed/cis_trans.sav', 'rb'))

    probs = svc.predict_proba(testdata)
    
    
    label=[]
    for i in range(len(probs)):
        if (probs[i,0]>0.45):
            label.append('cis&trans')         
        else:
            label.append('trans')
    
    with open(datapath+'output/'+outputname, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(zip(name,list(probs[:,0]),label))
    return