#########################################################################################
# CDPred is developed for predicting, desigining and scanning celiac disease causing    #
# peptides using sequence information. It is developed by Prof G. P. S. Raghava's group.#
# Please cite: https://webs.iiitd.edu.in/raghava/cdpred/                           #
#########################################################################################
import argparse
import warnings
import subprocess
import pkg_resources
import os
import sys
import numpy as np
import pandas as pd
import math
import itertools
from collections import Counter
from itertools import combinations
import pickle
import re
import glob
import time
import uuid
from time import sleep
from tqdm import tqdm
from sklearn.ensemble import ExtraTreesClassifier
import zipfile
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='Please provide following arguments') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence(s) in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
parser.add_argument("-m", "--method",type=int, choices = [1,2], help="Method Type: 1:Main (Disease Vs Random), 2:Alternate (Disease Vs Non-Disease),  by default 1")
parser.add_argument("-j", "--job",type=int, choices = [1,2,3,4,5], help="Job Type: 1:Predict, 2:Design, 3:Scan, 4:PQ Density, 5: Ensemble Method,  by default 1")
parser.add_argument("-t","--threshold", type=float, help="Threshold: Value between 0 to 1 by default 0.32")
parser.add_argument("-w","--winleng", type=int, choices =range(9, 21), help="Window Length: 8 to 20 (scan mode only), by default 9")
parser.add_argument("-p","--patleng", type=int, choices =range(3, 10), help="Pattern Length for PQ Density Module: 3 to 9, by default 5")
parser.add_argument("-d","--display", type=int, choices = [1,2], help="Display: 1:Only Disease-Causing Peptides, 2: All peptides, by default 1")
args = parser.parse_args()

# Function for generating all possible mutants
def mutants(file1,file2):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    cc = []
    dd = []
    ee = []
    df2 = file2
    df2.columns = ['Name']
    df1 = file1
    df1.columns = ['Seq']
    for k in range(len(df1)):
        cc.append(df1['Seq'][k])
        dd.append('Original_'+'Seq'+str(k+1))
        ee.append(df2['Name'][k])
        for i in range(0,len(df1['Seq'][k])):
            for j in std:
                if df1['Seq'][k][i]!=j:
                    dd.append('Mutant_'+df1['Seq'][k][i]+str(i+1)+j+'_Seq'+str(k+1))
                    cc.append(df1['Seq'][k][:i] + j + df1['Seq'][k][i + 1:])
                    ee.append(df2['Name'][k])
    xx = pd.concat([pd.DataFrame(ee),pd.DataFrame(dd),pd.DataFrame(cc)],axis=1)
    xx.columns = ['Seq_ID','Mutant_ID','Seq']
    return xx
# Function for generating pattern of a given length
def seq_pattern(file1,file2,num):
    df1 = file1
    df1.columns = ['Seq']
    df2 = file2
    df2.columns = ['Name']
    cc = []
    dd = []
    ee = []
    for i in range(len(df1)):
        for j in range(len(df1['Seq'][i])):
            xx = df1['Seq'][i][j:j+num]
            if len(xx) == num:
                cc.append(df2['Name'][i])
                dd.append('Pattern_'+str(j+1)+'_Seq'+str(i+1))
                ee.append(xx)
    df3 = pd.concat([pd.DataFrame(cc),pd.DataFrame(dd),pd.DataFrame(ee)],axis=1)
    df3.columns= ['Seq_ID','Pattern_ID','Seq']
    return df3
# Function to check the seqeunce
def readseq(file):
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    seqid = []
    seq = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '', ''.join(array[1:]).upper())
        seqid.append('>'+name)
        seq.append(sequence)
    if len(seqid) == 0:
        f=open(file,"r")
        data1 = f.readlines()
        for each in data1:
            seq.append(each.replace('\n',''))
        for i in range (1,len(seq)+1):
            seqid.append(">Seq_"+str(i))
    df1 = pd.DataFrame(seqid)
    df2 = pd.DataFrame(seq)
    return df1,df2
# Function to check the length of seqeunces
def lenchk(file1):
    cc = []
    df1 = file1
    df1.columns = ['seq']
    for i in range(len(df1)):
        if len(df1['seq'][i])>20:
            cc.append(df1['seq'][i][0:20])
        else:
            cc.append(df1['seq'][i])
    df2 = pd.DataFrame(cc)
    df2.columns = ['Seq']
    return df2
# Function to generate the features out of seqeunces
def feature_gen(file):
    std = list('ACDEFGHIKLMNPQRSTVWY')
    df1 = file
    dd = []
    for j in df1['Seq']:
        cc = []
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            cc.append(composition)
        dd.append(cc)
    df2 = pd.DataFrame(dd)
    head = []
    for mm in std:
        head.append('AAC_'+mm)
    df2.columns = head
    return df2
# Function to predict using PQ density
def PQ_abd(file,K):
    dd = []
    df1 = file
    for i in range(len(df1)):
        test_str = df1.Seq[i]
        res = [test_str[x:y] for x, y in combinations(range(len(test_str) + 1), r = 2) if len(test_str[x:y]) == K ]
        cc = []
        for j in res:
            cc.append((j.count('P')+j.count('Q'))/len(j))
        dd.append(max(cc))
    df2 = pd.DataFrame(dd)
    return df2
# Function to run ensemble method
def aac_comp_str(file):
    std = list('ACDEFGHIKLMNPQRSTVWY')
    cc = []
    for i in std:
        count = 0
        for k in file:
            temp1 = k
            if temp1 == i:
                count += 1
            composition = (count/len(file))*100
        cc.append(composition)
    df2 = pd.DataFrame(cc)
    df3 = df2.T
    head = []
    for mm in std:
        head.append('AAC_'+mm)
    df3.columns = head
    return df3
def ensmbl(file_1,file_2,thr,meth):
    file1 = pd.read_csv(file_1, header=None) #Motif file
    motifs = list(file1[0])
    df,df_ = readseq(file_2)
    df.columns = ['Seq_ID']
    df_.columns = ['Seq']
    file2 = df_
    aa = []
    cc = []
    dd = []
    ee = []
    ff = []
    if meth == 1:
        clfmain = pickle.load(open('model/main_model.pkl','rb'))
    if meth == 2:
        clfmain = pickle.load(open('model/alt_model.pkl','rb'))
    for s in file2.Seq:
        if any((c in s) for c in motifs):
            cc.append('Hit found')
            ee.append('-')
            dd.append('Disease Causing')
        else:
            cc.append('No hit found')
            xx = aac_comp_str(s)
            y_p_score1=clfmain.predict_proba(xx)
            ee.append(y_p_score1[0][1].round(2))
            if y_p_score1[0][1]>=float(thr):
                dd.append('Disease Causing')
            else:
                dd.append('Disease non-causing')
    df123 = pd.concat([df,df_,pd.DataFrame(cc),pd.DataFrame(ee),pd.DataFrame(dd)],axis=1)
    df123.Seq_ID = df123.Seq_ID.str.replace('>','')
    return df123
# Function to read and implement the model
def model_run(file1,file2):
    a = []
    data_test = file1
    clf = pickle.load(open(file2,'rb'))
    y_p_score1=clf.predict_proba(data_test)
    y_p_s1=y_p_score1.tolist()
    a.extend(y_p_s1)
    df = pd.DataFrame(a)
    df1 = df.iloc[:,-1].round(2)
    df2 = pd.DataFrame(df1)
    df2.columns = ['ML_score']
    return df2
# To Determine the label
def det_lab(file1,thr):
    df1 = file1
    df1['label'] = ['Disease Causing' if df1['ML_score'][i]>=float(thr) else 'Disease non-causing' for i in range(len(df1))]
    return df1
    
('############################################################################################')
print('# This program CDPred is developed for predicting, desigining and scanning celiac disease #')
print('# causing peptides, developed by Prof G. P. S. Raghava group.        #')
print('# Please cite: CDPred; available at https://webs.iiitd.edu.in/raghava/cdpred/  #')
print('############################################################################################')

# Parameter initialization or assigning variable for command level arguments

Sequence= args.input        # Input variable 
 
# Output file 
 
if args.output == None:
    result_filename= "outfile.csv" 
else:
    result_filename = args.output
         
# Method
if args.method == None:
        methd = 1
else:
        methd= int(args.method)
# Job Type 
if args.job == None:
        Job = int(1)
else:
        Job = int(args.job)
# Window Length 
if args.winleng == None:
        Win_len = int(9)
else:
        Win_len = int(args.winleng)
# Patter Length
if args.patleng == None:
        Pat_len = int(5)
else:
        Pat_len = int(args.patleng)
# Display
if args.display == None:
        dplay = int(1)
else:
        dplay = int(args.display)
# Threshold
if Job ==4:
    if args.threshold == None:
        Threshold = 0.41
    else:
        Threshold= float(args.threshold)
else:
    if args.threshold == None:
        Threshold = 0.32
    else:
        Threshold= float(args.threshold)

###########################################################################################
if Job==2:
    print("\n");
    print('##############################################################################')
    print('Summary of Parameters:')
    print('Input File: ',Sequence,'; Threshold: ', Threshold,'; Job Type: ',Job,'; Method Name: ',methd)
    print('Output File: ',result_filename,'; Window Length: ',Win_len,'; Display: ',dplay)
    print('##############################################################################')
if Job==4:
    print("\n");
    print('##############################################################################')
    print('Summary of Parameters:')
    print('Input File: ',Sequence,'; Threshold: ', Threshold,'; Job Type: ',Job,'; Method Name: ',methd)
    print('Output File: ',result_filename,'; Pattern Length: ',Pat_len,'; Display: ',dplay)
    print('##############################################################################')
else:
    print("\n");
    print('##############################################################################')
    print('Summary of Parameters:')
    print('Input File: ',Sequence,'; Threshold: ', Threshold,'; Job Type: ',Job)
    print('Output File: ',result_filename,'; Display: ',dplay,'; Method Name: ',methd)
    print('# ############################################################################')
#========================================Extracting Model====================================
if os.path.isdir('model') == False:
    with zipfile.ZipFile('./model.zip', 'r') as zip_ref:
        zip_ref.extractall('.')
else:
    pass
#======================= Prediction Module start from here =====================
if Job == 1:
    print('\n======= Thanks for using Predict module of CDPred. Your results will be stored in file :',result_filename,' =====\n')
    df_2,dfseq = readseq(Sequence)
    df1 = lenchk(dfseq)
    X = feature_gen(df1)
    if methd == 1:
        mlres = model_run(X,'model/main_model.pkl')
    else:
        mlres = model_run(X,'model/alt_model.pkl')
    mlres = mlres.round(3)
    df45 = det_lab(mlres,Threshold)
    df44 = pd.concat([df_2,dfseq,df45],axis=1)
    df44.columns = ['Seq_ID','Sequence','ML_Score','Prediction']
    df44.Seq_ID = df44.Seq_ID.str.replace('>','')
    if dplay == 1:
        df44 = df44.loc[df44.Prediction=="Disease Causing"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")    
#===================== Design Model Start from Here ======================
elif Job == 2:
    print('\n======= Thanks for using Design module of CDPred. Your results will be stored in file :',result_filename,' =====\n')
    print('==== Designing Peptides: Processing sequences please wait ...')
    df_2,dfseq = readseq(Sequence)
    df1 = lenchk(dfseq)
    df_1 = mutants(df1,df_2)
    dfseq = df_1[['Seq']]
    X = feature_gen(dfseq)
    if methd == 1:
        mlres = model_run(X,'model/main_model.pkl')
    else:
        mlres = model_run(X,'model/alt_model.pkl')
    mlres = mlres.round(3)
    df45 = det_lab(mlres,Threshold)
    df45.columns = ['ML_Score','Prediction']
    df44 = pd.concat([df_1,df45],axis=1)
    df44['Mutant_ID'] = ['_'.join(df44['Mutant_ID'][i].split('_')[:-1]) for i in range(len(df44))]
    df44.drop(columns=['Seq_ID'],inplace=True)
    df44['Seq_ID'] = [i.replace('>','') for i in df_1['Seq_ID']]
    df44['Sequence'] = df_1.Seq
    df44 = df44[['Seq_ID','Mutant_ID','Sequence','ML_Score','Prediction']]
    df44.Seq_ID = df44.Seq_ID.str.replace('>','')
    if dplay == 1:
        df44 = df44.loc[df44.Prediction=="Disease Causing"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")
#=============== Scan Model start from here ==================
elif Job==3:
    print('\n======= Thanks for using Scan module of CDPred. Your results will be stored in file :',result_filename,' =====\n')
    print('==== Scanning Peptides: Processing sequences please wait ...')
    df_2,dfseq = readseq(Sequence)
    df_1 = seq_pattern(dfseq,df_2,Win_len)
    dfseq = df_1[['Seq']]
    X = feature_gen(dfseq)
    if methd == 1:
        mlres = model_run(X,'model/main_model.pkl')
    else:
        mlres = model_run(X,'model/alt_model.pkl')
    mlres = mlres.round(3)
    df45 = det_lab(mlres,Threshold)
    df45.columns = ['ML_Score','Prediction']
    df44 = pd.concat([df_1,df45],axis=1)
    df44['Pattern_ID'] = ['_'.join(df44['Pattern_ID'][i].split('_')[:-1]) for i in range(len(df44))]
    df44.drop(columns=['Seq_ID'],inplace=True)
    df44['Seq_ID'] = [i.replace('>','') for i in df_1['Seq_ID']]
    df44['Sequence'] = df_1.Seq
    df44 = df44[['Seq_ID','Pattern_ID','Sequence','ML_Score','Prediction']]
    df44.Seq_ID = df44.Seq_ID.str.replace('>','')
    if dplay == 1:
        df44 = df44.loc[df44.Prediction=="Disease Causing"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")
#=============== PQ Desnsity Model start from here ==================
elif Job==4:
    print('\n======= Thanks for using PQ Density module of CDPred. Your results will be stored in file :',result_filename,' =====\n')
    print('==== Checking Density: Processing sequences please wait ...')
    df_2,dfseq = readseq(Sequence)
    df1 = lenchk(dfseq)
    X = PQ_abd(df1,Pat_len)
    X.columns = ['ML_score']
    df45 = det_lab(X,Threshold)
    df44 = pd.concat([df_2,dfseq,df45],axis=1)
    df44.columns = ['Seq_ID','Sequence','Abundance','Prediction']
    df44.Seq_ID = df44.Seq_ID.str.replace('>','')
    if dplay == 1:
        df44 = df44.loc[df44.Prediction=="Disease Causing"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")
#=============== Ensemble Model start from here ==================
elif Job==5:
    print('\n======= Thanks for using Ensemble Model module of CDPred. Your results will be stored in file :',result_filename,' =====\n')
    df44 = ensmbl('model/motif_ensemble',Sequence,Threshold,methd)
    df44.columns = ['Seq_ID','Sequence','Hit Found','ML_Score','Prediction']
    if dplay == 1:
        df44 = df44.loc[df44.Prediction=="Disease Causing"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")
print('\n======= Thanks for using CDPred. Your results are stored in file :',result_filename,' =====\n\n')
