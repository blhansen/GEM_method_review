
# coding: utf-8

# In[1]:


#!/usr/bin/env python
#usage: python corda_gdc_data_script -f FPKM_file_name -h High_Confidence_percentile -m Medium_Confidence_percentile -l Low_Confidence_percentile -id Run_id
#ex: python corda_gdc_data_script -f 0349f526-7816-4a7d-9967-1f75dd9ff00a.xml -h 90 -m 80 -l 70 -r 4

import sys
import cobra.io
import corda
import pandas as pd
import re
import numpy as np
import gzip
import pickle

#Used to load pre-existing dictionaries which should be the same between different samples
def load_obj(name ):
    with open('../obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

#data from tcga follows format: SBML_id\tFPKM_value
#converts to: HGNC_id\tFPKM_value
#necessary because the genes encoded in recon2_2 follow HGNC format
def translate_gene_exp_data(InFileName):
    #import dictionary of associations between SBML gene IDs and HGNC gene IDs
    sbml_ids=load_obj('sbml_ids')
    
    #create FPKM dict with structure: FPKM[HGNC]=FPKM_value
    Data=[]
    FPKM={}
    exp=0
    
    InFile=open(InFileName,'r')
    
    for Line in InFile:
        Line=str(Line)
        Line=re.sub("'",'',Line)
        Line=re.sub("b",'',Line)
        Line=Line.strip("\\n")
        Data=Line.split("\t")
        Data[0]=re.sub('\.\d+','',Data[0])
        exp=float(Data[1])
        Level=0
        try:
            FPKM[sbml_ids[Data[0]]]=exp
        except:
            pass
        
    InFile.close()
    
    return(FPKM)

#input desired percentiles and returns list of associated values for a sample
def find_percentile_values(FPKM,h,m,l):
    #Filter through FPKM values in file
    #Based on the distribution of values NOT EQUAL TO ZERO
    #Defines confidence groups based on arguments passed in beginning of function
    
    FPKM_values=list(filter((0.0.__ne__),list(FPKM.values())))
    
    high_threshold=np.percentile(FPKM_values,h)
    medium_threshold=np.percentile(FPKM_values,m)
    low_threshold=np.percentile(FPKM_values,l)
    
    #print(low_threshold,medium_threshold,high_threshold)
    return([high_threshold,medium_threshold,low_threshold])

#generate dictionary with dict[gene]=conf_level
def assign_confidence_scores(InFileName,h,m,l):
    
    FPKM=translate_gene_exp_data(InFileName)
    thresholds=find_percentile_values(FPKM,h,m,l)
    high=thresholds[0]
    medium=thresholds[1]
    low=thresholds[2]
    
    Levels={}
    #Loops through HGNC keys in FPKM dict and creates new dict with
    #Levels[HGNC]=confidence level
    #High=3,Medium=2,Low=1,Unsure=0,Not_Detected=-1
    for i in FPKM:
        Level=0
        value=FPKM[i]
        if value>high:
            Level=3
        elif value>medium:
            Level=2
        elif value>low:
            Level=1
        #assume that values equal to exactly 0.0 were not accounted for
        elif value>0:
            Level=-1
        #print(i,'\t',value,'\t',Level)
        Levels[i]=Level
    #Levels
    return(Levels)

#convert dict[gene]=conf_level into dict[reaction]=conf
#utilizes corda.reaction_confidence from CORDA module
def make_confidence_dict(GEM,Levels):
    #generate reaction_conf dict which is creaction_conf[reaction_id]=reaction_confidence
    #confidence levels are same as gene ones, but boolean gene reaction rules are used to filter reactions that are
    #under the control of several genes
    reaction_conf={}
    for r in GEM.reactions:
        reaction_conf[r.id]=corda.reaction_confidence(r.gene_reaction_rule,Levels)
    return(reaction_conf)

def make_OutFileName(InFileName):
    #InFileName=re.sub('\/data\/NCBI\/gene\_exp\_data\/','',InFileName)
    InFileName=re.sub('\/samples\/','',InFileName)
    InFileName=InFileName.strip(".FPKM.txt")
    InFileName=InFileName+".xml"
    return(InFileName)


# In[ ]:


#Initialize h,m,l variables for use later in defining gene confidence groups
h=0
m=0
l=0

#loop through list of arguments and assign variables accordingly
for i in range(len(sys.argv)):
    if sys.argv[i]=='-f':
        InFileName=sys.argv[i+1]
    elif sys.argv[i]=='-h':
        h=sys.argv[i+1]
    elif sys.argv[i]=='-m':
        m=sys.argv[i+1]
    elif sys.argv[i]=='-l':
        l=sys.argv[i+1]
    elif sys.argv[i]=='-id':
        run_id=sys.argv[i+1]


# In[3]:


#Assign Confidence Scores with levels defined by given parameters
Levels=assign_confidence_scores(InFileName,h,m,l)


# In[4]:


#import RECON2 Model
recon2 = cobra.io.read_sbml_model("../GEMs/recon2_2.xml")


# In[5]:


#Make confidence dict
reaction_conf = make_confidence_dict(recon2,Levels)


# In[6]:


#load list of required reactions
#This is the suspect region of code, Moon. Look here.
met_prod=load_obj('met_prod')

#print(met_prod)


# In[7]:


#initialize CORDA object
#model=recon2 'baseline' for your model, in this case recon2.2
#conf= dict of genes and confidence levels they are present in tissue of interest
# 3 = High Conf INCLUDE WHEN POSSIBLE through -1 = No Confidence EXCLUDE WHEN POSSIBLE
# 0 = unknown conf
model = corda.CORDA(recon2,reaction_conf,met_prod)
#use CORDA algorithm to construct tissue-specific model
#computationally intensive task, takes ~10 minutes
print("Building CORDA model...")
model.build()
print(model)


# In[8]:


#Manipulate InFileName for useful OutFileName for use at end of program
exp_file=make_OutFileName(InFileName)
#write COBRA model to computer
out_model=model.cobra_model()
OutFileName="../GEMs/"+exp_file
cobra.io.write_sbml_model(out_model,OutFileName)
print("GEM written to " +OutFileName)


# In[ ]:


#write csv that records which reactions are included in each model
included=pd.DataFrame(data=list(dict(model.included).items()),columns=['Reaction','Included'])
included.to_csv(path_or_buf="../data/results/set_"+run_id+"/"+ exp_file.strip(".xml"))

