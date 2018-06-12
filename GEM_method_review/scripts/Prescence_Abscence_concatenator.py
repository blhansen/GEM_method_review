
# coding: utf-8

# In[1]:


import cobra
import pandas as pd
import os


# In[2]:


#import recon2 GEM
recon2=cobra.io.read_sbml_model("../GEM_reconstructions/recon2_2.xml")


# In[3]:


#create list of metabolic reactions in recon2
recon2_reactions=[]
for r in recon2.reactions:
    recon2_reactions.append(r.id)
#recon2_reactions


# In[4]:


#define directory containing GEMs to loop through
data_dir='../GEM_reconstructions/GEM_reconstructions/single_cell/'

gem_files=[]
gem_filenames=[]
for subdir, dirs, files in os.walk(data_dir):
    for file in files:
        if ".xml" in file:
            #initialize dict to store results
            Included={}
            
            #remove .xml tag for readablity
            gem_filenames.append(file.strip(".xml"))
            
            #define infilename with absolute path
            InFileName=data_dir+file
            
            #load GEM to examine
            model=cobra.io.read_sbml_model(InFileName)
            
            #loop through reacions in recon2 and add to included dict
            # included[reaction]=1 if the reaction is in the gem, 0 otherwise
            for r in recon2_reactions:
                try: 
                    model.reactions.get_by_id(r)
                    Included[r]=1
                except:
                    Included[r]=0
            gem_files.append(Included)


# In[ ]:


#convert results to a data frame
GEMs_summary=pd.DataFrame(data=gem_files,columns=recon2_reactions,index=gem_filenames)


# In[ ]:


#write csv to file
GEMs_summary.to_csv(path_or_buf="~/results.csv")

