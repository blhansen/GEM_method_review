
import cobra
import re
import pandas
import os
import sys

model=cobra.io.read_sbml_model("../GEMs/recon2_2.xml")
model.objective=model.reactions.DM_atp_c_.id
model.objective.direction="max"

m_fva = pandas.DataFrame(cobra.flux_analysis.flux_variability_analysis(model))
m_differences=pandas.DataFrame(m_fva['maximum'])
m_differences.columns=['recon2']

path='../GEM_examples/'
data_dir=path

set_10=pandas.DataFrame()
set_10=m_differences

for subdir, dirs, files in os.walk(data_dir):
    for file in files:
        
        filename=data_dir+file
        
        test=cobra.io.read_sbml_model(filename)

        filename=re.sub(path,'',filename)
        filename=re.sub('.xml','',filename)

        test.add_reaction(model.reactions.DM_atp_c_)
        test.objective=test.reactions.DM_atp_c_.id
        test.objective.direction="max"
        
        fva=pandas.DataFrame()
        fva = pandas.DataFrame(cobra.flux_analysis.flux_variability_analysis(test))
        
        differences=pandas.DataFrame()
        differences=pandas.DataFrame(fva['maximum'])
        differences.columns=[filename]
        
        set_10=set_10.join(differences,how='left')

OutFileName="../flux_maxes_single_cell.csv"

set_10.to_csv(OutFileName)

