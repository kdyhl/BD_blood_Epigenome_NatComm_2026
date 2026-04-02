#!/usr/bin/env python
# coding: utf-8

# In[1]:
#source activate cMap

from cmapPy.pandasGEXpress.parse import parse
import pandas as pd
import csv
import os
import gseapy as gp
import sys
# Specify the path to your GCTX file



FN_IN_gctx = '/home/unix/leihou/lhou.compbio/data/expression/CMap/level5_beta_trt_cp_blood.gctx'
#FN_IN_cellMeta='/home/unix/leihou/lhou.compbio/data/expression/CMap/cellinfo_beta.txt'
#FN_IN_expMeta='/home/unix/leihou/lhou.compbio/data/expression/CMap/siginfo_beta.txt'
FN_IN_slctExp='./b_comp2CompoundSig_cMap_V1.1/GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.slct.sigID.txt'
# FNs_IN_degs={'dnDEG':"a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds/DEG.metaDP.DEGs.dn.gid_byDAVID.txt", 
#             'upDEG':"a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds/DEG.metaDP.DEGs.up.gid_byDAVID.txt"}

DIR_OU = "b_comp2CompoundSig_cMap_V1.1/b3_slctExp/"
os.makedirs(DIR_OU, exist_ok=True)
#FN_OU_gmt=DIR_OU + '/EpiDEGs.gmt'
FN_OU_tsv = DIR_OU + 'exp.slct.p0.001.tsv'



# In[3]:


#cell info
exp_id_df = pd.read_csv(FN_IN_slctExp, sep='\t')


# Use cmapPy to read the GCTX file
gctoo = parse(FN_IN_gctx)
# gctoo is a GCToo object containing the data from the GCTX file
# Now you can access the data frame of gene expression data
data_df = gctoo.data_df
# Print the first few rows of the data frame
#print(data_df.head())

#column_names = gctoo.col_metadata_df.index.tolist()

#data_df.head(2).to_csv("test.tsv", sep='\t', index=True)

# print(gctoo.col_metadata_df)

# In[18]:
exp_id_set=set(exp_id_df.iloc[:, 0])

# expIDs_overlap = list(set(exp_id_df.iloc[:, 0]).intersection(data_df.columns))
# len(expIDs_overlap)
# data_slct_df=data_df[expIDs_overlap]

data_slct_df = data_df.loc[:, data_df.columns.isin(exp_id_set)]

#data_slct_df.index = data_slct_df.index.map(lambda x: f'"{x}"')
#data_slct_df.to_csv(FN_OU_tsv, sep='\t', index=False)


results_df = pd.DataFrame(data_slct_df)
#results_df.set_index('experiment', inplace=True)
results_df.to_csv(FN_OU_tsv, sep="\t")


# In[ ]:




