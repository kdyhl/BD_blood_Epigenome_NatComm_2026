#!/usr/bin/env python
# coding: utf-8

# In[1]:

#V1.1
#b.testEpiDEGvsCMapSignatures_subGrpDP_V1.1.py

from cmapPy.pandasGEXpress.parse import parse
import pandas as pd
import csv
import os
import gseapy as gp
import sys
# Specify the path to your GCTX file

args = sys.argv
if len(args) > 1:
    DIR_OU = args[1]
    id_start = args[2]
    id_end = args[3]
# In[11]:

FN_IN_gctx = '/home/unix/leihou/lhou.compbio/data/expression/CMap/level5_beta_trt_cp_blood.gctx'
# FN_IN_cellMeta='/home/unix/leihou/lhou.compbio/data/expression/CMap/cellinfo_beta.txt'
# FN_IN_expMeta='/home/unix/leihou/lhou.compbio/data/expression/CMap/siginfo_beta.txt'

# FNs_IN_degs={'inflamSubGrps.dnDEG':"a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds_subGrpDP/inflamSubgrps.DEG.metaDP.DEGs.dn.gid_byDAVID.txt", 
#             'inflamSubGrps.upDEG':"a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds_subGrpDP/inflamSubgrps.DEG.metaDP.DEGs.up.gid_byDAVID.txt",
#             'nonInflamSubGrps.dnDEG':"a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds_subGrpDP/nonInflamSubgrps.DEG.metaDP.DEGs.dn.gid_byDAVID.txt", 
#             'nonInflamSubGrps.upDEG':"a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds_subGrpDP/nonInflamSubgrps.DEG.metaDP.DEGs.up.gid_byDAVID.txt"}

#DIR_OU = "b_comp2CompoundSig_cMap_inflamSubGrps/a3_V1.2.1_orderedDEGs/"
#os.makedirs(DIR_OU, exist_ok=True)
FN_IN_gmt=DIR_OU + '/EpiDEGs.gmt'
FN_OU_tsv = DIR_OU + '/GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood_'+id_start + '-' + id_end +'.1000perm.tsv'



# Use cmapPy to read the GCTX file
gctoo = parse(FN_IN_gctx)
# gctoo is a GCToo object containing the data from the GCTX file
# Now you can access the data frame of gene expression data
data_slct_df = gctoo.data_df
# Print the first few rows of the data frame


#Dictionary to store results
# Initialize a dictionary to hold your resultsl
results_dict = {'experiment': [], 
                'inflamSubGrps.upDEG.NES': [], 
                'inflamSubGrps.upDEG.NOM_pval': [], 
                'inflamSubGrps.dnDEG.NES': [], 
                'inflamSubGrps.dnDEG.NOM_pval': [],
                'nonInflamSubGrps.upDEG.NES': [], 
                'nonInflamSubGrps.upDEG.NOM_pval': [], 
                'nonInflamSubGrps.dnDEG.NES': [], 
                'nonInflamSubGrps.dnDEG.NOM_pval': []}


# Iterate through each column (experiment) in the DataFrame
id_start = int(id_start)
id_end = int(id_end)
# if(id_end>len(data_slct_df.columns)):
#     id_end=len(data_slct_df.columns)
for experiment in data_slct_df.columns[id_start:id_end]:
    print("id: " + experiment)
    # Sort the DataFrame based on the current column to get the ranked gene list
    ranked_genes = data_slct_df[experiment].sort_values(ascending=False)
    
    if(sum(ranked_genes!=0)/len(ranked_genes) < 0.05):
        continue
    #ranked_genes.index=[f'"{gene}"' for gene in genes]
    print(ranked_genes)
    # Write the ranked gene list to a temporary .rnk file
    FN_TMP_file = DIR_OU+"/"+ str(id_start)+ ".temp.rnk"
    ranked_genes.to_csv(FN_TMP_file, sep='\t', header=False)
    
    # Run gseapy.prerank with the .rnk file and the GMT file
    prerank_results = gp.prerank(rnk=FN_TMP_file, 
                                 gene_sets=FN_IN_gmt, 
                                 outdir=DIR_OU+f"/prerank_results/{experiment}", 
                                 permutation_num=1000, # Set to your desired permutation number
                                 #min_size=15, 
                                 max_size=2000, # Adjust based on your gene set sizes
                                 seed=666,
                                 processes=5) # Adjust based on your available CPU cores
    
    # Extract NES and nominal p-value for both 'up' and 'dn' gene sets from the results
    res2d = prerank_results.res2d  # The results DataFrame
    
    for i in range(0, len(res2d)):
        gset=res2d.loc[i, "Term"]
        nes = res2d.loc[i, 'NES']
        nom_pval = res2d.loc[i, 'NOM p-val']
        results_dict[f'{gset}.NES'].append(nes)
        results_dict[f'{gset}.NOM_pval'].append(nom_pval)
    
#     for gene_set in ['upDEG', 'dnDEG']:
#         # Extracting metrics and appending to results_dict
#         nes = res2d.loc[gene_set, 'NES'] if gene_set in res2d.index else None
#         nom_pval = res2d.loc[gene_set, 'pval'] if gene_set in res2d.index else None
        

    results_dict['experiment'].append(experiment)

# for experiment in data_slct_df.columns[0:10000]:
#     print(experiment)
#     # Sort the DataFrame based on the current column to get the ranked gene list
#     ranked_genes = data_slct_df[experiment].sort_values(ascending=False)
#     #ranked_genes.index=[f'"{gene}"' for gene in genes]
#     print(ranked_genes)
    
#     if(sum(ranked_genes!=0)/len(ranked_genes) < 0.05):
#         continue

#     # Write the ranked gene list to a temporary .rnk file
#     FN_TMP_file = DIR_OU+"/temp.rnk"
#     ranked_genes.to_csv(FN_TMP_file, sep='\t', header=False)
    
#     # Run gseapy.prerank with the .rnk file and the GMT file
#     prerank_results = gp.prerank(rnk=FN_TMP_file, 
#                                  gene_sets=FN_IN_gmt, 
#                                  outdir=DIR_OU+f"/prerank_results/{experiment}", 
#                                  permutation_num=100, # Set to your desired permutation number
#                                  #min_size=15, max_size=500, # Adjust based on your gene set sizes
#                                  seed=666,
#                                  processes=1) # Adjust based on your available CPU cores
    
#     # Extract NES and nominal p-value for both 'up' and 'dn' gene sets from the results
#     res2d = prerank_results.res2d  # The results DataFrame
#     for gene_set in ['upDEG', 'dnDEG']:
#         # Extracting metrics and appending to results_dict
#         nes = res2d.loc[gene_set, 'NES'] if gene_set in res2d.index else None
#         nom_pval = res2d.loc[gene_set, 'pval'] if gene_set in res2d.index else None
#         results_dict[f'{gene_set}.NES'].append(nes)
#         results_dict[f'{gene_set}.NOM_pval'].append(nom_pval)

#     results_dict['experiment'].append(experiment)

results_df = pd.DataFrame(results_dict)
#results_df.set_index('experiment', inplace=True)
results_df.to_csv(FN_OU_tsv, sep="\t")


# In[ ]:




