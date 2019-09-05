# coding: utf-8


#Import packages
from __future__ import print_function
import sys
import os
import pandas as pd
import numpy as np
import cmapPy
import h5py
from cmapPy.pandasGEXpress.parse import parse
from pandas.util.testing import assert_frame_equal
import time
from numpy import unravel_index,fill_diagonal,nanargmax,nanargmin,nanmax,nanmin
import random
import scipy
from scipy import stats
from scipy.stats import spearmanr, pearsonr
import math
from math import sqrt
from joblib import Parallel, delayed
import multiprocessing
from datetime import datetime
from optparse import OptionParser
startTime = datetime.now()


#Specify condition arguments for command line input
parser = OptionParser()
parser.add_option('--cell_line',default='MCF7',type=str,dest='Cell_line',help='Cell line for extraction (default MCF7)')
parser.add_option('--pert_time',default='24',type=str,dest='Perturbation_Time',help='Perturbation time for extraction (hours) (default 24h)')
parser.add_option('--pert_dose',default='10',type=str,dest='Perturbation_Dose',help='Perturbation dose for extraction (uM) (default 10)')
parser.add_option('--ncores',default='1',type=int,dest='cores',help='Number of cores (default 1)')
(options, args) = parser.parse_args()

Cell_line = options.Cell_line
Perturbation_Time = options.Perturbation_Time
Perturbation_Dose = options.Perturbation_Dose
cores = options.cores

Perturbation_Time=Perturbation_Time+" "+"h"
Perturbation_Dose_2 = Perturbation_Dose+".0"+" "+"um"
Perturbation_Dose=Perturbation_Dose+" "+"µM"



#Specify directories for files
data_dir = 'data/'

ds_path_Phase1=data_dir+"GSE92742_Broad_LINCS_Level5_*.gctx"
ds_path_Phase2=data_dir+"GSE70138_Broad_LINCS_Level5_*.gctx"
gene_info=pd.read_table(data_dir+"GSE92742_Broad_LINCS_gene_info.txt", sep='\t')
sig_info_1=pd.read_table(data_dir+"GSE92742_Broad_LINCS_sig_info.txt", sep="\t")
sig_info_2=pd.read_table(data_dir+"GSE70138_Broad_LINCS_sig_info.txt", sep="\t")

#Find compounds in phase 1
compounds_phase_1 = sig_info_1.loc[(sig_info_1['cell_id']==Cell_line) & (sig_info_1['pert_idose']==Perturbation_Dose) & (sig_info_1['pert_itime']==Perturbation_Time)]
compounds_phase_1 = compounds_phase_1['pert_iname'].tolist()

#find compounds in phase 2
compounds_phase_2 = sig_info_2.loc[(sig_info_2['cell_id']==Cell_line) & (sig_info_2['pert_idose']==Perturbation_Dose_2) & (sig_info_2['pert_itime']==Perturbation_Time)]
compounds_phase_2 = compounds_phase_2['pert_iname'].tolist()

compounds = list(set(compounds_phase_1) | set(compounds_phase_2))
compounds = [str(i) for i in compounds]
#compounds = compounds[0:100] subset compounds to test the code

#Initialise empty dataframe
empty_df=pd.DataFrame()

#Empty df of all the gene_IDs and compounds
row_names_lm = list((gene_info['pr_gene_id'][gene_info["pr_is_lm"] == 1]).sort_values())
column_names_lm = list(compounds)
matrix_test_lm = np.empty((len(row_names_lm),len(compounds)))
matrix_test_lm[:] = np.NaN
row_names_lm=[str(i) for i in row_names_lm]
df_lm = pd.DataFrame(matrix_test_lm, columns=column_names_lm, index=row_names_lm)
compound_geneExpression_lm=df_lm.sort_index()

#Create file to contain results
headers=['Compound_id','No. Replicates','CC','SS','TAS']
GE_results = open('TAS_Scores/'+'LINCS_'+Cell_line+'_'+Perturbation_Time+'_'+Perturbation_Dose+'_'+'TAS'+'.txt','w+')
GE_results.write('\t'.join(map(str,headers)) + '\n')

headers = row_names_lm
headers = ['Compound_id'] + headers

no_rep_mat = open('Compounds_With_No_Replicates/'+'LINCS_'+'all_'+Cell_line+'_'+Perturbation_Time+'_'+Perturbation_Dose+'_'+'No_Replicates'+'.txt','w+')
no_rep_mat.write('\t'.join(map(str,headers)) + '\n')
#no_rep_mat.close()
#no_rep_mat = open('Compounds_With_No_Replicates/'+'LINCS_'+'all_'+Cell_line+'_'+Perturbation_Time+'_'+Perturbation_Dose+'_'+'No_Replicates'+'.txt','a')

consensus_mat = open('Consensus_Signatures/'+'LINCS_'+'all_'+Cell_line+'_'+Perturbation_Time+'_'+Perturbation_Dose+'_'+'Consensus'+'.txt','w+')
consensus_mat.write('\t'.join(map(str,headers)) + '\n')

#Define function for extracting data and calculating TAS

def CalculateTAS(compound):
    print(str(compounds.index(compound)), 'out of', str((len(compounds))))
    sig_ids=[]
    multi_mat_1=pd.DataFrame()
    multi_mat_2=pd.DataFrame()
    
    sigs_phase_1=sig_info_1.loc[(sig_info_1['pert_iname'] == compound)&(sig_info_1['cell_id']==Cell_line)&(sig_info_1['pert_itime']==Perturbation_Time)&(sig_info_1['pert_idose']==Perturbation_Dose)]
    sig_ids=sigs_phase_1['sig_id']
    try:
        multi_mat_1= parse(ds_path_Phase1, rid = row_names_lm, cid=sig_ids)
        multi_mat_1=multi_mat_1.data_df.sort_index()
    except UnboundLocalError:
        multi_mat_1=pd.DataFrame()
        

    sigs_phase_2=sig_info_2.loc[(sig_info_2['pert_iname'] == compound)&(sig_info_2['cell_id']==Cell_line)&(sig_info_2['pert_itime']==Perturbation_Time)&(sig_info_2['pert_idose']==Perturbation_Dose_2)]
    sig_ids=sigs_phase_2['sig_id']
    try:
        multi_mat_2= parse(ds_path_Phase2,rid = row_names_lm, cid=sig_ids)
        multi_mat_2=multi_mat_2.data_df.sort_index()
    except UnboundLocalError:
        multi_mat_2=pd.DataFrame()
        
        
        
    if not multi_mat_2.equals(empty_df) and not multi_mat_1.equals(empty_df):
        multi_mat=pd.concat([multi_mat_1, multi_mat_2], axis=1)
    elif not multi_mat_1.equals(empty_df) and multi_mat_2.equals(empty_df): 
        multi_mat=multi_mat_1
    elif not multi_mat_2.equals(empty_df) and multi_mat_1.equals(empty_df):
        multi_mat=multi_mat_2
    elif multi_mat_2.equals(empty_df) and multi_mat_1.equals(empty_df):
        print("NA")
        multi_mat=pd.DataFrame()
        allstats = []
        single_sig = []
        consensus_sig = []
        return([allstats,consensus_sig,single_sig])
      
    if multi_mat.shape[1]==1:
        #Save data
        single_sig = multi_mat.median(axis=1)
        single_sig = single_sig.tolist()
        #print(len(consensus_sig))
        single_sig.insert(0,compound)
        single_sig = list(map(str, single_sig))
        allstats = []
        consensus_sig = []
        
        return([allstats, consensus_sig,single_sig])

        #NA result, cannot calculate statistics
        #allstats=[compound,multi_mat.shape[1],'NA','NA','NA']


        

    elif multi_mat.shape[1]>=2:
        #Calculate replicate correlation (CC)
        corr_mat = multi_mat.corr(method='spearman')
        upper_tri = corr_mat.values[np.triu_indices_from(corr_mat.values,1)] #get upper triangle
        CC = np.median(upper_tri)

        #Calculate signature strength (SS)
        greater_than_2 = multi_mat[multi_mat >= 2.0].stack().count()
        less_than_2 = multi_mat[multi_mat <= -2.0].stack().count()
        ss_notnorm = less_than_2 + greater_than_2
        ss_norm = ss_notnorm/len(multi_mat.columns) #divide by no. replicates


        #Calculate TAS
        TAS = sqrt((abs(CC) * ss_norm)/978)
        

        allstats=[compound,multi_mat.shape[1],CC,ss_norm,TAS]

        if multi_mat.shape[1]==2:
        
            #Save replicates
            fname = 'All_Replicates/'+str(compound)+"_"+str(Cell_line)+"_"+str(Perturbation_Time)+"_"+str(Perturbation_Dose)
            try:
                multi_mat.to_csv(fname+".csv")
            except:
                print("Error with saving replicates of compound "+compound)
        
            #Consensus signature is just median
            consensus_sig = multi_mat.median(axis=1)
            consensus_sig = consensus_sig.tolist()
            #print(len(consensus_sig))
            consensus_sig.insert(0,compound)
            consensus_sig = list(map(str, consensus_sig))
            single_sig = []

                        
        
        #GE_results.write('\t'.join(map(str,allstats)) + '\n')

        elif multi_mat.shape[1]>2:
                            
            #Save replicates
            fname = 'All_Replicates/'+str(compound)+"_"+str(Cell_line)+"_"+str(Perturbation_Time)+"_"+str(Perturbation_Dose)
            try:
                multi_mat.to_csv(fname+".csv")
            except:
                print("Error with saving replicates of compound "+compound)

            #Get consensus signature
            corr_reps = pd.DataFrame(multi_mat.corr(method='spearman'))
            a = np.array(corr_reps)
            np.fill_diagonal(a,'NaN')
            corr_reps_nan = pd.DataFrame(a)
            col_weights = []
            for col in corr_reps_nan:
                col_weight = corr_reps_nan[col].sum()
                col_weights.append(col_weight)
            col_weights_norm = [float(i)/sum(col_weights) for i in col_weights]
            linear_list = []
            for idx, weight in enumerate(col_weights_norm):
                linear = multi_mat.ix[:,idx]*col_weights_norm[idx]
                linear_list.append(linear)
            combined_sig = sum(linear_list)
            combined_sig = list(map(str, combined_sig))
            consensus_sig = combined_sig
            consensus_sig.insert(0,compound)
            
            single_sig =  []
            
        return([allstats,consensus_sig,single_sig])
        #Calculate modified signature strength (SS) - TO DO
        
        
         

results = Parallel(n_jobs=cores,backend="multiprocessing")(delayed(CalculateTAS)(compound) for compound in compounds)

for compound in results:
    
    stat = compound[0]
    consensus = compound[1]
    singlesig = compound[2]
    
    if len(stat)!=0:
        GE_results.write('\t'.join(map(str,stat)) + '\n')

    if len(consensus)!=0:
        consensus_mat.write('\t'.join(map(str,consensus)) + '\n')

    if len(singlesig)!=0
        no_rep_mat.write('\t'.join(map(str,singlesig)) + '\n')
        

                            
GE_results.close()
no_rep_mat.close()
consensus_mat.close()

print(datetime.now() - startTime)
print("Finished.")       
        
                
