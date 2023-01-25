# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:42:46 2019

@author: rmoge
"""

#%%
from __future__ import division
import os, math, pickle
#import mult_tools as mt
import numpy as np
import pandas as pd
import random
from decimal import Decimal
import time
from pathlib import Path
homedir = str(Path.home())
workingdir = homedir + '\\Documents\\GitHub\\XevoLong'
os.chdir(workingdir)

poptypes = ['rpl19a','pmr1','sch9','ypt','tor1','tif2','JD174','BY4742']
strains = ['BY4742']
output_to_keep=["SNP"]
my_anc_dirs=[81,82,83,84,85,86,87,88]


#%%

def get_sites_to_remove(ancestor_dirs):
    sites_to_remove = []
    for number in ancestor_dirs:
        ancestor_path = workingdir + '\\Pipeline_XL\\1_Carbonate\\breseq\\breseq_' + str(number) + '\\output\\output.gd'
        #C:\Users\rmoge\Documents\GitHub\Bifidobacterium\Pipeline_XL\1_Carbonate\breseq\breseq_25\output
        for i, line in enumerate(open(ancestor_path, 'r')):
            line_split = line.strip().split('\t')
            if line_split[0] in output_to_keep:
                sites_to_remove.append( line_split[3] + '_' + str(line_split[4]))
            else:
                pass
# =============================================================================
#     Below, you will add sites that did no show up in the ancestor, but are obvious problematic sites, either due to cross contamination , or due to a mutation in a founding population for one of the treatments specifically
# =============================================================================
#    sites_to_remove+=['CP001853_146344','CP001853_1305547','CP001853_1437318']
# =============================================================================
#     
# =============================================================================
    print(sites_to_remove)
    return sites_to_remove
#%%
                    
def get_gxp_matrix_N(output_name='',include_polymorphisms=True,output_to_keep=["SNP"]):

    gene_by_pop_dict = {}#initialize GxP matrix
    for strain in strains:
        sites_to_remove = get_sites_to_remove(my_anc_dirs)#note that if you were working with DIVERGENT ANCESTORS, the correct call for this argument would by strain, not my_anc_dirs. I.e., you would be calling get_sites_to_remove() once separately for each strain in strains .
        if strain == 'BY4742':
            dirs = []
            dirs = ['breseq_' + str(x) for x in range(1,81)]
            #dirs = ['breseq_11','breseq_17']
            ref_path = workingdir + '\\Pipeline_XL\\1_Carbonate\\S288C_NCBI.gb'
        

        #effective_gene_lengths, effective_gene_lengths_syn, Lsyn, Lnon, substitution_specific_synonymous_fraction = mt.calculate_synonymous_nonsynonymous_target_sizes(ref_path)
        for dir in dirs:
            print(dir)
            pop = dir.split('_')[1]
            gene_count_dict_pop = 0#initialize number of S mutations for THIS pop
            count = 0#initialize number of S mutations for THIS pop
            for i, line in enumerate(open(workingdir +'\\Pipeline_XL\\1_Carbonate\\breseq\\' + dir + '\\output\\evidence\\annotated.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] not in output_to_keep:
                    continue
                if line_split[3] + '_' + line_split[4] in sites_to_remove:
                    continue
                try: 
                    frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                except ValueError:
                    continue
                if frequency != 1 and include_polymorphisms==False:
                    continue
                if line_split[0] == 'SNP':
                    if [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsynonymous':# NOTE that nonsense mutations are thrown out when calculating dN/dS!!!!!! or [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsense':
                        r=random.random()
                        if frequency > r:
                            count += 1
            gene_by_pop_dict[int(pop)]=count  #does not have to be int(), can just be pop if youre not using numerical population assignments              
                        
                      
    
    return gene_by_pop_dict

#%%
def get_gxp_matrix_S(output_name='',include_polymorphisms=True,output_to_keep=["SNP"]):

    gene_by_pop_dict = {}#initialize GxP matrix
    for strain in strains:
        sites_to_remove = get_sites_to_remove(my_anc_dirs)#note that if you were working with DIVERGENT ANCESTORS, the correct call for this argument would by strain, not my_anc_dirs. I.e., you would be calling get_sites_to_remove() once separately for each strain in strains .
        if strain == 'BY4742':
            dirs = []
            dirs = ['breseq_' + str(x) for x in range(1,81)]
            #dirs = ['breseq_11','breseq_17']
            ref_path = workingdir + '\\Pipeline_XL\\1_Carbonate\\S288C_NCBI.gb'
        

        #effective_gene_lengths, effective_gene_lengths_syn, Lsyn, Lnon, substitution_specific_synonymous_fraction = mt.calculate_synonymous_nonsynonymous_target_sizes(ref_path)
        for dir in dirs:
            print(dir)
            pop = dir.split('_')[1]
            gene_count_dict_pop = 0#initialize number of S mutations for THIS pop
            count = 0#initialize number of S mutations for THIS pop
            for i, line in enumerate(open(workingdir +'\\Pipeline_XL\\1_Carbonate\\breseq\\' + dir + '\\output\\evidence\\annotated.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] not in output_to_keep:
                    continue
                if line_split[3] + '_' + line_split[4] in sites_to_remove:
                    continue
                try: 
                    frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                except ValueError:
                    continue
                if frequency != 1 and include_polymorphisms==False:
                    continue
                if line_split[0] == 'SNP':
                    if [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'synonymous':
                        r=random.random()
                        if frequency > r:
                            count += 1
            gene_by_pop_dict[int(pop)]=count  #does not have to be int(), can just be pop if youre not using numerical population assignments              
                        
                      
    
    return gene_by_pop_dict
#%%
start=time.time()
testN=get_gxp_matrix_N('XL_gxp',True)
testS=get_gxp_matrix_S('XL_gxp',True)
end=time.time()
print(end-start)
#%%
'''
def get_Dx(xdict):
    for evopop,mutations in xdict.items():
        #for gene,totalfreqs
        pass
    return
'''
#%%
'''
myc=0
for k,v in test1['1'].items():
    r=random.random()
    if v > r:
        myc+=1
        '''
#%%
def calc_targetsize(p=1/6,q=1/6,r=1/6,s=1/6,t=1/6,possATtoCG=0,possATtoGC=0,posstATtoTA=0,possCGtoGC=0,possCGtoTA=0,possGCtoTA=0):
    u=1-p-q-r-s-t
    NorS=6*(p*possATtoCG + q*possATtoGC + r*possCGtoGC + s*possCGtoTA + t*possGCtoTA + u*possGCtoTA)
    return NorS

def calc_dNdS_with_pseudocount(N=1,S=1,DN=0,DS=0):
    dN = -3/4*math.log(1-(4/3*DN/(N)))
    dS = -3/4*math.log(1-(4/3*(DS+1)/(S)))
    #print("dN = " + str(dN)); print("dS = " + str(dS))
    #print("dN/dS = " + str(dN/dS))
    return dN/dS

def calc_dNdS_no_zeroes(N=1,S=1,DN=0,DS=0):
    dN = -3/4*math.log(1-(4/3*DN/N))
    dS = -3/4*math.log(1-(4/3*DS/S))
    #print("dN = " + str(dN)); print("dS = " + str(dS))
    #print("dN/dS = " + str(dN/dS))
    return dN/dS

def calc_dNdS_fromtuple(intup):
    N=intup[0];S=intup[1];DN=intup[2];DS=intup[3]
    try:
        dN = -3/4*math.log(1-(4/3*DN/N))
    except ZeroDivisionError:
        dN = 0
    dS = -3/4*math.log(1-(4/3*DS/(S)))
    #print("dN = " + str(dN)); print("dS = " + str(dS))
    #print("dN/dS = " + str(dN/dS))
    return dN/dS

def calc_dNdS_fromtuple_with_pseudocount(intup):
    N=intup[0];S=intup[1];DN=intup[2];DS=intup[3]
    dN = -3/4*math.log(1-(4/3*DN/N))
    dS = -3/4*math.log(1-(4/3*(DS+1)/S))
    #print("dN = " + str(dN)); print("dS = " + str(dS))
    #print("dN/dS = " + str(dN/dS))
    return dN/dS




#calc_dNdS_fromtuple(mm9)
#%%
#values for N and S were calculated using gdtools count -b. See script: dNdS_gdtools.count.sh













#Spectrum-aware NYI










































yeastN=18535746
yeastS=5342553
dnds_dict=dict()
zerolist=[]
for k,v in testS.items():
    zerolist+=[v]
if min(zerolist)<=0:
    for i in range(1,81):
        eo = calc_dNdS_with_pseudocount(yeastN,yeastS,testN[i],testS[i])
        dnds_dict[i]=eo
    print("A pseudocount of 1 was added to every DS value.")
else:
    for i in range(1,81):
        eo = calc_dNdS_no_zeroes(yeastN,yeastS,testN[i],testS[i])
        dnds_dict[i]=eo    
    print("No populations had DS = 0, so no pseudocount of 1 was added.")

dnds_df = pd.DataFrame(dnds_dict,index=['dN/dS'])
dnds_df=dnds_df.T
dnds_df.fillna(0, inplace=True)

dnds_df_out = "~\\Documents\\GitHub\\XevoLong\\Pipeline_XL\\3_Python\\dnds_output.csv"
dnds_df.to_csv(dnds_df_out, sep = ',', index = True)