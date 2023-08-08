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
import statsmodels.api
import random
from decimal import Decimal
import time
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

homedir = str(Path.home())
workingdir = homedir + '\\Documents\\GitHub\\XevoLong'
os.chdir(workingdir)

poptypes = ['rpl19a','pmr1','sch9','ypt6','tor1','tif2','JD174','BY4742']
strains = ['BY4742']
output_to_keep=["SNP"]
my_anc_dirs=[81,82,83,84,85,86,87,88]

canddict=dict(); poptrtsL=list()
for g in poptypes:
    canddict[g+"_NL"]=dict()
    poptrtsL+=[g+"_NL"]
for g in poptypes:
    canddict[g+"_L"]=dict()
    poptrtsL+=[g+"_L"]
    
poptrtsD=dict()
CXL=1
for name in poptrtsL:
    poptrtsD[name]=range(CXL,CXL+5)
    CXL+=5
    
    

gxpcounts_fixed=pd.read_csv(workingdir + "\\Pipeline_XL\\3_Python\\XL_gxp_fixed.only.csv")
gxpcounts_p=pd.read_csv(workingdir + "\\Pipeline_XL\\3_Python\\XL_gxp.csv")    
    
    
'''
pseudocode for reading the spreadsheet
for row in spreadsheet with index that matches that pop trt:
    for each column:
        try:
            if the value of the cell > 0:
                dict entry += 1
        except KeyError:
            initialize a dict entry for name of that column with value = 0
            if the value of the cell > 0:
                dict entry += 1
                
                '''

#%%
def debugit(df):
    for ri,r in df.iterrows():
        if ri == 79:
            for ci,c in r.items():
                print(type(ci))
                print(c)

debugit(gxpcounts_fixed)

def iteration_test(df):
    for rowIndex, row in df.iterrows():
        print(rowIndex)
        for columnIndex,value in row.items():
            #print(type(columnIndex))
            #print(value)
            pass
    return

def iteration(rowgroupsD,df):
    multiplesD=dict();multiplesDL=dict()
    for k,v in rowgroupsD.items():
        multiplesD[k]=dict();multiplesDL[k]=list()
        tempDict=dict()
        for rowIndex, row in df.iterrows():
            if rowIndex+1 in v:
                for columnIndex,value in row.items():
                    if columnIndex == "Unnamed: 0":
                        continue
                    else:
                        try:
                            if value > 0:
                                tempDict[columnIndex]+=1
                        except KeyError:
                            tempDict[columnIndex]=0
                            tempDict[columnIndex]+=1            
            else:
                continue
            
        for k2,v2 in tempDict.items():
            if v2 > 1:
                multiplesD[k][k2]=v2;multiplesDL[k]+=[k2]
    return multiplesD, multiplesDL

#md,mdl=iteration(poptrtsD,gxpcounts_fixed)
#mdp,mdlp=iteration(poptrtsD,gxpcounts_p)
#%%
def output_candidates_from_gxp(rowgroupsD,df):
    multiplesD=dict();multiplesDL=dict();sumsD=dict()
    for k,v in rowgroupsD.items():
        multiplesD[k]=dict();multiplesDL[k]=list();sumsD[k]=dict()
        tempDict=dict();tempSums=dict()
        for rowIndex, row in df.iterrows():
            if rowIndex+1 in v:
                for columnIndex,value in row.items():
                    if columnIndex == "Unnamed: 0":
                        continue
                    else:
                        try:
                            if value > 0:
                                tempDict[columnIndex]+=1
                        except KeyError:
                            tempDict[columnIndex]=0
                            tempDict[columnIndex]+=1   
                        try:
                            tempSums[columnIndex]+=value
                        except KeyError:
                            tempSums[columnIndex]=0
                            tempSums[columnIndex]+=value
            else:
                continue
            
        for k2,v2 in tempDict.items():
            if v2 > 1:
                multiplesD[k][k2]=v2;multiplesDL[k]+=[k2]
        for ks2,vs2 in tempSums.items():
            sumsD[k][ks2]=vs2
    return multiplesD, multiplesDL, sumsD

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
    #print(sites_to_remove)
    return sites_to_remove

def get_freqs_list_N(output_name='',include_polymorphisms=True,output_to_keep=["SNP"],rowgroupsD=dict()):    
    for strain in strains:
        sites_to_remove = get_sites_to_remove(my_anc_dirs)#note that if you were working with DIVERGENT ANCESTORS, the correct call for this argument would be strain, not my_anc_dirs. I.e., you would be calling get_sites_to_remove() once separately for each strain in strains .
        if strain == 'BY4742':
            dirs = []
            dirs = ['breseq_' + str(x) for x in range(1,81)]
            #dirs = ['breseq_11','breseq_17']
            ref_path = workingdir + '\\Pipeline_XL\\1_Carbonate\\S288C_NCBI.gb'
        

        #effective_gene_lengths, effective_gene_lengths_syn, Lsyn, Lnon, substitution_specific_synonymous_fraction = mt.calculate_synonymous_nonsynonymous_target_sizes(ref_path)
        freqlistD=dict()
        for k,v in rowgroupsD.items():
            #print(k)
            freqlistD[k]=list()
            for dIndex,d in enumerate(dirs):
                if (dIndex+1) in v:
                    #print(dIndex)
                    #print(d)
                    pop = d.split('_')[1]
                    #gene_count_dict_pop = 0#initialize number of S mutations for THIS pop
                    #count = 0#initialize number of S mutations for THIS pop
                    for i, line in enumerate(open(workingdir +'\\Pipeline_XL\\1_Carbonate\\breseq\\' + d + '\\output\\evidence\\annotated.gd', 'r')):
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
                            if [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsynonymous' or [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsense':
                                freqlistD[k]+=[frequency]
                            else:
                                continue
                    #gene_by_pop_dict[int(pop)]=count  #does not have to be int(), can just be pop if youre not using numerical population assignments              
                else:
                    continue
                      
    
    return freqlistD

def get_freqs_list_S(output_name='',include_polymorphisms=True,output_to_keep=["SNP"],rowgroupsD=dict()):    
    for strain in strains:
        sites_to_remove = get_sites_to_remove(my_anc_dirs)#note that if you were working with DIVERGENT ANCESTORS, the correct call for this argument would be strain, not my_anc_dirs. I.e., you would be calling get_sites_to_remove() once separately for each strain in strains .
        if strain == 'BY4742':
            dirs = []
            dirs = ['breseq_' + str(x) for x in range(1,81)]
            #dirs = ['breseq_11','breseq_17']
            ref_path = workingdir + '\\Pipeline_XL\\1_Carbonate\\S288C_NCBI.gb'
        

        #effective_gene_lengths, effective_gene_lengths_syn, Lsyn, Lnon, substitution_specific_synonymous_fraction = mt.calculate_synonymous_nonsynonymous_target_sizes(ref_path)
        freqlistD=dict()
        for k,v in rowgroupsD.items():
            #print(k)
            freqlistD[k]=list()
            for dIndex,d in enumerate(dirs):
                if (dIndex+1) in v:
                    #print(dIndex)
                    #print(d)
                    pop = d.split('_')[1]
                    #gene_count_dict_pop = 0#initialize number of S mutations for THIS pop
                    #count = 0#initialize number of S mutations for THIS pop
                    for i, line in enumerate(open(workingdir +'\\Pipeline_XL\\1_Carbonate\\breseq\\' + d + '\\output\\evidence\\annotated.gd', 'r')):
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
                                freqlistD[k]+=[frequency]
                            else:
                                continue
                    #gene_by_pop_dict[int(pop)]=count  #does not have to be int(), can just be pop if youre not using numerical population assignments              
                else:
                    continue
                      
    
    return freqlistD
#%%
def index_genbank_features_multiplechroms(gb_records, feature_type, qualifier, qualifier_translation, AT_rel_rate, GC_rel_rate):
    
    maxperchromlist = list()
    for rec in gb_records:
        max_on_this_chrom = max([len(i) for i in rec.features[1:]])
        print(max_on_this_chrom)
        maxperchromlist.append(max_on_this_chrom)
    max_gene_length = max(maxperchromlist)
    rel_rate_MAX=max(AT_rel_rate,GC_rel_rate)
    print(str(max_gene_length) + " is max gene length")
    
    answer = dict()

    for rec in gb_records:
        for (index, feature) in enumerate(rec.features) :
            if feature.type==feature_type :
                if qualifier in feature.qualifiers :
                    #There should only be one locus_tag per feature, but there
                    #are usually several db_xref entries
                    for value in feature.qualifiers[qualifier] :
                        if value in answer :
                            print("WARNING - Duplicate key %s for %s features %i and %i" \
                               % (value, feature_type, answer[value], index))
                        else :
                            answer[value] = {}
                            answer[value]["index"] = index
                try:
                    mygg=feature.qualifiers[qualifier_translation]#attempt to access the gene name for this locus tag
                except KeyError:
                    mygg=feature.qualifiers[qualifier]#if there is no gene name, just call it by the locus tag
                if qualifier in feature.qualifiers:
                    for value in feature.qualifiers[qualifier] :
                        answer[value]['gene'] = mygg[0]
                if qualifier in feature.qualifiers:
                    for value in feature.qualifiers[qualifier] :
                        answer[value]['gene_length']=len(feature)
                if qualifier in feature.qualifiers:
                    for value in feature.qualifiers[qualifier] :        
                        answer[value]["gene_len_rel"]= (len(feature) / max_gene_length)
                if qualifier in feature.qualifiers:
                    for value in feature.qualifiers[qualifier]:
                        GCcount=0
                        for (index2,pos) in enumerate(rec.features[index]):
                            if rec[pos]=="C" or rec[pos]=="G":
                                GCcount+=1
                            else:
                                pass
                            #print(pos)
                        answer[value]["gene_GC_prop"]= (GCcount/len(feature))
                if qualifier in feature.qualifiers:
                    for value in feature.qualifiers[qualifier]:
                        answer[value]['mut_rel_rate']=((answer[value]['gene_GC_prop']*GC_rel_rate) + ((1 - answer[value]['gene_GC_prop'])*AT_rel_rate)) / rel_rate_MAX
                    
            
    return answer
#%%
def simulate_mutations_GC(mutation_freq_tup, indict, reps=10000):
    counter = 0; sim_master=dict(); sim_dict=dict()
    for k in indict.keys():
        sim_master[k] = list()
    while counter < reps:
        sim_temp=dict()
        for index, m in enumerate(mutation_freq_tup):
            while True:
                w = random.choice(list(indict));r = random.random();rGC=random.random()
                #print("w = " + str(w) + ", w len rel = " + str(indict[w]['gene_len_rel']) + " and r = " + str(r))
                if indict[w]['gene_len_rel'] >= r and indict[w]['mut_rel_rate'] >= rGC:
                    #do something###########################################
                    ##################################
                    #################################################
                    try:
                        sim_temp[w] += m
                    except KeyError:
                        sim_temp[w] = 0
                        sim_temp[w] += m
                    #print("value added = " + str(m))
                    break
                else:
                    continue
        sim_dict[counter] = sim_temp
        for k,v in sim_temp.items():
            sim_master[k].append(v)
        for km in sim_master.keys():
            if len(sim_master[km]) < counter+1:
                sim_master[km].append(0)
        counter += 1
        #pp.pprint(sim_temp)
    
    return sim_master

def simulate_mutations_noGC(mutation_freq_tup, indict, reps=10000):
    counter = 0; sim_master=dict(); sim_dict=dict()
    for k in indict.keys():
        sim_master[k] = list()
    while counter < reps:
        sim_temp=dict()
        for index, m in enumerate(mutation_freq_tup):
            while True:
                w = random.choice(list(indict));r = random.random()
                #print("w = " + str(w) + ", w len rel = " + str(indict[w]['gene_len_rel']) + " and r = " + str(r))
                if indict[w]['gene_len_rel'] >= r:
                    #do something###########################################
                    ##################################
                    #################################################
                    try:
                        sim_temp[w] += m
                    except KeyError:
                        sim_temp[w] = 0
                        sim_temp[w] += m
                    #print("value added = " + str(m))
                    break
                else:
                    continue
        sim_dict[counter] = sim_temp
        for k,v in sim_temp.items():
            sim_master[k].append(v)
        for km in sim_master.keys():
            if len(sim_master[km]) < counter+1:
                sim_master[km].append(0)
        counter += 1
        #pp.pprint(sim_temp)
    
    return sim_master

#%%
def generate_dict_of_simulation_for_each_poptrt(poptrtnames,freqlistsDict,index_genbank_features_Dict,reps_requested,filenamer=""):
    outdict=dict();osvar=str(time.time());os.mkdir(path=workingdir+"\\Pipeline_XL\\3_Python\\"+osvar)
    for ptI, poptrt in enumerate(poptrtnames):
        outdict[poptrt] = simulate_mutations_GC((freqlistsDict[poptrt]),DD,reps = reps_requested)
        pdexporter=pd.DataFrame.from_dict(outdict[poptrt])
        pdexporter.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline_XL\\3_Python\\"+osvar+"\\simout_"+filenamer+"_"+str(poptrt)+"_"+str(reps_requested)+".csv")
        
    return outdict


#%%
freqlistspN=get_freqs_list_N(output_name='XL_gxp',include_polymorphisms=True,rowgroupsD=poptrtsD)
freqlistspS=get_freqs_list_S(output_name='XL_gxp',include_polymorphisms=True,rowgroupsD=poptrtsD)
freqlistsfN=get_freqs_list_N(output_name='XL_gxp',include_polymorphisms=False,rowgroupsD=poptrtsD)
freqlistsfS=get_freqs_list_S(output_name='XL_gxp',include_polymorphisms=False,rowgroupsD=poptrtsD)

#%%
recs = [rec for rec in SeqIO.parse(workingdir + "\\Pipeline_XL\\1_Carbonate\\S288C_NCBI.gb", "genbank")]

# print the number of sequence records that were extracted
print(len(recs))

# print annotations for each sequence record
'''for rec in recs:
	print(rec.annotations)
    
    '''
#%%  
myrec=recs[0]#chrom I

# print the CDS sequence feature summary information for each feature in each
# sequence record
'''for rec in recs:
    feats = [feat for feat in rec.features if feat.type == "CDS"]
    for feat in feats:
        print(feat)
        
print(repr(myrec.seq))'''
#%%
atRelRateToInput=.063+.110+.144;gcRelRateToInput=.152+.182+.350
start=time.time()
DD = index_genbank_features_multiplechroms(recs, "gene", "locus_tag", 'gene', AT_rel_rate = atRelRateToInput, GC_rel_rate=gcRelRateToInput)
end=time.time();print(end-start)
#%%
 
start=time.time()
simoutspNDict=generate_dict_of_simulation_for_each_poptrt(poptrtsL,freqlistspN,DD,5357,"pN")
end=time.time();print(end-start)

#%%    
start=time.time();print(start)
simoutsfNDict=generate_dict_of_simulation_for_each_poptrt(poptrtsL,freqlistsfN,DD,5357,"fN")
end=time.time();print(end-start)

#%%
start=time.time();print(start)
simoutspSDict=generate_dict_of_simulation_for_each_poptrt(poptrtsL,freqlistspS,DD,5357,"pS")
end=time.time();print(end-start)    
    
#%%    
start=time.time();print(start)
simoutsfSDict=generate_dict_of_simulation_for_each_poptrt(poptrtsL,freqlistsfS,DD,5357,"fS")
end=time.time();print(end-start)    

#%%
mdf,mdfl,mdfs=output_candidates_from_gxp(poptrtsD,gxpcounts_fixed)
mdp,mdpl,mdps=output_candidates_from_gxp(poptrtsD,gxpcounts_p)
    
    

#%%
#Import simulation results from file
#%%
#You can import simulation results from a file.
start = time.time();print(start)
#recs = [rec for rec in SeqIO.parse(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\syn3A_genome\Synthetic.bacterium_JCVI-Syn3A.gb", "genbank")]
#myrec=recs[0]
#syn1recs = [rec for rec in SeqIO.parse(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\syn1.0_genome\Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb", "genbank")]  
#syn1rec=syn1recs[0]

#GD = gene_to_tag(myrec, "locus_tag", "gene")
#DD = index_genbank_features(myrec, "gene", "locus_tag", 'gene', AT_rel_rate = 35/0.76, GC_rel_rate=486/0.24)
#GD1 = gene_to_tag(syn1rec, "locus_tag", "gene")
#DD1 = index_genbank_features(syn1rec, "gene", "locus_tag", 'gene', AT_rel_rate=219/0.76,GC_rel_rate=1591/0.24)


testingBY4742_pd = pd.read_csv(workingdir+"\\Pipeline_XL\\3_Python\\fS\\simout_fS_BY4742_L_5357.csv")
testingBY4742_dic = pd.DataFrame.to_dict(testingBY4742_pd)
#%%

def get_pvalues_from_simulations(workingdir=workingdir,popTreatments=[],poptrtsDict=poptrtsD,gxpcounts={},dirname='',repsNumber=5357):
    import statsmodels
    candidatesMultiplyMutated,candidatesLists,candidatesCounts=output_candidates_from_gxp(poptrtsD,gxpcounts)
    candidatePvalsDict=dict()
    
    for trt in popTreatments:
        candidatePvalsDict[trt]=dict()
        for locusTag,v in candidatesMultiplyMutated[trt].items():
            candidatePvalsDict[trt][locusTag]=dict()
            candidatePvalsDict[trt][locusTag]["repsWithMutation"]=candidatesMultiplyMutated[trt][locusTag]
            candidatePvalsDict[trt][locusTag]['pvalue']=-9
            candidatePvalsDict[trt][locusTag]['padj']=-9
            
        
        simulationsDF=pd.read_csv(workingdir+"\\Pipeline_XL\\3_Python\\"+dirname+"\\simout_"+dirname+"_"+trt+"_"+str(repsNumber)+".csv")#dirName = 'fS'
        simsD=pd.DataFrame.to_dict(simulationsDF)
        for cand in candidatesLists[trt]:
            myMoreextremeCount=np.count_nonzero([i>=candidatesCounts[trt][cand] for i in simsD[cand].values()])
            candidatePvalsDict[trt][cand]['pvalue']=myMoreextremeCount/len(simsD[cand].values())
    
        pvalList=list();padjListIndexCounter=0
        for cand,candD in candidatePvalsDict[trt].items():
            pvalList.append(candidatePvalsDict[trt][cand]['pvalue'])
        fdrnparray=statsmodels.stats.multitest.fdrcorrection(pvalList)
        for cand,candD in candidatePvalsDict[trt].items():
            candidatePvalsDict[trt][cand]['padj']=fdrnparray[1][padjListIndexCounter]
            padjListIndexCounter+=1
    return candidatePvalsDict
#%%
start=time.time();print(start)
candidatePvalsDictFixedNonsyn = get_pvalues_from_simulations(workingdir=workingdir, popTreatments=poptrtsL, poptrtsDict=poptrtsD, gxpcounts=gxpcounts_fixed, dirname='fN',repsNumber=5357)
end=time.time();print(end-start)
#%%
start=time.time();print(start)
candidatePvalsDictPolymorphicNonsyn = get_pvalues_from_simulations(workingdir=workingdir, popTreatments=poptrtsL, poptrtsDict=poptrtsD, gxpcounts=gxpcounts_p, dirname='pN',repsNumber=5357)
end=time.time();print(end-start)
#%%
#os.mkdir(path=workingdir+"\\Pipeline_XL\\3_Python\\pvalues");os.mkdir(path=workingdir+"\\Pipeline_XL\\3_Python\\pvalues\\fN");os.mkdir(path=workingdir+"\\Pipeline_XL\\3_Python\\pvalues\\pN")
for trt in poptrtsL:
    pd.DataFrame.from_dict(candidatePvalsDictFixedNonsyn[trt]).to_csv(index=True,path_or_buf=workingdir+"\\Pipeline_XL\\3_Python\\pvalues\\fN\\"+str(trt)+"pvalues.csv")
    pd.DataFrame.from_dict(candidatePvalsDictPolymorphicNonsyn[trt]).to_csv(index=True,path_or_buf=workingdir+"\\Pipeline_XL\\3_Python\\pvalues\\pN\\"+str(trt)+"pvalues.csv")
#%%
#smo1pd = pd.read_csv(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\simulations_syn1.0_100000.csv")
smopd = pd.read_csv(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\simulations_syn3B_GC_only.synonymous_100000.csv")
smo1pd = pd.read_csv(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\simulations_syn1.0_GC_only.synonymous_100000.csv")


smo = pd.DataFrame.to_dict(smopd)
smo1 = pd.DataFrame.to_dict(smo1pd)
end = time.time()
print(end - start)
    
    
    
    #%%
#read in the file
#for each locus
#make a dict st
#k = locus tag
#v = gene name
#write out the dict
def output_gene_names(infile,inDict=DD):
    mypanda=pd.read_csv(infile);outD=dict();counter=0
    for locus_tag in mypanda["Gene"]:
        outD[counter]=inDict[locus_tag];counter+=1
        print(inDict[locus_tag])
        
    writeout=pd.DataFrame.from_dict(outD,orient="index")
    return writeout
    
fixedPanda=output_gene_names(r"C:\Users\rmoge\OneDrive\20230124_yeast-status\candidate.genes.pvalues_fixed_for.Python.csv",DD);
polyPanda=output_gene_names(r"C:\Users\rmoge\OneDrive\20230124_yeast-status\candidate.genes.pvalues_polymorphic_for.Python.csv",DD);
fixedPanda.to_csv(r"C:\Users\rmoge\OneDrive\20230124_yeast-status\candidate.genes.pvalues_fixed_from.Python.csv")
polyPanda.to_csv(r"C:\Users\rmoge\OneDrive\20230124_yeast-status\candidate.genes.pvalues_polymorphic_from.Python.csv")
    
    
#%%
def get_pvalues_from_simulations(workingdir=workingdir,popTreatments=[],poptrtsDict=poptrtsD,gxpcounts={},dirname='',repsNumber=5357):
    candidatesMultiplyMutated,candidatesLists,candidatesCounts=output_candidates_from_gxp(poptrtsD,gxpcounts)
    candidatePvalsDict=dict()
    
    for trt in popTreatments:
        candidatePvalsDict[trt]=dict()
        for locusTag,v in candidatesMultiplyMutated[trt].items():
            candidatePvalsDict[trt][locusTag]=dict()
            candidatePvalsDict[trt][locusTag]["repsWithMutation"]=candidatesMultiplyMutated[trt][locusTag]
            candidatePvalsDict[trt][locusTag]['pvalue']=-9
            candidatePvalsDict[trt][locusTag]['padj']=-9
            

        
        simulationsDF=pd.read_csv(workingdir+"\\Pipeline_XL\\3_Python\\"+dirname+"\\simout_"+dirname+"_"+trt+"_"+str(repsNumber)+".csv")#dirName = 'fS'
        simsD=pd.DataFrame.to_dict(simulationsDF)
        for cand in candidatesLists[trt]:
            myMoreextremeCount=np.count_nonzero([i>=candidatesCounts[trt][cand] for i in simsD[cand].values()])
            candidatePvalsDict[trt][cand]['pvalue']=myMoreextremeCount/len(simsD[cand].values())
        pass
    

        pvalList=list();padjListIndexCounter=0
        for cand,candD in candidatePvalsDict[trt].items():
            pvalList.append(candidatePvalsDict[trt][cand]['pvalue'])
        fdrnparray=statsmodels.stats.multitest.fdrcorrection(pvalList)
        for cand,candD in candidatePvalsDict.items():
            candidatePvalsDict[trt][cand]['padj']=fdrnparray[1][padjListIndexCounter]
            padjListIndexCounter+=1
    return candidatePvalsDict
def get_pvalues_from_simulations(workingdir=workingdir,popTreatments=[],poptrtsDict=poptrtsD,gxpcounts={},dirname='',repsNumber=5357):
    candidatesMultiplyMutated,candidatesLists,candidatesCounts=output_candidates_from_gxp(poptrtsD,gxpcounts)
    candidatePvalsDict=dict()
    
    for t in popTreatments:
        candidatePvalsDict[t]=dict()
        for locusTag,v in candidatesMultiplyMutated[t].items():
            candidatePvalsDict[t][locusTag]=dict()
            candidatePvalsDict[t][locusTag]["repsWithMutation"]=candidatesMultiplyMutated[t][locusTag]
            candidatePvalsDict[t][locusTag]['pvalue']=-9
            candidatePvalsDict[t][locusTag]['padj']=-9
            
    for trt in popTreatments:
        
        simulationsDF=pd.read_csv(workingdir+"\\Pipeline_XL\\3_Python\\"+dirname+"\\simout_"+dirname+"_"+t+"_"+str(repsNumber)+".csv")#dirName = 'fS'
        simsD=pd.DataFrame.to_dict(simulationsDF)
        for cand in candidatesLists[trt]:
            myMoreextremeCount=np.count_nonzero([i>=candidatesCounts[trt][cand] for i in simsD[cand].values()])
            candidatePvalsDict[trt][cand]['pvalue']=myMoreextremeCount/len(simsD[cand].values())
        pass
    
    for trt in popTreatments:
        pvalList=list();padjListIndexCounter=0
        for cand,candD in candidatePvalsDict[trt].items():
            pvalList.append(candidatePvalsDict[trt][cand]['pvalue'])
        fdrnparray=statsmodels.stats.multitest.fdrcorrection(pvalList)
        for cand,candD in candidatePvalsDict.items():
            candidatePvalsDict[trt][cand]['padj']=fdrnparray[1][padjListIndexCounter]
            padjListIndexCounter+=1
    return candidatePvalsDict
