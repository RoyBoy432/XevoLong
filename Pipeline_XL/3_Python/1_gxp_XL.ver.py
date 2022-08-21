##get_multiplicity_matrix() written originally by Will Shoemaker, modified by Roy Moger-Reischer

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

output_to_keep = ['INS', 'DEL', 'SNP', 'SUB']
poptypes = ['rpl19a','pmr1','sch9','ypt','tor1','tif2','JD174','BY4742']
strains = ['BY4742']
my_anc_dirs=[81,82,83,84,85,86,87,88]

#%%
#NYI
def get_sites_to_remove_divergent_ancestors(strain):
    if strain == 'rpl19a':
        ancestor_path = workingdir + '\\Pipeline_XL\\1_Carbonate\\breseq\\breseq_25\\output\\output.gd'
        #C:\Users\rmoge\Documents\GitHub\Bifidobacterium\Pipeline_XL\1_Carbonate\breseq\breseq_25\output
    sites_to_remove = []
    for i, line in enumerate(open(ancestor_path, 'r')):
        line_split = line.strip().split('\t')
        if line_split[0] in output_to_keep:
            sites_to_remove.append( line_split[3] + '_' + str(line_split[4]))
# =============================================================================
#     Below, you will add sites that did no show up in the ancestor, but are obvious problematic sites, either due to cross contamination , or due to a mutation in a founding population for one of the treatments specifically
# =============================================================================
    sites_to_remove+=['CP001853_146344','CP001853_1305547','CP001853_1437318']
# =============================================================================
#     
# =============================================================================
    return sites_to_remove

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
def get_gxp_matrix(output_name,include_polymorphisms=True):

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
            gene_count_dict_pop = {}#initialize the dict of genes for THIS pop
            gene_by_pop_dict[pop] = {}#GxP matrix add "rows", i.e., sub-dicts
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
                    if [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsynonymous' or [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsense':
                        locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                        frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                        if ';' in locus_tag:
                            for locus_tag_j in locus_tag.split(';'):
                                if locus_tag_j not in gene_count_dict_pop:
                                    gene_count_dict_pop[locus_tag_j] = 0
                                gene_count_dict_pop[locus_tag_j] += frequency
                        else:
                            if locus_tag not in gene_count_dict_pop:
                                gene_count_dict_pop[locus_tag] = 0
                            gene_count_dict_pop[locus_tag] += frequency
                    else:
                        continue                            
                ######################################
                              
                else:
                    if len([s for s in line_split if 'gene_position=coding' in s]) >= 1:
                        locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                        frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                        if frequency != 1 and include_polymorphisms == False:
                            continue
                        if ';' in locus_tag:
                            for locus_tag_j in locus_tag.split(';'):
                                if locus_tag_j not in gene_count_dict_pop:
                                    gene_count_dict_pop[locus_tag_j] = 0
                                gene_count_dict_pop[locus_tag_j] += frequency

                        else:
                            if locus_tag not in gene_count_dict_pop:
                                gene_count_dict_pop[locus_tag] = 0
                            gene_count_dict_pop[locus_tag] += frequency
                    else:        
                        if line_split[0] == 'DEL':
                            if len([s for s in line_split if 'gene_position' in s]) == 0:
                                
                                if len([s for s in line_split if 'locus_tags_inactivated=' in s]) > 0:
                                    locus_tags_comma_unseparated = [s for s in line_split if 'locus_tags_inactivated=' in s][0].split('=')[1]
                                    frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                                    #print(locus_tags_comma_unseparated)
                                    if frequency != 1 and include_polymorphisms==False:
                                        continue
                                    #print(type(locus_tags_comma_unseparated))
                                    if locus_tags_comma_unseparated[-1]=='"':
                                        #print('it is 17')
                                        for locus_tag_j in locus_tags_comma_unseparated[:-1].split(','):
                                            if locus_tag_j not in gene_count_dict_pop:
                                                gene_count_dict_pop[locus_tag_j] = 0
                                            gene_count_dict_pop[locus_tag_j] += frequency
                                    else:
                                        for locus_tag_j in locus_tags_comma_unseparated.split(','):
                                            if locus_tag_j not in gene_count_dict_pop:
                                                gene_count_dict_pop[locus_tag_j] = 0
                                            gene_count_dict_pop[locus_tag_j] += frequency                                        
                                
                                elif len([s for s in line_split if 'locus_tags_overlapping=' in s]) > 0:        
                                    locus_tags_comma_unseparated = [s for s in line_split if 'locus_tags_overlapping=' in s][0].split('=')[1]
                                    frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                                    if frequency != 1 and include_polymorphisms==False:
                                        continue
                                    #print(locus_tags_comma_unseparated)
                                    for locus_tag_j in locus_tags_comma_unseparated.split(','):
                                        if locus_tag_j not in gene_count_dict_pop:
                                            gene_count_dict_pop[locus_tag_j] = 0
                                        gene_count_dict_pop[locus_tag_j] += frequency
                            else:
                                continue
                        else:
                            continue
                #print(line_split[3] + "_" + line_split[4])        
    #return gene_count_dict_pop
#            for locus_tag_i in gene_count_dict_pop.keys():
#                gene_by_pop_dict[pop][locus_tag_i] = gene_count_dict_pop[locus_tag_i]

            for locus_tag_i in gene_count_dict_pop.keys():
#                 mult_i = gene_parallelism_statistics[locus_tag_i]['multiplicity']
#                 if mult_i > 0:
                #print(locus_tag_i)
                try:
                    locus_tag_i_num = locus_tag_i.split('_')[1]
                    gene_by_pop_dict[pop][locus_tag_i_num] = gene_count_dict_pop[locus_tag_i]
                except IndexError:
                    gene_by_pop_dict[pop][locus_tag_i] = gene_count_dict_pop[locus_tag_i]
    

    gene_by_pop_df = pd.DataFrame(gene_by_pop_dict)#outputs the GxP matrix
    gene_by_pop_df = gene_by_pop_df.T
    gene_by_pop_df.fillna(0, inplace=True)

    gene_by_pop_df_out = workingdir + '\\Pipeline_XL\\3_Python\\' + output_name + '.csv'
    #gene_by_pop_df.to_csv(gene_by_pop_df_out, sep = '\t', index = True)
    gene_by_pop_df.to_csv(gene_by_pop_df_out, sep = ',', index = True)
    
    
    return gene_by_pop_dict


#%%
start=time.time()
test1=get_gxp_matrix('XL_gxp',True)
test2=get_gxp_matrix('XL_gxp_fixed.only',False)
end=time.time()
print(end-start)
#%%
def get_DN_and_DS(Ndict,Ddict):
    
    return

myc=0
for k,v in test1['1'].items():
    r=random.random()
    if v > r:
        myc+=1