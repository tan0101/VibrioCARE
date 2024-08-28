#!/usr/bin/env python
# coding: utf-8


import sys
import os
import cobra 
import pandas as pd
import numpy as np
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
import pickle as pkl
from pathlib import Path


#Load Data
AllgeneIDs=pd.read_csv('Data/geneIDs.csv')

SampleID=sys.argv[1]
print(SampleID)

#Load Data
AllgeneIDs=pd.read_csv('Data/geneIDs.csv')
AllgeneIDs=AllgeneIDs[AllgeneIDs['Attribute'] == SampleID]

#Load Model
m = cobra.io.read_sbml_model('Data/Model_specifications/'+SampleID+'_model.xml')
# gene dictionaries for mapping between ID and names
gene_name_dic = {}
for g in m.genes:
    gene_name_dic[g.name] = g.id

gene_name_dic2 = {}
for g in m.genes:
    gene_name_dic2[g.id] = g.name
    
    
# Load Clinical Genes

Genes= pd.merge(pd.read_csv('Data/Clinical_inputs.csv'), AllgeneIDs, how='inner', on='Gene')


genes_clinical = []
for ind, g in enumerate(Genes['Value']):
    if str(g) != 'nan':
        for gene in m.genes:
            if gene.id in g:
                genes_clinical.append(gene.id)                    
            
genes_clincal= list(set(genes_clinical))

       
clinicalgenes = pd.DataFrame({'col':genes_clinical})

# save dataframe

clinicalgenes.to_pickle('GSMGenesClinical'+SampleID+'.pkl')
clinicalgenes.to_csv('GSMGenesClinical'+SampleID+'.csv')


#------------------------------ Gene Essentiality Clinical Genes- Rich Media-------------------------------------------#
#Load Model
m = cobra.io.read_sbml_model('Data/Model_specifications/'+SampleID+'_model.xml')
# gene dictionaries for mapping between ID and names
gene_name_dic = {}
for g in m.genes:
    gene_name_dic[g.name] = g.id

gene_name_dic2 = {}
for g in m.genes:
    gene_name_dic2[g.id] = g.name


genes_ALL = list(genes_clinical)

m.reactions.get_by_id('EX_o2_e').bounds = (-18.5, 0.0)
# get the carbon sources from the model
listSources=[]
for r in m.reactions:
    if 'EX_' in r.id:
        listSources.append(r.id)

m.reactions.get_by_id('EX_glc__D_e').lower_bound = 0.0
m.reactions.get_by_id('EX_o2_e').lower_bound = -18.5
medium = m.medium

    
# run in rich media
col = []
for source in listSources:
    m.reactions.get_by_id(source).lower_bound=-1.0

# get values for WT

try:
    sol_bio = m.optimize()
    if sol_bio.status != 'infeasible':
        col.append(sol_bio['Growth'])
        WT_growth = sol_bio['Growth']
    else:
        WT_growth = 0.0
        col.append(0.0)

except:
    WT_growth = 0.0
    col.append(0.0)



# run for all knockouts
count = 1
KO_growth_rich = {}
for gene in genes_ALL:
    with m: 
        col = []
        m.genes.get_by_id(gene).knock_out() 

        try:
            sol_bio = m.optimize()
            if sol_bio.status != 'infeasible':
                KO_growth_rich[gene] = sol_bio['Growth']
                col.append(sol_bio['Growth'])
            else:
                KO_growth_rich[gene] = 0.0
                col.append(0.0)

        except:
            KO_growth_rich[gene.id] = 0.0
            col.append(0.0)
            
growth_limiting_rich = []
count = 0
for i, j in KO_growth_rich.items():
    if j < 0.0001:
        
        growth_limiting_rich.append(i)


Path("GrowthlimitingRich_results").mkdir(parents=True, exist_ok=True)
growthlimitingrich = pd.DataFrame({'col':growth_limiting_rich})
growthlimitingrich.to_csv('GrowthlimitingRich_results/clinicalessentialinrich_'+SampleID+'.csv')

#------------------------------ Gene Essentiality Clinical Genes- Rich Media-------------------------------------------#
m = cobra.io.read_sbml_model('Data/Model_specifications/'+SampleID+'_model.xml')

#In silico M9 minimal media + glucose aerobic;

m.reactions.get_by_id('EX_glc__D_e').lower_bound = -10.0

m.reactions.get_by_id('EX_o2_e').bounds = (-18.5, 0.0)
medium = m.medium


# get values for WT

try:
    sol_bio = m.optimize()
    if sol_bio.status != 'infeasible':
        WT_growth = sol_bio['Growth']
    else:
        WT_growth = 0.0

except:
    WT_growth = 0.0


# run for all knockouts

KO_growth = {}
for gene in genes_ALL:
    
    with m: 
        col = []
        m.genes.get_by_id(gene).knock_out() 

        try:
            sol_bio = m.optimize()
            if sol_bio.status != 'infeasible':
                KO_growth[gene] = sol_bio['Growth']
            else:
                KO_growth[gene] = 0.0

        except:
            KO_growth[gene] = 0.0
growth_limiting_minimal = []
count = 0
for i, j in KO_growth.items():
    if j < 0.0001:
        if i not in growth_limiting_rich:
            
            growth_limiting_minimal.append(i)
            count += 1

Path("GrowthlimitingMinimal_results").mkdir(parents=True, exist_ok=True)
growthlimitingminimal = pd.DataFrame({'col':growth_limiting_minimal})
growthlimitingminimal.to_csv('GrowthlimitingMinimal_results/clinicalessentialinminimal_'+SampleID+'.csv')



# --------------------------Auxotrophy Clinical----------------------------#
m = cobra.io.read_sbml_model('Data/Model_specifications/'+SampleID+'_model.xml')

# aerobic
m.reactions.get_by_id('EX_o2_e').bounds = (-18.5, 0.0)

# gene dictionaries for mapping between ID and names
gene_name_dic = {}
for g in m.genes:
    gene_name_dic[g.name] = g.id
    
gene_name_dic2 = {}
for g in m.genes:
    gene_name_dic2[g.id] = g.name
    
listSources=[]
for r in m.reactions:
    if 'EX_' in r.id:
        listSources.append(r.id)

#m.reactions.get_by_id('EX_glc-D(e)').bounds = (-10.0, 1000.0)
m.reactions.get_by_id('EX_o2_e').lower_bound = -18.5
medium = m.medium

CgrowthCapabilities_auxo = pd.DataFrame(index=listSources)


# run for all knockouts
count = 0
for gene in growth_limiting_rich:
    
    
    with m: 
        col = []
        #gene_id = gene_name_dic[gene]
        m.genes.get_by_id(gene).knock_out() 
        
        for source in listSources:
            m.reactions.get_by_id(source).bounds= (-10.0, 1000.0)


            try:
                sol_bio = m.optimize()
                if sol_bio.status != 'infeasible':
                    col.append(sol_bio['Growth'])
                else:
                    col.append(0.0)

            except:
                col.append(0.0)
    
    

            # turn transporter off
            if source not in medium.keys():
                m.reactions.get_by_id(source).bounds = (0.0, 0.0)
        
        CgrowthCapabilities_auxo.insert(count, column = gene, value = col)

Path("Auxotrophy").mkdir(parents=True, exist_ok=True)
CgrowthCapabilities_auxo.to_csv('Auxotrophy/'+SampleID+'_auxotrophy_clinical.csv')
# --------------------------Alternative Carbon Sources Clinical----------------------------#
m.reactions.get_by_id('EX_o2_e').bounds = (-18.5, 0.0)

# get the carbon sources from the model
listPotCarbonSources=[]
for r in m.reactions:
    if 'EX_' in r.id:
        for met in r.metabolites:
            if 'C' in met.formula:
                listPotCarbonSources.append(r.id)
                
m.reactions.get_by_id('EX_glc__D_e').lower_bound = 0.0

m.reactions.get_by_id('EX_o2_e').lower_bound = -18.5
medium = m.medium

CgrowthCapabilities_carbon = pd.DataFrame(index=listPotCarbonSources)

# get values for WT
 
col = []
for source in listPotCarbonSources:
    m.reactions.get_by_id(source).lower_bound=-10.0
    
    try:
        sol_bio = m.optimize()
        if sol_bio.status != 'infeasible':
            col.append(sol_bio['Growth']) 
        else:
            col.append(0.0)

    except:
        col.append(0.0)
        

# run for all knockouts
count = 0
for gene in genes_clincal:
    
    with m: 
        col = []
        m.genes.get_by_id(gene).knock_out() 
        
        for source in listPotCarbonSources:
            
            m.reactions.get_by_id(source).lower_bound=-10.0


            try:
                sol_bio = m.optimize()
                if sol_bio.status != 'infeasible':
                    col.append(sol_bio['Growth'])
                else:
                    col.append(0.0)

            except:
                col.append(0.0)
    
    
        
        # turn transporter off
        if source not in medium.keys():
            m.reactions.get_by_id(source).lower_bound = 0.0
        
        CgrowthCapabilities_carbon.insert(count, column = gene_name_dic2[gene], value = col)

Path("Alternativecarbon").mkdir(parents=True, exist_ok=True)
CgrowthCapabilities_carbon.to_csv('Alternativecarbon/'+SampleID+'_alternativecarbon_clinical.csv')

# --------------------------FVA analysis - Clinical----------------------------#

m = cobra.io.read_sbml_model('Data/Model_specifications/'+SampleID+'_model.xml')
# gene dictionaries for mapping between ID and names
gene_name_dic = {}
for g in m.genes:
    gene_name_dic[g.name] = g.id

gene_name_dic2 = {}
for g in m.genes:
    gene_name_dic2[g.id] = g.name
    
genes_ALL = pd.read_pickle('GSMGenesClinical'+SampleID+'.pkl')

## FVA on glucose only
m.reactions.EX_glc__D_e.bounds = (-10.0, 1000.0)
m.reactions.EX_o2_e.bounds = (-18.5, 1000.0)


    
for reaction in m.reactions:
    if reaction.bounds != (0.0, 0.0):
        if reaction.reversibility:
            reaction.bounds = (-1, 1)
        else:
            reaction.bounds = (0, 1)

fva_res = {}

# compute for WT first
res = cobra.flux_analysis.flux_variability_analysis(m, fraction_of_optimum = 0.0, processes = 12)
fva_res['WT'] = res.to_dict()

for index, gene in genes_ALL.iterrows():

    with m:
        m.genes.get_by_id(gene['col']).knock_out()
        res = cobra.flux_analysis.flux_variability_analysis(m, fraction_of_optimum = 0.0, processes = 12)
    fva_res[gene['col']] = res.to_dict()



reactions = []
for r in m.reactions:
    reactions.append(r.id)
FVA_genes_clinical = pd.DataFrame(index = reactions)

count = 0
for strain, res in fva_res.items():
    results = pd.DataFrame.from_dict(res)
    
    flux_span = pd.DataFrame(results['maximum']-results['minimum'])
    if strain != 'WT':
        FVA_genes_clinical.insert(count, column = strain, value = flux_span.values)
    else:
        FVA_genes_clinical.insert(count, column = strain, value = flux_span.values)
        count += 1 
        
SignificantFVAReactions = []
SignificantFVAGenes = []
for gene in FVA_genes_clinical.columns:
    if gene != 's0001':
        for ind, i in enumerate(FVA_genes_clinical[gene]):
            if FVA_genes_clinical['WT'][ind] > FVA_genes_clinical[gene][ind]:
                    if FVA_genes_clinical['WT'][ind] > 0.01:
                        if FVA_genes_clinical[gene][ind]/FVA_genes_clinical['WT'][ind] < 0.9:
                            SignificantFVAReactions.append((m.genes.get_by_id(gene).name, FVA_genes_clinical['WT'].index[ind]))
                            SignificantFVAGenes.append(m.genes.get_by_id(gene).name)

SignificantFVAreacts = pd.DataFrame({'Result':SignificantFVAReactions})
SignificantFVAgenes = pd.DataFrame({'Result':SignificantFVAGenes})
SignificantFVAgenes=SignificantFVAgenes.drop_duplicates()
SignificantFVAgenes = pd.merge(left=SignificantFVAgenes, right=AllgeneIDs, how='left',right_on='Value',left_on='Result')
Path("FVAresults_clinical").mkdir(parents=True, exist_ok=True)
SignificantFVAreacts.to_csv('FVAresults_clinical/SignificantFVAreactions_'+SampleID+'.csv')
SignificantFVAgenes.to_csv('FVAresults_clinical/SignificantFVAgenes_'+SampleID+'.csv')

#-------------------------------FBA analysis Clinical-----------------------------------------------#
m = cobra.io.read_sbml_model('Data/Model_specifications/'+SampleID+'_model.xml')

gene_name_dic = {}
for g in m.genes:
    gene_name_dic[g.name] = g.id
    
gene_name_dic2 = {}
for g in m.genes:
    gene_name_dic2[g.id] = g.name

genes_temp = pd.read_pickle('GSMGenesClinical'+SampleID+'.pkl')
genes_ALL=genes_temp['col'].to_list()

#####Calculate on HPC
# turn off growth as objective function
m.reactions.Growth.objective_coefficient = 0.0

metabolites = []
for met in m.metabolites:
    metabolites.append(met.id)
    
# calculate for wild type
wt_yields = {}
with m:
    for met in metabolites:
        
        dm_id = str(met) + '_dummy_drain'
        DM_met = cobra.core.Reaction(dm_id)
        DM_met.name = met + ' dummy drain'
        DM_met.lower_bound = 0.0
        DM_met.upper_bound = 0.0
        DM_met.add_metabolites({m.metabolites.get_by_id(met) : -1.0})
        m.add_reactions([DM_met])

        m.reactions.get_by_id(dm_id).bounds = (0.0, 1000.0)
        m.objective = m.reactions.get_by_id(dm_id)           
        # set direction
        m.objective_direction = 'max'
        try:
            sol = m.optimize()
            
            if sol.status == 'optimal':
                if sol[dm_id] > 0.00001:
                    wt_yields[dm_id] = sol[dm_id]
                else:
                    wt_yields[dm_id] = 0.0

            else:
                wt_yields[dm_id] = 0.0
        except:
            wt_yields[dm_id] = 0.0



# repeat for all the genetic determinants
MetYields_all = pd.DataFrame(index = metabolites)
MetYields_all.insert(0, column = 'wt', value = wt_yields.values())

MetYields_changes_all = pd.DataFrame(index = metabolites)
MetYields_changes_all.insert(0, column = 'wt', value = [1]*len(wt_yields.keys()))

count = 1
for gene in genes_ALL:

    ko_yields = []
    ko_changes = []
    with m:
        m.genes.get_by_id(gene).knock_out()
        for met in metabolites:
            
            # add drain reaction for metabolite
            dm_id = str(met) + '_dummy_drain'
            DM_met = cobra.core.Reaction(dm_id)
            DM_met.name = met + ' dummy drain'
            DM_met.lower_bound = 0.0
            DM_met.upper_bound = 0.0
            DM_met.add_metabolites({m.metabolites.get_by_id(met) : -1.0})
            m.add_reactions([DM_met])

            m.reactions.get_by_id(dm_id).bounds = (0.0, 1000.0)
            m.objective = m.reactions.get_by_id(dm_id)   # set the maximisation of metabolite production flux as objective        
            # set direction
            m.objective_direction = 'max'
            try:
                sol = m.optimize()
                if sol.status == 'optimal':
                    if sol[dm_id] > 0.00001:
                        ko_yields.append(sol[dm_id])
                        ko_changes.append(sol[dm_id]/wt_yields[dm_id])
                    else:
                        ko_yields.append(0.0)
                        ko_changes.append(0.0)
                else:
                    ko_yields.append(0.0)
                    ko_changes.append(0.0)
            except:
                ko_yields.append(0.0)
                ko_changes.append(0.0)
    MetYields_all.insert(count, column = gene, value = ko_yields) # add the yields of metabolites to dataframe
    MetYields_changes_all.insert(count, column = gene, value = ko_changes) # add change in metabolite yield compared to WT metabolite yield
    count += 1

# get the metabolite yields for just the genes in our list
MetYieldsClinical = pd.DataFrame(index = MetYields_all.index)
count = 0
MetYieldsClinical.insert(0, column = 'wt', value = list(MetYields_all['wt'].values))
for ind, i in enumerate(MetYields_all.columns):
    if i in genes_ALL:
        count += 1
        MetYieldsClinical.insert(count, column = i, value = list(MetYields_all[i].values))

SignificantMetabolites = []
SignificantGenes = []
for ind, g in enumerate(MetYields_all.columns):
    if g != 'wt':
        for ind_met, met in enumerate(MetYields_all.index):
            if MetYields_all['wt'][ind_met] > 0.0001:
                if MetYields_all[g][ind_met]/MetYields_all['wt'][ind_met] < 0.00001:
                   SignificantMetabolites.append((m.genes.get_by_id(g).name, met))
                   SignificantGenes.append((m.genes.get_by_id(g).name))


SignificantFBAgenes=pd.DataFrame.from_dict(SignificantGenes)
SignificantFBAmets = pd.DataFrame.from_dict(SignificantMetabolites) 
SignificantFBAgenes=SignificantFBAgenes.drop_duplicates()

SignificantFBAgenes.rename({0:'gene', 1:'metabolite'},axis='columns',inplace=True)
SignificantFBAgenes = pd.merge(left=SignificantFBAgenes, right=AllgeneIDs, how='left',right_on='Value',left_on='gene')
Path("FBAresults_clinical").mkdir(parents=True, exist_ok=True)
SignificantFBAmets.to_csv('FBAresults_clinical/SignificantFBAmetabolites_'+SampleID+'.csv')
SignificantFBAgenes.to_csv('FBAresults_clinical/SignificantFBAgenes_'+SampleID+'.csv')

#---------------------------------Load Data Linage--------------------------------#

m = cobra.io.read_sbml_model('Data/Model_specifications/'+SampleID+'_model.xml')
# gene dictionaries for mapping between ID and names
gene_name_dic = {}
for g in m.genes:
    gene_name_dic[g.name] = g.id

gene_name_dic2 = {}
for g in m.genes:
    gene_name_dic2[g.id] = g.name
    
    
# Load Genes
Genes= pd.merge(pd.read_csv('Data/Lineage_inputs.csv'), AllgeneIDs, how='inner', on='Gene')  
genes_lineage = []
for ind, g in enumerate(Genes['Value']):
    if str(g) != 'nan':
        for gene in m.genes:
            if gene.id in g:
                genes_lineage.append(gene.id)                    
            
genes_lineage= list(set(genes_lineage))



lineagegenes = pd.DataFrame({'col':genes_lineage})

# save dataframe

lineagegenes.to_pickle('GSMGenesLineage'+SampleID+'.pkl')
lineagegenes.to_csv('GSMGenesLineage'+SampleID+'.csv')


# ---------------------------------- Gene Essentiality Linage - Rich Media --------------------------------------------#
genes_ALL = list(genes_lineage)

m.reactions.get_by_id('EX_o2_e').bounds = (-18.5, 0.0)
# get the carbon sources from the model
listSources=[]
for r in m.reactions:
    if 'EX_' in r.id:
        listSources.append(r.id)

m.reactions.get_by_id('EX_glc__D_e').lower_bound = 0.0
m.reactions.get_by_id('EX_o2_e').lower_bound = -18.5
medium = m.medium


# run in rich media
col = []
for source in listSources:
    m.reactions.get_by_id(source).lower_bound=-1.0

# get values for WT

try:
    sol_bio = m.optimize()
    if sol_bio.status != 'infeasible':
        col.append(sol_bio['Growth'])
        WT_growth = sol_bio['Growth']
    else:
        WT_growth = 0.0
        col.append(0.0)

except:
    WT_growth = 0.0
    col.append(0.0)

    
# run for all knockouts
count = 1
KO_growth_rich = {}
for gene in genes_ALL:
    
    with m: 
        col = []
        m.genes.get_by_id(gene).knock_out() 

        try:
            sol_bio = m.optimize()
            if sol_bio.status != 'infeasible':
                KO_growth_rich[gene] = sol_bio['Growth']
                col.append(sol_bio['Growth'])
            else:
                KO_growth_rich[gene] = 0.0
                col.append(0.0)

        except:
            KO_growth_rich[gene.id] = 0.0
            col.append(0.0)
            
growth_limiting_rich = []
count = 0
for i, j in KO_growth_rich.items():
    if j < 0.0001:
        
        growth_limiting_rich.append(i)

Path("GrowthlimitingRich_results").mkdir(parents=True, exist_ok=True)
growthlimitingrich = pd.DataFrame({'col':growth_limiting_rich})
growthlimitingrich.to_csv('GrowthlimitingRich_results/lineageessentialinrich_'+SampleID+'.csv')

# ---------------------------------- Gene Essentiality Linage - Minimal Media --------------------------------------------#
m = cobra.io.read_sbml_model('Data/Model_specifications/'+SampleID+'_model.xml')

#In silico M9 minimal media + glucose aerobic;

m.reactions.get_by_id('EX_glc__D_e').lower_bound = -10.0

m.reactions.get_by_id('EX_o2_e').bounds = (-18.5, 0.0)
medium = m.medium


# get values for WT

try:
    sol_bio = m.optimize()
    if sol_bio.status != 'infeasible':
        WT_growth = sol_bio['Growth']
    else:
        WT_growth = 0.0

except:
    WT_growth = 0.0


# run for all knockouts

KO_growth = {}
for gene in genes_ALL:
    
    with m: 
        col = []
        m.genes.get_by_id(gene).knock_out() 

        try:
            sol_bio = m.optimize()
            if sol_bio.status != 'infeasible':
                KO_growth[gene] = sol_bio['Growth']
            else:
                KO_growth[gene] = 0.0

        except:
            KO_growth[gene] = 0.0
growth_limiting_minimal = []
count = 0
for i, j in KO_growth.items():
    if j < 0.0001:
        if i not in growth_limiting_rich:
            
            growth_limiting_minimal.append(i)
            count += 1
Path("GrowthlimitingMinimal_results").mkdir(parents=True, exist_ok=True)
growthlimitingminimal = pd.DataFrame({'col':growth_limiting_minimal})
growthlimitingminimal.to_csv('GrowthlimitingMinimal_results/lineageessentialinminimal_'+SampleID+'.csv')

# --------------------------Auxotrophy lineage----------------------------#
m = cobra.io.read_sbml_model('Data/Model_specifications/'+SampleID+'_model.xml')

# aerobic
m.reactions.get_by_id('EX_o2_e').bounds = (-18.5, 0.0)

# gene dictionaries for mapping between ID and names
gene_name_dic = {}
for g in m.genes:
    gene_name_dic[g.name] = g.id
    
gene_name_dic2 = {}
for g in m.genes:
    gene_name_dic2[g.id] = g.name
    
listSources=[]
for r in m.reactions:
    if 'EX_' in r.id:
        listSources.append(r.id)

#m.reactions.get_by_id('EX_glc-D(e)').bounds = (-10.0, 1000.0)
m.reactions.get_by_id('EX_o2_e').lower_bound = -18.5
medium = m.medium

CgrowthCapabilities_auxo = pd.DataFrame(index=listSources)


# run for all knockouts
count = 0
for gene in growth_limiting_rich:
    
    
    with m: 
        col = []
        #gene_id = gene_name_dic[gene]
        m.genes.get_by_id(gene).knock_out() 
        
        for source in listSources:
            m.reactions.get_by_id(source).bounds= (-10.0, 1000.0)


            try:
                sol_bio = m.optimize()
                if sol_bio.status != 'infeasible':
                    col.append(sol_bio['Growth'])
                else:
                    col.append(0.0)

            except:
                col.append(0.0)
    
    

            # turn transporter off
            if source not in medium.keys():
                m.reactions.get_by_id(source).bounds = (0.0, 0.0)
        
        CgrowthCapabilities_auxo.insert(count, column = gene, value = col)

Path("Auxotrophy").mkdir(parents=True, exist_ok=True)
CgrowthCapabilities_auxo.to_csv('Auxotrophy/'+SampleID+'_auxotrophy_lineage.csv')

# --------------------------Alternative Carbon Sources Lineage----------------------------#
m.reactions.get_by_id('EX_o2_e').bounds = (-18.5, 0.0)

# get the carbon sources from the model
listPotCarbonSources=[]
for r in m.reactions:
    if 'EX_' in r.id:
        for met in r.metabolites:
            if 'C' in met.formula:
                listPotCarbonSources.append(r.id)
                
m.reactions.get_by_id('EX_glc__D_e').lower_bound = 0.0

m.reactions.get_by_id('EX_o2_e').lower_bound = -18.5
medium = m.medium

CgrowthCapabilities_carbon = pd.DataFrame(index=listPotCarbonSources)

# get values for WT
 
col = []
for source in listPotCarbonSources:
    m.reactions.get_by_id(source).lower_bound=-10.0
    
    try:
        sol_bio = m.optimize()
        if sol_bio.status != 'infeasible':
            col.append(sol_bio['Growth']) 
        else:
            col.append(0.0)

    except:
        col.append(0.0)
        
  

# run for all knockouts
count = 0
for gene in genes_lineage:
    
    with m: 
        col = []
        m.genes.get_by_id(gene).knock_out() 
        
        for source in listPotCarbonSources:
            
            m.reactions.get_by_id(source).lower_bound=-10.0


            try:
                sol_bio = m.optimize()
                if sol_bio.status != 'infeasible':
                    col.append(sol_bio['Growth'])
                else:
                    col.append(0.0)

            except:
                col.append(0.0)
    
    
        
        # turn transporter off
        if source not in medium.keys():
            m.reactions.get_by_id(source).lower_bound = 0.0
        
        CgrowthCapabilities_carbon.insert(count, column = gene_name_dic2[gene], value = col)
        
Path("Alternativecarbon").mkdir(parents=True, exist_ok=True)
CgrowthCapabilities_carbon.to_csv('Alternativecarbon/'+SampleID+'_alternativecarbon_lineage.csv')




# ---------------------------------- FVA Analysis Linage  --------------------------------------------#

m = cobra.io.read_sbml_model('Data/Model_specifications/'+SampleID+'_model.xml')
# gene dictionaries for mapping between ID and names
gene_name_dic = {}
for g in m.genes:
    gene_name_dic[g.name] = g.id

gene_name_dic2 = {}
for g in m.genes:
    gene_name_dic2[g.id] = g.name
    
genes_ALL = pd.read_pickle('GSMGenesLineage'+SampleID+'.pkl')

## FVA on glucose only
m.reactions.EX_glc__D_e.bounds = (-10.0, 1000.0)
m.reactions.EX_o2_e.bounds = (-18.5, 1000.0)


    
for reaction in m.reactions:
    if reaction.bounds != (0.0, 0.0):
        if reaction.reversibility:
            reaction.bounds = (-1, 1)
        else:
            reaction.bounds = (0, 1)

fva_res = {}

# compute for WT first
res = cobra.flux_analysis.flux_variability_analysis(m, fraction_of_optimum = 0.0, processes = 12)
fva_res['WT'] = res.to_dict()

for index, gene in genes_ALL.iterrows():

    with m:
        m.genes.get_by_id(gene['col']).knock_out()
        res = cobra.flux_analysis.flux_variability_analysis(m, fraction_of_optimum = 0.0, processes = 12)
    fva_res[gene['col']] = res.to_dict()



reactions = []
for r in m.reactions:
    reactions.append(r.id)
FVA_genes_lineage = pd.DataFrame(index = reactions)

count = 0
for strain, res in fva_res.items():
    results = pd.DataFrame.from_dict(res)
    
    flux_span = pd.DataFrame(results['maximum']-results['minimum'])
    if strain != 'WT':
        FVA_genes_lineage.insert(count, column = strain, value = flux_span.values)
    else:
        FVA_genes_lineage.insert(count, column = strain, value = flux_span.values)
        count += 1 
        
SignificantFVAReactions = []
SignificantFVAGenes = []
for gene in FVA_genes_lineage.columns:
    if gene != 's0001':
        for ind, i in enumerate(FVA_genes_lineage[gene]):
            if FVA_genes_lineage['WT'][ind] > FVA_genes_lineage[gene][ind]:
                    if FVA_genes_lineage['WT'][ind] > 0.01:
                        if FVA_genes_lineage[gene][ind]/FVA_genes_lineage['WT'][ind] < 0.9:
                            SignificantFVAReactions.append((m.genes.get_by_id(gene).name, FVA_genes_lineage['WT'].index[ind]))
                            SignificantFVAGenes.append(m.genes.get_by_id(gene).name)

SignificantFVAreacts = pd.DataFrame({'Result':SignificantFVAReactions})
SignificantFVAgenes = pd.DataFrame({'Result':SignificantFVAGenes})
SignificantFVAgenes=SignificantFVAgenes.drop_duplicates()
SignificantFVAgenes = pd.merge(left=SignificantFVAgenes, right=AllgeneIDs, how='left',right_on='Value',left_on='Result')
Path("FVAresults_lineage").mkdir(parents=True, exist_ok=True)
SignificantFVAreacts.to_csv('FVAresults_lineage/SignificantFVAreactions_'+SampleID+'.csv')
SignificantFVAgenes.to_csv('FVAresults_lineage/SignificantFVAgenes_'+SampleID+'.csv')



# ---------------------------------- FBA Analysis Linage  --------------------------------------------#
m = cobra.io.read_sbml_model('Data/Model_specifications/'+SampleID+'_model.xml')

gene_name_dic = {}
for g in m.genes:
    gene_name_dic[g.name] = g.id
    
gene_name_dic2 = {}
for g in m.genes:
    gene_name_dic2[g.id] = g.name

genes_temp = pd.read_pickle('GSMGenesLineage'+SampleID+'.pkl')
genes_ALL=genes_temp['col'].to_list()

#####Calculate on HPC
# turn off growth as objective function
m.reactions.Growth.objective_coefficient = 0.0

metabolites = []
for met in m.metabolites:
    metabolites.append(met.id)
    
# calculate for wild type
wt_yields = {}
with m:
    for met in metabolites:
        
        dm_id = str(met) + '_dummy_drain'
        DM_met = cobra.core.Reaction(dm_id)
        DM_met.name = met + ' dummy drain'
        DM_met.lower_bound = 0.0
        DM_met.upper_bound = 0.0
        DM_met.add_metabolites({m.metabolites.get_by_id(met) : -1.0})
        m.add_reactions([DM_met])

        m.reactions.get_by_id(dm_id).bounds = (0.0, 1000.0)
        m.objective = m.reactions.get_by_id(dm_id)           
        # set direction
        m.objective_direction = 'max'
        try:
            sol = m.optimize()
            
            if sol.status == 'optimal':
                if sol[dm_id] > 0.00001:
                    wt_yields[dm_id] = sol[dm_id]
                else:
                    wt_yields[dm_id] = 0.0

            else:
                wt_yields[dm_id] = 0.0
        except:
            wt_yields[dm_id] = 0.0



# repeat for all the genetic determinants
MetYields_all = pd.DataFrame(index = metabolites)
MetYields_all.insert(0, column = 'wt', value = wt_yields.values())

MetYields_changes_all = pd.DataFrame(index = metabolites)
MetYields_changes_all.insert(0, column = 'wt', value = [1]*len(wt_yields.keys()))

count = 1
for gene in genes_ALL:

    ko_yields = []
    ko_changes = []
    with m:
        m.genes.get_by_id(gene).knock_out()
        for met in metabolites:
            
            # add drain reaction for metabolite
            dm_id = str(met) + '_dummy_drain'
            DM_met = cobra.core.Reaction(dm_id)
            DM_met.name = met + ' dummy drain'
            DM_met.lower_bound = 0.0
            DM_met.upper_bound = 0.0
            DM_met.add_metabolites({m.metabolites.get_by_id(met) : -1.0})
            m.add_reactions([DM_met])

            m.reactions.get_by_id(dm_id).bounds = (0.0, 1000.0)
            m.objective = m.reactions.get_by_id(dm_id)   # set the maximisation of metabolite production flux as objective        
            # set direction
            m.objective_direction = 'max'
            try:
                sol = m.optimize()
                if sol.status == 'optimal':
                    if sol[dm_id] > 0.00001:
                        ko_yields.append(sol[dm_id])
                        ko_changes.append(sol[dm_id]/wt_yields[dm_id])
                    else:
                        ko_yields.append(0.0)
                        ko_changes.append(0.0)
                else:
                    ko_yields.append(0.0)
                    ko_changes.append(0.0)
            except:
                ko_yields.append(0.0)
                ko_changes.append(0.0)
    MetYields_all.insert(count, column = gene, value = ko_yields) # add the yields of metabolites to dataframe
    MetYields_changes_all.insert(count, column = gene, value = ko_changes) # add change in metabolite yield compared to WT metabolite yield
    count += 1

# get the metabolite yields for just the genes in our list
MetYieldsClinical = pd.DataFrame(index = MetYields_all.index)
count = 0
MetYieldsClinical.insert(0, column = 'wt', value = list(MetYields_all['wt'].values))
for ind, i in enumerate(MetYields_all.columns):
    if i in genes_ALL:
        count += 1
        MetYieldsClinical.insert(count, column = i, value = list(MetYields_all[i].values))

SignificantMetabolites = []
SignificantGenes = []
for ind, g in enumerate(MetYields_all.columns):
    if g != 'wt':
        for ind_met, met in enumerate(MetYields_all.index):
            if MetYields_all['wt'][ind_met] > 0.0001:
                if MetYields_all[g][ind_met]/MetYields_all['wt'][ind_met] < 0.00001:
                   SignificantMetabolites.append((m.genes.get_by_id(g).name, met))
                   SignificantGenes.append((m.genes.get_by_id(g).name))


SignificantFBAgenes=pd.DataFrame.from_dict(SignificantGenes)
SignificantFBAmets = pd.DataFrame.from_dict(SignificantMetabolites) 
SignificantFBAgenes=SignificantFBAgenes.drop_duplicates()


SignificantFBAgenes.rename({0:'gene', 1:'metabolite'},axis='columns',inplace=True)
SignificantFBAgenes = pd.merge(left=SignificantFBAgenes, right=AllgeneIDs, how='left',right_on='Value',left_on='gene')
Path("FBAresults_lineage").mkdir(parents=True, exist_ok=True)
SignificantFBAmets.to_csv('FBAresults_lineage/SignificantFBAmetabolites_'+SampleID+'.csv')
SignificantFBAgenes.to_csv('FBAresults_lineage/SignificantFBAgenes_'+SampleID+'.csv')

