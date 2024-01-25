from multiprocessing import freeze_support
import cobra
import numpy as np
import pandas as pd
from cobra import Reaction, Metabolite
from cobra.flux_analysis import flux_variability_analysis
import re


def from_gene_to_rxns(model, gene):
    rxns = model.genes.get_by_id(gene).reactions
    ids = [r.id for r in rxns]
    return ids


def from_rxn_to_genes(model, rxn):
    genes = model.reactions.get_by_id(rxn).genes
    ids = [g.id for g in genes]
    return ids



def prepare_model():
    model = cobra.io.load_matlab_model('../model/iML1515.mat')
    # add C4 rxns
    ala_e = Metabolite('5aop_e',
                       compartment='e')
    model.add_metabolites(ala_e)
    ala_syth = Reaction('ala_C4')
    ala_syth.add_metabolites({
        model.metabolites.get_by_id('succoa_c'): -1,
        model.metabolites.get_by_id('gly_c'): -1,
        model.metabolites.get_by_id('h_c'): -1,
        model.metabolites.get_by_id('5aop_c'): 1,
        model.metabolites.get_by_id('co2_c'): 1,
        model.metabolites.get_by_id('coa_c'): 1
    })
    ala_syth.gene_reaction_rule = 'ALAS'
    model.add_reactions([ala_syth])

    ala_exchange = Reaction('ala_ex')
    ala_exchange.bounds = (-1000, 1000)
    ala_exchange.add_metabolites({
        model.metabolites.get_by_id('5aop_c'): -1,
        model.metabolites.get_by_id('5aop_e'): 1,
    })
    ala_exchange.gene_reaction_rule = 'eamA'
    model.add_reactions([ala_exchange])

    ala_transport = Reaction('EX_ALA_e')
    ala_transport.bounds = (0, 1000)
    ala_transport.add_metabolites({
        model.metabolites.get_by_id('5aop_e'): -1,
    })
    model.add_reactions([ala_transport])
    model.objective = model.reactions.get_by_id(ala_transport.id)
    return model


def get_rxn_score(r, targets):
    gene = [g.id for g in r.genes]
    score = [targets.loc[gid, 0] for gid in gene if gid in targets.index.tolist()]
    if len(score) > 1:
        assert np.min(score) > 1 or np.max(score) < 1
    return score


targets = pd.read_excel('../output/First_round_targets.xlsx', index_col=0)
genes = targets.index.values.tolist()
model = prepare_model()

ids = []
for g in genes:
    ids = ids + from_gene_to_rxns(model, g)


rxnids = list(set(ids))
rxns = [model.reactions.get_by_id(id) for id in rxnids]

# if __name__ == '__main__':
#     freeze_support()
#     fva = flux_variability_analysis(model, rxns, fraction_of_optimum=0.95, loopless=True)
lable = pd.Series(index=rxnids)
lable.loc[:] = 'ideal target'
a = re.compile('and')
b = re.compile('or')
for r in rxns:
    score = get_rxn_score(r, targets)
    # Only "OR" in GPR for upregulation reactions.
    if np.min(score) > 1 and len(a.findall(r.gene_reaction_rule)) != 0:
        lable.loc[r.id] = 'Not ideal'
    # Only "AND" in GPR for downregulation reactions.
    if np.max(score) < 1 and len(b.findall(r.gene_reaction_rule)) != 0:
        lable.loc[r.id] = 'Not ideal'
    # exclude reactions related to formate, lactate and acetate
    byproduct = ['for_e', 'for_p', 'for_c', 'lac__D_e', 'lac__L_e', 'lac__L_p', 'ac_c', 'ac_e', 'ac_p']
    mets = list(r.metabolites.keys())
    metsid = [mets[i].id for i in range(len(mets))]
    bypdt = [p for p in metsid if p in byproduct]
    if len(bypdt) != 0:
        lable.loc[r.id] = 'Not ideal'

ideal_rxn = [ir for ir in rxnids if lable.loc[ir] == 'ideal target']
ideal_genes_tmp = []
for rr in ideal_rxn:
    ideal_genes_tmp = ideal_genes_tmp + from_rxn_to_genes(model, rr)

ideal_genes = [g for g in ideal_genes_tmp if g in genes]
ideal_genes = list(set(ideal_genes))

def zz_sub(str):
    import re
    sub_ = re.findall(re.compile(r'\[\'.*\'\]'), str)
    sub__ = re.findall(re.compile(r'\w.*\w'), sub_[0])
    return sub__

def get_sub(model, gene):
    reactions = [i for i in model.genes.get_by_id(gene).reactions]
    sub_all = [reactions[j].subsystem for j in range(len(reactions))]
    sub = [zz_sub(str)[0] for str in sub_all]
    sub = list(set(sub))
    return sub

sub = pd.DataFrame(index=ideal_genes, columns=['subsystem'])
for g in ideal_genes:
    if g != 'eamA' and g != 'ALAS':
        sub.loc[g, 'subsystem'] = get_sub(model, g)
    if g == 'eamA':
        sub.loc[g, 'subsystem'] = 'Transport'
    if g == 'ALAS' or g == 'b0369' or g == 'b3805' or g == 'b3804' or g == 'b0475' or g == 'b2436' or g == 'b3867' or g == 'b0475'or g == 'b3802' or g == 'b3850' or g == 'b3997':
        sub.loc[g, 'subsystem'] = 'Porphyrin metabolism'

sub.to_excel('../output/Final_targets.xlsx')
