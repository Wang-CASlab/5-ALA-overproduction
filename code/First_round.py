import cobra
from cobra.io import read_sbml_model
from cobra import Reaction, Metabolite
import numpy as np
import pandas as pd
import copy


def simulateGrowth(model, alpha):
    tmpmodel = model.copy()

    # max growth
    tmpmodel.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
    sol = tmpmodel.optimize()

    # fix growth and max product
    tmpgrow = sol.objective_value * 0.999 * alpha
    tmpmodel.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').bounds = (tmpgrow, tmpgrow)
    tmpmodel.objective = 'EX_ALA_e'
    heme_max = tmpmodel.optimize()
    print(heme_max.objective_value)

    # fix product and perform pFBA
    tmpmodel.reactions.get_by_id('EX_ALA_e').bounds = (heme_max.objective_value * 0.999,
                                                       heme_max.objective_value * 0.999)
    sol_pfba = cobra.flux_analysis.pfba(tmpmodel)
    return sol_pfba


def div_(start_, end_, num):
    l = []
    interv = (end_ - start_) / num
    co = 0
    while co < num + 1:
        l.append(start_ + co * interv)
        co += 1
    return l


def compare_substrate(model):
    flux_WT = simulateGrowth(model, 1)
    aex = flux_WT.fluxes.loc['BIOMASS_Ec_iML1515_core_75p37M']
    alpha = div_(0.2 * aex, aex * 0.8, 15)
    v_matrix = []
    k_matrix = []
    for a in alpha:
        tmpflux = simulateGrowth(model, a)
        try:
            v_matrix.append(tmpflux.fluxes)
            k_matrix.append(np.asarray(tmpflux.fluxes) /
                            np.asarray(flux_WT.fluxes))
        except AttributeError:
            v_matrix.append(tmpflux)
            k_matrix.append(tmpflux)
        print(a)
        print(alpha.index(a))
    return v_matrix, k_matrix, alpha


model = cobra.io.load_matlab_model('../model/iML1515.mat')

# add C4 rxns and mets
ala_e = Metabolite('5aop_e', compartment='e')
model.add_metabolites(ala_e)

ala_syth = Reaction('ala_C4')
ala_syth.add_metabolites({
    model.metabolites.get_by_id('succoa_c'): -1,
    model.metabolites.get_by_id('gly_c'): -1,
    model.metabolites.get_by_id('h_c'): -1,
    model.metabolites.get_by_id('5aop_c'): 1,
    model.metabolites.get_by_id('co2_c'): 1,
    model.metabolites.get_by_id('coa_c'): 1})
ala_syth.gene_reaction_rule = 'ALAS'
model.add_reactions([ala_syth])

ala_exchange = Reaction('ala_ex')
ala_exchange.bounds = (-1000, 1000)
ala_exchange.add_metabolites({
    model.metabolites.get_by_id('5aop_c'): -1,
    model.metabolites.get_by_id('5aop_e'): 1})
ala_exchange.gene_reaction_rule = 'eamA'
model.add_reactions([ala_exchange])

ala_transport = Reaction('EX_ALA_e')
ala_transport.bounds = (0, 1000)
ala_transport.add_metabolites({
    model.metabolites.get_by_id('5aop_e'): -1})
model.add_reactions([ala_transport])

model.objective = model.reactions.get_by_id(ala_transport.id)
solu = model.optimize()

medium = model.medium
medium['EX_gly_e'] = 1.0
tmpmodel = model.copy()
tmpmodel.medium = medium
model.objective = model.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M')
solu2 = model.optimize()

# Follow Olena P. Ishchuk's paper (https://doi.org/10.1073/pnas.2108245119)
v_matrix, k_matrix, alpha = compare_substrate(model)
all_flux = pd.DataFrame(v_matrix)
all_flux = all_flux.T
all_k = pd.DataFrame(k_matrix)
all_k = all_k.T
all_k.index = all_flux.index

f1 = list()
for h in all_flux.index:
    if bool(model.reactions.get_by_id(h).gene_reaction_rule == ''):
        f1.append(h)

all_flux.drop(index=f1, inplace=True)
all_k.drop(index=f1, inplace=True)
print()
print('delete all nan')
all_flux.dropna(axis=0, how='all', inplace=True)
all_k.dropna(axis=0, how='all', inplace=True)
print(len(all_k.index))

print('nan --> 1')
all_flux.fillna(1, inplace=True)
all_k.fillna(1, inplace=True)
print(len(all_k.index))
# inf --> 1000
all_k.replace([np.inf, -np.inf], 1000, inplace=True)

print('delete inconsistent rxn')
f5 = []
for i in all_k.index:
    big = any(all_k.loc[i, :] > 1)
    small = any(all_k.loc[i, :] < 1)
    if big == small == True:
        f5.append(i)
all_k.drop(index=f5, inplace=True)
print(len(all_k.index))

# order k(descend)
orderd_k = copy.deepcopy(all_k)
orderd_k['meank'] = orderd_k.mean(axis=1)
orderd_k = orderd_k.sort_values(by=["meank"], ascending=False)

print('get gene')
cons_g = {}
single_g = {}
for ge in model.genes:
    # 取基因的k值
    ge_rxn = [r.id for r in ge.reactions if r.id in orderd_k.index]
    ge_k = [orderd_k.loc[r.id, 'meank'] for r in ge.reactions if r.id in orderd_k.index]
    # delete inconsistent gene
    if ge_k != []:
        gk = np.asarray(ge_k)
        if all(gk > 1) == True or all(gk < 1) == True:
            cons_g[ge.id] = np.mean(ge_k)
print(len(cons_g))
print(len(single_g))

cons_g_f = {}
for k, v in cons_g.items():
    if (v < 0.002) or (v > 500):
        cons_g_f[k] = v
    rea = model.genes.get_by_id(k).reactions
    for r in rea:
        for mets in r.products + r.reactants:
            if mets.id == '5aop_c' or mets.id == '5aop_e':
                cons_g_f[k] = v
# 排序gene的k值
cons_g_f = pd.Series(cons_g_f)
cons_g_f.sort_values(ascending=False, inplace=True)
cons_g_f.to_excel('../output/First_round_targets.xlsx')
