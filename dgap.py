#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:27:09 2022

@author: omidard
"""

#diverse template gapfilling: dgap

#imports


import argparse
import json
import re
import zipfile
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import multiprocessing as mp
import numpy as np
import os
from functools import partial
from pathlib import Path


LEGACY_GAP_DIR = '/home/omidard/allgems/Limosilactobacillus_fermentum/biomassed'
LEGACY_TEMPLATE = 'GCF_009556455.1.json'


def resolve_model_file(model_ref, models_dir=None):
    candidate = Path(model_ref)
    if candidate.is_file():
        return str(candidate)

    if models_dir is None:
        raise FileNotFoundError(f"Could not resolve model path for {model_ref}")

    base = Path(models_dir)
    guesses = [
        base / model_ref,
        base / f"{model_ref}.json",
    ]
    for guess in guesses:
        if guess.is_file():
            return str(guess)

    matches = list(base.glob(f"{model_ref}*.json"))
    if len(matches) == 1:
        return str(matches[0])

    raise FileNotFoundError(f"Could not resolve model path for {model_ref} in {models_dir}")


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


BIOMASS2_STOICHIOMETRY = {
    # Matches the packaged post-gapfilled JSON models in this repository.
    "CPS_LBR_c": -0.078,
    "DNA_LBR_c": -0.205,
    "LIP_LBR_c": -0.106,
    "LTAtotal_LBR_c": -0.006,
    "PGlac2_c": -0.009,
    "PROT_LBR_c": -3.311,
    "RNA_LBR_c": -0.926,
    "atp_c": -27.2,
    "coa_c": -0.002,
    "h2o_c": -27.2,
    "nad_c": -0.002,
    "pydx5p_c": -0.000001,
    "udcpdp_c": -0.002,
    "adp_c": 27.2,
    "h_c": 27.2,
    "pi_c": 27.2,
}


M9_LOWER_BOUNDS = {
    "EX_arg_L_e": -2,
    "EX_cys_L_e": -2,
    "EX_glu_L_e": -2,
    "EX_ile_L_e": -2,
    "EX_leu_L_e": -2,
    "EX_met_L_e": -2,
    "EX_tyr_L_e": -2,
    "EX_phe_L_e": -2,
    "EX_thr_L_e": -2,
    "EX_val_L_e": -2,
    "EX_gly_e": -2,
    "EX_ala_L_e": -2,
    "EX_asp_L_e": -2,
    "EX_his_L_e": -2,
    "EX_lys_L_e": -2,
    "EX_pro_L_e": -2,
    "EX_ser_L_e": -2,
    "EX_trp_L_e": -2,
    "EX_glc_D_e": -15,
    "EX_ac_e": -1,
    "EX_cit_e": -1,
    "EX_thymd_e": -0.1,
    "EX_ura_e": -1,
    "EX_gua_e": -1,
    "EX_ins_e": -1,
    "EX_ade_e": -1,
    "EX_xan_e": -1,
    "EX_orot_e": -1,
    "EX_btn_e": -0.1,
    "EX_pnto_R_e": -0.5,
    "EX_thm_e": -0.1,
    "EX_pydam_e": -0.1,
    "EX_pydxn_e": -0.1,
    "EX_ribflv_e": -0.1,
    "EX_fol_e": -0.1,
    "EX_ascb_L_e": -0.5,
    "EX_4abz_e": -1,
    "EX_nac_e": -1,
    "EX_cl_e": -10,
    "EX_h_e": -1000,
    "EX_h2o_e": -10,
    "EX_nh4_e": -10,
    "EX_ca2_e": -10,
    "EX_co_e": -10,
    "EX_co2_e": -10,
    "EX_pi_e": -100,
    "EX_cobalt2_e": -10,
    "EX_cu2_e": -10,
    "EX_fe3_e": -10,
    "EX_k_e": -10,
    "EX_mn2_e": -10,
    "EX_so4_e": -10,
    "EX_na1_e": -10,
    "EX_mg2_e": -10,
    "EX_zn2_e": -10,
}


def has_id(items, item_id):
    return hasattr(items, "has_id") and items.has_id(item_id)


def clean_model_id(model, model_path=None):
    model_id = str(model.id or (Path(model_path).stem if model_path else "model"))
    while model_id.endswith(".json"):
        model_id = Path(model_id).stem
    model_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", model_id).strip("._")
    return model_id or "model"


def ensure_biomass2(model, replace=False):
    if has_id(model.reactions, "BIOMASS2"):
        reaction = model.reactions.get_by_id("BIOMASS2")
        if not replace:
            model.objective = "BIOMASS2"
            return reaction
        reaction.subtract_metabolites(reaction.metabolites.copy())
    else:
        reaction = Reaction("BIOMASS2")
        model.add_reactions([reaction])

    missing = [met_id for met_id in BIOMASS2_STOICHIOMETRY if not has_id(model.metabolites, met_id)]
    if missing:
        raise KeyError(f"Cannot build BIOMASS2; missing metabolites: {', '.join(missing)}")

    reaction.name = "Biomass production"
    reaction.subsystem = "Biomass"
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites(
        {
            model.metabolites.get_by_id(met_id): coefficient
            for met_id, coefficient in BIOMASS2_STOICHIOMETRY.items()
        }
    )
    reaction.gene_reaction_rule = "BIOMASS"
    model.objective = "BIOMASS2"
    model.repair()
    return reaction


def set_growth_objective(model, objective="BIOMASS2", replace_biomass2=False):
    if objective == "BIOMASS2":
        ensure_biomass2(model, replace=replace_biomass2)
    elif has_id(model.reactions, objective):
        model.objective = objective
    elif has_id(model.reactions, "BIOMASS2"):
        model.objective = "BIOMASS2"
    elif has_id(model.reactions, "BIOMASS"):
        model.objective = "BIOMASS"
    else:
        raise KeyError(f"Model {model.id} has no usable biomass objective.")
    return model


def optimize_growth(model):
    solution = model.optimize()
    if solution.status != "optimal" or solution.objective_value is None:
        return 0.0, solution.status
    return float(solution.objective_value), solution.status


def apply_open_exchanges(model, uptake=-1000, solver=None):
    if solver:
        model.solver = solver
    set_growth_objective(model)
    for reaction in model.reactions:
        if reaction.id.startswith("EX_"):
            reaction.lower_bound = min(reaction.lower_bound, uptake)
            reaction.upper_bound = max(reaction.upper_bound, 1000)
    if has_id(model.reactions, "ATPM"):
        model.reactions.get_by_id("ATPM").lower_bound = 1
        model.reactions.get_by_id("ATPM").upper_bound = 1
    return model


def gap_reaction_count(model):
    if has_id(model.genes, "GAP"):
        return len(model.genes.get_by_id("GAP").reactions)
    return 0




#formulation of a chemically defined media
def m9(model, solver=None, replace_biomass2=False):
    if solver:
        model.solver = solver
    set_growth_objective(model, "BIOMASS2", replace_biomass2=replace_biomass2)
    for reaction in model.reactions:
        if reaction.id.startswith("EX_"):
            reaction.lower_bound = 0

    for reaction_id, lower_bound in M9_LOWER_BOUNDS.items():
        if has_id(model.reactions, reaction_id):
            model.reactions.get_by_id(reaction_id).lower_bound = lower_bound

    for reaction in model.reactions:
        if "HCLTr" in reaction.id:
            reaction.lower_bound = -10
            reaction.upper_bound = 10
        if "PTRCt2" in reaction.id:
            reaction.lower_bound = -10

    if has_id(model.reactions, "ATPM"):
        model.reactions.get_by_id("ATPM").lower_bound = 1
        model.reactions.get_by_id("ATPM").upper_bound = 1
    if has_id(model.reactions, "EX_ptrc_e"):
        model.reactions.get_by_id("EX_ptrc_e").lower_bound = 0
        model.reactions.get_by_id("EX_ptrc_e").upper_bound = 0

    model.objective = "BIOMASS2"
    return model
    


#scanning models and retrive important information
def scan(directory, solver=None, growth_threshold=0.01):
    models =glob('%s/*.json'%directory)
    gr=[]
    status=[]
    mid=[]
    reactions_total = []
    reactions_gap = []
    failed=[] #failed models
    temps=[] #gapfilled mopdels
    for mod in models:
        model=load_json_model(mod)
        m9(model, solver=solver)
        growth, solution_status = optimize_growth(model)
        if growth >= growth_threshold:
            temps.append(model.id)
        else:
            failed.append(model.id)
        gr.append(growth)
        status.append(solution_status)
        mid.append(model.id)
        reactions_gap.append(gap_reaction_count(model))
        reactions_total.append(len(model.reactions))
    all_gems = pd.DataFrame()
    all_gems['id']=mid
    all_gems['growth']=gr
    all_gems['status']=status
    all_gems['total_reactions']=reactions_total
    all_gems['gapfilled_reactions']=reactions_gap
    return all_gems,failed,temps
 
        
 
#find best template
def tempfind(all_gems,failed,temps):
    tr=[]
    gapr=[]
    for i in range(len(all_gems)):
        if all_gems.id[i] in temps:
            tr.append(all_gems.total_reactions[i])
            gapr.append(all_gems.gapfilled_reactions[i])
    temps_inf = pd.DataFrame()
    temps_inf['id'] = temps
    temps_inf['total_reactions'] = tr
    temps_inf['gapfilled_reactions']=gapr
    temps_inf.sort_values(by='total_reactions', axis=0, ascending=False, inplace=True, kind='quicksort', na_position='last', ignore_index=True, key=None)
    return temps_inf



#find candidate reactions
def gaps(failed,temp_name, models_dir=LEGACY_GAP_DIR, template_path=None):
    temp_re_id=[] #list1
    if template_path is None:
        template_path = resolve_model_file(temp_name, models_dir)
    template = load_json_model(template_path)
    for reaction in template.reactions:
        temp_re_id.append(reaction.id)
    allmissing=[]
    for i in failed:
        missing=pd.DataFrame()
        model = load_json_model(resolve_model_file(i, models_dir))
        failed_re_id=[] #list two
        for reaction in model.reactions:
            failed_re_id.append(reaction.id)
        missing_reactions = list(set(temp_re_id).difference(failed_re_id))
        missing[model.id] = missing_reactions
        allmissing.append(missing)
    return allmissing



#alternate function development
def zx(allmissing):
    xlist=[] #potential gaps
    modelsz=[] #models
    for i in allmissing:
        for x in i:
            modelsz.append(x)
            zlist=[]
            for z in i[x]:
                zlist.append(z)
            xlist.append(zlist)
    return xlist,modelsz



#add candidate reactions and save models based on gapfilling status
def addgaps(
    allgapz,
    modelsz,
    temp_name,
    models_dir=LEGACY_GAP_DIR,
    feasible_dir=None,
    failed_dir=None,
    template_path=None,
    solver=None,
    growth_threshold=0.01,
):
    if feasible_dir is None:
        feasible_dir = os.path.join(models_dir, 'feasible')
    if failed_dir is None:
        failed_dir = os.path.join(models_dir, 'failed2')
    ensure_dir(feasible_dir)
    ensure_dir(failed_dir)
    if template_path is None:
        template_path = resolve_model_file(temp_name, models_dir)

    for i in range(len(modelsz)):
        model = load_json_model(resolve_model_file(modelsz[i], models_dir))
        template = load_json_model(template_path)
        for x in allgapz[i]:
            if x != 'no':
                reaction = template.reactions.get_by_id(x)
                reaction2 = Reaction(reaction.id)
                reaction2.name = reaction.name
                reaction2.subsystem = reaction.subsystem
                reaction2.lower_bound = reaction.lower_bound
                reaction2.upper_bound = reaction.upper_bound
                reaction2.add_metabolites(reaction.metabolites)
                reaction2.gene_reaction_rule = '(GAP)'
                model.add_reactions([reaction2])
                model.repair()
                print(reaction2.id)
                print('-----done')
                m9(model, solver=solver)
                print(model.optimize().fluxes.BIOMASS2)
        m9(model, solver=solver)
        out_path = os.path.join(
            feasible_dir if model.optimize().fluxes.BIOMASS2 > growth_threshold else failed_dir,
            f"{Path(model.id).stem}.json",
        )
        cobra.io.json.save_json_model(model, out_path)



def fluxanalyze(rea, template_path=None, solver=None, growth_threshold=0.1):
    if template_path is None:
        template_path = resolve_model_file(LEGACY_TEMPLATE, LEGACY_GAP_DIR)
    template = load_json_model(template_path)
    m9(template, solver=solver)
    template.reactions.get_by_id(rea).lower_bound = 0
    template.reactions.get_by_id(rea).upper_bound = 0
    if template.optimize().fluxes.BIOMASS2 < growth_threshold:
        tar = rea
    else:
        tar = 'no'
    return tar




def gap_extract(mod):
    model=load_json_model(mod)
    total_gaps =pd.DataFrame({model.id:[re.id for re in model.genes.GAP.reactions]})
    tg = total_gaps.T
    return tg

#find non essential gaps
def none_ess_gap(total_gaps,directory1,directory2):
    model = load_json_model(directory1+total_gaps.columns[0])
    print('in process',model.id)
    m9(model)
    for re in total_gaps[total_gaps.columns[0]]:
        try:
            m9(model)
            low = model.reactions.get_by_id(re).lower_bound
            up = model.reactions.get_by_id(re).upper_bound
            model.reactions.get_by_id(re).lower_bound = 0
            model.reactions.get_by_id(re).upper_bound = 0
            if model.optimize().fluxes.BIOMASS2 > 0.1:
                reaction = model.reactions.get_by_id(re)
                reaction.remove_from_model()
                print(re,'removed')
            else:
                model.reactions.get_by_id(re).lower_bound = low
                model.reactions.get_by_id(re).upper_bound = up
        except AttributeError:
            pass
    cobra.io.json.save_json_model(model,directory2+model.id)
    print(directory2,'finished')



def total_dir(rootdir):
    dirs = []
    for file in os.listdir(rootdir):
        local = os.path.join(rootdir, file)
        if os.path.isdir(local):
            dirs.append(local)
    return dirs




def exchange_reactions(directory):
    models =glob('%s/*.json'%directory)
    model = load_json_model(models[0])
    model = m9(model)
    exchanges = []
    for reaction in model.reactions:
        if 'EX_' in reaction.id and model.reactions.get_by_id(reaction.id).lower_bound < 0:
            exchanges.append(reaction.id)
    df = pd.DataFrame({'metabolites': exchanges})
    df2 = df.set_index(df['metabolites'])
    substrates = df2.drop('metabolites', axis=1)
    return substrates,exchanges



def essential_substrates(models,substrates):
    for mod in models:
        model = load_json_model(mod)
        growth_rates =[]
        for i in substrates[1]:
            m9(model)
            model.reactions.get_by_id(i).lower_bound = 0
            try:
                growth_rates.append(model.optimize().fluxes.BIOMASS2)
            except:
                growth_rates.append(0)
        col=model.id
        name=col.replace('.json','')
        substrates[0][name]=growth_rates
    return substrates



def exchanges_fluxes(model):
    ex_flux=[]
    ex_reaction=[]
    m9(model)
    all_flux = model.optimize().fluxes
    for i in range(len(all_flux)):
        if 'EX_' in all_flux.index[i]:
            ex_flux.append(all_flux[i])
            ex_reaction.append(all_flux.index[i])
    flux_frame = pd.DataFrame()
    flux_frame['exchanges']=ex_reaction
    flux_frame[model.id]=ex_flux
    flux_frame.sort_values(by=['exchanges'],inplace=True, ignore_index=True)
    flux_frame.set_index('exchanges',inplace=True)
    return flux_frame


def knockout_fluxes(model):
    ex_flux=[]
    ex_reaction=[]
    m9(model)
    all_flux = model.optimize().fluxes
    for i in range(len(all_flux)):
        if 'EX_' in all_flux.index[i]:
            ex_flux.append(all_flux[i])
            ex_reaction.append(all_flux.index[i])
    flux_frame = pd.DataFrame()
    flux_frame['exchanges']=ex_reaction
    flux_frame[model.id]=ex_flux
    flux_frame.sort_values(by=['exchanges'],inplace=True, ignore_index=True)
    flux_frame.set_index('exchanges',inplace=True)
    return flux_frame


def ess_re(mod):
    ess_collector=[]
    model=load_json_model(mod)
    for re in model.reactions:
        if 'EX_' not in re.id:
            model = load_json_model(mod)
            m9(model)
            model.reactions.get_by_id(re.id).lower_bound = 0
            model.reactions.get_by_id(re.id).upper_bound = 0
            try:
                if model.optimize().fluxes.BIOMASS2 <= 0:
                    ess_collector.append(re.id)
            except AttributeError:
                ess_collector.append(re.id)
    ess_collector2=pd.DataFrame()
    ess_collector2[model.id]=ess_collector
    return ess_collector2



def allinf (modelid):
    aux= ['EX_orot_e', 'EX_co_e', 'EX_pydam_e', 'EX_pydxn_e', 'EX_h_e', 'EX_h2o_e', 'EX_4abz_e', 'EX_ac_e', 'EX_ade_e', 'EX_ala_L_e', 'EX_arg_L_e', 'EX_asp_L_e', 'EX_btn_e', 'EX_ca2_e', 'EX_cit_e', 'EX_cl_e', 'EX_co2_e', 'EX_cys_L_e', 'EX_fol_e', 'EX_glc_D_e', 'EX_glu_L_e', 'EX_gly_e', 'EX_gua_e', 'EX_his_L_e', 'EX_ile_L_e', 'EX_ascb_L_e', 'EX_ins_e', 'EX_k_e', 'EX_leu_L_e', 'EX_lys_L_e', 'EX_met_L_e', 'EX_nh4_e', 'EX_phe_L_e', 'EX_pi_e', 'EX_pnto_R_e', 'EX_pro_L_e', 'EX_ribflv_e', 'EX_ser_L_e', 'EX_thm_e', 'EX_thr_L_e', 'EX_thymd_e', 'EX_trp_L_e', 'EX_tyr_L_e', 'EX_ura_e', 'EX_val_L_e', 'EX_xan_e']
    carbs=['EX_acgam_e', 'EX_gal_e', 'EX_2ddglcn_e', 'EX_cellb_e', 'EX_fru_e', 'EX_gam_e', 'EX_glc_D_e', 'EX_glcn_e', 'EX_lcts_e', 'EX_malt_e', 'EX_man_e', 'EX_raffin_e', 'EX_rib_D_e', 'EX_sbt_D_e', 'EX_sucr_e']
    exchanges = ['EX_ac_e', 'EX_actn_S_e', 'EX_ade_e', 'EX_adn_e', 'EX_akg_e', 'EX_ala_D_e', 'EX_ala_L_e', 'EX_amp_e', 'EX_arg_L_e', 'EX_ascb_L_e', 'EX_asp_L_e', 'EX_btn_e', 'EX_cit_e', 'EX_co2_e', 'EX_cys_L_e', 'EX_fol_e', 'EX_for_e', 'EX_fum_e', 'EX_gam_e', 'EX_glc_D_e', 'EX_glu_L_e', 'EX_gly_e', 'EX_glyc_e', 'EX_gua_e', 'EX_h2o_e', 'EX_h2s_e', 'EX_h_e', 'EX_hdca_e', 'EX_his_L_e', 'EX_hxan_e', 'EX_ile_L_e', 'EX_ins_e', 'EX_leu_L_e', 'EX_lys_L_e', 'EX_mal_L_e', 'EX_man_e', 'EX_met_L_e', 'EX_mnl_e', 'EX_nh4_e', 'EX_orn_e', 'EX_orot_e', 'EX_oxa_e', 'EX_phe_L_e', 'EX_pi_e', 'EX_pnto_R_e', 'EX_pro_L_e', 'EX_ptrc_e', 'EX_pydx_e', 'EX_pydxn_e', 'EX_pyr_e', 'EX_ribflv_e', 'EX_ser_L_e', 'EX_succ_e', 'EX_thm_e', 'EX_thr_L_e', 'EX_thymd_e', 'EX_trp_L_e', 'EX_tyr_L_e', 'EX_ura_e', 'EX_uri_e', 'EX_val_L_e', 'EX_xan_e', 'EX_lipoate_e', 'EX_co_e', 'EX_dad_5_e', 'EX_pydam_e', 'EX_12ppd_S_e', 'EX_acgam_e', 'EX_acnam_e', 'EX_fuc_L_e', 'EX_gal_e', 'EX_lac_L_e', 'EX_lald_L_e', 'EX_12dgr180_e', 'EX_13ppd_e', 'EX_2ddglcn_e', 'EX_2dmmq8_e', 'EX_2obut_e', 'EX_34dhpha_e', 'EX_34dhphe_e', 'EX_3mop_e', 'EX_4abut_e', 'EX_4abz_e', 'EX_4hbz_e', 'EX_5htrp_e', 'EX_5mthf_e', 'EX_acald_e', 'EX_adocbl_e', 'EX_ahcys_e', 'EX_alltn_e', 'EX_arab_L_e', 'EX_arbt_e', 'EX_arsenb_e', 'EX_asn_L_e', 'EX_btd_RR_e', 'EX_butso3_e', 'EX_C02528_e', 'EX_ca2_e', 'EX_cd2_e', 'EX_cellb_e', 'EX_cgly_e', 'EX_chol_e', 'EX_cholate_e', 'EX_chols_e', 'EX_cl_e', 'EX_cobalt2_e', 'EX_csn_e', 'EX_ctbt_e', 'EX_cu2_e', 'EX_cytd_e', 'EX_dad_2_e', 'EX_dcyt_e', 'EX_ddca_e', 'EX_dextrin_e', 'EX_dgsn_e', 'EX_diact_e', 'EX_din_e', 'EX_dopa_e', 'EX_dpcoa_e', 'EX_drib_e', 'EX_etha_e', 'EX_ethso3_e', 'EX_etoh_e', 'EX_fe3_e', 'EX_fecrm_e', 'EX_fru_e', 'EX_galt_e', 'EX_galur_e', 'EX_gcald_e', 'EX_gchola_e', 'EX_glcn_e', 'EX_glcur_e', 'EX_gln_L_e', 'EX_glyb_e', 'EX_glyclt_e', 'EX_gsn_e', 'EX_gthox_e', 'EX_gthrd_e', 'EX_h2_e', 'EX_hexs_e', 'EX_hg2_e', 'EX_hista_e', 'EX_ind3ac_e', 'EX_indole_e', 'EX_inost_e', 'EX_k_e', 'EX_lac_D_e', 'EX_lcts_e', 'EX_Lcyst_e', 'EX_malt_e', 'EX_malthx_e', 'EX_malttr_e', 'EX_melib_e', 'EX_met_D_e', 'EX_metsox_S_L_e', 'EX_mn2_e', 'EX_mops_e', 'EX_mqn8_e', 'EX_mso3_e', 'EX_n2o_e', 'EX_nac_e', 'EX_ncam_e', 'EX_ni2_e', 'EX_nmn_e', 'EX_no_e', 'EX_no2_e', 'EX_no3_e', 'EX_o2_e', 'EX_ocdcea_e', 'EX_pb_e', 'EX_pheme_e', 'EX_ppa_e', 'EX_ppi_e', 'EX_q8_e', 'EX_raffin_e', 'EX_rib_D_e', 'EX_rmn_e', 'EX_salcn_e', 'EX_sbt_D_e', 'EX_ser_D_e', 'EX_sheme_e', 'EX_so4_e', 'EX_spmd_e', 'EX_srtn_e', 'EX_sucr_e', 'EX_sulfac_e', 'EX_taur_e', 'EX_tchola_e', 'EX_tdchola_e', 'EX_thf_e', 'EX_tre_e', 'EX_trypta_e', 'EX_tsul_e', 'EX_ttdca_e', 'EX_tym_e', 'EX_urea_e', 'EX_xyl_D_e', 'EX_na1_e', 'EX_alaasp_e', 'EX_alagln_e', 'EX_alaglu_e', 'EX_alagly_e', 'EX_alahis_e', 'EX_alaleu_e', 'EX_alathr_e', 'EX_crn_e', 'EX_glyasn_e', 'EX_glyasp_e', 'EX_glycys_e', 'EX_glygln_e', 'EX_glyglu_e', 'EX_glyleu_e', 'EX_glymet_e', 'EX_glyphe_e', 'EX_glypro_e', 'EX_glytyr_e', 'EX_isetac_e', 'EX_mantr_e', 'EX_metsox_R_L_e', 'EX_mg2_e', 'EX_Ser_Thr_e', 'EX_stys_e', 'EX_zn2_e', 'EX_dha_e', 'EX_12ppd_R_e', 'EX_4ahmmp_e']
    auxf = []
    carbsf=[]
    carbsp = []
    exf =[]
    wgr = []
    names = []
    ids = modelid.replace('/home/omidard/gems/allgems/','')
    model = load_json_model(modelid)
    m9(model)
    fx = model.optimize().fluxes
    wgr.append(model.optimize().fluxes.BIOMASS2)
    for e in exchanges:
        exf.append(fx[e])
        
    for c in carbs:
        model = load_json_model(modelid)
        m9(model)
        model.reactions.EX_glc_D_e.lower_bound = 0
        model.reactions.get_by_id(c).lower_bound = -15
        fxc = model.optimize().fluxes
        carbsf.append(model.optimize().fluxes.BIOMASS2)
        for e in exchanges:
            carbsp.append(fxc[e])
            name = e+'_'+'growth_on_'+c
            names.append(name)
            
    for a in aux:
        model = load_json_model(modelid)
        m9(model)
        model.reactions.get_by_id(a).lower_bound = 0
        auxf.append(model.optimize().fluxes.BIOMASS2)   
    
    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    df3 = pd.DataFrame()
    df4 = pd.DataFrame()
    #df1 setup
    df1.index = ['growth_on_'+c for c in carbs]
    df1[ids.replace('.json','')] = carbsf
    df1 = df1.T
    #df2
    df2.index = ['eliminated_'+a for a in aux]
    df2[ids.replace('.json','')] = auxf
    df2 = df2.T
    #df3
    df3.index = [n for n in names]
    df3[ids.replace('.json','')] = carbsp
    df3 = df3.T
    #df4
    df4.index = [e for e in exchanges]
    df4[ids.replace('.json','')] = exf
    df4 = df4.T
    #sumup
    df = pd.concat([df1,df2,df3,df4],axis = 1)
    df['wild_type_growth_rate'] = wgr
    return df


def all_reactions(modelid):
    reactions_list=[]
    model = load_json_model(modelid)
    for re in model.reactions:
        reactions_list.append(re.id)
    df = pd.DataFrame()
    df['all_reactions'] = reactions_list
    return df




def gene_associated_reactions (modelid):
    model = load_json_model(modelid)
    nonmetgenes=[]
    metgenes=[]
    model = load_json_model(modelid)
    for re in model.genes.GAP.reactions:
        nonmetgenes.append(re.id)
    for re in model.genes.EXCHANGE.reactions:
        nonmetgenes.append(re.id)
    for re in model.genes.ORPHAN.reactions:
        nonmetgenes.append(re.id)
    for re in model.genes.DEMAND.reactions:
        nonmetgenes.append(re.id)
    for re in model.genes.BIOMASS.reactions:
        nonmetgenes.append(re.id)
    for re in model.genes.spontaneous.reactions:
        nonmetgenes.append(re.id)
    for re in model.reactions:
        if re.id not in nonmetgenes:
            if 'SINK' not in re.id:
                metgenes.append(re.id)
    
    df = pd.DataFrame()
    df[model.id]=metgenes
    df.set_index(model.id,inplace=True)
    df[model.id]=[1 for n in range(len(df.index))]
    
    df2 = pd.DataFrame()
    df2[model.id]=nonmetgenes
    df2.set_index(model.id,inplace=True)
    df2[model.id]=[0 for n in range(len(df2.index))]
    
    df3=pd.concat([df,df2],axis=0)
    
    return df3

    

def basic_counts(modelid):
    reactions=[]
    genes=[]
    gaps=[]
    model = load_json_model(modelid)
    reactions.append(len(model.reactions))
    genes.append(len(model.genes))
    gaps.append(len(model.genes.GAP.reactions))
    df = pd.DataFrame()
    df['reactions']=reactions
    df['genes']=genes
    df['gaps']=gaps
    df['id']=[model.id]
    df.set_index('id',inplace=True)
    return df


def all_fluxes(modelid):
    model = load_json_model(modelid)
    m9(model)
    reaction = []
    fx = model.optimize().fluxes
    for i in range(len(fx)):
        if fx[i] != 0:
            reaction.append(fx.index[i])
            
    dfx = pd.DataFrame()
    dfx['reactions']=reaction
    dfx[model.id] = [1 for v in range(len(dfx))]
    dfx.set_index('reactions',inplace=True)
    return dfx





def product_associated_genes(targetmet):
    allflux=pd.read_csv('/home/omidard/allflux.csv')
    allflux.set_index('reactions',inplace=True)
    producers=[]
    for c in allflux.columns:
        if allflux[c][targetmet] == 1:
            producers.append(allflux[c])
    allflux_producers_df = pd.concat(producers,axis=1)
    allflux_producers_df.fillna(0, inplace=True)
    return allflux_producers_df
            


def product_essential_genes(mod):
    allflux_producers_df=pd.read_csv('/home/omidard/allflux_producers_df.csv')
    allflux_producers_df.set_index('reactions',inplace=True)
    inv_genes_collect = []
    for i in allflux_producers_df.index:
        if allflux_producers_df[mod][i] == 1:
            df = pd.DataFrame()
            model = load_json_model('/home/omidard/gems/allgems/'+mod)
            m9(model)
            df['id'] = [model.id]
            df.set_index('id',inplace=True)
            df['gene']=[i]
            df['presence'] = [model.optimize().fluxes.EX_mnl_e]
            model.reactions.get_by_id(i).lower_bound = 0
            model.reactions.get_by_id(i).upper_bound = 0
            if model.optimize().fluxes.BIOMASS2 != 0:
                df['absence'] = [model.optimize().fluxes.EX_mnl_e]
            if model.optimize().fluxes.BIOMASS2 == 0:
                df['absence'] = ['essential gene']
            inv_genes_collect.append(df)
    gene_product_effect = pd.concat(inv_genes_collect,axis = 0)
    return gene_product_effect



            
def product_based_essential_genes(mod):
    targetmet = 'EX_4abut_e'  #change target met
    model = load_json_model(mod)
    m9(model)
    fluxed_reactions = []
    fx = model.optimize().fluxes
    if fx['EX_4abut_e'] !=0: #change target met
        for i in range(len(fx)):
            if fx[i] != 0:
                fluxed_reactions.append(fx.index[i])
    coll=[]
    if len(fluxed_reactions) > 0:
        for i in fluxed_reactions:
            model = load_json_model(mod)
            m9(model)
            df = pd.DataFrame()
            df['gene'] = [i]
            df.set_index('gene',inplace=True)
            df['wildtype'+model.id] = [model.optimize().fluxes.EX_4abut_e] #change target met
            model.reactions.get_by_id(i).lower_bound = 0
            model.reactions.get_by_id(i).upper_bound = 0
            if model.optimize().fluxes.BIOMASS2 !=0:
                df['knocked_out'+model.id] = [model.optimize().fluxes.EX_4abut_e] #change target met
            if model.optimize().fluxes.BIOMASS2 ==0:
                df['knocked_out'+model.id] = ['essential']
            coll.append(df)
    if len(fluxed_reactions) == 0:
        df = pd.DataFrame()
        coll.append(df)
    genes_effect = pd.concat(coll,axis=0)
    
    return genes_effect
                
                

                
def coreflux(mod):
    model = load_json_model(mod)
    m9(model)
    reaction=[]
    flux=[]
    flx = model.optimize().fluxes
    for i in flx.index:
        if flx[i] !=0:
            reaction.append(i)
            flux.append(flx[i])
    df = pd.DataFrame()
    df['reactions'] = reaction
    df[model.id] = flux
    df.set_index('reactions',inplace=True)
    return df       
    
        
        
def biomass2(model, out_path=None, replace=True):
    if isinstance(model, (str, Path)):
        model = load_json_model(str(model))
    ensure_biomass2(model, replace=replace)
    if out_path:
        cobra.io.json.save_json_model(model, str(out_path))
    return model
        
        
        
        
        
def producers(mod):
    producers=[]
    model = load_json_model(mod)
    m9(model)
    if model.optimize().fluxes.EX_mnl_e != 0: #targetmet
        producers.append(model.id)
    return producers

def pro_act_ge(mod):
    reactions=[]
    model = load_json_model('/home/omidard/gems/allgems/'+mod)
    m9(model)
    fx = model.optimize().fluxes
    for i in fx.index:
        if fx[i] !=0:
            reactions.append(i)
    df = pd.DataFrame()
    df[model.id] = reactions
    return df

def pro_ess_ge(df):
    non_ess = []
    for i in df[df.columns[0]]:
        model = load_json_model('/home/omidard/gems/allgems/'+df.columns[0])
        m9(model)
        model.reactions.get_by_id(i).lower_bound=0
        model.reactions.get_by_id(i).upper_bound=0
        if model.optimize().fluxes.BIOMASS2 != 0:
            non_ess.append(i)
    noess = pd.DataFrame()
    noess['index'] = non_ess
    noess.set_index('index',inplace=True)
    noess[df.columns[0]] = [1000 for v in range(len(noess))]
    return noess

def pro_ess_eff(df):
    flxc=[]
    for i in df.index:
        model = load_json_model('/home/omidard/gems/allgems/'+df.columns[0])
        m9(model)
        model.reactions.get_by_id(i).lower_bound=0
        model.reactions.get_by_id(i).upper_bound=0
        flx = model.optimize().fluxes.EX_mnl_e #targetmet
        flxc.append(flx)
    df[df.columns[0]] = flxc
    return df


def resolve_model_inputs(model_refs):
    if isinstance(model_refs, (str, Path)):
        model_refs = [model_refs]

    resolved = []
    for model_ref in model_refs:
        ref = str(model_ref)
        path = Path(ref)
        if path.is_dir():
            matches = sorted(path.glob("*.json"))
        elif path.is_file():
            matches = [path]
        else:
            matches = [Path(match) for match in sorted(glob(ref))]

        for match in matches:
            if match.is_file() and match.suffix == ".json":
                resolved.append(match)

    unique = []
    seen = set()
    for path in resolved:
        key = str(path.resolve())
        if key not in seen:
            unique.append(path)
            seen.add(key)

    if not unique:
        raise FileNotFoundError(f"No JSON model files matched: {model_refs}")
    return unique


def is_gap_gene_rule(gene_reaction_rule):
    tokens = re.findall(r"[A-Za-z0-9_.-]+", str(gene_reaction_rule))
    return "GAP" in tokens


def collect_packaged_gap_ids(zip_patterns=None, base_dir="."):
    zip_patterns = zip_patterns or ["*_json.zip"]
    base_dir = Path(base_dir)
    zip_paths = []
    for pattern in zip_patterns:
        pattern_path = Path(pattern)
        if pattern_path.is_absolute():
            matches = glob(str(pattern_path))
        else:
            matches = glob(str(base_dir / pattern))
        zip_paths.extend(Path(match) for match in sorted(matches))

    gap_ids = set()
    member_count = 0
    for zip_path in sorted(set(zip_paths)):
        if not zip_path.is_file():
            continue
        with zipfile.ZipFile(zip_path) as archive:
            for member in archive.namelist():
                if not member.endswith(".json"):
                    continue
                member_count += 1
                data = json.loads(archive.read(member))
                for reaction in data.get("reactions", []):
                    if is_gap_gene_rule(reaction.get("gene_reaction_rule", "")):
                        gap_ids.add(reaction["id"])

    return gap_ids, {"zip_count": len(set(zip_paths)), "json_model_count": member_count}


def eligible_template_reaction(reaction, include_exchanges=False):
    if reaction.id in {"BIOMASS", "BIOMASS2"}:
        return False
    if not include_exchanges and reaction.id.startswith("EX_"):
        return False
    return True


def copy_template_reaction_as_gap(template, reaction_id, source):
    reaction = template.reactions.get_by_id(reaction_id).copy()
    reaction.gene_reaction_rule = "GAP"
    reaction.notes = dict(reaction.notes)
    reaction.notes["gapfill_source"] = source
    return reaction


def add_candidate_reactions(model, template, candidate_ids, source, added_sources):
    added = []
    for reaction in template.reactions:
        reaction_id = reaction.id
        if reaction_id not in candidate_ids or has_id(model.reactions, reaction_id):
            continue
        model.add_reactions([copy_template_reaction_as_gap(template, reaction_id, source)])
        added_sources[reaction_id] = source
        added.append(reaction_id)
    if added:
        model.repair()
    return added


def prune_added_reactions(model, template, added_sources, growth_threshold):
    kept = []
    pruned = []
    template_order = [reaction.id for reaction in template.reactions if reaction.id in added_sources]

    for reaction_id in template_order:
        if not has_id(model.reactions, reaction_id):
            continue
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.remove_from_model(remove_orphans=False)
        model.repair()
        growth, status = optimize_growth(model)
        if growth >= growth_threshold:
            pruned.append(
                {
                    "reaction_id": reaction_id,
                    "reaction_source": added_sources[reaction_id],
                    "growth_after_removal": growth,
                    "status_after_removal": status,
                }
            )
            continue

        model.add_reactions(
            [copy_template_reaction_as_gap(template, reaction_id, added_sources[reaction_id])]
        )
        model.repair()
        kept.append(reaction_id)

    return kept, pruned


def evaluate_model_media(model, solver=None, open_exchange_uptake=-1000):
    m9_model = model.copy()
    m9(m9_model, solver=solver)
    m9_growth, m9_status = optimize_growth(m9_model)

    open_model = model.copy()
    apply_open_exchanges(open_model, uptake=open_exchange_uptake, solver=solver)
    open_growth, open_status = optimize_growth(open_model)

    return m9_growth, m9_status, open_growth, open_status


def complete_single_model(
    model_path,
    template,
    package_gap_ids,
    out_model_path,
    candidate_mode="auto",
    prune=True,
    solver=None,
    growth_threshold=0.01,
    open_exchange_uptake=-1000,
    include_exchanges=False,
):
    model = load_json_model(str(model_path))
    had_biomass = has_id(model.reactions, "BIOMASS")
    had_biomass2 = has_id(model.reactions, "BIOMASS2")
    original_reactions = len(model.reactions)
    original_gap_reactions = gap_reaction_count(model)
    model.id = clean_model_id(model, model_path)

    m9(model, solver=solver, replace_biomass2=True)
    initial_m9_growth, initial_m9_status = optimize_growth(model)
    _, _, initial_open_growth, initial_open_status = evaluate_model_media(
        model, solver=solver, open_exchange_uptake=open_exchange_uptake
    )

    template_ids = {
        reaction.id
        for reaction in template.reactions
        if eligible_template_reaction(reaction, include_exchanges=include_exchanges)
    }
    package_candidates = {
        reaction_id
        for reaction_id in package_gap_ids
        if reaction_id in template_ids and not has_id(model.reactions, reaction_id)
    }

    candidate_rows = []
    selected_rows = []
    pruned_rows = []
    added_sources = {}

    growth = initial_m9_growth
    package_added = []
    reference_added = []
    used_reference_fallback = False

    if growth < growth_threshold and candidate_mode in {"auto", "packaged-gaps"}:
        package_added = add_candidate_reactions(
            model, template, package_candidates, "packaged_gap", added_sources
        )
        candidate_rows.extend(
            {
                "reaction_id": reaction_id,
                "reaction_source": "packaged_gap",
                "added": True,
            }
            for reaction_id in package_added
        )
        m9(model, solver=solver)
        growth, _ = optimize_growth(model)

    if growth < growth_threshold and candidate_mode in {"auto", "reference-missing"}:
        used_reference_fallback = True
        reference_candidates = {
            reaction.id
            for reaction in template.reactions
            if reaction.id in template_ids and not has_id(model.reactions, reaction.id)
        }
        reference_added = add_candidate_reactions(
            model, template, reference_candidates, "reference_missing", added_sources
        )
        candidate_rows.extend(
            {
                "reaction_id": reaction_id,
                "reaction_source": "reference_missing",
                "added": True,
            }
            for reaction_id in reference_added
        )
        m9(model, solver=solver)
        growth, _ = optimize_growth(model)

    if growth >= growth_threshold and prune and added_sources:
        kept_ids, pruned_rows = prune_added_reactions(model, template, added_sources, growth_threshold)
        m9(model, solver=solver)
        growth, _ = optimize_growth(model)
    else:
        kept_ids = [reaction.id for reaction in template.reactions if reaction.id in added_sources]

    selected_rows.extend(
        {
            "reaction_id": reaction_id,
            "reaction_source": added_sources[reaction_id],
        }
        for reaction_id in kept_ids
        if has_id(model.reactions, reaction_id)
    )

    final_m9_growth, final_m9_status, final_open_growth, final_open_status = evaluate_model_media(
        model, solver=solver, open_exchange_uptake=open_exchange_uptake
    )
    m9(model, solver=solver)
    cobra.io.json.save_json_model(model, str(out_model_path), pretty=False)

    summary = {
        "source_path": str(model_path),
        "model_id": model.id,
        "output_path": str(out_model_path),
        "had_biomass": had_biomass,
        "had_biomass2": had_biomass2,
        "original_reactions": original_reactions,
        "original_gap_reactions": original_gap_reactions,
        "initial_m9_growth": initial_m9_growth,
        "initial_m9_status": initial_m9_status,
        "initial_open_growth": initial_open_growth,
        "initial_open_status": initial_open_status,
        "package_candidates_added": len(package_added),
        "reference_candidates_added": len(reference_added),
        "used_reference_fallback": used_reference_fallback,
        "pruned_reactions": len(pruned_rows),
        "selected_gap_reactions": len(selected_rows),
        "final_reactions": len(model.reactions),
        "final_gap_reactions": gap_reaction_count(model),
        "final_m9_growth": final_m9_growth,
        "final_m9_status": final_m9_status,
        "final_open_growth": final_open_growth,
        "final_open_status": final_open_status,
        "completed": final_m9_growth >= growth_threshold,
    }

    return summary, candidate_rows, selected_rows, pruned_rows


def run_gapfill_pipeline(
    models,
    out_dir,
    template_path="LBReactome.json",
    gap_zip_globs=None,
    candidate_mode="auto",
    prune=True,
    solver=None,
    threads=None,
    growth_threshold=0.01,
    open_exchange_uptake=-1000,
    include_exchanges=False,
):
    if threads:
        print("Note: --threads is accepted for compatibility but this completion workflow is serial.")

    model_paths = resolve_model_inputs(models)
    out_dir = Path(out_dir)
    ensure_dir(out_dir)
    completed_dir = Path(ensure_dir(out_dir / "completed_models"))
    failed_dir = Path(ensure_dir(out_dir / "failed_models"))

    resolved_template = resolve_model_file(template_path)
    template = load_json_model(resolved_template)
    m9(template, solver=solver, replace_biomass2=True)
    template_growth, template_status = optimize_growth(template)
    if template_growth < growth_threshold:
        raise ValueError(
            f"Template {resolved_template} does not grow under m9 "
            f"(growth={template_growth}, status={template_status})."
        )

    package_gap_ids, package_stats = collect_packaged_gap_ids(
        zip_patterns=gap_zip_globs,
        base_dir=Path.cwd(),
    )
    if candidate_mode == "packaged-gaps" and not package_gap_ids:
        raise ValueError("No packaged GAP reactions found. Use --candidate-mode reference-missing.")

    template_reaction_ids = {reaction.id for reaction in template.reactions}
    candidate_df = pd.DataFrame(
        {
            "reaction_id": sorted(package_gap_ids),
            "reaction_source": "packaged_gap",
            "present_in_template": [reaction_id in template_reaction_ids for reaction_id in sorted(package_gap_ids)],
        }
    )
    candidate_df.to_csv(out_dir / "packaged_gap_candidates.csv", index=False)

    summaries = []
    all_candidates = []
    all_selected = []
    all_pruned = []
    used_output_names = set()

    for model_path in model_paths:
        preview = load_json_model(str(model_path))
        model_id = clean_model_id(preview, model_path)
        output_name = f"{model_id}.json"
        index = 2
        while output_name in used_output_names:
            output_name = f"{model_id}_{index}.json"
            index += 1
        used_output_names.add(output_name)
        out_model_path = completed_dir / output_name

        print(f"Completing {model_path} -> {out_model_path}")
        summary, candidate_rows, selected_rows, pruned_rows = complete_single_model(
            model_path=model_path,
            template=template,
            package_gap_ids=package_gap_ids,
            out_model_path=out_model_path,
            candidate_mode=candidate_mode,
            prune=prune,
            solver=solver,
            growth_threshold=growth_threshold,
            open_exchange_uptake=open_exchange_uptake,
            include_exchanges=include_exchanges,
        )

        if not summary["completed"]:
            failed_path = failed_dir / output_name
            Path(summary["output_path"]).replace(failed_path)
            summary["output_path"] = str(failed_path)

        summaries.append(summary)
        for row in candidate_rows:
            row.update({"model_id": summary["model_id"], "source_path": str(model_path)})
            all_candidates.append(row)
        for row in selected_rows:
            row.update({"model_id": summary["model_id"], "source_path": str(model_path)})
            all_selected.append(row)
        for row in pruned_rows:
            row.update({"model_id": summary["model_id"], "source_path": str(model_path)})
            all_pruned.append(row)

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(out_dir / "completion_summary.csv", index=False)
    summary_df[
        [
            "model_id",
            "source_path",
            "output_path",
            "initial_m9_growth",
            "initial_open_growth",
            "final_m9_growth",
            "final_open_growth",
            "completed",
        ]
    ].to_csv(out_dir / "scan_summary.csv", index=False)
    pd.DataFrame(all_candidates).to_csv(out_dir / "candidate_gap_reactions.csv", index=False)
    pd.DataFrame(all_selected).to_csv(out_dir / "selected_gap_reactions.csv", index=False)
    pd.DataFrame(all_pruned).to_csv(out_dir / "pruned_gap_reactions.csv", index=False)
    summary_df[
        [
            "model_id",
            "output_path",
            "final_m9_growth",
            "final_m9_status",
            "final_open_growth",
            "final_open_status",
            "completed",
        ]
    ].to_csv(out_dir / "verification_summary.csv", index=False)

    return {
        "summary": summary_df,
        "template_path": resolved_template,
        "template_growth": template_growth,
        "package_stats": package_stats,
        "completed_dir": str(completed_dir),
        "failed_dir": str(failed_dir),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Complete draft GEM JSON models by adding BIOMASS2 and pruned GAP reactions."
    )
    parser.add_argument(
        "--models",
        nargs="+",
        required=True,
        help="Input JSON files, directories containing JSON models, or glob patterns.",
    )
    parser.add_argument("--out", required=True, help="Output directory for completed models and reports.")
    parser.add_argument(
        "--template",
        default="LBReactome.json",
        help="Reference/template model to copy candidate reactions from.",
    )
    parser.add_argument(
        "--gap-zip-glob",
        action="append",
        default=None,
        help="Glob for packaged JSON ZIPs used to seed GAP candidates. Defaults to *_json.zip.",
    )
    parser.add_argument(
        "--candidate-mode",
        choices=["auto", "packaged-gaps", "reference-missing"],
        default="auto",
        help="Candidate source strategy. auto tries packaged GAP reactions, then missing template reactions if needed.",
    )
    parser.add_argument("--solver", default=None, help="Optional COBRA solver name, for example glpk.")
    parser.add_argument("--threads", type=int, default=None, help="Accepted for compatibility; workflow is serial.")
    parser.add_argument(
        "--growth-threshold",
        type=float,
        default=0.01,
        help="Minimum BIOMASS2 flux for a model to count as completed.",
    )
    parser.add_argument(
        "--open-exchange-uptake",
        type=float,
        default=-1000,
        help="Lower bound applied to exchange reactions for open-exchange verification.",
    )
    parser.add_argument("--no-prune", action="store_true", help="Keep all added candidates without pruning.")
    parser.add_argument(
        "--include-exchanges",
        action="store_true",
        help="Allow exchange reactions to be copied from the template as gapfill candidates.",
    )
    args = parser.parse_args()

    results = run_gapfill_pipeline(
        models=args.models,
        out_dir=args.out,
        template_path=args.template,
        gap_zip_globs=args.gap_zip_glob,
        candidate_mode=args.candidate_mode,
        prune=not args.no_prune,
        solver=args.solver,
        threads=args.threads,
        growth_threshold=args.growth_threshold,
        open_exchange_uptake=args.open_exchange_uptake,
        include_exchanges=args.include_exchanges,
    )
    completed = int(results["summary"]["completed"].sum())
    total = len(results["summary"])
    print(f"Template path: {results['template_path']}")
    print(f"Template m9 growth: {results['template_growth']}")
    print(
        "Packaged GAP scan: "
        f"{results['package_stats']['json_model_count']} JSON models in "
        f"{results['package_stats']['zip_count']} ZIP files"
    )
    print(f"Completed models: {completed}/{total}")
    print(f"Completed model directory: {results['completed_dir']}")
    print(f"Reports written to: {Path(args.out)}")


if __name__ == "__main__":
    main()
