import numpy as np


def read_metabolites(metabolites_path):
    met_dicts = []
    with open(metabolites_path) as f:
        f_lines = f.readlines()
        # num_metabolites = f_lines[0]
        for met_info in f_lines[1:]:
            met_splt = met_info.split('[')
            met_id = met_splt[0][1:]
            compartment = met_splt[1][0]
            met_dicts.append({'full_id': met_info[1:-2], 'met_id': met_id, 'compartment': compartment})
    return met_dicts


def read_reactions(reactions_path):
    rxn_dicts = []
    with open(reactions_path) as f:
        f_lines = f.readlines()
        for rxn_info in f_lines[1:]:
            rxn_str = rxn_info[:-1]
            if rxn_str == '** manually added':
                pass
            elif not rxn_str.isspace():
                rxn_dicts.append({'rxn_id': rxn_str[1:-1], 'type': 0, 'lower_bound': 0, 'upper_bound': 0})
    return rxn_dicts


def read_off_reactions(off_rxns_path):
    off_rxn_ids = []
    with open(off_rxns_path) as f:
        f_lines = f.readlines()
        for rxn_info in f_lines[1:]:
            rxn_str = rxn_info[:-1]
            if rxn_str[0] == '*':
                off_rxn_ids.append(rxn_str[2:-1])
            else:
                off_rxn_ids.append(rxn_str[1:-1])
    return off_rxn_ids


def read_uptake_reactions(uptake_rxns_path):
    uptake_rxn_ids = []
    with open(uptake_rxns_path) as f:
        f_lines = f.readlines()
        for rxn_info in f_lines[1:]:
            rxn_str = rxn_info[:-1]
            uptake_rxn_ids.append(rxn_str[1:-1])
    return uptake_rxn_ids


def read_stoichimetry_matrix(stoichimetry_matrix_path):
    stoic_elements = []
    with open(stoichimetry_matrix_path) as f:
        f_lines = f.readlines()
        for stoch_info in f_lines[1:]:
            if stoch_info[:16] == '* Manually added':
                pass
            elif not stoch_info.isspace():
                stoch_splt = stoch_info.split(' ')
                met_rxn_info = stoch_splt[0].split('.')
                stoch_coeff = float(stoch_splt[1])
                met_id = met_rxn_info[0][1:-1]
                rxn_id = met_rxn_info[1][1:-1]
                stoic_elements.append({'met_full_id': met_id, 'rxn_id': rxn_id, 'coeff': stoch_coeff})
    return stoic_elements


def find_reaction(rxn_id, reaction_dicts):
    for rxn_dict in reaction_dicts:
        if rxn_dict['rxn_id'] == rxn_id:
            return rxn_dict
    print("missed id is: " + rxn_id)
    raise Exception


def update_rxn_types(reaction_dicts, rxn_types_path):
    with open(rxn_types_path) as f:
        f_lines = f.readlines()
        for type_info in f_lines[1:]:
            type_str = type_info[:-1]
            if type_str == '** manually added':
                pass
            elif type_str and not type_str.isspace():
                type_splt = type_str.split(' ')
                rxn_id = type_splt[0][1:-1]
                rxn_dict = find_reaction(rxn_id, reaction_dicts)
                rxn_dict['type'] = int(type_splt[1])
    return reaction_dicts


def modify_bounds(reaction_dicts, updake_reactions_ids, off_reactions_ids):
    for rxn_dict in reaction_dicts:
        rxn_id = rxn_dict['rxn_id']
        rxn_type = rxn_dict['type']
        if rxn_type == 0:
            rxn_dict['lower_bound'] = 0
            rxn_dict['upper_bound'] = 1000
        elif rxn_type == 1:
            rxn_dict['lower_bound'] = -1000
            rxn_dict['upper_bound'] = 1000
        elif rxn_type == 2:
            rxn_dict['lower_bound'] = 0
            rxn_dict['upper_bound'] = 0
        elif rxn_type == 4:
            rxn_dict['lower_bound'] = 0
            rxn_dict['upper_bound'] = 1000
            
        if rxn_id in updake_reactions_ids:
            rxn_dict['lower_bound'] = -1000
        elif rxn_id == 'Ec_biomass_iAF1260_core_59p81M':
            rxn_dict['lower_bound'] = 0
            rxn_dict['upper_bound'] = 0
        elif rxn_id in off_reactions_ids:
            rxn_dict['lower_bound'] = 0
            rxn_dict['upper_bound'] = 0
        elif rxn_id == 'ATPM':
            rxn_dict['lower_bound'] = 8.39
            rxn_dict['upper_bound'] = 8.39
        elif rxn_id == 'EX_glc(e)':
            rxn_dict['lower_bound'] = -10;
        elif rxn_id == 'EX_o2(e)':
            rxn_dict['lower_bound'] = -20;
        elif rxn_id == 'LDH_D_f':
            rxn_dict['lower_bound'] = 0
            rxn_dict['upper_bound'] = 0
        elif rxn_id == 'LDH_D2':
            rxn_dict['lower_bound'] = 0
            rxn_dict['upper_bound'] = 0
        elif rxn_id == 'L-LACD2':
            rxn_dict['lower_bound'] = 0
            rxn_dict['upper_bound'] = 0
        elif rxn_id == 'L-LACD3':
            rxn_dict['lower_bound'] = 0
            rxn_dict['upper_bound'] = 0
        elif rxn_id == 'PFL':
            rxn_dict['lower_bound'] = 0
            rxn_dict['upper_bound'] = 0
    return reaction_dicts


def get_met_index_in_order(met_id, metabolites_dicts):
    index = 0
    for met_dict in metabolites_dicts:
        if met_dict['full_id'] == met_id:
            return index
        index += 1
    print("missed id: " + met_id)
    raise Exception


def get_rxn_index_in_order(rxn_id, reaction_dicts):
    index = 0
    for rxn_dict in reaction_dicts:
        if rxn_dict['rxn_id'] == rxn_id:
            return index
        index += 1
    print("missed id: " + rxn_id)
    raise Exception


def make_stoichiometry_matrix(stoic_elements, reaction_dicts, metabolites_dicts, _m, _n):
    stoic = np.zeros((_m, _n))
    for the_element in stoic_elements:
        met_id = get_met_index_in_order(the_element['met_full_id'], metabolites_dicts)
        rxn_id = get_rxn_index_in_order(the_element['rxn_id'], reaction_dicts)
        stoic[met_id, rxn_id] = the_element['coeff']
    return stoic


def make_bounds_vector(reaction_dicts, _n):
    lb = np.zeros(_n)
    ub = np.zeros(_n)
    index = 0
    for rxn_dict in reaction_dicts:
        # print(rxn_dict)
        lb[index] = rxn_dict['lower_bound']
        ub[index] = rxn_dict['upper_bound']
        index += 1
    return lb, ub


def get_model_metabolite(rxn_vec, mdl_metabolites, _m):
    rxn_metabolites = {}
    for i in range(_m):
        if rxn_vec[i] != 0:
            the_metabolite = mdl_metabolites[i]
            rxn_metabolites[the_metabolite] = rxn_vec[i]
    return rxn_metabolites
