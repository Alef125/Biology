import numpy as np
import pandas as pd
from scipy.optimize import linprog


class Reaction:
    def __init__(self, description):
        self.description = description
        self.substrates, self.products, self.enzyme = parse_reaction(self.description)  # Todo: number of metabolites
        if self.substrates == ['']:
            self.substrates = []
        if self.products == ['']:
            self.products = []


def get_metabolites():
    metabolites_file = open('Metabolites.txt', 'r')
    metabolites_list = metabolites_file.read().splitlines()
    return metabolites_list


def parse_reaction(the_reaction):
    reaction_eq, enzyme_text = the_reaction.rsplit(" | ")
    left, right = reaction_eq.rsplit(" —> ")
    substrates = left.rsplit(" + ")
    products = right.rsplit(" + ")
    enzyme = enzyme_text[8:]  # Omits "Enzyme: "
    return substrates, products, enzyme


def get_reactions():
    all_reactions = []
    reactions_file = open('Reactions.txt', 'r')
    reactions_list = reactions_file.readlines()
    for reaction_description in reactions_list:
        all_reactions.append(Reaction(reaction_description.rstrip()))
    return all_reactions


def get_enzymes_name(all_reactions):
    enzymes_name = []
    for a_reaction in all_reactions:
        enzymes_name.append(a_reaction.enzyme)
    return enzymes_name


def make_stoichiometry_matrix(all_metabolites, all_reactions):
    m = len(all_metabolites)
    n = len(all_reactions)
    stoichiometry_matrix = np.zeros(shape=(m, n))
    for j in range(n):
        a_reaction = all_reactions[j]
        for substrate in a_reaction.substrates:
            i = all_metabolites.index(substrate)
            stoichiometry_matrix[i][j] -= 1
        for product in a_reaction.products:
            i = all_metabolites.index(product)
            stoichiometry_matrix[i][j] += 1
    stoichiometry_dataframe = pd.DataFrame(data=stoichiometry_matrix,
                                           columns=get_enzymes_name(all_reactions))
    stoichiometry_dataframe.index = all_metabolites
    return stoichiometry_dataframe


def get_objective_function(num_reactions, atp_idx):
    # for "minimize c.x" problem
    c = [0] * num_reactions
    c[atp_idx] = -1
    return c


def get_bounds(num_reactions, glucose_entrance_idx: int):
    bounds_pre = [(0.0, None) for _ in range(glucose_entrance_idx-1)]
    bounds_glucose = [(0.9, 1.0)]
    bounds_post = [(0.0, None) for _ in range(num_reactions-glucose_entrance_idx)]
    bounds = bounds_pre + bounds_glucose + bounds_post
    return bounds


def optimize_stoichiometry(stoichiometry_dataframe, glc_idx, atp_idx):
    m, n = stoichiometry_dataframe.shape  # m= num_metabolites, n= num_reactions
    b = [0] * m
    c = get_objective_function(n, atp_idx=atp_idx)
    bounds = get_bounds(n, glucose_entrance_idx=glc_idx+1)
    res = linprog(c=c,
                  A_eq=stoichiometry_dataframe.to_numpy(),
                  b_eq=b,
                  bounds=bounds)
    print('ATP produced:  ' + str(-round(res['fun'], 1)))
    return res['x']


def print_optimal_path(reactions_vector, all_reactions):
    m = len(all_reactions)
    for j in range(m):
        if reactions_vector[j] > 0.1:
            a_reaction = all_reactions[j]
            reaction_flux = round(reactions_vector[j], 1)
            print('(' + str(reaction_flux) + ') ' + a_reaction.description)


def get_glucose_and_atp_reactions_idx(all_reactions):
    num_reactions = len(all_reactions)
    glc_idx = -1
    atp_idx = -1
    for j in range(num_reactions):
        a_reaction = all_reactions[j]
        if a_reaction.description == ' —> Glucose | Enzyme: Exchange Reaction':
            glc_idx = j
        elif a_reaction.description == 'ATP —>  | Enzyme: Exchange Reaction':
            atp_idx = j
    if glc_idx == -1 or atp_idx == -1:
        raise Exception
    return glc_idx, atp_idx


metabolites = get_metabolites()
reactions = get_reactions()
S = make_stoichiometry_matrix(metabolites, reactions)
glucose_input_idx, ATP_output_idx = get_glucose_and_atp_reactions_idx(reactions)
optimal_reactions_vector = optimize_stoichiometry(S, glucose_input_idx, ATP_output_idx)
print_optimal_path(optimal_reactions_vector, reactions)
