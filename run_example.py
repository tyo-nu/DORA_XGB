"""
Example script on how to initialize and deploy DORA_XGB reaction feasibility classification models
"""
from DORA_XGB import DORA_XGB

cofactor_positioning = 'by_descending_MW' # either 'by_ascending_MW', 'by_descending_MW', 'add_concat',or 'add_subtract'
max_species = 4 # maximum number of participating species allowed on either side of a reaction (always 4 by default)
fp_type = 'ecfp4' # fingerprinting method - our consolidated models were trained with ecfp4 and 2048 bits
model_type = 'main' # choose from 'main' (cleaner models used to report results in paper) or 'spare' (noisier)
cofactors_filepath = './cofactors/expanded_cofactors_no_stereochem.tsv'

# initialize a model with selected hyperparameters
model = DORA_XGB.feasibility_classifier(cofactor_positioning = cofactor_positioning,
                                        max_species = max_species,
                                        fp_type = fp_type,
                                        nBits = 2048,
                                        model_type = model_type,
                                        cofactors_filepath = cofactors_filepath)

if __name__ == "__main__":

    # Following are examples of using each type of model to predict feasibility score and a label

    ### Example 1: positive (downhill) rule0003 reaction from BRENDA (EC: 1.1.1.1_083)
    # https://www.brenda-enzymes.org/structure.php?show=reaction&id=434984&type=S&displayType=marvin
    substrate = "CC(=O)C(=O)c1ccccc1"
    NADH = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    NAD_plus = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    original_product = "CC(O)C(=O)c1ccccc1"
    alternate_product = "CC(=O)C(O)c1ccccc1"

    # reaction strings representing the feasible, reported reaction and the alternate reaction
    # these reaction strings are written in the form "A + B = C + D"
    original_rxn_str = f"{substrate} + {NADH} = {original_product} + {NAD_plus}"
    alternate_rxn_str = f"{substrate} + {NADH} = {alternate_product} + {NAD_plus}"

    # predict the feasibility score and label of this reported reaction using selected model and hyperparameters
    feasibility_score = model.predict_proba(original_rxn_str)
    feasibility_label = model.predict_label(original_rxn_str)

    print('')
    print('Example 1 original reaction:')
    print('----------------------------')
    print(f'Predicted feasibility score with {cofactor_positioning} model: {feasibility_score:.3f} (written as "A + B = C + D")')
    print(f'Predicted feasibility label with {cofactor_positioning} model: {feasibility_label:.3f} (written as "A + B = C + D")')

    # reaction strings can also be written in the form "A.B>>C.D"
    original_rxn_str_format2 = f"{substrate}.{NADH}>>{original_product}.{NAD_plus}"
    alternate_rxn_str_format2 = f"{substrate}.{NADH}>>{alternate_product}.{NAD_plus}"

    # predict the feasibility score and label of this reported reaction using selected model and hyperparameters
    feasibility_score = model.predict_proba(original_rxn_str_format2)
    feasibility_label = model.predict_label(original_rxn_str_format2)

    print(f'Predicted feasibility score with {cofactor_positioning} model: {feasibility_score:.3f} (written as "A.B>>C.D")')
    print(f'Predicted feasibility label with {cofactor_positioning} model: {feasibility_label:.3f} (written as "A.B>>C.D")')

    # predict the feasibility score and label of the alternate reaction using selected model and hyperparameters
    feasibility_score = model.predict_proba(alternate_rxn_str)
    feasibility_label = model.predict_label(alternate_rxn_str)

    print('')
    print('Example 1 alternate reaction:')
    print('----------------------------')
    print(f'Predicted feasibility score with {cofactor_positioning} model: {feasibility_score:.3f} (written as "A + B = C + D")')
    print(f'Predicted feasibility label with {cofactor_positioning} model: {feasibility_label:.3f} (written as "A + B = C + D")')

    feasibility_score = model.predict_proba(alternate_rxn_str_format2)
    feasibility_label = model.predict_label(alternate_rxn_str_format2)

    print(f'Predicted feasibility score with {cofactor_positioning} model: {feasibility_score:.3f} (written as "A.B>>C.D")')
    print(f'Predicted feasibility label with {cofactor_positioning} model: {feasibility_label:.3f} (written as "A.B>>C.D")')

    ### Example 2: positive (downhill) rule0003 reaction from BRENDA (reversible)
    # https://www.brenda-enzymes.org/structure.php?show=reaction&id=186894&type=S&displayType=marvin (EC: 1.1.1.1._093)
    substrate = "O=CCCCCCCCCCCC(=O)O"
    NADH = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    NAD_plus = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    original_product = "O=C(O)CCCCCCCCCCCO"
    alternate_product = "O=CCCCCCCCCCCC(O)O"

    # reaction strings representing the feasible, reported reaction and the infeasible, alternate reaction
    # these reaction strings are written in the form "A + B = C + D"
    original_rxn_str = f"{substrate} + {NADH} = {original_product} + {NAD_plus}"
    alternate_rxn_str = f"{substrate} + {NADH} = {alternate_product} + {NAD_plus}"

    # predict the feasibility score and label of this reported reaction using selected model and hyperparameters
    feasibility_score = model.predict_proba(original_rxn_str)
    feasibility_label = model.predict_label(original_rxn_str)

    print('')
    print('Example 2 original reaction:')
    print('----------------------------')
    print(f'Predicted feasibility score with {cofactor_positioning} model: {feasibility_score:.3f} (written as "A + B = C + D")')
    print(f'Predicted feasibility label with {cofactor_positioning} model: {feasibility_label:.3f} (written as "A + B = C + D")')

    # reaction strings can also be written in the form "A.B>>C.D"
    original_rxn_str_format2 = f"{substrate}.{NADH}>>{original_product}.{NAD_plus}"
    alternate_rxn_str_format2 = f"{substrate}.{NADH}>>{alternate_product}.{NAD_plus}"

    print(f'Predicted feasibility score with {cofactor_positioning} model: {feasibility_score:.3f} (written as "A.B>>C.D")')
    print(f'Predicted feasibility label with {cofactor_positioning} model: {feasibility_label:.3f} (written as "A.B>>C.D")')

    # predict the feasibility scores of this alternate reaction using selected model and hyperparameters
    feasibility_score = model.predict_proba(alternate_rxn_str)
    feasibility_label = model.predict_label(alternate_rxn_str)

    print('')
    print('Example 2 alternate reaction:')
    print('----------------------------')
    print(f'Predicted feasibility score with {cofactor_positioning} model: {feasibility_score:.3f} (written as "A + B = C + D")')
    print(f'Predicted feasibility label with {cofactor_positioning} model: {feasibility_label:.3f} (written as "A + B = C + D")')

    feasibility_score = model.predict_proba(alternate_rxn_str_format2)
    feasibility_label = model.predict_label(alternate_rxn_str_format2)

    print(f'Predicted feasibility score with {cofactor_positioning} model: {feasibility_score:.3f} (written as "A.B>>C.D")')
    print(f'Predicted feasibility label with {cofactor_positioning} model: {feasibility_label:.3f} (written as "A.B>>C.D")')