import os
import pickle
import numpy as np
import pandas as pd
from .featurizations import featurizations

dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'

class feasibility_classifier:
    """
    Convenience wrapper for enzymatic reaction feasibility prediction via an XGBoost classification model.
    Our DORA-XGB models can predict the feasibility score and feasibility label of a given reaction string.
    This reaction string can be in the form "A + B = C" or "A.B>>C". Multisubstrate reactions are supported.
    Stereochemical information is presently not supported by our DORA-XGB models.
    """

    def __init__(
        self,
        cofactor_positioning: str = 'by_descending_MW',
        max_species: int = 4,
        fp_type: str = 'ecfp4',
        nBits: int = 2048,
        model_type: str = 'main',
        cofactors_filepath: str = dir_path + 'cofactors/expanded_cofactors_no_stereochem.tsv'):

        """
        Initialize feasibility classifier

        Parameters:
        ----------
        cofactor_positioning: str
            Configuration in which cofactors should be arranged along a reaction fingerprint for feasibility prediction.
            Choose from: (1) 'by_descending_MW', (2) 'by_ascending_MW', (3) 'add_subtract', and (4) 'add_concat'.
            Default: 'by_descending_MW'

        max_species: int
            The maximum number of species a reaction can have on either side.
            All models were trained with by allowing a maximum of 4 species on either side.
            Default: 4

        fp_type: str
            The fingerprinting method to first construct molecular fingerprints and subsequently, reaction fingerprints.
            Smaller alcohol dehydrogenase models were trained with a variety of fingerprinting techniques.
            Our consolidated DORA-XGB models, however, were trained with ecfp4 fingerprints only with 2048 bits.
            Default: 'ecfp4'

        nBits: int
            Number of bits to use when converting a molecular structure into a molecular fingerprint.
            Since our consolidated DORA-XGB models were trained with ecfp4 fingerprints only, this is typically 2048.
            Default: 2048

        model_type: str
            Choose from either 'main' or 'spare'.
            Our 'main' models are most updated and trained on cleanest data. These were used for results in the paper.
            Our 'spare' models are less updated and trained on noisier data. These are early iterations of our models.
            Default: 'main'

        cofactors_filepath: str
            Filepath to the list of cofactors that we currently support when constructing a reaction fingerprint.
            It is important that any input reaction contain cofactors from this list of 40+ cofactors.
            This will enable the chosen cofactor configuration to be applied correctly.
            If a cofactor not in this list is part of an input reaction string, it will be treated as a substrate.
            Default: './expanded_cofactors_no_stereochem.tsv'

        Attributes
        ----------
        cofactor_positioning: str
            Configuration in which to arrange supported cofactors along a reaction fingerprint.

        max_species: int
            Maximum number of molecular species allowed on either side of an input reaction string.

        fp_type: str
            Fingerprinting method with which to fingerprint molecules and eventually construct reaction fingerprints.

        nBits: int
            Number of bits used to construct the molecular fingerprint of each participating species.

        model_type: str
            The models in the 'main' directory are our updated and cleanest models while others are noisier.

        cofactors_filepath: str
            Filepath to a tsv file comprising 40+ supported cofactor SMILES without their stereochemical information.

        all_cofactors_wo_stereo: set
            Set of unique cofactor SMILES without their stereochemical information.

        self.DORA_XGB_model: any
            Unpickled XGBoost model that is loaded into disk depending on user's selected cofactor configuration.

        self.DORA_XGB_threshold: float
            Feasibility threshold corresponding to the loaded feasibility DORA-XGB model.

        Methods
        ----------
        predict_proba:
            Predict the feasibility score of the input reaction (between 0 and 1, the higher the better)

        predict_label:
            Predict the feasibility label of the input reaction (either 0 or 1, 0 is infeasible and 1 is feasible)
        """

        self.cofactor_positioning = cofactor_positioning
        self.max_species = max_species
        self.fp_type = fp_type
        self.nBits = nBits
        self.model_type = model_type
        self.cofactors_filepath = cofactors_filepath
        self.all_cofactors_wo_stereo = set(pd.read_csv(self.cofactors_filepath, delimiter = ',')['SMILES'])

        # Select DORA-XGB model and its corresponding feasibility threshold based on user-inputs
        if self.model_type == 'main':
            model_filepath = dir_path + f'models/{self.model_type}/all_BKM_rxns_{self.fp_type}_XGBoost_{self.max_species}_{self.cofactor_positioning}.pkl'
            threshold_filepath = dir_path + f'models/{self.model_type}/all_BKM_rxns_{self.fp_type}_XGBoost_{self.max_species}_{self.cofactor_positioning}_feasibility_threshold.txt'
            self.DORA_XGB_model = pickle.load(open(model_filepath, 'rb'))
            self.DORA_XGB_threshold = np.loadtxt(threshold_filepath).item()

        if self.model_type == 'spare':
            model_filepath = dir_path + f'models/{self.model_type}/xgboost_{self.fp_type}_{self.nBits}_{self.max_species}_{self.cofactor_positioning}.pkl'
            threshold_filepath = dir_path + f'models/{self.model_type}/xgboost_{self.fp_type}_{self.nBits}_{self.max_species}_{self.cofactor_positioning}_feasibility_threshold.txt'
            self.DORA_XGB_model = pickle.load(open(model_filepath, 'rb'))
            self.DORA_XGB_threshold = np.loadtxt(threshold_filepath).item()

    def predict_proba(self, rxn_str: str) -> float:

        # initialize a reaction object with our custom-built featurizations package
        rxn_object = featurizations.reaction(rxn_str)

        # convert reaction string to a reaction fingerprint
        rxn_fp = rxn_object.rxn_2_fp_w_positioning(fp_type = self.fp_type,
                                                   cofactor_positioning = self.cofactor_positioning,
                                                   all_cofactors_wo_stereo = self.all_cofactors_wo_stereo)

        # reshape since only single sample
        rxn_fp = rxn_fp.reshape(1, -1)

        rxn_feasibility_score = self.DORA_XGB_model.predict_proba(rxn_fp)[:,1][0]
        return rxn_feasibility_score

    def predict_label(self, rxn_str) -> int:

        rxn_feasibility_score = self.predict_proba(rxn_str)

        if rxn_feasibility_score >= self.DORA_XGB_threshold:
            return 1
        else:
            return 0

if __name__ == "__main__":

    # We instante 4 DORA-XGB feasibility models by selecting each of the 4 cofactor configurations
    by_asc_MW_model = feasibility_classifier(cofactor_positioning = 'by_ascending_MW')
    by_desc_MW_model = feasibility_classifier(cofactor_positioning = 'by_descending_MW')
    add_concat_model = feasibility_classifier(cofactor_positioning = 'add_concat')
    add_subtract_model = feasibility_classifier(cofactor_positioning = 'add_subtract')

    # Following are examples of using each type of model to predict feasibility score and a label

    ### Example 1: positive (downhill) rule0003 reaction from BRENDA (EC: 1.1.1.1_083)
    # https://www.brenda-enzymes.org/structure.php?show=reaction&id=434984&type=S&displayType=marvin
    substrate = "CC(=O)C(=O)c1ccccc1"
    NADH = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    NAD_plus = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    original_product = "CC(O)C(=O)c1ccccc1"
    alternate_product = "CC(=O)C(O)c1ccccc1"

    # reaction strings representing the feasible, reported reaction and the alternate reaction
    original_rxn_str = f"{substrate} + {NADH} = {original_product} + {NAD_plus}"
    alternate_rxn_str = f"{substrate} + {NADH} = {alternate_product} + {NAD_plus}"

    # predict the feasibility scores of this reported reaction using all models
    by_asc_MW_model_score = by_asc_MW_model.predict_proba(original_rxn_str)
    by_desc_MW_model_score = by_desc_MW_model.predict_proba(original_rxn_str)
    add_concat_model_score = add_concat_model.predict_proba(original_rxn_str)
    add_subtract_model_score = add_subtract_model.predict_proba(original_rxn_str)

    print('')
    print('Example 1 original reaction:')
    print('----------------------------')
    print(f'Feasibility score with by ascending MW model: {by_asc_MW_model_score:.3f}')
    print(f'Feasibility score with by descending MW model: {by_desc_MW_model_score:.3f}')
    print(f'Feasibility score with add-concatenate model: {add_concat_model_score:.3f}')
    print(f'Feasibility score with add-subtract model: {add_subtract_model_score:.3f}')

    # predict the feasibility labels of this reported reaction using all models
    by_asc_MW_model_label = by_asc_MW_model.predict_label(original_rxn_str)
    by_desc_MW_model_label = by_desc_MW_model.predict_label(original_rxn_str)
    add_concat_model_label = add_concat_model.predict_label(original_rxn_str)
    add_subtract_model_label = add_subtract_model.predict_label(original_rxn_str)

    print('')
    print(f'Feasibility label with by ascending MW model: {by_asc_MW_model_label}')
    print(f'Feasibility label with by descending MW model: {by_desc_MW_model_label}')
    print(f'Feasibility label with add-concatenate model: {add_concat_model_label}')
    print(f'Feasibility label with add-subtract model: {add_subtract_model_label}')

    # predict the feasibility scores of this alternate reaction using all models
    by_asc_MW_model_score = by_asc_MW_model.predict_proba(alternate_rxn_str)
    by_desc_MW_model_score = by_desc_MW_model.predict_proba(alternate_rxn_str)
    add_concat_model_score = add_concat_model.predict_proba(alternate_rxn_str)
    add_subtract_model_score = add_subtract_model.predict_proba(alternate_rxn_str)

    print('')
    print('Example 1 alternate reaction:')
    print('----------------------------')
    print(f'Feasibility score with by ascending MW model: {by_asc_MW_model_score:.3f}')
    print(f'Feasibility score with by descending MW model: {by_desc_MW_model_score:.3f}')
    print(f'Feasibility score with add-concatenate model: {add_concat_model_score:.3f}')
    print(f'Feasibility score with add-subtract model: {add_subtract_model_score:.3f}')

    # predict the feasibility labels of this reported reaction using all models
    by_asc_MW_model_label = by_asc_MW_model.predict_label(alternate_rxn_str)
    by_desc_MW_model_label = by_desc_MW_model.predict_label(alternate_rxn_str)
    add_concat_model_label = add_concat_model.predict_label(alternate_rxn_str)
    add_subtract_model_label = add_subtract_model.predict_label(alternate_rxn_str)

    print('')
    print(f'Feasibility label with by ascending MW model: {by_asc_MW_model_label}')
    print(f'Feasibility label with by descending MW model: {by_desc_MW_model_label}')
    print(f'Feasibility label with add-concatenate model: {add_concat_model_label}')
    print(f'Feasibility label with add-subtract model: {add_subtract_model_label}')

    ### Example 2: positive (downhill) rule0003 reaction from BRENDA (reversible)
    # https://www.brenda-enzymes.org/structure.php?show=reaction&id=186894&type=S&displayType=marvin (EC: 1.1.1.1._093)
    substrate = "O=CCCCCCCCCCCC(=O)O"
    NADH = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    original_product = "O=C(O)CCCCCCCCCCCO"
    NAD_plus = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"

    alternate_product = "O=CCCCCCCCCCCC(O)O"

    # original reaction string representing the feasible, reported reaction
    original_rxn_str = f"{substrate} + {NADH} = {original_product} + {NAD_plus}"

    # predict the feasibility scores of this reported reaction using all models
    by_asc_MW_model_score = by_asc_MW_model.predict_proba(original_rxn_str)
    by_desc_MW_model_score = by_desc_MW_model.predict_proba(original_rxn_str)
    add_concat_model_score = add_concat_model.predict_proba(original_rxn_str)
    add_subtract_model_score = add_subtract_model.predict_proba(original_rxn_str)

    print('')
    print('Example 2 original reaction:')
    print('----------------------------')
    print(f'Feasibility score with by ascending MW model: {by_asc_MW_model_score:.3f}')
    print(f'Feasibility score with by descending MW model: {by_desc_MW_model_score:.3f}')
    print(f'Feasibility score with add-concatenate model: {add_concat_model_score:.3f}')
    print(f'Feasibility score with add-subtract model: {add_subtract_model_score:.3f}')

    by_asc_MW_model_label = by_asc_MW_model.predict_label(original_rxn_str)
    by_desc_MW_model_label = by_desc_MW_model.predict_label(original_rxn_str)
    add_concat_model_label = add_concat_model.predict_label(original_rxn_str)
    add_subtract_model_label = add_subtract_model.predict_label(original_rxn_str)

    print('')
    print(f'Feasibility label with by ascending MW model: {by_asc_MW_model_label}')
    print(f'Feasibility label with by descending MW model: {by_desc_MW_model_label}')
    print(f'Feasibility label with add-concatenate model: {add_concat_model_label}')
    print(f'Feasibility label with add-subtract model: {add_subtract_model_label}')

    # alternate reaction string representing the infeasible, reported reaction
    alt_rxn_str = f"{substrate} + {NADH} = {alternate_product} + {NAD_plus}"

    # predict the feasibility scores of this alternate reaction using all models
    by_asc_MW_model_score = by_asc_MW_model.predict_proba(alternate_rxn_str)
    by_desc_MW_model_score = by_desc_MW_model.predict_proba(alternate_rxn_str)
    add_concat_model_score = add_concat_model.predict_proba(alternate_rxn_str)
    add_subtract_model_score = add_subtract_model.predict_proba(alternate_rxn_str)

    print('')
    print('Example 2 alternate reaction:')
    print('----------------------------')
    print(f'Feasibility score with by ascending MW model: {by_asc_MW_model_score:.3f}')
    print(f'Feasibility score with by descending MW model: {by_desc_MW_model_score:.3f}')
    print(f'Feasibility score with add-concatenate model: {add_concat_model_score:.3f}')
    print(f'Feasibility score with add-subtract model: {add_subtract_model_score:.3f}')

    # predict the feasibility labels of this reported reaction using all models
    by_asc_MW_model_label = by_asc_MW_model.predict_label(alternate_rxn_str)
    by_desc_MW_model_label = by_desc_MW_model.predict_label(alternate_rxn_str)
    add_concat_model_label = add_concat_model.predict_label(alternate_rxn_str)
    add_subtract_model_label = add_subtract_model.predict_label(alternate_rxn_str)

    print('')
    print(f'Feasibility label with by ascending MW model: {by_asc_MW_model_label}')
    print(f'Feasibility label with by descending MW model: {by_desc_MW_model_label}')
    print(f'Feasibility label with add-concatenate model: {add_concat_model_label}')
    print(f'Feasibility label with add-subtract model: {add_subtract_model_label}')



