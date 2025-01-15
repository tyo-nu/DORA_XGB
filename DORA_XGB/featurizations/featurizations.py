"""
Author: Yash Chainani

This is a custom-designed package built atop of several rdkit functionalities
The goal of this package is to provide various fingerprinting methods to featurize both compounds and reactions
Compounds can be featurized by converting their SMILES into either hashed fingerprints (ecfp4, atom-pair)
or descriptor based fingerprints (Mordred and MACCS)
Reactions can be featurized by concatenating the fingerprints of constituent compounds together into a single vector
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import DataStructs

try:
    from mordred import Calculator, descriptors
except Exception as e:
    pass
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
from typing import Tuple, Union, Optional

# MHFP fingerprints were installed from https://github.com/reymond-group/mhfp
# There is also a nice blog post on using them from https://xinhaoli74.github.io/posts/2020/05/TMAP/
# Journal publication with these fingerprints: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0321-8
from mhfp.encoder import MHFPEncoder

# Silence non-critical RDKit warnings to minimize unnecessary outputs
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

class compound:
    """
    Compound class for various fingerprinting methods for a single compound
    """

    def __init__(self, smiles: str):
        self.smiles = smiles

    # checked and tested
    def _canonicalize_smiles(self) -> str:
        """
    Canonicalize the SMILES string of a molecule.

    This method uses RDKit's Chem module to convert an original SMILES string to its canonical form.
    The canonical form ensures a consistent representation of the molecule.
    If the canonicalization fails (e.g., due to invalid input SMILES),
    the method will catch the exception and return the original SMILES string.

    :return: A string of the canonicalized SMILES or the original SMILES if the canonicalization fails.
    """
        uncanon_smi = self.smiles  # assume original smiles not canonical

        try:
            canon_smi = Chem.MolToSmiles(
                Chem.MolFromSmiles(uncanon_smi)
            )  # try to canonicalize
        except:
            canon_smi = uncanon_smi  # return original if canonicalization failed

        return canon_smi

    # checked and tested
    def remove_stereo(self) -> str:
        """
        Removes stereochemistry if any is present after first canonicalizing smiles string
        :return: Smiles string without stereochemistry information
        """

        # canonicalize input smiles first then assume original smiles have stereochemistry
        smiles_w_stereo = self._canonicalize_smiles()

        try:
            mol = Chem.MolFromSmiles(smiles_w_stereo)
            Chem.RemoveStereochemistry(mol)
            smiles_wo_stereo = Chem.MolToSmiles(mol)
        except:
            smiles_wo_stereo = smiles_w_stereo

        return smiles_wo_stereo

    # checked and tested
    def _smiles_2_morganfp(self, radius: int, nBits: int) -> Union[np.ndarray, None]:
        """
        Generate Morgan fingerprints for a compound after canonicalizing smiles and removing stereochemistry
        :param radius: Radius for fragmentation
        :param nBits: Output dimensions of fingerprint
        :return: Morgan/ ecfp4 fingerprint of specified radius and dimension
        """
        try:
            canon_smi_wo_stereo = self.remove_stereo()
            mol = Chem.MolFromSmiles(canon_smi_wo_stereo)
        except Exception as E:
            return None

        if mol:
            fp = Chem.AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits)
            return np.array(fp)

        return None

    # checked and tested
    def _smiles_2_MACCS(self) -> Union[np.ndarray, None]:
        """
        Generate MACCS keys of a compound after canonicalizing smiles and removing stereochemistry
        Some code is inspired from the following link:
        https://chem.libretexts.org/Courses/Intercollegiate_Courses/Cheminformatics/06%3A_Molecular_Similarity/6.04%3A_Python_Assignment
        :return: MACCS fingerprint
        """
        try:
            canon_smi_wo_stereo = self.remove_stereo()
            mol = Chem.MolFromSmiles(canon_smi_wo_stereo)
        except Exception as E:
            return None

        if mol:
            fp = MACCSkeys.GenMACCSKeys(mol)
            return np.array(fp)

        return None

    # checked and tested
    def _smiles_2_atompair(self, nBits:int) -> Union[np.ndarray, None]:
        """
        Calculate Atom Pair fingerprints for a list of SMILES strings.

        :param smiles_list: List of SMILES strings representing the molecules.
        :param nBits: The length of the fingerprint.
        :return: Numpy array where each row corresponds to the Atom Pair fingerprint of a molecule.
        """
        canon_smi_wo_stereo = self.remove_stereo()
        mol = Chem.MolFromSmiles(canon_smi_wo_stereo)
        if mol is not None:  # Check if the molecule is correctly parsed
            fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=nBits)
            arr = np.zeros((1,))
            # Convert the RDKit explicit bit vector to a numpy array
            DataStructs.ConvertToNumpyArray(fp, arr)
            return np.array(arr)
        else:
            return None

    # checked and tested
    def _smiles_2_modred(self) -> Union[np.ndarray, None]:
        """
        Generate modred fingerprints and replace None/ Nan values with zero
        :return:
        """
        try:
            canon_smi_wo_stereo = self.remove_stereo()
            mol = Chem.MolFromSmiles(canon_smi_wo_stereo)
        except Exception as E:
            return None

        if mol:
            calc = Calculator(descriptors, ignore_3D=True)
            fp = calc(mol)
            fp = np.array(fp).reshape(-1, 1)
            fp_no_nan = np.nan_to_num(
                fp, copy=True
            )  # converts any nan values in modred descriptor to zero

            # Scale (Z-score) the vector since descriptors can vary significantly in magnitudes
            scaler = StandardScaler()
            scaler.fit(fp_no_nan)
            fp_scaled = scaler.transform(fp_no_nan)

            return fp_scaled.reshape(len(fp_scaled))

        return None

    # checked and tested
    def _smiles_2_MHFP(self, radius:int) -> Union[np.ndarray,None]:
        mhfp_encoder = MHFPEncoder()

        try:
            canon_smi_wo_stereo = self.remove_stereo()
            mol = Chem.MolFromSmiles(canon_smi_wo_stereo)
        except Exception as E:
            return None

        if mol:
            fp = mhfp_encoder.encode(canon_smi_wo_stereo, radius = radius)
            return np.array(fp)

        return None

    # checked and tested
    def smiles_2_fp(self,
                    fp_type,
                    radius: int = 2,
                    nBits:int = 2048,
                    is_folded: bool = True,
                    dim: int = 2048) -> np.ndarray:

        if fp_type == "morgan" or fp_type == "ecfp4":
            reactant_fp = self._smiles_2_morganfp(radius = radius, nBits = nBits)
            return reactant_fp

        if fp_type == "MACCS":
            reactant_fp = self._smiles_2_MACCS()
            return reactant_fp

        if fp_type == "atom_pair":
            reactant_fp = self._smiles_2_atompair(nBits = nBits)
            return reactant_fp

        if fp_type == "mordred":
            reactant_fp = self._smiles_2_modred()
            return reactant_fp

        if fp_type == "MHFP4": # the default radius here is actually 6
            reactant_fp = self._smiles_2_MHFP(radius = radius)
            return reactant_fp

        else:
            raise ValueError("Please enter a valid fingerprinting method, such as: morgan/ecfp4, MACCS, atom_pair, modred, or MHFP4")

class reaction:
    """
    Reaction class for various fingerprinting methods for a full reaction
    """
    def __init__(self, rxn_str: str):
        self.rxn_str = rxn_str

    # checked and tested
    def _rxn_2_cpds(self) -> Tuple[str, str]:
        """
        Parse a reaction string to return two lists: a reactants list and a products list
        :return reactants_str: string of reactants on the LHS of the reaction (index 0)
        :return products_str: string of products on the RHS of the reaction (index 1)
        """
        rxn_str = self.rxn_str
        if " = " in rxn_str:
            reactants_str = rxn_str.split(" = ")[0]
            products_str = rxn_str.split(" = ")[1]
            return reactants_str, products_str
        if ">>" in rxn_str:
            reactants_str = rxn_str.split(">>")[0]
            products_str = rxn_str.split(">>")[1]
            return reactants_str, products_str

    # checked and tested
    def _load_cofactors_set_w_stereo(self, cofactors_filepath: str) -> set:
        """
        Get the full list of cofactors available on memory
        Stereochemistry of cofactors is NOT removed here
        :param cofactors_filepath: relative path to local cofactors list
        :return: set of cofactors WITH stereochemistry
        """
        try:
            cofactors_df = pd.read_csv(cofactors_filepath)
            cofactors_set_w_stereo = set(cofactors_df["SMILES"])
        except KeyError:
            cofactors_df = pd.read_csv(cofactors_filepath,delimiter = '\t')
            cofactors_set_w_stereo = set(cofactors_df["SMILES"])
        return cofactors_set_w_stereo

    # checked and tested
    def print_full_cofactors_set_wo_stereo(self, cofactors_filepath: str) -> set:
        """
        Get the full list of cofactors available on memory
        Stereochemistry of cofactors IS removed here (even if the initial file did not contain stereochemistry)
        :param cofactors_filepath: relative path to local cofactors list
        :return: set of cofactors WITHOUT any stereochemistry
        """
        cofactors_set_w_stereo = self._load_cofactors_set_w_stereo(cofactors_filepath)
        cofactors_set_wo_stereo = []
        for cofactor_smiles in cofactors_set_w_stereo:
            cpd_w_stereo = compound(cofactor_smiles)
            cpd_wo_stereo = cpd_w_stereo.remove_stereo()
            cofactors_set_wo_stereo.append(cpd_wo_stereo)
        return set(cofactors_set_wo_stereo)

    # checked and tested
    def get_substrates(self, all_cofactors_wo_stereo: set) -> list:
        """
        Extract the substrate/s from a reaction string of the form "A + B = C + D" or "A.B>>C.D"
        :param all_cofactors_wo_stereo: set of cofactors
        :return:
        """
        # separate rxn str into two lists of strings - reactants and products
        reactants_str, products_str = self._rxn_2_cpds()
        substrates = []

        # first, determine the delimiter that separates the substrates
        if " + " in reactants_str:
            delimiter = " + "
        elif "." in reactants_str:
            delimiter = "."
        else:
            delimiter = None

        # then, split the reactants based on their corresponding delimiter
        reactant_smiles_list = reactants_str.split(delimiter) if delimiter else [reactants_str]

        # for each reactant's SMILES string, create a compound object
        # then canonicalize and remove any stereochemical information
        for reactant_smiles in reactant_smiles_list:
            reactant = compound(reactant_smiles)
            canon_smi_wo_stereo = reactant.remove_stereo()

            # if canonicalized SMILES are not in the set of cofactors
            if canon_smi_wo_stereo not in all_cofactors_wo_stereo:
                # then this is a substrate
                substrates.append(canon_smi_wo_stereo)

        return substrates

    # checked and tested
    def get_products(self, all_cofactors_wo_stereo: set) -> list:
        """
        Extract the product/s from a reaction string of the form "A + B = C + D" or "A.B>>C.D"
        :param all_cofactors_wo_stereo:
        :return:
        """
        # separate rxn str into two lists of strings - reactants and products
        reactants_str, products_str = self._rxn_2_cpds()
        products = []

        # first, determine the delimiter that separates the products
        if " + " in products_str:
            delimiter = " + "
        elif "." in products_str:
            delimiter = "."
        else:
            delimiter = None

        # then, split the products based on their corresponding delimiter
        product_smiles_list = products_str.split(delimiter) if delimiter else [products_str]

        # for each product's SMILES string, create a compound object
        # then canonicalize and remove any stereochemical information
        for product_smiles in product_smiles_list:
            product = compound(product_smiles)
            canon_smi_wo_stereo = product.remove_stereo()

            # if canonicalized SMILES are not in the set of cofactors
            if canon_smi_wo_stereo not in all_cofactors_wo_stereo:
                # then this is a product
                products.append(canon_smi_wo_stereo)

        return products

    # checked and tested
    def get_lhs_cofactors(self, all_cofactors_wo_stereo: set) -> list:

        # separate rxn str into two lists of strings - reactants and products
        reactants_str, products_str = self._rxn_2_cpds()
        lhs_cofactors = []

        # first, determine the delimiter that separates the reactants (since searching for cofactors on the LHS)
        if " + " in reactants_str:
            delimiter = " + "
        elif "." in reactants_str:
            delimiter = "."
        else:
            delimiter = None

        # then, split the reactants based on their corresponding delimiter
        reactants_smiles_list = reactants_str.split(delimiter) if delimiter else [reactants_str]

        # for each reactant's SMILES string, create a compound object
        # then canonicalize and remove any stereochemical information
        for reactant_smiles in reactants_smiles_list:
            reactant = compound(reactant_smiles)
            canon_smi_wo_stereo = reactant.remove_stereo()

            # if canonicalized SMILES are in cofactors list
            if canon_smi_wo_stereo in all_cofactors_wo_stereo:
                lhs_cofactors.append(canon_smi_wo_stereo)  # then this is a cofactor on the lhs

        return lhs_cofactors

    # checked and tested
    def get_rhs_cofactors(self, all_cofactors_wo_stereo: set) -> list:

        # separate rxn str into two lists of strings - reactants and products
        reactants_str, products_str = self._rxn_2_cpds()
        rhs_cofactors = []

        # first, determine the delimiter that separates the products (since searching for cofactors on the RHS)
        if " + " in products_str:
            delimiter = " + "
        elif "." in products_str:
            delimiter = "."
        else:
            delimiter = None

        # then, split the products based on their corresponding delimiter
        products_smiles_list = products_str.split(delimiter) if delimiter else [products_str]

        # for each product's SMILES string, create a compound object
        # then canonicalize and remove any stereochemical information
        for product_smiles in products_smiles_list:
            product = compound(product_smiles)
            canon_smi_wo_stereo = product.remove_stereo()

            # if canonicalized SMILES are in cofactors list
            if canon_smi_wo_stereo in all_cofactors_wo_stereo:
                rhs_cofactors.append(canon_smi_wo_stereo) # then this is a cofactor on the rhs

        return rhs_cofactors

    # checked and tested
    def _reorder_cofactors_by_MW(self, cofactors_list: list, ascending: bool) -> list:
        """
        Rearrange cofactors in ascending or descending molecular weights
        :param cofactors_list: list of cofactor SMILES
        :param ascending: cofactors will be arranged from lowest to highest molecular weight if true
        :return: list of rearranged cofactors
        """
        MW_list = []
        for cofactor_smiles in cofactors_list:
            mol = Chem.MolFromSmiles(cofactor_smiles)
            mw = Descriptors.ExactMolWt(mol)
            MW_list.append(mw)

        sorted_cofactors_list = [val for (_, val) in sorted(zip(MW_list, cofactors_list), key=lambda x: x[0])]

        if ascending:
            return sorted_cofactors_list

        else:
            return sorted_cofactors_list[::-1]

    # checked and tested
    def _reorder_substrates_by_MW(self, substrates_list: list, ascending: bool) -> list:
        """
        Rearrange substrates in ascending or descending molecular weights
        :param substrates_list: list of substrate SMILES
        :param ascending: substrates will be arranged from lowest to highest molecular weight if true
        :return: list of rearranged substrates
        """
        MW_list = []
        for substrate_SMILES in substrates_list:
            mol = Chem.MolFromSmiles(substrate_SMILES)
            mw = Descriptors.ExactMolWt(mol)
            MW_list.append(mw)

        sorted_substrates_list = [val for (_, val) in sorted(zip(MW_list, substrates_list), key=lambda x: x[0])]

        if ascending:
            return sorted_substrates_list

        else:
            return sorted_substrates_list[::-1]

    # checked and tested
    def _reorder_products_by_MW(self, products_list: list, ascending: bool) -> list:
        """
        Rearrange products in ascending or descending molecular weights
        :param products_list: list of product SMILES
        :param ascending: products will be arranged from lowest to highest molecular weight if true
        :return: list of rearranged products
        """
        MW_list = []
        for product_SMILES in products_list:
            mol = Chem.MolFromSmiles(product_SMILES)
            mw = Descriptors.ExactMolWt(mol)
            MW_list.append(mw)

        sorted_products_list = [val for (_, val) in sorted(zip(MW_list, products_list), key=lambda x: x[0])]

        if ascending:
            return sorted_products_list

        else:
            return sorted_products_list[::-1]

    # checked and tested
    def rxn_2_fp(self, type: str, max_species: int) -> np.ndarray:
        """
        Fingerprint a reaction string of form "substrate_smiles + cofactor_smiles = product_smiles + cofactor_smiles"


        :param type: Type of fingerprint to generate (eg: morgan, ecfp4, MACCS)
        :param max_species: Number of species on either side of reaction
                            If the number of species specified if less than the actual number, pad with extra zeroes
        :return:
        """
        # first, the input reaction string is broken into its reactions (LHS) and products (RHS)
        reactants_str, products_str = self._rxn_2_cpds()

        # initialize empty numpy arrays to stack fingerprints in a horizontal fashion (i.e. side by side)
        all_reactants_fp = np.array([])
        all_products_fp = np.array([])

        # initialize a counter to track the number of reactants that have been featurized
        reactant_counter = 0

        # initialize a counter to track the number of products that have been featurized
        product_counter = 0

        # featurize all reactants on the LHS of the reaction string (including cofactors)
        for reactant_smiles in reactants_str.split(" + "):
            reactant_counter += 1
            reactant_object = compound(reactant_smiles)

            if type == "morgan" or type == "ecfp4":
                reactant_fp = reactant_object._smiles_2_morganfp(radius=2, nBits=2048)

            if type == "MACCS":
                reactant_fp = reactant_object._smiles_2_MACCS()

            if type == "atom_pair":
                reactant_fp = reactant_object._smiles_2_atompair(nBits=2048)

            if type == "mordred":
                reactant_fp = reactant_object._smiles_2_modred()

            if type == "MHFP4":
                reactant_fp = reactant_object._smiles_2_MHFP(radius = 2)

            num_features = len(reactant_fp)  # length of each reactant's fingerprint
            all_reactants_fp = np.concatenate((all_reactants_fp, reactant_fp), axis=None)

        # if number of reactants featurized is less than the maximum number of species
        # then pad the LHS fingerprint with extra zeroes (dummy reactants)
        if reactant_counter < max_species:
            lhs_diff = max_species - reactant_counter
            dummy_fp_for_lhs = np.zeros(num_features)  # create a dummy fingerprint of the same length

            # add as many dummy fingerprints as the difference in number of reactants and max species
            for i in range(0, lhs_diff):
                all_reactants_fp = np.concatenate((all_reactants_fp, dummy_fp_for_lhs), axis=None)

        ### Featurize all products on the RHS of the products string (includes cofactors)
        for product_smiles in products_str.split(" + "):
            product_counter += 1
            product_object = compound(product_smiles)

            if type == "morgan" or type == "ecfp4":
                product_fp = product_object._smiles_2_morganfp(radius=2, nBits=2048)

            if type == "MACCS":
                product_fp = product_object._smiles_2_MACCS()

            if type == "atom_pair":
                product_fp = product_object._smiles_2_atompair(nBits=2048)

            if type == "mordred":
                product_fp = product_object._smiles_2_modred()

            if type == "MHFP4":
                product_fp = product_object._smiles_2_MHFP(radius=2)

            all_products_fp = np.concatenate((all_products_fp, product_fp), axis=None)
            num_features = len(product_fp)  # length of each product's fingerprint

        # if number of products featurized is less than the maximum number of species
        # then pad the RHS fingerprint with extra zeroes (dummy products)
        if product_counter < max_species:
            rhs_diff = max_species - product_counter
            dummy_fp_for_rhs = np.zeros(num_features)  # create a dummy fingerprint of the same length

            # add as many dummy fingerprints as the difference in number of reactants and max species
            for i in range(0, rhs_diff):
                all_products_fp = np.concatenate((all_products_fp, dummy_fp_for_rhs), axis=None)

        ### Finally, concatenate fingerprints on both sides of the reaction for a full reaction fingerprint
        reaction_fp = np.concatenate((all_reactants_fp, all_products_fp), axis=None)

        return reaction_fp

    # checked and tested and ascending MW and descending MW configurations
    # checked and tested add_concat
    # now checking and tested add_random
    # checked half_random but need to find a way to test it also
    # checked full_random but need to find a way to test it also
    def rxn_2_fp_w_positioning(self,
                               fp_type: str,
                               radius: int = 2,
                               nBits: int = 2048,
                               is_folded: bool = True,
                               dim: int = 2048,
                               max_species: int = 4,
                               cofactor_positioning: str = None,
                               all_cofactors_wo_stereo: set = None) -> np.ndarray:
        """
        :param fp_type: Type of reaction fingerprint to generate ( 'morgan/ecfp4', 'MACCS', 'mordred','atom_pair', 'MHFP4')
        :param radius: Radius of fragmentation if using morgan or ecfp4
        :param nBits: Number of bits if using morgan or ecfp4
        :param is_folded: If fingerprint should be folded or not when using MAP4 (MAP4 now removed)
        :param dim: Number of bits if using MAP4 (MAP4 now removed)
        :param max_species: Maximum number of species on each side (pad zeroes if less species than max are present)
        :param cofactor_positioning: Arrangement of cofactors: in increasing MW, decreasing MW, or as per reaction rule
        :param all_cofactors_wo_stereo: set of all cofactors without stereochemistry
        :return: reaction fingerprint
        """

        # extract the substrates and cofactors on the LHS of the reaction from the reaction string
        substrates_list = self.get_substrates(all_cofactors_wo_stereo)
        lhs_cofactors_list = self.get_lhs_cofactors(all_cofactors_wo_stereo)

        # extract the products and cofactors on the RHS of the reaction from the reaction string
        products_list = self.get_products(all_cofactors_wo_stereo)
        rhs_cofactors_list = self.get_rhs_cofactors(all_cofactors_wo_stereo)

        # initialize empty arrays to store fingerprints for both the LHS and RHS of reaction
        all_lhs_fp = np.array([])
        all_rhs_fp = np.array([])

        # initialize counter to track reactants that are featurized
        reactant_counter = 0

        # initialize counter to track products that are featurized
        product_counter = 0

        # for reaction fingerprints created using the "by_ascending_MW" and "by_descending_MW" cofactor positioning,

        if cofactor_positioning == "by_ascending_MW":

            # If cofactors are present, rearrange from lightest to heaviest (ascending molecular weights)
            if lhs_cofactors_list:
                lhs_cofactors = self._reorder_cofactors_by_MW(lhs_cofactors_list, ascending = True)

            if rhs_cofactors_list:
                rhs_cofactors = self._reorder_cofactors_by_MW(rhs_cofactors_list, ascending = True)

            if substrates_list:
                reordered_substrates = self._reorder_substrates_by_MW(substrates_list, ascending = True)

            if products_list:
                reordered_products = self._reorder_products_by_MW(products_list, ascending = True)

            # Featurize all reactants on the LHS of the reaction string - begin with substrates then do cofactors
            if substrates_list != []:
                for substrate_smiles in reordered_substrates:
                    reactant_counter += 1
                    reactant_object = compound(substrate_smiles)
                    reactant_fp = reactant_object.smiles_2_fp(fp_type = fp_type,
                                                              radius = radius,
                                                              nBits = nBits,
                                                              is_folded = is_folded,
                                                              dim = dim)

                    try:
                        num_features = len(reactant_fp)  # length of each reactant's fingerprint
                    except:
                        # if error occurs because no reactants are present, assume num_features = 2048
                        num_features = 2048
                    all_lhs_fp = np.concatenate((all_lhs_fp, reactant_fp), axis=None)

            # then repeat for cofactors on the LHS
            if lhs_cofactors_list != []:
                for lhs_cofactor_smiles in lhs_cofactors:
                    reactant_counter += 1
                    lhs_cofactor_object = compound(lhs_cofactor_smiles)
                    lhs_cofactor_fp = lhs_cofactor_object.smiles_2_fp(fp_type = fp_type,
                                                                      radius = radius,
                                                                      nBits = nBits,
                                                                      is_folded = is_folded,
                                                                      dim = dim)

                    try:
                        num_features = len(lhs_cofactor_fp)  # length of each cofactor's fingerprint
                    except:
                        # if error occurs because no cofactors are present on the lhs, assume num_features = 2048
                        num_features = 2048
                    all_lhs_fp = np.concatenate((all_lhs_fp, lhs_cofactor_fp), axis=None)

            # if number of reactants featurized is less than the maximum number of species
            # then pad the LHS fingerprint with extra zeroes (dummy reactants)
            if reactant_counter < max_species:
                lhs_diff = max_species - reactant_counter
                dummy_fp_for_lhs = np.zeros(num_features)  # create a dummy fingerprint of the same length

                # add as many dummy fingerprints as the difference in number of reactants and max species
                for i in range(0, lhs_diff):
                    all_lhs_fp = np.concatenate((all_lhs_fp, dummy_fp_for_lhs), axis=None)

            # Featurize all products on the RHS of the reaction string - begin with products then do cofactors
            if products_list != []:
                for product_smiles in reordered_products:
                    product_counter += 1
                    product_object = compound(product_smiles)
                    product_fp = product_object.smiles_2_fp(fp_type = fp_type,
                                                            radius = radius,
                                                            nBits = nBits,
                                                            is_folded = is_folded,
                                                            dim = dim)

                    num_features = len(product_fp)  # length of each reactant's fingerprint
                    all_rhs_fp = np.concatenate((all_rhs_fp, product_fp), axis=None)

            # then repeat for cofactors on the RHS
            if rhs_cofactors_list != []:
                for rhs_cofactor_smiles in rhs_cofactors:
                    product_counter += 1
                    rhs_cofactor_object = compound(rhs_cofactor_smiles)
                    rhs_cofactor_fp = rhs_cofactor_object.smiles_2_fp(fp_type = fp_type,
                                                            radius = radius,
                                                            nBits = nBits,
                                                            is_folded = is_folded,
                                                            dim = dim)

                    num_features = len(rhs_cofactor_fp)  # length of each reactant's fingerprint
                    all_rhs_fp = np.concatenate((all_rhs_fp, rhs_cofactor_fp), axis=None)

            # if number of products featurized is less than the maximum number of species
            # then pad the RHS fingerprint with extra zeroes (dummy products)
            if product_counter < max_species:
                rhs_diff = max_species - product_counter
                dummy_fp_for_rhs = np.zeros(num_features)

                # add as many dummy fingerprints as the difference in number of products and max species
                for i in range(0, rhs_diff):
                    all_rhs_fp = np.concatenate((all_rhs_fp, dummy_fp_for_rhs), axis=None)

            ### Finally, concatenate fingerprints on both sides of the reaction for a full reaction fingerprint
            reaction_fp = np.concatenate((all_lhs_fp, all_rhs_fp), axis=None)

            return reaction_fp

        if cofactor_positioning == "by_descending_MW":
            
            # If cofactors present, rearrange from lightest to heaviest (ascending molecular weights)
            if lhs_cofactors_list:
                lhs_cofactors = self._reorder_cofactors_by_MW(lhs_cofactors_list, ascending=False)

            if rhs_cofactors_list:
                rhs_cofactors = self._reorder_cofactors_by_MW(rhs_cofactors_list, ascending=False)

            if substrates_list:
                reordered_substrates = self._reorder_substrates_by_MW(substrates_list, ascending=False)

            if products_list:
                reordered_products = self._reorder_products_by_MW(products_list, ascending=False)
                
            # Featurize all reactants on the LHS of the reaction string - begin with substrates then do cofactors
            if substrates_list:
                
                for substrate_smiles in reordered_substrates:
                    reactant_counter += 1
                    reactant_object = compound(substrate_smiles)
                    
                    reactant_fp = reactant_object.smiles_2_fp(fp_type = fp_type,
                                                              radius = radius,
                                                              nBits = nBits,
                                                              is_folded = is_folded,
                                                              dim = dim)

                    num_features = len(reactant_fp)  # length of each reactant's fingerprint
                    all_lhs_fp = np.concatenate((all_lhs_fp, reactant_fp), axis=None)

            # then repeat for cofactors on the LHS
            if lhs_cofactors_list:
                
                for lhs_cofactor_smiles in lhs_cofactors:
                    reactant_counter += 1
                    lhs_cofactor_object = compound(lhs_cofactor_smiles)
                    lhs_cofactor_fp = lhs_cofactor_object.smiles_2_fp(fp_type = fp_type,
                                                                      radius = radius,
                                                                      nBits = nBits,
                                                                      is_folded = is_folded,
                                                                      dim = dim)

                    num_features = len(lhs_cofactor_fp)  # length of each reactant's fingerprint
                    all_lhs_fp = np.concatenate((all_lhs_fp, lhs_cofactor_fp), axis=None)

            # if number of reactants featurized is less than the maximum number of species
            # then pad the LHS fingerprint with extra zeroes (dummy reactants)
            if reactant_counter < max_species:
                lhs_diff = max_species - reactant_counter
                dummy_fp_for_lhs = np.zeros(num_features)  # create a dummy fingerprint of the same length

                # add as many dummy fingerprints as the difference in number of reactants and max species
                for i in range(0, lhs_diff):
                    all_lhs_fp = np.concatenate((all_lhs_fp, dummy_fp_for_lhs), axis=None)

            # Featurize all products on the RHS of the reaction string - begin with products then do cofactors
            if products_list:
                
                for product_smiles in reordered_products:
                    product_counter += 1
                    product_object = compound(product_smiles)
                    product_fp = product_object.smiles_2_fp(fp_type = fp_type,
                                                            radius = radius,
                                                            nBits = nBits,
                                                            is_folded = is_folded,
                                                            dim = dim)

                    num_features = len(product_fp)  # length of each reactant's fingerprint
                    all_rhs_fp = np.concatenate((all_rhs_fp, product_fp), axis=None)

            # then repeat for cofactors on the RHS
            if rhs_cofactors_list:
                
                for rhs_cofactor_smiles in rhs_cofactors:
                    product_counter += 1
                    rhs_cofactor_object = compound(rhs_cofactor_smiles)
                    rhs_cofactor_fp = rhs_cofactor_object.smiles_2_fp(fp_type = fp_type,
                                                                      radius = radius,
                                                                      nBits = nBits,
                                                                      is_folded = is_folded,
                                                                      dim = dim)

                    num_features = len(rhs_cofactor_fp)  # length of each reactant's fingerprint
                    all_rhs_fp = np.concatenate((all_rhs_fp, rhs_cofactor_fp), axis=None)

            # if number of products featurized is less than the maximum number of species
            # then pad the RHS fingerprint with extra zeroes (dummy products)
            if product_counter < max_species:
                rhs_diff = max_species - product_counter
                dummy_fp_for_rhs = np.zeros(num_features)

                # add as many dummy fingerprints as the difference in number of products and max species
                for i in range(0, rhs_diff):
                    all_rhs_fp = np.concatenate((all_rhs_fp, dummy_fp_for_rhs), axis=None)
            
            ### Finally, concatenate fingerprints on both sides of the reaction for a full reaction fingerprint
            reaction_fp = np.concatenate((all_lhs_fp, all_rhs_fp), axis=None)
            
            return reaction_fp

        if cofactor_positioning == "add_concat":
            # generate ecfp4 fingerprints for all species
            # add up all reactants' fingerprints elementwise into a single vector
            # add up all products' fingerprints elementwise into a single vector
            # concatenate the resultant reactants' and products' fingerprints side by side
            # this approach does not require any padding

            if fp_type == "morgan" or fp_type == "ecfp4":
                fp_length = 2048

            if fp_type == "atom_pair":
                fp_length = 2048

            if fp_type == "MACCS":
                fp_length = 167

            if fp_type == "mordred":
                fp_length = 1613

            if fp_type == "MHFP4":
                fp_length = 2048

            # initialize a fingerprint for all reactants on the LHS (incl. cofactors)
            lhs_fp = np.zeros(fp_length)

            # initialize a fingerprint for all products on the RHS (incl. cofactors)
            rhs_fp = np.zeros(fp_length)

            for substrate_smiles in substrates_list:
                compound_object = compound(substrate_smiles)
                fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)
                lhs_fp = lhs_fp + fp

            for lhs_cofactor_smiles in lhs_cofactors_list:
                compound_object = compound(lhs_cofactor_smiles)
                fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)
                lhs_fp = lhs_fp + fp

            for product_smiles in products_list:
                compound_object = compound(product_smiles)
                fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)
                rhs_fp = rhs_fp + fp

            for rhs_cofactor_smiles in rhs_cofactors_list:
                compound_object = compound(rhs_cofactor_smiles)
                fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)
                rhs_fp = rhs_fp + fp

            ### Finally, concatenate fingerprints on both sides of the reaction for a full reaction fingerprint
            reaction_fp = np.concatenate((lhs_fp, rhs_fp), axis=None)

            return reaction_fp

        if cofactor_positioning == "half_random":
            # randomly rearrange reactants anywhere along positions 1 to 4 of the final reaction fingerprint
            # randomly rearrange products anywhere along positions 5 to 8 for the final reaction fingerprint
            # concatenate the two
            # this approach will still require padding

            # start by initializing a vector with 16384 zeros
            rxn_fp = np.zeros(16384)

            # randomly select the "slots" at which to slot in all reactants (i.e. substrates as well as LHS cofactors)
            # any species on the LHS of a reaction will occupy slots 1 through 4 (so either slot 1,2,3, or 4)
            random_positions_for_lhs_species = np.random.permutation(4)

            # randomly, select the "slots" at which to slot in all products (i.e. products as well as RHS cofactors)
            # any species on the RHS of a reaction will occupy slots 5 through 8 (so either slot 5,6,7, or 8)
            random_positions_for_rhs_species = 4 + np.random.permutation(4)

            for substrate_smiles in substrates_list:
                reactant_counter += 1 # keep track of substrate
                compound_object = compound(substrate_smiles)
                substrate_fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)

                # grab a position for the substrate from the previously initialized random positions for reactants
                substrate_position = random_positions_for_lhs_species[0]

                # remove this position from the list of positions since it is no longer available
                random_positions_for_lhs_species = np.delete(random_positions_for_lhs_species,0)

                # slot this substrate fingerprint into its allocated position
                # if this substrate fingerprint was allocated the very first slot, then just add it to the left
                if substrate_position == 0:
                    rxn_fp[:2048] = substrate_fp
                    assert np.array_equal(rxn_fp[:2048], substrate_fp)

                # but if this substrate fingerprint was allocated either slots 2, 3, or 4
                # then a bit more work needs to be done to insert this fingerprint at the correct position
                else:
                    rxn_fp[ (substrate_position-1) * 2048 : substrate_position * 2048 ] = substrate_fp
                    assert np.array_equal(rxn_fp[ (substrate_position-1) * 2048 : substrate_position * 2048 ],
                                          substrate_fp)

            for lhs_cofactor_smiles in lhs_cofactors_list:
                reactant_counter += 1 # keep track of lhs_cofactor
                compound_object = compound(lhs_cofactor_smiles)
                lhs_cofactor_fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)

                # grab a position for the lhs_cofactor from the previously initialized random positions for reactants
                lhs_cofactor_position = random_positions_for_lhs_species[0]

                # remove this position from the list of positions since it is no longer available
                random_positions_for_lhs_species = np.delete(random_positions_for_lhs_species,0)

                # slot this lhs_cofactor fingerprint into its allocated position
                # if this cofactor fingerprint was allocated the very first slot, then just add it to the left
                if lhs_cofactor_position == 0:
                    rxn_fp[:2048] = lhs_cofactor_fp
                    assert np.array_equal(rxn_fp[:2048], lhs_cofactor_fp)

                # but if this cofactor fingerprint was allocated either slots 2, 3, or 4
                # then a bit more work needs to be done to insert this fingerprint at the correct position
                else:
                    rxn_fp[ (lhs_cofactor_position-1) * 2048 : lhs_cofactor_position * 2048 ] = lhs_cofactor_fp
                    assert np.array_equal( rxn_fp[ (lhs_cofactor_position-1) * 2048 : lhs_cofactor_position * 2048 ],
                                           lhs_cofactor_fp)

            for product_smiles in products_list:
                product_counter += 1 # keep track of product
                compound_object = compound(product_smiles)
                product_fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)

                # grab a position for the product
                product_position = random_positions_for_rhs_species[0]

                # remove this position from the list of positions since it is no longer available
                random_positions_for_rhs_species = np.delete(random_positions_for_rhs_species,0)

                # slot the product fingerprint into its allocated position
                rxn_fp[(product_position - 1) * 2048: product_position * 2048] = product_fp
                assert np.array_equal(rxn_fp[(product_position - 1) * 2048: product_position * 2048],
                                      product_fp)

            for rhs_cofactor_smiles in rhs_cofactors_list:
                product_counter += 1 # keep track of rhs_cofactor
                compound_object = compound(rhs_cofactor_smiles)
                rhs_cofactor_fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)

                # grab a position for the rhs_cofactor
                rhs_cofactor_position = random_positions_for_rhs_species[0]

                # remove this position from the list of positions since it is no longer available
                random_positions_for_rhs_species = np.delete(random_positions_for_rhs_species,0)

                # slot the rhs fingerprint into its allocated position
                rxn_fp[(rhs_cofactor_position - 1) * 2048: rhs_cofactor_position * 2048] = rhs_cofactor_fp
                assert np.array_equal(rxn_fp[(rhs_cofactor_position - 1) * 2048: rhs_cofactor_position * 2048],
                                      rhs_cofactor_fp)

            return rxn_fp

        if cofactor_positioning == "full_random":
            # Here, we randomly arrange fingerprints at any position along a reaction feature vector

            rxn_fp = np.zeros(16384)

            random_positions_for_all_species = np.random.permutation(8)

            for substrate_smiles in substrates_list:
                reactant_counter += 1 # keep track of substrate
                compound_object = compound(substrate_smiles)
                substrate_fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                         radius = radius,
                                                         nBits = nBits,
                                                         is_folded = is_folded,
                                                         dim = dim)

                # grab a position for the substrate
                substrate_position = random_positions_for_all_species[0]

                # remove this position from the list of positions since it is no longer available
                random_positions_for_all_species = np.delete(random_positions_for_all_species,0)

                # slot the substrate fingerprint into its allocated position
                if substrate_position == 0:
                    rxn_fp[:2048] = substrate_fp
                    assert np.array_equal(rxn_fp[:2048], substrate_fp)

                else:
                    rxn_fp[ (substrate_position-1) * 2048 : substrate_position * 2048 ] = substrate_fp
                    assert np.array_equal(rxn_fp[ (substrate_position-1) * 2048 : substrate_position * 2048 ],
                                          substrate_fp)

            for lhs_cofactor_smiles in lhs_cofactors_list:
                reactant_counter += 1 # keep track of lhs_cofactor
                compound_object = compound(lhs_cofactor_smiles)
                lhs_cofactor_fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)

                # grab a position for the lhs_cofactor
                lhs_cofactor_position = random_positions_for_all_species[0]

                # remove this position from the list of positions since it is no longer available
                random_positions_for_all_species = np.delete(random_positions_for_all_species,0)

                # slot the lhs_cofactor fingerprint into its allocated position
                if lhs_cofactor_position == 0:
                    rxn_fp[:2048] = lhs_cofactor_fp
                    assert np.array_equal(rxn_fp[:2048], lhs_cofactor_fp)

                else:
                    rxn_fp[ (lhs_cofactor_position-1) * 2048 : lhs_cofactor_position * 2048 ] = lhs_cofactor_fp
                    assert np.array_equal( rxn_fp[ (lhs_cofactor_position-1) * 2048 : lhs_cofactor_position * 2048 ],
                                           lhs_cofactor_fp)

            for product_smiles in products_list:
                product_counter += 1 # keep track of rhs_cofactor
                compound_object = compound(product_smiles)
                product_fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)

                # grab a position for the lhs_cofactor
                product_position = random_positions_for_all_species[0]

                # remove this position from the list of positions since it is no longer available
                random_positions_for_all_species = np.delete(random_positions_for_all_species,0)

                # slot the lhs_cofactor fingerprint into its allocated position
                if product_position == 0:
                    rxn_fp[:2048] = product_fp
                    assert np.array_equal(rxn_fp[:2048], product_fp)

                else:
                    rxn_fp[ (product_position-1) * 2048 : product_position * 2048 ] = product_fp
                    assert np.array_equal( rxn_fp[ (product_position-1) * 2048 : product_position * 2048 ], product_fp)

            for rhs_cofactor_smiles in rhs_cofactors_list:
                product_counter += 1 # keep track of rhs_cofactor
                compound_object = compound(rhs_cofactor_smiles)
                rhs_cofactor_fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)

                # grab a position for the lhs_cofactor
                rhs_cofactor_position = random_positions_for_all_species[0]

                # remove this position from the list of positions since it is no longer available
                random_positions_for_all_species = np.delete(random_positions_for_all_species,0)

                # slot the lhs_cofactor fingerprint into its allocated position
                if rhs_cofactor_position == 0:
                    rxn_fp[:2048] = rhs_cofactor_fp
                    assert np.array_equal(rxn_fp[:2048], rhs_cofactor_fp)

                else:
                    rxn_fp[ (rhs_cofactor_position-1) * 2048 : rhs_cofactor_position * 2048 ] = rhs_cofactor_fp
                    assert np.array_equal( rxn_fp[ (rhs_cofactor_position-1) * 2048 : rhs_cofactor_position * 2048 ],rhs_cofactor_fp)

            return rxn_fp

        if cofactor_positioning == "add_subtract":
            # generate ecfp4 fingerprints for all species
            # add up all reactants' fingerprints elementwise into a single vector
            # add up all products' fingerprints elementwise into a single vector
            # subtract the resultant reactants' fingerprint from the products' fingerprint
            # this approach does not require any padding

            # initialize a fingerprint for all reactants on the LHS (incl. cofactors)
            lhs_fp = np.zeros(2048)

            # initialize a fingerprint for all products on the RHS (incl. cofactors)
            rhs_fp = np.zeros(2048)

            for substrate_smiles in substrates_list:
                compound_object = compound(substrate_smiles)
                fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)
                lhs_fp = lhs_fp + fp

            for lhs_cofactor_smiles in lhs_cofactors_list:
                compound_object = compound(lhs_cofactor_smiles)
                fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)
                lhs_fp = lhs_fp + fp

            for product_smiles in products_list:
                compound_object = compound(product_smiles)
                fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)
                rhs_fp = rhs_fp + fp

            for rhs_cofactor_smiles in rhs_cofactors_list:
                compound_object = compound(rhs_cofactor_smiles)
                fp = compound_object.smiles_2_fp(fp_type = fp_type,
                                                 radius = radius,
                                                 nBits = nBits,
                                                 is_folded = is_folded,
                                                 dim = dim)
                rhs_fp = rhs_fp + fp

            ### Finally, subtract the substrates' fingerprints (lhs_fp) from the products' fingerprints (rhs_fp)
            reaction_fp = rhs_fp - lhs_fp

            return reaction_fp
