""""
A series of test to confirm the functionality of our custom package DORA_XGB/featurizations/featurizations.py
I have excluded tests here on using mordred descriptors because older version of mordred and networkx conflict
This conflict arises due to networkx 2.2's use of 'np.int', which has been deprecated by numpy and is now just 'int'
If you want to use mordred descriptors but run into conflicts with networkx, consider manually changing np.int to int
"""

from DORA_XGB.featurizations import featurizations
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, DataStructs, MACCSkeys
from sklearn.preprocessing import StandardScaler
from mhfp.encoder import MHFPEncoder
import pandas as pd

#TODO - import and write tests for mordred descriptors

# --------- Tests with the compound class of our custom featurizations package  ---------

### Test if the original fingerprinting functions are working independently of this custom package

def test_if_atom_pair_fp_working():
    # Define a known molecule, for example, benzene
    benzene_smiles = 'c1ccccc1'
    benzene_mol = Chem.MolFromSmiles(benzene_smiles)

    # Generate the Atom Pair fingerprint for benzene as a hashed bit vector
    fingerprint_size = 2048  # typically a power of 2, e.g. 1024, 2048, etc.
    benzene_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(benzene_mol, nBits=fingerprint_size)

    # Convert fingerprint to numpy array for easy comparison
    benzene_fp_array = np.zeros((1,), dtype=int)
    DataStructs.ConvertToNumpyArray(benzene_fp, benzene_fp_array)

    # As a simple test, confirm fingerprint is not all zeros
    assert benzene_fp_array.any()

def test_if_ecfp_working():
    benzene_smiles = 'c1ccccc1'
    mol = Chem.MolFromSmiles(benzene_smiles)
    benzene_fp = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048))
    assert len(benzene_fp) == 2048, "The length of the ECFP4 fingerprint should be 2048 bits."

def test_if_MACCS_descriptors_working():
    benzene_smiles = 'c1ccccc1'

    # Generate the MACCS fingerprint for benzene
    mol = Chem.MolFromSmiles(benzene_smiles)
    benzene_fp = MACCSkeys.GenMACCSKeys(mol)

    # MACCS keys are a set of 166 bits, but the RDKit implementation adds an extra zero bit at the beginning making the length 167
    assert len(benzene_fp) == 167, "The length of the MACCS fingerprint should be 167 bits."

# def test_if_MAP4_descriptors_working():
#     # SMILES for water and ethanol
#     smiles_water = "O"
#     smiles_ethanol = "CCO"
#
#     # Generate the MAP4 fingerprint
#     MAP4_f = MAP4Calculator(dimensions=2048, is_folded=True) #folded
#
#     fp_water = MAP4_f.calculate(Chem.MolFromSmiles(smiles_water))
#     fp_ethanol = MAP4_f.calculate(Chem.MolFromSmiles(smiles_ethanol))
#     fp_water_duplicate = MAP4_f.calculate(Chem.MolFromSmiles(smiles_water))
#
#     # Assert that the same molecule gives the same fingerprint
#     assert np.array_equal(fp_water,fp_water_duplicate)
#
#     # Assert that different molecules give different fingerprints
#     assert (np.array_equal(fp_water,fp_ethanol) == False)
#
#     # Optionally, you can add more checks, e.g., the length of the fingerprint, if it's consistent with your settings.
#     expected_length = 2048
#     assert len(fp_water) == expected_length

### Begin testing our custom featurizations package
def test_smiles_canonicalization():
    """
    Tests for the _canonicalize_smiles method in the compound class of the featurization module
    :return:
    """

    cpd_object = featurizations.compound("OCC")
    assert cpd_object._canonicalize_smiles() == "CCO"

    cpd_object = featurizations.compound("OCCC")
    assert cpd_object._canonicalize_smiles() == "CCCO"

def test_remove_stereochemistry():
    """
    Tests for the remove_stereo method in the compound class of the featurization module
    We want to ensure this method removes stereochemistry while canonicalizing a smiles string
    :return:
    """

    # Check if a canonical SMILES string with stereochemistry can have its stereochemistry successfully removed
    cpd_object = featurizations.compound("CC[C@@H](C)C")
    assert cpd_object.remove_stereo() == "CCC(C)C"

    # Check if a canonical SMILES string with stereochemistry can have its stereochemistry successfully removed
    cpd_object = featurizations.compound("CC[C@@H](C)[C@@H](C)C[C@@H](C)C")
    assert cpd_object.remove_stereo() == "CCC(C)C(C)CC(C)C"

    # Check if a non-canonical SMILES with stereochemistry can be canonicalized along with removing its stereochemistry
    cpd_object = featurizations.compound("C[C@](F)(Cl)N")
    assert cpd_object.remove_stereo() == "CC(N)(F)Cl"

    # Check if a non-canonical SMILES with stereochemistry can be canonicalized along with removing its stereochemistry
    cpd_object = featurizations.compound("OC[C@H](O)C(O)=O")
    assert cpd_object.remove_stereo() == "O=C(O)C(O)CO"

# --------- Tests to confirm computing various fingerprints for compounds  ---------
def test_smiles_2_morgan_fp_non_canon_smi_w_stereo():
    """
    Tests for the smiles_2_morgan_fp method in the compound class of the featurization module
    the canonical and non-canonical, stereochemical smiles of a compound should give the same fingerprint
    :param smi_1: input smiles string (non-canonical) for compound
    :param smi_2: input smiles string (canonical) for compound
    :param radius: radius of fragmentation
    :param nBits: output dimensions of fingerprint
    :return:
    """
    smi_1 = "C[C@](F)(Cl)N"
    smi_2 = "CC(N)(F)Cl"
    radius = 2
    nBits = 2048

    cpd_1 = featurizations.compound(smi_1)
    cpd_2 = featurizations.compound(smi_2)
    fp_1 = cpd_1._smiles_2_morganfp(radius=radius, nBits=nBits)
    fp_2 = cpd_2._smiles_2_morganfp(radius=radius, nBits=nBits)

    # check if morgan/ ecfp4 fingerprints generated using the custom featurizations module are in a numpy array
    assert type(fp_1) is np.ndarray
    assert type(fp_2) is np.ndarray

    # Check if morgan fingerprints are of the right length
    assert len(fp_1) == nBits
    assert len(fp_2) == nBits

    # Check if morgan fingerprints from both compounds (canon and non-canon SMILES) are equal
    assert np.array_equal(fp_1, fp_2)

def test_smiles_2_MACCS_fp_for_canonical_smiles():
    """
    Tests for the smiles_2_MACCS method in the compound class of the featurization module

    :param smi: input smiles string (non-canonical) for compound
    :return:
    """
    smi = "CCO"
    cpd = featurizations.compound(smi)
    fp_from_featurizations_module = cpd._smiles_2_MACCS()

    mol = Chem.MolFromSmiles(smi)
    direct_fp = MACCSkeys.GenMACCSKeys(mol)

    # Check if MACCS fingerprints are numpy arrays
    assert type(fp_from_featurizations_module) is np.ndarray

    # Check if MACCS fingerprints are of the right length
    assert len(fp_from_featurizations_module) == 167
    assert len(direct_fp) == 167

    # Check if MACCS fingerprints from both compounds (canon and non-canon SMILES) are equal
    assert np.array_equal(fp_from_featurizations_module, direct_fp)

def test_smiles_2_MACCS_fp_for_non_canon_smiles():
    """
    Tests for the smiles_2_MACCS_fp method in the compound class of the featurization module
    The canonical and non-canonical smiles string of the same compound should give the same fingerprint
    :param smi_1: input smiles string (non-canonical) for compound
    :param smi_2: input smiles string (canonical) for compound
    :param radius: radius of fragmentation
    :param nBits: output dimensions of fingerprint
    :return:
    """
    smi_1 = "OCC"
    smi_2 = "CCO"
    cpd_1 = featurizations.compound(smi_1)
    cpd_2 = featurizations.compound(smi_2)
    fp_1 = cpd_1._smiles_2_MACCS()
    fp_2 = cpd_2._smiles_2_MACCS()

    # check if MACCS fingerprints generated using the custom featurizations module are in a numpy array
    assert type(fp_1) is np.ndarray
    assert type(fp_2) is np.ndarray

    # Check if MACCS fingerprints are of the right length
    assert len(fp_1) == 167
    assert len(fp_2) == 167

    # Check if MACCS fingerprints from both compounds (canon and non-canon SMILES) are equal
    assert np.array_equal(fp_1,fp_2)

def test_smiles_2_MACCS_fp_for_non_canonical_smiles_w_stereochem():
    """
    Tests for the smiles_2_MACCS method in the compound class of the featurization module
    Want to ensure that the canonical and non-canonical, stereochemical smiles of compound give the same fingerprint
    :param smi_1: input smiles string (non-canonical) for compound
    :param smi_2: input smiles string (canonical) for compound
    :param radius: radius of fragmentation
    :param nBits: output dimensions of fingerprint
    :return:
    """
    smi_1 = "C[C@](F)(Cl)N"
    smi_2 = "CC(N)(F)Cl"
    cpd_1 = featurizations.compound(smi_1)
    cpd_2 = featurizations.compound(smi_2)
    fp_1 = cpd_1._smiles_2_MACCS()
    fp_2 = cpd_2._smiles_2_MACCS()

    # check if MACCS fingerprints generated using the custom featurizations module are in a numpy array
    assert type(fp_1) is np.ndarray
    assert type(fp_2) is np.ndarray

    # Check if MACCS fingerprints are of the right length
    assert len(fp_1) == 167
    assert len(fp_2) == 167

    # Check if MACCS fingerprints from both compounds (canon and non-canon SMILES) are equal
    assert np.array_equal(fp_1, fp_2)

def test_smiles_2_atompair_fp_for_canonical_smiles():
    """
    Tests for the smiles_2_atompair method in the compound class of the featurization module

    :param smi: input smiles string (non-canonical) for compound
    :return:
    """
    smi = "CCO"
    cpd = featurizations.compound(smi)
    fp_from_featurizations_module = cpd._smiles_2_atompair(nBits=2048)

    mol = Chem.MolFromSmiles(smi)
    direct_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol,nBits=2048)

    # Check if atom-pair fingerprints are numpy arrays
    assert type(fp_from_featurizations_module) is np.ndarray

    # Check if atom-pair fingerprints are of the right length
    assert len(fp_from_featurizations_module) == 2048
    assert len(direct_fp) == 2048

    # Check if atom-pair fingerprints from both compounds (canon and non-canon SMILES) are equal
    assert np.array_equal(fp_from_featurizations_module, direct_fp)

def test_smiles_2_atompair_fp_for_non_canonical_smiles():
    """
    Tests for the smiles_2_atompair_fp method in the compound class of the featurization module
    Want to ensure that the canonical and non-canonical smiles string of the same compound give the same fingerprint
    :param smi_1: input smiles string (non-canonical) for compound
    :param smi_2: input smiles string (canonical) for compound
    :return:
    """
    smi_1 = "OCC"
    smi_2 = "CCO"
    cpd_1 = featurizations.compound(smi_1)
    cpd_2 = featurizations.compound(smi_2)
    fp_1 = cpd_1._smiles_2_atompair(nBits=2048)
    fp_2 = cpd_2._smiles_2_atompair(nBits=2048)

    # Check if atom pair fingerprints are numpy arrays
    assert type(fp_1) is np.ndarray
    assert type(fp_2) is np.ndarray

    # Check if atom-pair fingerprints are of the right length
    assert len(fp_1) == 2048
    assert len(fp_1) == 2048

    # Check if atom pair fingerprints from both compounds (canon and non-canon SMILES) are equal
    assert np.array_equal(fp_1, fp_2)

def test_smiles_2_atompair_fp_for_non_canon_smi_w_stereo( ):
    """
    Tests for the smiles_2_atompair_fp method in the compound class of the featurization module
    Want to ensure that the canonical and non-canonical, stereochemical smiles of compound give the same fingerprint
    :param smi_1: input smiles string (non-canonical) for compound
    :param smi_2: input smiles string (canonical) for compound
    :param radius: radius of fragmentation
    :param nBits: output dimensions of fingerprint
    :return:
    """
    smi_1 = "C[C@](F)(Cl)N"
    smi_2 = "CC(N)(F)Cl"
    cpd_1 = featurizations.compound(smi_1)
    cpd_2 = featurizations.compound(smi_2)
    fp_1 = cpd_1._smiles_2_atompair(nBits = 2048)
    fp_2 = cpd_2._smiles_2_atompair(nBits = 2048)

    # check if MACCS fingerprints generated using the custom featurizations module are in a numpy array
    assert type(fp_1) is np.ndarray
    assert type(fp_2) is np.ndarray

    # Check if MACCS fingerprints are of the right length
    assert len(fp_1) == 2048
    assert len(fp_2) == 2048

    # Check if MACCS fingerprints from both compounds (canon and non-canon SMILES) are equal
    assert np.array_equal(fp_1,fp_2)

# def test_smiles_2_MAP4_for_canonical_smiles():
#     """
#     Tests for the smiles_2_MAP4_fp method in the compound class of the featurization module
#     Want to ensure that the canonical and non-canonical smiles string of the same compound give the same fingerprint
#     :param smi: input smiles string (non-canonical) for compound
#     :return:
#     """
#     smi = "CCO"
#     is_folded = True
#     dim = 2048
#     cpd = featurizations.compound(smi)
#     fp_from_featurizations_module = cpd._smiles_2_MAP4(is_folded=is_folded, dim=dim)
#
#     # Check if MAP4 fingerprints are numpy arrays
#     assert type(fp_from_featurizations_module) is np.ndarray
#
#     # Check if MAP4 fingerprints are of correct lengths
#     assert len(fp_from_featurizations_module) == dim
#
#     # Generate fingerprint directly
#     mol = Chem.MolFromSmiles("CCO")
#     MAP4_folded = MAP4Calculator(dimensions=dim, is_folded=is_folded)
#     direct_fp = MAP4_folded.calculate(mol)
#
#     # Check if MAP4 fingerprints from both compounds (canon and non-canon SMILES) are equal
#     assert np.array_equal(fp_from_featurizations_module, direct_fp)

# def test_smiles_2_MAP4_for_non_canon_smiles():
#     """
#     Tests for the smiles_2_MAP4_fp method in the compound class of the featurization module
#     Want to ensure that the canonical and non-canonical smiles string of the same compound give the same fingerprint
#     :param smi_1: input smiles string (non-canonical) for compound
#     :param smi_2: input smiles string (canonical) for compound
#     :return:
#     """
#     smi_1 = "OCC"
#     smi_2 = "CCO"
#     is_folded = True
#     dim = 2048
#     cpd_1 = featurizations.compound(smi_1)
#     cpd_2 = featurizations.compound(smi_2)
#     fp_1 = cpd_1._smiles_2_MAP4(is_folded=is_folded, dim=dim)
#     fp_2 = cpd_2._smiles_2_MAP4(is_folded=is_folded, dim=dim)
#
#     # Check if MAP4 fingerprints are numpy arrays
#     assert type(fp_1) is np.ndarray
#     assert type(fp_2) is np.ndarray
#
#     # Check if MAP4 fingerprints are of correct lengths
#     assert len(fp_1) == dim
#     assert len(fp_2) == dim
#
#     # Check if MAP4 fingerprints from both compounds (canon and non-canon SMILES) are equal
#     assert np.array_equal(fp_1, fp_2)

# def test_smiles_2_MAP4_for_non_canonical_smi_w_stereo():
#     """
#     Tests for the smiles_2_atompair_fp method in the compound class of the featurization module
#     Want to ensure that the canonical and non-canonical, stereochemical smiles of compound give the same fingerprint
#     :param smi_1: input smiles string (non-canonical) for compound
#     :param smi_2: input smiles string (canonical) for compound
#     :param radius: radius of fragmentation
#     :param nBits: output dimensions of fingerprint
#     :return:
#     """
#     smi_1 = "C[C@](F)(Cl)N"
#     smi_2 = "CC(N)(F)Cl"
#     is_folded = True
#     dim = 2048
#     cpd_1 = featurizations.compound(smi_1)
#     cpd_2 = featurizations.compound(smi_2)
#     fp_1 = cpd_1._smiles_2_MAP4(is_folded=is_folded, dim=dim)
#     fp_2 = cpd_2._smiles_2_MAP4(is_folded=is_folded, dim=dim)
#
#     # check if MACCS fingerprints generated using the custom featurizations module are in a numpy array
#     assert type(fp_1) is np.ndarray
#     assert type(fp_2) is np.ndarray
#
#     # Check if MACCS fingerprints are of the right length
#     assert len(fp_1) == dim
#     assert len(fp_2) == dim
#
#     fp_1 = np.nan_to_num(fp_1, copy=True)
#     fp_2 = np.nan_to_num(fp_2, copy=True)
#
#     # Check if MAP4 fingerprints from both compounds (canon and non-canon SMILES) are equal
#     assert np.array_equal(fp_1, fp_2)

def test_smiles_2_mhfp_for_canonical_smiles():
    smi = "CCO"
    cpd = featurizations.compound(smi)
    fp_from_featurizations_module = cpd._smiles_2_MHFP(radius = 2)
    assert len(fp_from_featurizations_module) == 2048

    mhfp_encoder = MHFPEncoder()
    direct_fp = mhfp_encoder.encode(in_smiles = smi, # these are the default parameters
                                    radius = 2,
                                    rings = True,
                                    kekulize = True,
                                    sanitize = True)

    assert np.array_equal(fp_from_featurizations_module,direct_fp)

# --------- Tests to confirm that all fingerprinting technique produce unique vectors  ---------
def test_all_fp_methods_and_combined_fp_function():
    smiles = 'NC1=CC=C(C(O)=O)C=C1'
    cpd = featurizations.compound(smiles)

    # first, we calculate fingerprints with the individual fingerprinting methods
    ecfp4_fp_a = cpd._smiles_2_morganfp(radius = 2, nBits = 2048)
    MACCS_fp_a = cpd._smiles_2_MACCS()
    atompair_fp_a = cpd._smiles_2_atompair(nBits = 2048)
    # MAP4_fp_a = cpd._smiles_2_MAP4(is_folded = True, dim = 2048)
    mhfp_a = cpd._smiles_2_MHFP(radius = 2)

    def are_arrays_distinct(arrays):
        n = len(arrays)
        for i in range(n):
            for j in range(i + 1, n):
                if np.array_equal(arrays[i], arrays[j]):
                    return False
        return True

    # now check if fingerprints are unique
    assert are_arrays_distinct([ecfp4_fp_a,
                                MACCS_fp_a,
                                atompair_fp_a,
                                #MAP4_fp_a,
                                mhfp_a])

    # next, we calculate fingerprints with the combined fingerprinting methods
    ecpf4_fp_b = cpd.smiles_2_fp(fp_type = 'ecfp4')
    MACCS_fp_b = cpd.smiles_2_fp(fp_type = 'MACCS')
    atompair_fp_b = cpd.smiles_2_fp(fp_type = 'atom_pair')
    #MAP4_fp_b = cpd.smiles_2_fp(fp_type = 'MAP4')
    mhfp_b = cpd.smiles_2_fp(fp_type = 'MHFP4')

    # again check if fingerprints are unique
    assert are_arrays_distinct([ecpf4_fp_b,
                                MACCS_fp_b,
                                atompair_fp_b,
                                #MAP4_fp_b,
                                mhfp_b])

    # finally, check that fingerprints calculated either through the individual function
    # or the combined function are equal
    assert np.array_equal(ecfp4_fp_a,ecpf4_fp_b)
    assert np.array_equal(MACCS_fp_a,MACCS_fp_b)
    assert np.array_equal(atompair_fp_a,atompair_fp_b)
    #assert np.array_equal(MAP4_fp_a,MAP4_fp_b)
    assert np.array_equal(mhfp_a,mhfp_b)

# --------- Tests with the reaction class of our custom featurizations package  ---------
def test_rxn_2_cpds():
    rxn_str = "O=P(O)(O)O + O=P(O)(O)O + O=c1ccn(C2OC(COP(=O)(O)O)C(O)C2O)c(=O)[nH]1 = O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1 + O + O"
    rxn = featurizations.reaction(rxn_str)
    reactants_str, products_str = rxn._rxn_2_cpds()
    assert reactants_str == "O=P(O)(O)O + O=P(O)(O)O + O=c1ccn(C2OC(COP(=O)(O)O)C(O)C2O)c(=O)[nH]1"
    assert products_str == "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1 + O + O"

    rxn_str = "COc1cc(CC(O)C(=O)O)ccc1O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = COc1cc(CC(=O)C(=O)O)ccc1O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn = featurizations.reaction(rxn_str)
    reactants_str, products_str = rxn._rxn_2_cpds()
    assert reactants_str == "COc1cc(CC(O)C(=O)O)ccc1O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    assert products_str == "COc1cc(CC(=O)C(=O)O)ccc1O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

    rxn_str = "CC(C)=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 + O=O = CC(=O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + O"
    rxn = featurizations.reaction(rxn_str)
    reactants_str, products_str = rxn._rxn_2_cpds()
    assert reactants_str == "CC(C)=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 + O=O"
    assert products_str == "CC(=O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + O"

def test_load_cofactors_set_wo_stereo():
    rxn_str = ""
    rxn = featurizations.reaction(rxn_str)
    cofactors_filepath = '../cofactors/expanded_cofactors_no_stereochem.tsv'
    cofactors_set = rxn.print_full_cofactors_set_wo_stereo(cofactors_filepath=cofactors_filepath)

    # Each of the following stereochemistries were manually removed to check if we got the same result as RDKit
    assert "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O" in cofactors_set
    assert "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)O)C(O)C1O" in cofactors_set
    assert "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)c2cc1C" in cofactors_set
    assert "Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc4c(N)ncnc43)C(O)C1O)c1[nH]c(=O)[nH]c(=O)c1N2" in cofactors_set
    assert "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O" in cofactors_set
    assert "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)O)C(O)C1O" in cofactors_set
    assert "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1" in cofactors_set
    assert "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1" in cofactors_set
    assert "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OS(=O)(=O)O)C(OP(=O)(O)O)C1O" in cofactors_set
    assert "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)O)C(OP(=O)(O)O)C1O" in cofactors_set
    assert "C[S+](CCC(N)C(=O)O)CC1OC(n2cnc3c(N)ncnc32)C(O)C1O" in cofactors_set
    assert "Nc1ncnc2c1ncn2C1OC(CSCCC(N)C(=O)O)C(O)C1O" in cofactors_set
    assert "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)C(O)C2O)c(=O)[nH]1" in cofactors_set
    assert "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1" in cofactors_set
    assert "*c1c(*)c(O)c(*)c(*)c1O" in cofactors_set
    assert "*C1=C(*)C(=O)C(*)=C(*)C1=O" in cofactors_set
    assert "CC(C)=CCOP(=O)(O)OP(=O)(O)O" in cofactors_set
    assert "O=P(O)(O)OP(=O)(O)O" in cofactors_set
    assert "O=C(O)CCC(=O)C(=O)O" in cofactors_set
    assert "NC(CCC(=O)O)C(=O)O" in cofactors_set
    assert "*C(=O)CCC(NC(=O)c1ccc(N(C=O)CC2CNc3nc(N)[nH]c(=O)c3N2)cc1)C(=O)O" in cofactors_set
    assert "*C(=O)CCC(NC(=O)c1ccc(NCC2CNc3nc(N)[nH]c(=O)c3N2)cc1)C(=O)O" in cofactors_set
    assert "[O]C1=C(O)C(=O)OC1C(O)CO" in cofactors_set
    assert "O=C1OC(C(O)CO)C(O)=C1O" in cofactors_set
    assert "*C(=O)C(C)OP(=O)(O)OCC(O)C(O)C(O)Cn1c2nc(=O)[nH]c(=O)c-2cc2ccc(O)cc21" in cofactors_set
    assert "*C(=O)C(C)OP(=O)(O)OCC(O)C(O)C(O)CN1c2cc(O)ccc2Cc2c1[nH]c(=O)[nH]c2=O" in cofactors_set
    assert "CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O" in cofactors_set
    assert "C#[O+]" in cofactors_set
    assert "O=C=O" in cofactors_set
    assert "O=C(O)O" in cofactors_set
    assert  "CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCS" in cofactors_set
    assert "[H+]" in cofactors_set
    assert "OO" in cofactors_set
    assert "Br" in cofactors_set
    assert "C#N" in cofactors_set
    assert "Cl" in cofactors_set
    assert "F" in cofactors_set
    assert "I" in cofactors_set
    assert "N" in cofactors_set
    assert "O=O" in cofactors_set
    assert "O=P(O)(O)OP(=O)(O)O" in cofactors_set
    assert "O=P(O)(O)O" in cofactors_set
    assert "O=S(=O)(O)O" in cofactors_set
    assert "O=S(O)O" in cofactors_set
    assert "O" in cofactors_set
    assert "O=O" in cofactors_set

all_cofactors_wo_stereo_filepath = '../cofactors/expanded_cofactors_no_stereochem.tsv'
cofactors_df = pd.read_csv(all_cofactors_wo_stereo_filepath, delimiter=',')
all_cofactors_wo_stereo = set(cofactors_df['SMILES'])

# test extracting substrates from a reaction string
def test_get_substrates_frm_rxn_1():
    rxn_str = "O=P(O)(O)O + O=P(O)(O)O + O=c1ccn(C2OC(COP(=O)(O)O)C(O)C2O)c(=O)[nH]1 = O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1 + O + O"
    rxn = featurizations.reaction(rxn_str)
    substrates = rxn.get_substrates(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert substrates == ["O=c1ccn(C2OC(COP(=O)(O)O)C(O)C2O)c(=O)[nH]1"]

def test_get_substrates_frm_rxn_2():
    rxn_str = "COc1cc(CC(O)C(=O)O)ccc1O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = COc1cc(CC(=O)C(=O)O)ccc1O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn = featurizations.reaction(rxn_str)
    substrates = rxn.get_substrates(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert substrates == ["COc1cc(CC(O)C(=O)O)ccc1O"]

def test_get_substrates_frm_rxn_3():
    rxn_str = "CC(C)=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 + O=O = CC(=O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + O"
    rxn = featurizations.reaction(rxn_str)
    substrates = rxn.get_substrates(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert substrates == ["CC(C)=O"]

def test_get_substrates_frm_rxn_4():
    rxn_str = "CCC + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 + O=O = CC(C)O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + O"
    rxn = featurizations.reaction(rxn_str)
    substrates = rxn.get_substrates(all_cofactors_wo_stereo=all_cofactors_wo_stereo)
    assert substrates == ["CCC"]

def test_get_substrates_frm_rxn_5():
    rxn_str = "O=P(O)(O)OC1OC(CO)C(O)C(O)C1O + *OC1OC(CO)C(O)C(O)C1O = *OC1OC(COP(=O)(O)OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O + O"
    rxn = featurizations.reaction(rxn_str)
    substrates = rxn.get_substrates(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert substrates == ["O=P(O)(O)OC1OC(CO)C(O)C(O)C1O","*OC1OC(CO)C(O)C(O)C1O"]

def test_get_substrates_frm_rxn_6():
    rxn_str = "CCOP(=O)(OCC)Oc1ccc([N+](=O)O)cc1 + O = CCOP(=O)(O)OCC + O=[N+](O)c1ccc(O)cc1"
    rxn = featurizations.reaction(rxn_str)
    substrates = rxn.get_substrates(all_cofactors_wo_stereo=all_cofactors_wo_stereo)
    assert substrates == ["CCOP(=O)(OCC)Oc1ccc([N+](=O)O)cc1"]

# test extracting products from a reaction string
def test_get_products_frm_rxn_1():
    rxn_str = "O=P(O)(O)O + O=P(O)(O)O + O=c1ccn(C2OC(COP(=O)(O)O)C(O)C2O)c(=O)[nH]1 = O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1 + O + O"
    rxn = featurizations.reaction(rxn_str)
    products = rxn.get_products(all_cofactors_wo_stereo=all_cofactors_wo_stereo)
    assert products == ["O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"]

def test_get_products_frm_rxn_2():
    rxn_str = "COc1cc(CC(O)C(=O)O)ccc1O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = COc1cc(CC(=O)C(=O)O)ccc1O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn = featurizations.reaction(rxn_str)
    products = rxn.get_products(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert products == ["COc1cc(CC(=O)C(=O)O)ccc1O"]

def test_get_products_frm_rxn_3():
    rxn_str = "CC(C)=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 + O=O = CC(=O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + O"
    rxn = featurizations.reaction(rxn_str)
    products = rxn.get_products(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert products == ["CC(=O)CO"]

def test_get_products_frm_rxn_4():
    rxn_str = "CCC + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 + O=O = CC(C)O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + O"
    rxn = featurizations.reaction(rxn_str)
    products = rxn.get_products(all_cofactors_wo_stereo=all_cofactors_wo_stereo)
    assert products == ["CC(C)O"]

def test_get_products_frm_rxn_5():
    rxn_str = "O=P(O)(O)OC1OC(CO)C(O)C(O)C1O + *OC1OC(CO)C(O)C(O)C1O = *OC1OC(COP(=O)(O)OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O + O"
    rxn = featurizations.reaction(rxn_str)
    products = rxn.get_products(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert products == ["*OC1OC(COP(=O)(O)OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O"]

def test_get_products_frm_rxn_6():
    rxn_str = "CCOP(=O)(OCC)Oc1ccc([N+](=O)O)cc1 + O = CCOP(=O)(O)OCC + O=[N+](O)c1ccc(O)cc1"
    rxn = featurizations.reaction(rxn_str)
    products = rxn.get_products(all_cofactors_wo_stereo=all_cofactors_wo_stereo)
    assert products == ["CCOP(=O)(O)OCC","O=[N+](O)c1ccc(O)cc1"]

# test extracting lhs cofactors from a reaction string
def test_get_lhs_cofactors_frm_rxn_1():
    rxn_str = "O=P(O)(O)O + O=P(O)(O)O + O=c1ccn(C2OC(COP(=O)(O)O)C(O)C2O)c(=O)[nH]1 = O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1 + O + O"
    rxn = featurizations.reaction(rxn_str)
    lhs_cofactors = rxn.get_lhs_cofactors(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert lhs_cofactors == ["O=P(O)(O)O","O=P(O)(O)O"]

def test_get_lhs_cofactors_frm_rxn_2():
    rxn_str = "COc1cc(CC(O)C(=O)O)ccc1O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = COc1cc(CC(=O)C(=O)O)ccc1O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn = featurizations.reaction(rxn_str)
    lhs_cofactors = rxn.get_lhs_cofactors(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert lhs_cofactors == ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"]

def test_get_lhs_cofactors_frm_rxn_3():
    rxn_str = "CC(C)=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 + O=O = CC(=O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + O"
    rxn = featurizations.reaction(rxn_str)
    lhs_cofactors = rxn.get_lhs_cofactors(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert lhs_cofactors == ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1","O=O"]

def test_get_lhs_cofactors_frm_rxn_4():
    rxn_str = "CCC + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 + O=O = CC(C)O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + O"
    rxn = featurizations.reaction(rxn_str)
    lhs_cofactors = rxn.get_lhs_cofactors(all_cofactors_wo_stereo=all_cofactors_wo_stereo)
    assert lhs_cofactors == ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1","O=O"]

def test_get_lhs_cofactors_frm_rxn_5():
    rxn_str = "O=P(O)(O)OC1OC(CO)C(O)C(O)C1O + *OC1OC(CO)C(O)C(O)C1O = *OC1OC(COP(=O)(O)OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O + O"
    rxn = featurizations.reaction(rxn_str)
    lhs_cofactors = rxn.get_lhs_cofactors(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert lhs_cofactors == []

def test_get_lhs_cofactors_frm_rxn_6():
    rxn_str = "CCOP(=O)(OCC)Oc1ccc([N+](=O)O)cc1 + O = CCOP(=O)(O)OCC + O=[N+](O)c1ccc(O)cc1"
    rxn = featurizations.reaction(rxn_str)
    lhs_cofactors = rxn.get_lhs_cofactors(all_cofactors_wo_stereo=all_cofactors_wo_stereo)
    assert lhs_cofactors == ["O"]

# test extracting rhs cofactors from a reaction string
def test_get_rhs_cofactors_frm_rxn_1():
    rxn_str = "O=P(O)(O)O + O=P(O)(O)O + O=c1ccn(C2OC(COP(=O)(O)O)C(O)C2O)c(=O)[nH]1 = O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1 + O + O"
    rxn = featurizations.reaction(rxn_str)
    rhs_cofactors = rxn.get_rhs_cofactors(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert rhs_cofactors == ["O","O"]

def test_get_rhs_cofactors_frm_rxn_2():
    rxn_str = "COc1cc(CC(O)C(=O)O)ccc1O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = COc1cc(CC(=O)C(=O)O)ccc1O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn = featurizations.reaction(rxn_str)
    rhs_cofactors = rxn.get_rhs_cofactors(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert rhs_cofactors == ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"]

def test_get_rhs_cofactors_frm_rxn_3():
    rxn_str = "CC(C)=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 + O=O = CC(=O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + O"
    rxn = featurizations.reaction(rxn_str)
    rhs_cofactors = rxn.get_rhs_cofactors(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert rhs_cofactors == ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1","O"]

def test_get_rhs_cofactors_frm_rxn_4():
    rxn_str = "CCC + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 + O=O = CC(C)O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + O"
    rxn = featurizations.reaction(rxn_str)
    lhs_cofactors = rxn.get_rhs_cofactors(all_cofactors_wo_stereo=all_cofactors_wo_stereo)
    assert lhs_cofactors == ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1","O"]

def test_get_rhs_cofactors_frm_rxn_5():
    rxn_str = "O=P(O)(O)OC1OC(CO)C(O)C(O)C1O + *OC1OC(CO)C(O)C(O)C1O = *OC1OC(COP(=O)(O)OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O + O"
    rxn = featurizations.reaction(rxn_str)
    rhs_cofactors = rxn.get_rhs_cofactors(all_cofactors_wo_stereo = all_cofactors_wo_stereo)
    assert rhs_cofactors == ["O"]

def test_get_rhs_cofactors_frm_rxn_6():
    rxn_str = "CCOP(=O)(OCC)Oc1ccc([N+](=O)O)cc1 + O = CCOP(=O)(O)OCC + O=[N+](O)c1ccc(O)cc1"
    rxn = featurizations.reaction(rxn_str)
    rhs_cofactors = rxn.get_rhs_cofactors(all_cofactors_wo_stereo=all_cofactors_wo_stereo)
    assert rhs_cofactors == []

# test reordering cofactors, substrates, and products by their molecular weights
def test_reorder_cofactors_by_MW_1():
    rxn_str = ""
    rxn_object = featurizations.reaction(rxn_str)
    cofactors_list = ["O","O=O"]
    assert rxn_object._reorder_cofactors_by_MW(cofactors_list = cofactors_list, ascending = True) == ["O","O=O"]

def test_reorder_cofactors_by_MW_2():
    rxn_str = ""
    rxn_object = featurizations.reaction(rxn_str)
    cofactors_list = ["O", "O=O"]
    assert rxn_object._reorder_cofactors_by_MW(cofactors_list = cofactors_list, ascending = False) == ["O=O","O"]

def test_reorder_cofactors_by_MW_3():
    rxn_str = ""
    rxn_object = featurizations.reaction(rxn_str)
    NADH_smi = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    NADPH_smi = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1"
    cofactors_list = [NADH_smi, NADPH_smi]
    assert rxn_object._reorder_cofactors_by_MW(cofactors_list = cofactors_list, ascending = False) == [NADPH_smi, NADH_smi]

def test_reorder_cofactors_by_MW_4():
    rxn_str = ""
    rxn_object = featurizations.reaction(rxn_str)
    NADH_smi = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    NADPH_smi = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1"
    cofactors_list = [NADH_smi,NADPH_smi]
    assert rxn_object._reorder_cofactors_by_MW(cofactors_list=cofactors_list, ascending = True) == [NADH_smi, NADPH_smi]

def test_reorder_cofactors_by_MW_5():
    rxn_str = ""
    rxn_object = featurizations.reaction(rxn_str)
    water_smi = "O"
    NADH_smi = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    NADPH_smi = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1"
    cofactors_list = [water_smi, NADH_smi, NADPH_smi]
    assert rxn_object._reorder_cofactors_by_MW(cofactors_list = cofactors_list, ascending=False) == [NADPH_smi, NADH_smi,water_smi]

def test_reorder_cofactors_by_MW_6():
    rxn_str = ""
    rxn_object = featurizations.reaction(rxn_str)
    water_smi = "O"
    NADH_smi = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    NADPH_smi = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1"
    cofactors_list = [water_smi, NADH_smi, NADPH_smi]
    assert rxn_object._reorder_cofactors_by_MW(cofactors_list = cofactors_list, ascending=True) == [water_smi, NADH_smi, NADPH_smi]

def test_reorder_substrates_and_products_by_MW_1():
    rxn_str = ""
    rxn_object = featurizations.reaction(rxn_str)
    smiles_list = ["C","CC","CCC"]
    assert rxn_object._reorder_substrates_by_MW(substrates_list = smiles_list, ascending=True) == ["C","CC","CCC"]
    assert rxn_object._reorder_products_by_MW(products_list = smiles_list, ascending=True) == ["C","CC","CCC"]

def test_reorder_substrates_and_products_by_MW_2():
    rxn_str = ""
    rxn_object = featurizations.reaction(rxn_str)
    smiles_list = ["C","CC","CCC"]
    assert rxn_object._reorder_substrates_by_MW(substrates_list=smiles_list, ascending=False) == ["CCC", "CC", "C"]
    assert rxn_object._reorder_products_by_MW(products_list=smiles_list, ascending=False) == ["CCC", "CC", "C"]

# --------- Tests to create reaction fingerprints with our custom featurizations package's reaction class  ---------
def test_rxn_2_fp_an_AdH_rxn_w_ecfp4_max_species_2():
    rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn_object = featurizations.reaction(rxn_str)

    # first, we automatically generate a reaction fingerprint using our featurizations package's rxn_2_fp method
    rxn_fp = rxn_object.rxn_2_fp(type = "ecfp4", max_species = 2)

    # then, we manually create the same reaction fingerprint and check if our rxn_2_fp method got it right
    substrate_fp = featurizations.compound("OCC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = "ecfp4")
    lhs_cofactor_fp = featurizations.compound("NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1").smiles_2_fp(fp_type = "ecfp4")
    product_fp = featurizations.compound("O=CC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = "ecfp4")
    rhs_cofactor_fp = featurizations.compound("NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1").smiles_2_fp(fp_type = "ecfp4")
    assert type(substrate_fp) == np.ndarray
    assert type(lhs_cofactor_fp) == np.ndarray
    assert type(product_fp) == np.ndarray
    assert type(rhs_cofactor_fp) == np.ndarray

    assert len(substrate_fp) == 2048
    assert len(lhs_cofactor_fp) == 2048
    assert len(product_fp) == 2048
    assert len(rhs_cofactor_fp) == 2048

    manual_rxn_fp = np.concatenate((substrate_fp,lhs_cofactor_fp,product_fp,rhs_cofactor_fp),axis=None)

    # finally, we check that the reaction fingerprint generated by our featurizations package equals the manual one
    assert np.array_equal(rxn_fp,manual_rxn_fp)

def test_rxn_2_fp_an_AdH_rxn_w_ecfp4_max_species_4():
    rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn_object = featurizations.reaction(rxn_str)

    # first, we automatically generate a reaction fingerprint using our featurizations package's rxn_2_fp method
    rxn_fp = rxn_object.rxn_2_fp(type = "ecfp4", max_species = 4)

    # then, we manually create the same reaction fingerprint and check if our rxn_2_fp method got it right
    substrate_fp = featurizations.compound("OCC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = "ecfp4")
    lhs_cofactor_fp = featurizations.compound("NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1").smiles_2_fp(fp_type = "ecfp4")
    product_fp = featurizations.compound("O=CC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = "ecfp4")
    rhs_cofactor_fp = featurizations.compound("NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1").smiles_2_fp(fp_type = "ecfp4")

    assert type(substrate_fp) == np.ndarray
    assert type(lhs_cofactor_fp) == np.ndarray
    assert type(product_fp) == np.ndarray
    assert type(rhs_cofactor_fp) == np.ndarray

    assert len(substrate_fp) == 2048
    assert len(lhs_cofactor_fp) == 2048
    assert len(product_fp) == 2048
    assert len(rhs_cofactor_fp) == 2048

    manual_rxn_fp = np.concatenate((substrate_fp,lhs_cofactor_fp,np.zeros(2048),np.zeros(2048),
                                    product_fp,rhs_cofactor_fp,np.zeros(2048),np.zeros(2048)),
                                   axis=None)

    # finally, we check that the reaction fingerprint generated by our featurizations package equals the manual one
    assert np.array_equal(rxn_fp,manual_rxn_fp)

def test_rxn_2_fp_an_AdH_rxn_w_atompair_max_species_2():
    rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn_object = featurizations.reaction(rxn_str)

    fp_type = "atom_pair"
    max_species = 2

    # first, we automatically generate a reaction fingerprint using our featurizations package's rxn_2_fp method
    rxn_fp = rxn_object.rxn_2_fp(type = fp_type, max_species = max_species)

    # then, we manually create the same reaction fingerprint and check if our rxn_2_fp method got it right
    substrate_fp = featurizations.compound("OCC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    lhs_cofactor_fp = featurizations.compound("NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1").smiles_2_fp(fp_type = fp_type)
    product_fp = featurizations.compound("O=CC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    rhs_cofactor_fp = featurizations.compound("NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1").smiles_2_fp(fp_type = fp_type)
    assert type(substrate_fp) == np.ndarray
    assert type(lhs_cofactor_fp) == np.ndarray
    assert type(product_fp) == np.ndarray
    assert type(rhs_cofactor_fp) == np.ndarray

    assert len(substrate_fp) == 2048
    assert len(lhs_cofactor_fp) == 2048
    assert len(product_fp) == 2048
    assert len(rhs_cofactor_fp) == 2048

    manual_rxn_fp = np.concatenate((substrate_fp,lhs_cofactor_fp,product_fp,rhs_cofactor_fp),axis=None)

    # finally, we check that the reaction fingerprint generated by our featurizations package equals the manual one
    assert np.array_equal(rxn_fp,manual_rxn_fp)

def test_rxn_2_fp_an_AdH_rxn_w_atompair_max_species_4():
    rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn_object = featurizations.reaction(rxn_str)

    fp_type = "atom_pair"
    max_species = 4

    # first, we automatically generate a reaction fingerprint using our featurizations package's rxn_2_fp method
    rxn_fp = rxn_object.rxn_2_fp(type = fp_type, max_species = max_species)

    # then, we manually create the same reaction fingerprint and check if our rxn_2_fp method got it right
    substrate_fp = featurizations.compound("OCC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    lhs_cofactor_fp = featurizations.compound("NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1").smiles_2_fp(fp_type = fp_type)
    product_fp = featurizations.compound("O=CC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    rhs_cofactor_fp = featurizations.compound("NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1").smiles_2_fp(fp_type = fp_type)
    assert type(substrate_fp) == np.ndarray
    assert type(lhs_cofactor_fp) == np.ndarray
    assert type(product_fp) == np.ndarray
    assert type(rhs_cofactor_fp) == np.ndarray

    assert len(substrate_fp) == 2048
    assert len(lhs_cofactor_fp) == 2048
    assert len(product_fp) == 2048
    assert len(rhs_cofactor_fp) == 2048

    manual_rxn_fp = np.concatenate((substrate_fp,lhs_cofactor_fp,np.zeros(2048),np.zeros(2048),
                                    product_fp,rhs_cofactor_fp,np.zeros(2048),np.zeros(2048)),axis=None)

    # finally, we check that the reaction fingerprint generated by our featurizations package equals the manual one
    assert np.array_equal(rxn_fp,manual_rxn_fp)

def test_rxn_2_fp_an_AdH_rxn_w_MACCS_max_species_2():
    rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn_object = featurizations.reaction(rxn_str)

    fp_type = "MACCS"
    max_species = 2

    # first, we automatically generate a reaction fingerprint using our featurizations package's rxn_2_fp method
    rxn_fp = rxn_object.rxn_2_fp(type = fp_type, max_species = max_species)

    # then, we manually create the same reaction fingerprint and check if our rxn_2_fp method got it right
    substrate_fp = featurizations.compound("OCC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    lhs_cofactor_fp = featurizations.compound("NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1").smiles_2_fp(fp_type = fp_type)
    product_fp = featurizations.compound("O=CC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    rhs_cofactor_fp = featurizations.compound("NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1").smiles_2_fp(fp_type = fp_type)
    assert type(substrate_fp) == np.ndarray
    assert type(lhs_cofactor_fp) == np.ndarray
    assert type(product_fp) == np.ndarray
    assert type(rhs_cofactor_fp) == np.ndarray

    assert len(substrate_fp) == 167
    assert len(lhs_cofactor_fp) == 167
    assert len(product_fp) == 167
    assert len(rhs_cofactor_fp) == 167

    manual_rxn_fp = np.concatenate((substrate_fp,lhs_cofactor_fp,product_fp,rhs_cofactor_fp),axis=None)

    # finally, we check that the reaction fingerprint generated by our featurizations package equals the manual one
    assert np.array_equal(rxn_fp,manual_rxn_fp)

def test_rxn_2_fp_an_AdH_rxn_w_MACCS_max_species_4():
    rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn_object = featurizations.reaction(rxn_str)

    fp_type = "MACCS"
    max_species = 4

    # first, we automatically generate a reaction fingerprint using our featurizations package's rxn_2_fp method
    rxn_fp = rxn_object.rxn_2_fp(type = fp_type, max_species = max_species)

    # then, we manually create the same reaction fingerprint and check if our rxn_2_fp method got it right
    substrate_fp = featurizations.compound("OCC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    lhs_cofactor_fp = featurizations.compound("NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1").smiles_2_fp(fp_type = fp_type)
    product_fp = featurizations.compound("O=CC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    rhs_cofactor_fp = featurizations.compound("NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1").smiles_2_fp(fp_type = fp_type)
    assert type(substrate_fp) == np.ndarray
    assert type(lhs_cofactor_fp) == np.ndarray
    assert type(product_fp) == np.ndarray
    assert type(rhs_cofactor_fp) == np.ndarray

    assert len(substrate_fp) == 167
    assert len(lhs_cofactor_fp) == 167
    assert len(product_fp) == 167
    assert len(rhs_cofactor_fp) == 167

    manual_rxn_fp = np.concatenate((substrate_fp,lhs_cofactor_fp,np.zeros(167),np.zeros(167),
                                    product_fp,rhs_cofactor_fp,np.zeros(167),np.zeros(167)),axis=None)

    # finally, we check that the reaction fingerprint generated by our featurizations package equals the manual one
    assert np.array_equal(rxn_fp,manual_rxn_fp)

# def test_rxn_2_fp_an_AdH_rxn_w_MAP4_max_species_2():
#     rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
#     rxn_object = featurizations.reaction(rxn_str)
#
#     fp_type = "MAP4"
#     max_species = 2
#
#     # first, we automatically generate a reaction fingerprint using our featurizations package's rxn_2_fp method
#     rxn_fp = rxn_object.rxn_2_fp(type = fp_type, max_species = max_species)
#
#     # then, we manually create the same reaction fingerprint and check if our rxn_2_fp method got it right
#     substrate_fp = featurizations.compound("OCC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
#     lhs_cofactor_fp = featurizations.compound("NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1").smiles_2_fp(fp_type = fp_type)
#     product_fp = featurizations.compound("O=CC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
#     rhs_cofactor_fp = featurizations.compound("NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1").smiles_2_fp(fp_type = fp_type)
#     assert type(substrate_fp) == np.ndarray
#     assert type(lhs_cofactor_fp) == np.ndarray
#     assert type(product_fp) == np.ndarray
#     assert type(rhs_cofactor_fp) == np.ndarray
#
#     assert len(substrate_fp) == 2048
#     assert len(lhs_cofactor_fp) == 2048
#     assert len(product_fp) == 2048
#     assert len(rhs_cofactor_fp) == 2048
#
#     manual_rxn_fp = np.concatenate((substrate_fp,lhs_cofactor_fp,product_fp,rhs_cofactor_fp),axis=None)
#
#     # finally, we check that the reaction fingerprint generated by our featurizations package equals the manual one
#     assert np.array_equal(rxn_fp,manual_rxn_fp)

# def test_rxn_2_fp_an_AdH_rxn_w_MAP4_max_species_4():
#     rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
#     rxn_object = featurizations.reaction(rxn_str)
#
#     fp_type = "MAP4"
#     max_species = 4
#
#     # first, we automatically generate a reaction fingerprint using our featurizations package's rxn_2_fp method
#     rxn_fp = rxn_object.rxn_2_fp(type = fp_type, max_species = max_species)
#
#     # then, we manually create the same reaction fingerprint and check if our rxn_2_fp method got it right
#     substrate_fp = featurizations.compound("OCC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
#     lhs_cofactor_fp = featurizations.compound("NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1").smiles_2_fp(fp_type = fp_type)
#     product_fp = featurizations.compound("O=CC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
#     rhs_cofactor_fp = featurizations.compound("NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1").smiles_2_fp(fp_type = fp_type)
#     assert type(substrate_fp) == np.ndarray
#     assert type(lhs_cofactor_fp) == np.ndarray
#     assert type(product_fp) == np.ndarray
#     assert type(rhs_cofactor_fp) == np.ndarray
#
#     assert len(substrate_fp) == 2048
#     assert len(lhs_cofactor_fp) == 2048
#     assert len(product_fp) == 2048
#     assert len(rhs_cofactor_fp) == 2048
#
#     manual_rxn_fp = np.concatenate((substrate_fp,lhs_cofactor_fp,np.zeros(2048),np.zeros(2048),
#                                     product_fp,rhs_cofactor_fp,np.zeros(2048),np.zeros(2048)),axis=None)
#
#     # finally, we check that the reaction fingerprint generated by our featurizations package equals the manual one
#     assert np.array_equal(rxn_fp,manual_rxn_fp)

def test_rxn_2_fp_an_AdH_rxn_w_MHFP4_max_species_2():
    rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn_object = featurizations.reaction(rxn_str)

    fp_type = "MHFP4"
    max_species = 2

    # first, we automatically generate a reaction fingerprint using our featurizations package's rxn_2_fp method
    rxn_fp = rxn_object.rxn_2_fp(type = fp_type, max_species = max_species)

    # then, we manually create the same reaction fingerprint and check if our rxn_2_fp method got it right
    substrate_fp = featurizations.compound("OCC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    lhs_cofactor_fp = featurizations.compound("NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1").smiles_2_fp(fp_type = fp_type)
    product_fp = featurizations.compound("O=CC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    rhs_cofactor_fp = featurizations.compound("NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1").smiles_2_fp(fp_type = fp_type)
    assert type(substrate_fp) == np.ndarray
    assert type(lhs_cofactor_fp) == np.ndarray
    assert type(product_fp) == np.ndarray
    assert type(rhs_cofactor_fp) == np.ndarray

    assert len(substrate_fp) == 2048
    assert len(lhs_cofactor_fp) == 2048
    assert len(product_fp) == 2048
    assert len(rhs_cofactor_fp) == 2048

    manual_rxn_fp = np.concatenate((substrate_fp,lhs_cofactor_fp,product_fp,rhs_cofactor_fp),axis=None)

    # finally, we check that the reaction fingerprint generated by our featurizations package equals the manual one
    assert np.array_equal(rxn_fp,manual_rxn_fp)

def test_rxn_2_fp_an_AdH_rxn_w_MHFP4_max_species_4():
    rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rxn_object = featurizations.reaction(rxn_str)

    fp_type = "MHFP4"
    max_species = 4

    # first, we automatically generate a reaction fingerprint using our featurizations package's rxn_2_fp method
    rxn_fp = rxn_object.rxn_2_fp(type = fp_type, max_species = max_species)

    # then, we manually create the same reaction fingerprint and check if our rxn_2_fp method got it right
    substrate_fp = featurizations.compound("OCC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    lhs_cofactor_fp = featurizations.compound("NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1").smiles_2_fp(fp_type = fp_type)
    product_fp = featurizations.compound("O=CC(O)C(O)C(O)C(O)CO").smiles_2_fp(fp_type = fp_type)
    rhs_cofactor_fp = featurizations.compound("NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1").smiles_2_fp(fp_type = fp_type)
    assert type(substrate_fp) == np.ndarray
    assert type(lhs_cofactor_fp) == np.ndarray
    assert type(product_fp) == np.ndarray
    assert type(rhs_cofactor_fp) == np.ndarray

    assert len(substrate_fp) == 2048
    assert len(lhs_cofactor_fp) == 2048
    assert len(product_fp) == 2048
    assert len(rhs_cofactor_fp) == 2048

    manual_rxn_fp = np.concatenate((substrate_fp,lhs_cofactor_fp,np.zeros(2048),np.zeros(2048),
                                    product_fp,rhs_cofactor_fp,np.zeros(2048),np.zeros(2048)),axis=None)

    # finally, we check that the reaction fingerprint generated by our featurizations package equals the manual one
    assert np.array_equal(rxn_fp,manual_rxn_fp)

# finally, we check if all the reaction fingerprinting methods are generating unique fingerprints
# def test_if_all_rxn_fingerprinting_methods_are_unique():
#     rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
#     rxn_object = featurizations.reaction(rxn_str)
#
#     # first, we automatically generate a reaction fingerprint using our featurizations package's rxn_2_fp method
#     ecfp4_fp = rxn_object.rxn_2_fp(type="ecfp4", max_species=2)
#     atompair_fp = rxn_object.rxn_2_fp(type="atom_pair",max_species=2)
#     MACCS_fp = rxn_object.rxn_2_fp(type="MACCS",max_species=2)
#     MAP4_fp = rxn_object.rxn_2_fp(type="MAP4",max_species=2)
#     MHFP4_fp = rxn_object.rxn_2_fp(type="MHFP4",max_species=2)
#
#     def are_arrays_distinct(arrays):
#         n = len(arrays)
#         for i in range(n):
#             for j in range(i + 1, n):
#                 if np.array_equal(arrays[i], arrays[j]):
#                     return False
#         return True
#
#     # now check if all of these arrays are unique
#     # now check if fingerprints are unique
#     assert are_arrays_distinct([ecfp4_fp,
#                                 atompair_fp,
#                                 MACCS_fp,
#                                 MAP4_fp,
#                                 MHFP4_fp])

# --------- Tests to create reaction fingerprints with ascending MWs on an alcohol dehydrogenase reaction ---------

# the alcohol dehydrogenase reaction that we are using for the following tests
AdH_rxn_str = "OCC(O)C(O)C(O)C(O)CO + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = O=CC(O)C(O)C(O)C(O)CO + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

# Test featurizing an alcohol dehydrogenase reaction to ecfp4 both with and without padding
def test_rxn_2_ecfp4_without_padding_for_an_AdH_rxn():

    # first, use the featurizations module to create a reaction object using this reaction string
    reaction_object = featurizations.reaction(AdH_rxn_str)

    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type = "ecfp4",
        radius = 2,
        nBits = 2048,
        is_folded = True,
        dim = 2048,
        max_species = 2,
        cofactor_positioning = "by_ascending_MW",
        all_cofactors_wo_stereo = all_cofactors_wo_stereo)

    # next, manually generate a morgan/ ecfp4 reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(substrate_smiles), radius=2, nBits=2048)
    lhs_cof_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(lhs_cof_smiles), radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(product_smiles), radius=2, nBits=2048)
    rhs_cof_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(rhs_cof_smiles), radius=2, nBits=2048)

    # All fingerprints should always be assembled in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_without_padding = np.concatenate((substrate_fp, lhs_cof_fp, product_fp, rhs_cof_fp),axis=None)

    # Check the length of the manually assembled fingerprint
    assert (len(manual_rxn_fp_without_padding) == 4 * 2048)  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 4 * 2048  # Since 4 species total without padding

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp_without_padding)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp_without_padding[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    assert np.array_equal(reaction_fp, manual_rxn_fp_without_padding)

def test_rxn_2_ecfp4_with_padding_of_one_for_an_AdH_rxn():

    # first, use the featurizations module to create a reaction object using this reaction string
    reaction_object = featurizations.reaction(AdH_rxn_str)

    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=3, # since we are testing for a padding of one
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    # next, manually generate a morgan/ ecfp4 reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_cof_mol = Chem.MolFromSmiles(lhs_cof_smiles)
    rhs_cof_mol = Chem.MolFromSmiles(rhs_cof_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    lhs_cof_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_cof_mol, radius=2, nBits=2048)
    rhs_cof_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_cof_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)

    # All fingerprints should always be assembled in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_with_padding = np.concatenate((substrate_fp, lhs_cof_fp, np.zeros(2048),
                                                product_fp, rhs_cof_fp, np.zeros(2048)),axis=None)

    # Check the length of the manually assembled fingerprint
    assert (len(manual_rxn_fp_with_padding) == 6 * 2048)  # since 6 species total with padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 6 * 2048  # since 6 species total with padding

    # finally, check if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp_with_padding[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp_with_padding)

def test_rxn_2_ecfp4_with_padding_of_two_for_an_AdH_rxn():

    # first, use the featurizations module to create a reaction object using this reaction string
    reaction_object = featurizations.reaction(AdH_rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,  # since we are testing for a padding of one
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    # next, manually generate a morgan/ ecfp4 reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_cof_mol = Chem.MolFromSmiles(lhs_cof_smiles)
    rhs_cof_mol = Chem.MolFromSmiles(rhs_cof_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    lhs_cof_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_cof_mol, radius=2, nBits=2048)
    rhs_cof_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_cof_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)

    # All fingerprints should always be assembled in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_with_padding = np.concatenate((substrate_fp,lhs_cof_fp,np.zeros(2048),np.zeros(2048),
                                                 product_fp,rhs_cof_fp,np.zeros(2048),np.zeros(2048)),axis=None)

    # Check the length of the manually assembled fingerprint
    assert (len(manual_rxn_fp_with_padding) == 8 * 2048)  # since 8 species total with padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048  # since 8 species total with padding

    # finally, check if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp_with_padding[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp_with_padding)

# Test featurizing an alcohol dehydrogenase reaction to MACCS fingerprints both with and without padding
def test_rxn_2_MACCS_without_padding_for_an_AdH_rxn():

    # first, use the featurizations module to generate a MACCS reaction fingerprint
    reaction_object = featurizations.reaction(AdH_rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="MACCS",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=2,  # since we are testing for a padding of one
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    # next, manually generate a MACCS reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_cof_mol = Chem.MolFromSmiles(lhs_cof_smiles)
    rhs_cof_mol = Chem.MolFromSmiles(rhs_cof_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    substrate_fp = MACCSkeys.GenMACCSKeys(substrate_mol)
    product_fp = MACCSkeys.GenMACCSKeys(product_mol)
    lhs_cof_fp = MACCSkeys.GenMACCSKeys(lhs_cof_mol)
    rhs_cof_fp = MACCSkeys.GenMACCSKeys(rhs_cof_mol)

    # All fingerprints should always be assembled in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_without_padding = np.concatenate(
        (substrate_fp, lhs_cof_fp, product_fp, rhs_cof_fp)
    )

    # finally, check if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp_without_padding[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp_without_padding)

def test_rxn_2_MACCS_with_padding_of_one():

    # first, use the featurizations module to generate a MACCS reaction fingerprint
    reaction_object = featurizations.reaction(AdH_rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="MACCS",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=3,  # since we are testing for a padding of one
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    # next, manually generate a MACCS reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_cof_mol = Chem.MolFromSmiles(lhs_cof_smiles)
    rhs_cof_mol = Chem.MolFromSmiles(rhs_cof_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    substrate_fp = MACCSkeys.GenMACCSKeys(substrate_mol)
    product_fp = MACCSkeys.GenMACCSKeys(product_mol)
    lhs_cof_fp = MACCSkeys.GenMACCSKeys(lhs_cof_mol)
    rhs_cof_fp = MACCSkeys.GenMACCSKeys(rhs_cof_mol)

    # All fingerprints should always be assembled in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_without_padding = np.concatenate(
        (
            substrate_fp,
            lhs_cof_fp,
            np.zeros(len(substrate_fp)),
            product_fp,
            rhs_cof_fp,
            np.zeros(len(product_fp)),
        )
    )

    # finally, check if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp_without_padding[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp_without_padding)

def test_rxn_2_MACCS_with_padding_of_two():
    # first, use the featurizations module to generate a MACCS reaction fingerprint
    reaction_object = featurizations.reaction(AdH_rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="MACCS",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,  # since we are testing for a padding of one
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    # next, manually generate a MACCS reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_cof_mol = Chem.MolFromSmiles(lhs_cof_smiles)
    rhs_cof_mol = Chem.MolFromSmiles(rhs_cof_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    substrate_fp = MACCSkeys.GenMACCSKeys(substrate_mol)
    product_fp = MACCSkeys.GenMACCSKeys(product_mol)
    lhs_cof_fp = MACCSkeys.GenMACCSKeys(lhs_cof_mol)
    rhs_cof_fp = MACCSkeys.GenMACCSKeys(rhs_cof_mol)

    # All fingerprints should always be assembled in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_without_padding = np.concatenate(
        (substrate_fp, lhs_cof_fp, np.zeros(len(substrate_fp)), np.zeros(len(substrate_fp)),
         product_fp, rhs_cof_fp, np.zeros(len(product_fp)), np.zeros(len(product_fp))),axis=None)

    # finally, check if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp_without_padding[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp_without_padding)

# Test featurizing an alcohol dehydrogenase reaction to atom_pair fingerprints both with and without padding
def test_rxn_2_atom_pair_fp_without_padding():

    reaction_object = featurizations.reaction(AdH_rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="atom_pair",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=2,  # since we are testing for a padding of one
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    # next, manually generate a atom-pair reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_cof_mol = Chem.MolFromSmiles(lhs_cof_smiles)
    rhs_cof_mol = Chem.MolFromSmiles(rhs_cof_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    substrate_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(substrate_mol,nBits=2048)
    product_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(product_mol,nBits=2048)
    lhs_cof_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(lhs_cof_mol,nBits=2048)
    rhs_cof_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(rhs_cof_mol,nBits=2048)

    # All fingerprints should always be assembled in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_without_padding = np.concatenate(
        (substrate_fp, lhs_cof_fp, product_fp, rhs_cof_fp)
    )

    # finally, check if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp_without_padding[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp_without_padding)

def test_rxn_2_atom_pair_with_padding_of_one():
    # first, use the featurizations module to generate a MACCS reaction fingerprint
    reaction_object = featurizations.reaction(AdH_rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="atom_pair",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=3,  # since we are testing for a padding of one
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    # next, manually generate a atom-pair reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_cof_mol = Chem.MolFromSmiles(lhs_cof_smiles)
    rhs_cof_mol = Chem.MolFromSmiles(rhs_cof_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    substrate_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(substrate_mol, nBits=2048)
    product_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(product_mol, nBits=2048)
    lhs_cof_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(lhs_cof_mol, nBits=2048)
    rhs_cof_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(rhs_cof_mol, nBits=2048)

    # All fingerprints should always be assembled in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_without_padding = np.concatenate(
        (
            substrate_fp,
            lhs_cof_fp,
            np.zeros(len(substrate_fp)),
            product_fp,
            rhs_cof_fp,
            np.zeros(len(product_fp)),
        )
    )

    # finally, check if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp_without_padding[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp_without_padding)

def test_rxn_2_atom_pair_with_padding_of_two():
    # first, use the featurizations module to generate a MACCS reaction fingerprint
    reaction_object = featurizations.reaction(AdH_rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="atom_pair",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,  # since we are testing for a padding of one
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    # next, manually generate a atom-pair reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_cof_mol = Chem.MolFromSmiles(lhs_cof_smiles)
    rhs_cof_mol = Chem.MolFromSmiles(rhs_cof_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    substrate_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(substrate_mol, nBits=2048)
    product_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(product_mol, nBits=2048)
    lhs_cof_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(lhs_cof_mol, nBits=2048)
    rhs_cof_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(rhs_cof_mol, nBits=2048)

    # All fingerprints should always be assembled in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_without_padding = np.concatenate(
        (
            substrate_fp,
            lhs_cof_fp,
            np.zeros(len(substrate_fp)),
            np.zeros(len(substrate_fp)),
            product_fp,
            rhs_cof_fp,
            np.zeros(len(product_fp)),
            np.zeros(len(product_fp)),
        )
    )

    # finally, check if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp_without_padding[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp_without_padding)

# Test featurizing an alcohol dehydrogenase reaction to MAP4 fingerprints both with and without padding
# def test_rxn_2_MAP4_without_padding():
#     # first, use the featurizations module to generate a MACCS reaction fingerprint
#     reaction_object = featurizations.reaction(AdH_rxn_str)
#     reaction_fp = reaction_object.rxn_2_fp_w_positioning(
#         fp_type="MAP4",
#         radius=2,
#         nBits=2048,
#         is_folded=True,
#         dim=2048,
#         max_species=2,  # since we are testing for a padding of two
#         cofactor_positioning="by_ascending_MW",
#         all_cofactors_wo_stereo=all_cofactors_wo_stereo,
#     )
#
#     # next, manually generate a MAP4 reaction fingerprint
#     substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
#     substrate_object = featurizations.compound(substrate_smiles)
#     substrate_modred_fp = substrate_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     substrate_modred_fp = np.nan_to_num(
#         substrate_modred_fp, copy=True
#     )  # remove nan values from modred descriptors
#
#     lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
#     lhs_cof_object = featurizations.compound(lhs_cof_smiles)
#     lhs_cof_modred_fp = lhs_cof_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     lhs_cof_modred_fp = np.nan_to_num(lhs_cof_modred_fp, copy=True)
#
#     product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
#     product_object = featurizations.compound(product_smiles)
#     product_modred_fp = product_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     product_modred_fp = np.nan_to_num(product_modred_fp, copy=True)
#
#     rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
#     rhs_cof_object = featurizations.compound(rhs_cof_smiles)
#     rhs_cof_modred_fp = rhs_cof_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     rhs_cof_modred_fp = np.nan_to_num(rhs_cof_modred_fp, copy=True)
#
#     # All fingerprints should always be assembled in the form of "substrate + cofactor = product + cofactor"
#     manual_rxn_fp_without_padding = np.concatenate(
#         (substrate_modred_fp, lhs_cof_modred_fp, product_modred_fp, rhs_cof_modred_fp)
#     )
#
#     # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
#     assert np.array_equal(reaction_fp, manual_rxn_fp_without_padding)

# def test_rxn_2_MAP4_with_padding_of_one():
#     # first, use the featurizations module to generate a MACCS reaction fingerprint
#     reaction_object = featurizations.reaction(AdH_rxn_str)
#     reaction_fp = reaction_object.rxn_2_fp_w_positioning(
#         fp_type="MAP4",
#         radius=2,
#         nBits=2048,
#         is_folded=True,
#         dim=2048,
#         max_species=3,  # since we are testing for a padding of two
#         cofactor_positioning="by_ascending_MW",
#         all_cofactors_wo_stereo=all_cofactors_wo_stereo)
#
#     # next, manually generate a modred reaction fingerprint
#     substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
#     substrate_object = featurizations.compound(substrate_smiles)
#     substrate_modred_fp = substrate_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     substrate_modred_fp = np.nan_to_num(
#         substrate_modred_fp, copy=True
#     )  # remove nan values from modred descriptors
#
#     lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
#     lhs_cof_object = featurizations.compound(lhs_cof_smiles)
#     lhs_cof_modred_fp = lhs_cof_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     lhs_cof_modred_fp = np.nan_to_num(
#         lhs_cof_modred_fp, copy=True
#     )  # remove nan values from modred descriptors
#
#     product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
#     product_object = featurizations.compound(product_smiles)
#     product_modred_fp = product_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     product_modred_fp = np.nan_to_num(
#         product_modred_fp, copy=True
#     )  # remove nan values from modred descriptors
#
#     rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
#     rhs_cof_object = featurizations.compound(rhs_cof_smiles)
#     rhs_cof_modred_fp = rhs_cof_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     rhs_cof_modred_fp = np.nan_to_num(
#         rhs_cof_modred_fp, copy=True
#     )  # remove nan values from modred descriptors
#
#     # All fingerprints should be in the form of "substrate + cofactor = product + cofactor"
#     manual_rxn_fp_with_padding = np.concatenate(
#         (
#             substrate_modred_fp,
#             lhs_cof_modred_fp,
#             np.zeros(len(substrate_modred_fp)),
#             product_modred_fp,
#             rhs_cof_modred_fp,
#             np.zeros(len(product_modred_fp)),
#         )
#     )
#
#     # finally, check if both methods actually yield the same fingerprint
#     for i in range(0, len(reaction_fp)):
#         element_frm_rxn_fp = int(reaction_fp[i])
#         element_frm_manually_constructed_fp = int(manual_rxn_fp_with_padding[i])
#         assert element_frm_rxn_fp == element_frm_manually_constructed_fp
#
#     # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
#     assert np.array_equal(reaction_fp, manual_rxn_fp_with_padding)

# def test_rxn_2_MAP4_with_padding_of_two():
#     # first, use the featurizations module to generate a MACCS reaction fingerprint
#     reaction_object = featurizations.reaction(AdH_rxn_str)
#     reaction_fp = reaction_object.rxn_2_fp_w_positioning(
#         fp_type="MAP4",
#         radius=2,
#         nBits=2048,
#         is_folded=True,
#         dim=2048,
#         max_species=4,  # since we are testing for a padding of two
#         cofactor_positioning="by_ascending_MW",
#         all_cofactors_wo_stereo=all_cofactors_wo_stereo)
#
#     # next, manually generate a modred reaction fingerprint
#     substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
#     substrate_object = featurizations.compound(substrate_smiles)
#     substrate_modred_fp = substrate_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     substrate_modred_fp = np.nan_to_num(
#         substrate_modred_fp, copy=True
#     )  # remove nan values from modred descriptors
#
#     lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
#     lhs_cof_object = featurizations.compound(lhs_cof_smiles)
#     lhs_cof_modred_fp = lhs_cof_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     lhs_cof_modred_fp = np.nan_to_num(
#         lhs_cof_modred_fp, copy=True
#     )  # remove nan values from modred descriptors
#
#     product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
#     product_object = featurizations.compound(product_smiles)
#     product_modred_fp = product_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     product_modred_fp = np.nan_to_num(
#         product_modred_fp, copy=True
#     )  # remove nan values from modred descriptors
#
#     rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
#     rhs_cof_object = featurizations.compound(rhs_cof_smiles)
#     rhs_cof_modred_fp = rhs_cof_object._smiles_2_MAP4(is_folded=True, dim=2048)
#     rhs_cof_modred_fp = np.nan_to_num(
#         rhs_cof_modred_fp, copy=True
#     )  # remove nan values from modred descriptors
#
#     # All fingerprints should be in the form of "substrate + cofactor = product + cofactor"
#     manual_rxn_fp_with_padding = np.concatenate(
#         (
#             substrate_modred_fp,
#             lhs_cof_modred_fp,
#             np.zeros(len(substrate_modred_fp)),
#             np.zeros(len(substrate_modred_fp)),
#             product_modred_fp,
#             rhs_cof_modred_fp,
#             np.zeros(len(product_modred_fp)),
#             np.zeros(len(product_modred_fp)),
#         )
#     )
#
#     # finally, check if both methods actually yield the same fingerprint
#     for i in range(0, len(reaction_fp)):
#         element_frm_rxn_fp = int(reaction_fp[i])
#         element_frm_manually_constructed_fp = int(manual_rxn_fp_with_padding[i])
#         assert element_frm_rxn_fp == element_frm_manually_constructed_fp
#
#     # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
#     assert np.array_equal(reaction_fp, manual_rxn_fp_with_padding)

# Test featurizing an alcohol dehydrogenase reaction to MHFP4 fingerprints both with and without padding
def test_rxn_2_MHFP4_without_padding():
    # first, use the featurizations module to generate a MACCS reaction fingerprint
    reaction_object = featurizations.reaction(AdH_rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="MHFP4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=2,  # since we are testing for a padding of two
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    # next, manually generate a MHFP reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    substrate_object = featurizations.compound(substrate_smiles)
    substrate_modred_fp = substrate_object._smiles_2_MHFP(radius=2)
    substrate_modred_fp = np.nan_to_num(
        substrate_modred_fp, copy=True
    )  # remove nan values from modred descriptors

    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    lhs_cof_object = featurizations.compound(lhs_cof_smiles)
    lhs_cof_modred_fp = lhs_cof_object._smiles_2_MHFP(radius=2)
    lhs_cof_modred_fp = np.nan_to_num(lhs_cof_modred_fp, copy=True)

    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    product_object = featurizations.compound(product_smiles)
    product_modred_fp = product_object._smiles_2_MHFP(radius=2)
    product_modred_fp = np.nan_to_num(product_modred_fp, copy=True)

    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rhs_cof_object = featurizations.compound(rhs_cof_smiles)
    rhs_cof_modred_fp = rhs_cof_object._smiles_2_MHFP(radius=2)
    rhs_cof_modred_fp = np.nan_to_num(rhs_cof_modred_fp, copy=True)

    # All fingerprints should always be assembled in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_without_padding = np.concatenate(
        (substrate_modred_fp, lhs_cof_modred_fp, product_modred_fp, rhs_cof_modred_fp)
    )

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp_without_padding)

def test_rxn_2_MHFP4_with_padding_of_one():
    # first, use the featurizations module to generate a MACCS reaction fingerprint
    reaction_object = featurizations.reaction(AdH_rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="MHFP4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=3,  # since we are testing for a padding of two
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    # next, manually generate a modred reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    substrate_object = featurizations.compound(substrate_smiles)
    substrate_modred_fp = substrate_object._smiles_2_MHFP(radius=2)
    substrate_modred_fp = np.nan_to_num(
        substrate_modred_fp, copy=True
    )  # remove nan values from modred descriptors

    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    lhs_cof_object = featurizations.compound(lhs_cof_smiles)
    lhs_cof_modred_fp = lhs_cof_object._smiles_2_MHFP(radius=2)
    lhs_cof_modred_fp = np.nan_to_num(
        lhs_cof_modred_fp, copy=True
    )  # remove nan values from modred descriptors

    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    product_object = featurizations.compound(product_smiles)
    product_modred_fp = product_object._smiles_2_MHFP(radius=2)
    product_modred_fp = np.nan_to_num(
        product_modred_fp, copy=True
    )  # remove nan values from modred descriptors

    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rhs_cof_object = featurizations.compound(rhs_cof_smiles)
    rhs_cof_modred_fp = rhs_cof_object._smiles_2_MHFP(radius=2)
    rhs_cof_modred_fp = np.nan_to_num(
        rhs_cof_modred_fp, copy=True
    )  # remove nan values from modred descriptors

    # All fingerprints should be in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_with_padding = np.concatenate(
        (
            substrate_modred_fp,
            lhs_cof_modred_fp,
            np.zeros(len(substrate_modred_fp)),
            product_modred_fp,
            rhs_cof_modred_fp,
            np.zeros(len(product_modred_fp)),
        )
    )

    # finally, check if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp_with_padding[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp_with_padding)

def test_rxn_2_MHFP4_with_padding_of_two():
    # first, use the featurizations module to generate a MACCS reaction fingerprint
    reaction_object = featurizations.reaction(AdH_rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="MHFP4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,  # since we are testing for a padding of two
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    # next, manually generate a modred reaction fingerprint
    substrate_smiles = "OCC(O)C(O)C(O)C(O)CO"
    substrate_object = featurizations.compound(substrate_smiles)
    substrate_modred_fp = substrate_object._smiles_2_MHFP(radius=2)
    substrate_modred_fp = np.nan_to_num(
        substrate_modred_fp, copy=True
    )  # remove nan values from modred descriptors

    lhs_cof_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    lhs_cof_object = featurizations.compound(lhs_cof_smiles)
    lhs_cof_modred_fp = lhs_cof_object._smiles_2_MHFP(radius=2)
    lhs_cof_modred_fp = np.nan_to_num(
        lhs_cof_modred_fp, copy=True
    )  # remove nan values from modred descriptors

    product_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    product_object = featurizations.compound(product_smiles)
    product_modred_fp = product_object._smiles_2_MHFP(radius=2)
    product_modred_fp = np.nan_to_num(
        product_modred_fp, copy=True
    )  # remove nan values from modred descriptors

    rhs_cof_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    rhs_cof_object = featurizations.compound(rhs_cof_smiles)
    rhs_cof_modred_fp = rhs_cof_object._smiles_2_MHFP(radius=2)
    rhs_cof_modred_fp = np.nan_to_num(
        rhs_cof_modred_fp, copy=True
    )  # remove nan values from modred descriptors

    # All fingerprints should be in the form of "substrate + cofactor = product + cofactor"
    manual_rxn_fp_with_padding = np.concatenate(
        (
            substrate_modred_fp,
            lhs_cof_modred_fp,
            np.zeros(len(substrate_modred_fp)),
            np.zeros(len(substrate_modred_fp)),
            product_modred_fp,
            rhs_cof_modred_fp,
            np.zeros(len(product_modred_fp)),
            np.zeros(len(product_modred_fp)),
        )
    )

    # finally, check if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp_with_padding[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp_with_padding)

# --------- Tests to create reaction fingerprints on a monooxygenase reaction in various configurations ---------

# test featurizing a rule0004 reaction by ordering cofactors with ascending molecular weights
# so this reaction should be fingerprinted as: [substrate][O2][NADH][padding][product][Water][NAD+][padding]
def test_rxn_fp_rule_0004_ascending_cofactor_MW_01(rxn_str="Oc1ccccc1 + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 + O=O = Oc1ccccc1O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    substrate_smiles = 'Oc1ccccc1'
    lhs_NADH_smiles = 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'
    lhs_O2_smiles = 'O=O'
    product_smiles = 'Oc1ccccc1O'
    rhs_NAD_plus_smiles = 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'
    rhs_H2O_smiles = 'O'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_NADH_mol = Chem.MolFromSmiles(lhs_NADH_smiles)
    lhs_O2_mol = Chem.MolFromSmiles(lhs_O2_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NAD_plus_mol = Chem.MolFromSmiles(rhs_NAD_plus_smiles)
    rhs_H2O_mol = Chem.MolFromSmiles(rhs_H2O_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    lhs_NADH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NADH_mol, radius=2, nBits=2048)
    lhs_O2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_O2_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NAD_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NAD_plus_mol, radius=2, nBits=2048)
    rhs_H2O_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_H2O_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                     lhs_O2_fp,
                     lhs_NADH_fp,
                     np.zeros(2048),
                     product_fp,
                     rhs_H2O_fp,
                     rhs_NAD_plus_fp,
                     np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# test featurizing a rule0004 reaction by ordering cofactors with ascending molecular weights but reshuffled str
# so despite this reaction being the exact same as the reaction above albeit with species reshuffled along the str
# we would still want this reaction to be fingerprinted as: [substrate][O2][NADH][padding][product][Water][NAD+][padding]
def test_rxn_fp_rule_0004_ascending_cofactor_MW_02(rxn_str="NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 + O=O + Oc1ccccc1 = O + Oc1ccccc1O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    substrate_smiles = 'Oc1ccccc1'
    lhs_NADH_smiles = 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'
    lhs_O2_smiles = 'O=O'
    product_smiles = 'Oc1ccccc1O'
    rhs_NAD_plus_smiles = 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'
    rhs_H2O_smiles = 'O'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_NADH_mol = Chem.MolFromSmiles(lhs_NADH_smiles)
    lhs_O2_mol = Chem.MolFromSmiles(lhs_O2_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NAD_plus_mol = Chem.MolFromSmiles(rhs_NAD_plus_smiles)
    rhs_H2O_mol = Chem.MolFromSmiles(rhs_H2O_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    lhs_NADH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NADH_mol, radius=2, nBits=2048)
    lhs_O2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_O2_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NAD_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NAD_plus_mol, radius=2, nBits=2048)
    rhs_H2O_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_H2O_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                     lhs_O2_fp,
                     lhs_NADH_fp,
                     np.zeros(2048),
                     product_fp,
                     rhs_H2O_fp,
                     rhs_NAD_plus_fp,
                     np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# test featurizing a rule0004 reaction by ordering cofactors with descending molecular weights
# so this reaction should be fingerprinted as: [substrate][NADH][O2][padding][product][NAD+][Water][padding]
def test_rxn_fp_rule_0004_descending_cofactor_MW_01(
        rxn_str="Oc1ccccc1 + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 + O=O = Oc1ccccc1O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_descending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    substrate_smiles = 'Oc1ccccc1'
    lhs_O2_smiles = 'O=O'
    lhs_NADH_smiles = 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'
    product_smiles = 'Oc1ccccc1O'
    rhs_H2O_smiles = 'O'
    rhs_NAD_plus_smiles = 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_NADH_mol = Chem.MolFromSmiles(lhs_NADH_smiles)
    lhs_O2_mol = Chem.MolFromSmiles(lhs_O2_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NAD_plus_mol = Chem.MolFromSmiles(rhs_NAD_plus_smiles)
    rhs_H2O_mol = Chem.MolFromSmiles(rhs_H2O_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    lhs_NADH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NADH_mol, radius=2, nBits=2048)
    lhs_O2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_O2_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NAD_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NAD_plus_mol, radius=2, nBits=2048)
    rhs_H2O_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_H2O_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                     lhs_NADH_fp,
                     lhs_O2_fp,
                     np.zeros(2048),
                     product_fp,
                     rhs_NAD_plus_fp,
                     rhs_H2O_fp,
                     np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# test featurizing a rule0004 reaction by ordering cofactors with descending molecular weights but reshuffled str
# so despite this reaction being the exact same as the reaction above albeit with species reshuffled along the str
# we would still want this reaction fingerprinted as: [substrate][NADH][O2][padding][product][NAD+][Water][padding]
def test_rxn_fp_rule_0004_descending_cofactor_MW_02(
        rxn_str="NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 + Oc1ccccc1 + O=O = Oc1ccccc1O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_descending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    substrate_smiles = 'Oc1ccccc1'
    lhs_O2_smiles = 'O=O'
    lhs_NADH_smiles = 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'
    product_smiles = 'Oc1ccccc1O'
    rhs_H2O_smiles = 'O'
    rhs_NAD_plus_smiles = 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_NADH_mol = Chem.MolFromSmiles(lhs_NADH_smiles)
    lhs_O2_mol = Chem.MolFromSmiles(lhs_O2_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NAD_plus_mol = Chem.MolFromSmiles(rhs_NAD_plus_smiles)
    rhs_H2O_mol = Chem.MolFromSmiles(rhs_H2O_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    lhs_NADH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NADH_mol, radius=2, nBits=2048)
    lhs_O2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_O2_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NAD_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NAD_plus_mol, radius=2, nBits=2048)
    rhs_H2O_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_H2O_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                     lhs_NADH_fp,
                     lhs_O2_fp,
                     np.zeros(2048),
                     product_fp,
                     rhs_NAD_plus_fp,
                     rhs_H2O_fp,
                     np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# --------- Tests to create reaction fingerprints on another alcohol dehydrogenase reaction ---------

# test featurizing a rule0002 reaction by ordering cofactors with ascending molecular weights
# so this reaction should be fingerprinted as: [substrate][NADP+][padding][padding][product][NADPH][padding][padding]
# notice that the reported reaction uses NADPH and NADP+ as cofactors
# rather than NADH and NAD+
# sterochemistry has also been left in here (because stereochemistry is actually not captured in ecfp fingerprints)
def test_rxn_fp_rule_0002_ascending_cofactor_MW_01(
    rxn_str="OCCO + NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1 = O=CCO + NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    substrate_smiles = 'OCCO'
    lhs_NADPH_smiles = 'NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1'
    product_smiles = 'O=CCO'
    rhs_NADP_plus_smiles = 'NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_NADPH_mol = Chem.MolFromSmiles(lhs_NADPH_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NADP_plus_mol = Chem.MolFromSmiles(rhs_NADP_plus_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    lhs_NADPH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NADPH_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NADP_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NADP_plus_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                     lhs_NADPH_fp,
                     np.zeros(2048),
                     np.zeros(2048),
                     product_fp,
                     rhs_NADP_plus_fp,
                     np.zeros(2048),
                     np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# repeat the above but with a jumbled up reaction string to ensure we get the same reaction fingerprint
def test_rxn_fp_rule_0002_ascending_cofactor_MW_02(
    rxn_str="NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1 + OCCO = O=CCO + NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    substrate_smiles = 'OCCO'
    lhs_NADPH_smiles = 'NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1'
    product_smiles = 'O=CCO'
    rhs_NADP_plus_smiles = 'NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_NADPH_mol = Chem.MolFromSmiles(lhs_NADPH_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NADP_plus_mol = Chem.MolFromSmiles(rhs_NADP_plus_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    lhs_NADPH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NADPH_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NADP_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NADP_plus_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                     lhs_NADPH_fp,
                     np.zeros(2048),
                     np.zeros(2048),
                     product_fp,
                     rhs_NADP_plus_fp,
                     np.zeros(2048),
                     np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# test featurizing a rule0002 reaction by ordering cofactors with descending molecular weights
# so this reaction should be fingerprinted as: [substrate][NADP+][padding][padding][product][NADPH][padding][padding]
def test_rxn_fp_rule_0002_descending_cofactor_MW_01(
    rxn_str="OCCO + NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1 = O=CCO + NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_descending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    substrate_smiles = 'OCCO'
    lhs_NADPH_smiles = 'NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1'
    product_smiles = 'O=CCO'
    rhs_NADP_plus_smiles = 'NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_NADPH_mol = Chem.MolFromSmiles(lhs_NADPH_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NADP_plus_mol = Chem.MolFromSmiles(rhs_NADP_plus_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    lhs_NADPH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NADPH_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NADP_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NADP_plus_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                     lhs_NADPH_fp,
                     np.zeros(2048),
                     np.zeros(2048),
                     product_fp,
                     rhs_NADP_plus_fp,
                     np.zeros(2048),
                     np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# repeat the above but with a jumbled up reaction string
def test_rxn_fp_rule_0002_descending_cofactor_MW_02(
    rxn_str="NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1 + OCCO = NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1 + O=CCO"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_descending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    substrate_smiles = 'OCCO'
    lhs_NADPH_smiles = 'NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1'
    product_smiles = 'O=CCO'
    rhs_NADP_plus_smiles = 'NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    lhs_NADPH_mol = Chem.MolFromSmiles(lhs_NADPH_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NADP_plus_mol = Chem.MolFromSmiles(rhs_NADP_plus_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    lhs_NADPH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NADPH_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NADP_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NADP_plus_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                     lhs_NADPH_fp,
                     np.zeros(2048),
                     np.zeros(2048),
                     product_fp,
                     rhs_NADP_plus_fp,
                     np.zeros(2048),
                     np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# --------- Tests to create reaction fingerprints on a rule0026 reaction ---------
def test_rxn_fp_rule_0026_ascending_cofactor_MW_01(
    rxn_str = "CC(=O)C=O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = CC(=O)C(=O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    substrate_smiles = 'CC(=O)C=O'
    water_smiles = 'O'
    lhs_NAD_plus_smiles = 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'
    product_smiles = 'CC(=O)C(=O)O'
    rhs_NADH_smiles = 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    lhs_NAD_plus_mol = Chem.MolFromSmiles(lhs_NAD_plus_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NADH_mol = Chem.MolFromSmiles(rhs_NADH_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    lhs_NAD_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NAD_plus_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NADH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NADH_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                                    water_fp,
                                    lhs_NAD_plus_fp,
                                    np.zeros(2048),
                                    product_fp,
                                    rhs_NADH_fp,
                                    np.zeros(2048),
                                    np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# repeat the above but with a jumbled up reaction string to make sure we get the same fingerprint
def test_rxn_fp_rule_0026_ascending_cofactor_MW_02(
    rxn_str = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + CC(=O)C=O + O = CC(=O)C(=O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    substrate_smiles = 'CC(=O)C=O'
    water_smiles = 'O'
    lhs_NAD_plus_smiles = 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'
    product_smiles = 'CC(=O)C(=O)O'
    rhs_NADH_smiles = 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    lhs_NAD_plus_mol = Chem.MolFromSmiles(lhs_NAD_plus_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NADH_mol = Chem.MolFromSmiles(rhs_NADH_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    lhs_NAD_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NAD_plus_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NADH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NADH_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                                    water_fp,
                                    lhs_NAD_plus_fp,
                                    np.zeros(2048),
                                    product_fp,
                                    rhs_NADH_fp,
                                    np.zeros(2048),
                                    np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# test featurizing a rule0026 reaction by ordering cofactors with descending molecular weights
# so this reaction should be fingerprinted as: [substrate][NAD+][Water][padding][product][NADH][padding][padding]
def test_rxn_fp_rule_0026_descending_cofactor_MW_01(rxn_str = "CC(=O)C=O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = CC(=O)C(=O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_descending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    substrate_smiles = 'CC(=O)C=O'
    water_smiles = 'O'
    lhs_NAD_plus_smiles = 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'
    product_smiles = 'CC(=O)C(=O)O'
    rhs_NADH_smiles = 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    lhs_NAD_plus_mol = Chem.MolFromSmiles(lhs_NAD_plus_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NADH_mol = Chem.MolFromSmiles(rhs_NADH_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    lhs_NAD_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NAD_plus_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NADH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NADH_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                                    lhs_NAD_plus_fp,
                                    water_fp,
                                    np.zeros(2048),
                                    product_fp,
                                    rhs_NADH_fp,
                                    np.zeros(2048),
                                    np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

def test_rxn_fp_rule_0026_descending_cofactor_MW_02(rxn_str = "O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 + CC(=O)C=O = CC(=O)C(=O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_descending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    substrate_smiles = 'CC(=O)C=O'
    water_smiles = 'O'
    lhs_NAD_plus_smiles = 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'
    product_smiles = 'CC(=O)C(=O)O'
    rhs_NADH_smiles = 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    lhs_NAD_plus_mol = Chem.MolFromSmiles(lhs_NAD_plus_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NADH_mol = Chem.MolFromSmiles(rhs_NADH_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    lhs_NAD_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NAD_plus_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NADH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NADH_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                                    lhs_NAD_plus_fp,
                                    water_fp,
                                    np.zeros(2048),
                                    product_fp,
                                    rhs_NADH_fp,
                                    np.zeros(2048),
                                    np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# --------- Tests to create reaction fingerprints a rule0007 reaction ---------
def test_rxn_fp_rule_0007_ascending_cofactor_MW_01(
    rxn_str = "OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O + O = OCC1OC(O)C(O)C(O)C1O + OCC1OC(O)C(O)C(O)C1O"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    substrate_smiles = 'OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O'
    water_smiles = 'O'
    product_1_smiles = 'OCC1OC(O)C(O)C(O)C1O'
    product_2_smiles = 'OCC1OC(O)C(O)C(O)C1O'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    product_1_mol = Chem.MolFromSmiles(product_1_smiles)
    product_2_mol = Chem.MolFromSmiles(product_2_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    product_1_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_1_mol, radius=2, nBits=2048)
    product_2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_2_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                                    water_fp,
                                    np.zeros(2048),
                                    np.zeros(2048),
                                    product_1_fp,
                                    product_2_fp,
                                    np.zeros(2048),
                                    np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# descending MW for the same reaction as above
def test_rxn_fp_rule_0007_descending_cofactor_MW_01(
    rxn_str = "OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O + O = OCC1OC(O)C(O)C(O)C1O + OCC1OC(O)C(O)C(O)C1O"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_descending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    substrate_smiles = 'OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O'
    water_smiles = 'O'
    product_1_smiles = 'OCC1OC(O)C(O)C(O)C1O'
    product_2_smiles = 'OCC1OC(O)C(O)C(O)C1O'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    product_1_mol = Chem.MolFromSmiles(product_1_smiles)
    product_2_mol = Chem.MolFromSmiles(product_2_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    product_1_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_1_mol, radius=2, nBits=2048)
    product_2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_2_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                                    water_fp,
                                    np.zeros(2048),
                                    np.zeros(2048),
                                    product_1_fp,
                                    product_2_fp,
                                    np.zeros(2048),
                                    np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# ascending MW
def test_rxn_fp_rule_0007_ascending_cofactor_MW_02(
    rxn_str = "CC(=O)NC1C(OP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)OC(CO)C(O)C1O + O = O=c1ccn(C2OC(CO)C(O)C2O)c(=O)[nH]1 + CC(=O)NC1C(OP(=O)(O)OP(=O)(O)O)OC(CO)C(O)C1O"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    substrate_smiles = 'CC(=O)NC1C(OP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)OC(CO)C(O)C1O'
    water_smiles = 'O'
    product_1_smiles = 'O=c1ccn(C2OC(CO)C(O)C2O)c(=O)[nH]1'
    product_2_smiles = 'CC(=O)NC1C(OP(=O)(O)OP(=O)(O)O)OC(CO)C(O)C1O'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    product_1_mol = Chem.MolFromSmiles(product_1_smiles)
    product_2_mol = Chem.MolFromSmiles(product_2_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    product_1_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_1_mol, radius=2, nBits=2048)
    product_2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_2_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                                    water_fp,
                                    np.zeros(2048),
                                    np.zeros(2048),
                                    product_1_fp,
                                    product_2_fp,
                                    np.zeros(2048),
                                    np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# descending MW for the same reaction as above
def test_rxn_fp_rule_0007_descending_cofactor_MW_02(
    rxn_str = "CC(=O)NC1C(OP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)OC(CO)C(O)C1O + O = O=c1ccn(C2OC(CO)C(O)C2O)c(=O)[nH]1 + CC(=O)NC1C(OP(=O)(O)OP(=O)(O)O)OC(CO)C(O)C1O"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_descending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    substrate_smiles = 'CC(=O)NC1C(OP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)OC(CO)C(O)C1O'
    water_smiles = 'O'
    product_1_smiles = 'O=c1ccn(C2OC(CO)C(O)C2O)c(=O)[nH]1'
    product_2_smiles = 'CC(=O)NC1C(OP(=O)(O)OP(=O)(O)O)OC(CO)C(O)C1O'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    product_1_mol = Chem.MolFromSmiles(product_1_smiles)
    product_2_mol = Chem.MolFromSmiles(product_2_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    product_1_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_1_mol, radius=2, nBits=2048)
    product_2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_2_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_fp,
                                    water_fp,
                                    np.zeros(2048),
                                    np.zeros(2048),
                                    product_2_fp,
                                    product_1_fp,
                                    np.zeros(2048),
                                    np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# ------------------- Test cofactor positioning for a rule0006 reaction -------------------
def test_rxn_fp_rule_0006_ascending_cofactor_MW_01(
    rxn_str = "OCC1OC(O)C(O)C(O)C1O + OCC1OC(O)C(O)C(O)C1O = OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O + O"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo)

    substrate_1_smiles = 'OCC1OC(O)C(O)C(O)C1O'
    substrate_2_smiles = 'OCC1OC(O)C(O)C(O)C1O'
    product_smiles = 'OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O'
    water_smiles = 'O'

    substrate_1_mol = Chem.MolFromSmiles(substrate_1_smiles)
    substrate_2_mol = Chem.MolFromSmiles(substrate_2_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)

    substrate_1_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_1_mol, radius=2, nBits=2048)
    substrate_2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_2_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_1_fp,
                                    substrate_2_fp,
                                    np.zeros(2048),
                                    np.zeros(2048),
                                    product_fp,
                                    water_fp,
                                    np.zeros(2048),
                                    np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

def test_rxn_fp_rule_0006_descending_cofactor_MW_01(
    rxn_str = "OCC1OC(O)C(O)C(O)C1O + OCC1OC(O)C(O)C(O)C1O = OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O + O"):
    reaction_object = featurizations.reaction(rxn_str)
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type="ecfp4",
        radius=2,
        nBits=2048,
        is_folded=True,
        dim=2048,
        max_species=4,
        cofactor_positioning="by_ascending_MW",
        all_cofactors_wo_stereo=all_cofactors_wo_stereo,
    )

    substrate_1_smiles = 'OCC1OC(O)C(O)C(O)C1O'
    substrate_2_smiles = 'OCC1OC(O)C(O)C(O)C1O'
    product_smiles = 'OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O'
    water_smiles = 'O'

    substrate_1_mol = Chem.MolFromSmiles(substrate_1_smiles)
    substrate_2_mol = Chem.MolFromSmiles(substrate_2_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)

    substrate_1_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_1_mol, radius=2, nBits=2048)
    substrate_2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_2_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)

    manual_rxn_fp = np.concatenate((substrate_1_fp,
                                    substrate_2_fp,
                                    np.zeros(2048),
                                    np.zeros(2048),
                                    product_fp,
                                    water_fp,
                                    np.zeros(2048),
                                    np.zeros(2048)))

    # Check the length of the manually assembled fingerprint
    assert (
            len(manual_rxn_fp) == 8 * 2048
    )  # since 4 species total without padding

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 8 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# ------------------- Test cofactor positioning with add_concat configuration and ecfp fingerprints -------------------
def test_rxn_fp_rule_0026_rxn_with_add_concat_configuration_ecfp_01(
    rxn_str = "CC(=O)C=O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = CC(=O)C(=O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"):

    reaction_object = featurizations.reaction(rxn_str)

    # we use our reaction fingerprinting function to automatically generate a reaction fingerprint
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type = "ecfp4",
        radius = 2,
        nBits = 2048,
        is_folded = True,
        dim = 2048,
        max_species = 4,
        cofactor_positioning = "add_concat", # the add-concat configuration is used here
        all_cofactors_wo_stereo = all_cofactors_wo_stereo)

    # then, we manually create a reaction fingerprint using the add-concat configuration
    substrate_smiles = 'CC(=O)C=O'
    water_smiles = 'O'
    lhs_NAD_plus_smiles = 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'
    product_smiles = 'CC(=O)C(=O)O'
    rhs_NADH_smiles = 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    lhs_NAD_plus_mol = Chem.MolFromSmiles(lhs_NAD_plus_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NADH_mol = Chem.MolFromSmiles(rhs_NADH_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    lhs_NAD_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NAD_plus_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NADH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NADH_mol, radius=2, nBits=2048)

    # we initialize the lhs_fp with 2048 zeros and do the same for the rhs_fp
    lhs_fp = np.zeros(2048)
    rhs_fp = np.zeros(2048)

    # all manually fingerprinted compound fingerprints are summed up amongst reactants and products
    lhs_fp = lhs_fp + substrate_fp + water_fp + lhs_NAD_plus_fp
    rhs_fp = rhs_fp + product_fp + rhs_NADH_fp

    # finally, to create the reaction fingerprint, the sum of reactant fingerprints is concatenated with ...
    # ... the sum of product fingerprints
    manual_rxn_fp = np.concatenate((lhs_fp,rhs_fp),axis=None)

    # Check the length of the manually assembled fingerprint
    assert len(manual_rxn_fp) == 2 * 2048

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 2 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

def test_rxn_fp_rule_0006_rxn_with_add_concat_configuration_ecfp_01(
    rxn_str = "OCC1OC(O)C(O)C(O)C1O + OCC1OC(O)C(O)C(O)C1O = OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O + O"):

    reaction_object = featurizations.reaction(rxn_str)

    # we use our reaction fingerprinting function to automatically generate a reaction fingerprint
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type = "ecfp4",
        radius = 2,
        nBits = 2048,
        is_folded = True,
        dim = 2048,
        max_species = 4,
        cofactor_positioning = "add_concat", # the add-concat configuration is used here
        all_cofactors_wo_stereo = all_cofactors_wo_stereo)

    # then, we manually create a reaction fingerprint using the add-concat configuration
    substrate_1_smiles = 'OCC1OC(O)C(O)C(O)C1O'
    substrate_2_smiles = 'OCC1OC(O)C(O)C(O)C1O'
    product_smiles = 'OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O'
    water_smiles = 'O'

    substrate_1_mol = Chem.MolFromSmiles(substrate_1_smiles)
    substrate_2_mol = Chem.MolFromSmiles(substrate_2_smiles)

    product_mol = Chem.MolFromSmiles(product_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)

    substrate_1_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_1_mol, radius=2, nBits=2048)
    substrate_2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_2_mol, radius=2, nBits=2048)

    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)

    # we initialize the lhs_fp with 2048 zeros and do the same for the rhs_fp
    lhs_fp = np.zeros(2048)
    rhs_fp = np.zeros(2048)

    # all manually fingerprinted compound fingerprints are summed up amongst reactants and products
    lhs_fp = lhs_fp + substrate_1_fp + substrate_2_fp
    rhs_fp = rhs_fp + water_fp + product_fp

    # finally, to create the reaction fingerprint, the sum of reactant fingerprints is concatenated with ...
    # ... the sum of product fingerprints
    manual_rxn_fp = np.concatenate((lhs_fp,rhs_fp),axis=None)

    # Check the length of the manually assembled fingerprint
    assert len(manual_rxn_fp) == 2 * 2048

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 2 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

def test_rxn_fp_rule_0007_rxn_with_add_concat_configuration_ecfp_01(
    rxn_str = "CC(=O)NC1C(OP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)OC(CO)C(O)C1O + O = O=c1ccn(C2OC(CO)C(O)C2O)c(=O)[nH]1 + CC(=O)NC1C(OP(=O)(O)OP(=O)(O)O)OC(CO)C(O)C1O"):

    reaction_object = featurizations.reaction(rxn_str)

    # we use our reaction fingerprinting function to automatically generate a reaction fingerprint
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type = "ecfp4",
        radius = 2,
        nBits = 2048,
        is_folded = True,
        dim = 2048,
        max_species = 4,
        cofactor_positioning = "add_concat", # the add-concat configuration is used here
        all_cofactors_wo_stereo = all_cofactors_wo_stereo)

    # then, we manually create a reaction fingerprint using the add-concat configuration
    substrate_smiles = 'CC(=O)NC1C(OP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)OC(CO)C(O)C1O'
    water_smiles = 'O'
    product_1_smiles = 'O=c1ccn(C2OC(CO)C(O)C2O)c(=O)[nH]1'
    product_2_smiles = 'CC(=O)NC1C(OP(=O)(O)OP(=O)(O)O)OC(CO)C(O)C1O'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    product_1_mol = Chem.MolFromSmiles(product_1_smiles)
    product_2_mol= Chem.MolFromSmiles(product_2_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)

    product_1_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_1_mol, radius=2, nBits=2048)
    product_2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_2_mol, radius=2, nBits=2048)

    # we initialize the lhs_fp with 2048 zeros and do the same for the rhs_fp
    lhs_fp = np.zeros(2048)
    rhs_fp = np.zeros(2048)

    # all manually fingerprinted compound fingerprints are summed up amongst reactants and products
    lhs_fp = lhs_fp + substrate_fp + water_fp
    rhs_fp = rhs_fp + product_1_fp + product_2_fp

    # finally, to create the reaction fingerprint, the sum of reactant fingerprints is concatenated with ...
    # ... the sum of product fingerprints
    manual_rxn_fp = np.concatenate((lhs_fp,rhs_fp),axis=None)

    # Check the length of the manually assembled fingerprint
    assert len(manual_rxn_fp) == 2 * 2048

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 2 * 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

# ------------------- Test cofactor positioning with add_subtract configuration and ecfp fingerprints -------------------
def test_rxn_fp_rule_0026_rxn_with_add_subtract_configuration_ecfp_01(
    rxn_str = "CC(=O)C=O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = CC(=O)C(=O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"):

    reaction_object = featurizations.reaction(rxn_str)

    # we use our reaction fingerprinting function to automatically generate a reaction fingerprint
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type = "ecfp4",
        radius = 2,
        nBits = 2048,
        is_folded = True,
        dim = 2048,
        max_species = 4,
        cofactor_positioning = "add_subtract", # the add-subtract configuration is used here
        all_cofactors_wo_stereo = all_cofactors_wo_stereo)

    # then, we manually create a reaction fingerprint using the add-concat configuration
    substrate_smiles = 'CC(=O)C=O'
    water_smiles = 'O'
    lhs_NAD_plus_smiles = 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'
    product_smiles = 'CC(=O)C(=O)O'
    rhs_NADH_smiles = 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    lhs_NAD_plus_mol = Chem.MolFromSmiles(lhs_NAD_plus_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    rhs_NADH_mol = Chem.MolFromSmiles(rhs_NADH_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    lhs_NAD_plus_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lhs_NAD_plus_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)
    rhs_NADH_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(rhs_NADH_mol, radius=2, nBits=2048)

    # we initialize the lhs_fp with 2048 zeros and do the same for the rhs_fp
    lhs_fp = np.zeros(2048)
    rhs_fp = np.zeros(2048)

    # all manually fingerprinted compound fingerprints are summed up amongst reactants and products
    lhs_fp = lhs_fp + substrate_fp + water_fp + lhs_NAD_plus_fp
    rhs_fp = rhs_fp + product_fp + rhs_NADH_fp

    # finally, to create the reaction fingerprint, the sum of reactant fingerprints is subtracted from ...
    # ... the sum of product fingerprints
    manual_rxn_fp = rhs_fp - lhs_fp

    # Check the length of the manually assembled fingerprint
    # the add-subtract configuration only has 2048 elements since the sum of substrate fingerprints is ..
    # ... subtracted from the sum of product fingerprints
    assert len(manual_rxn_fp) == 2048

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

def test_rxn_fp_rule_0006_rxn_with_add_subtract_configuration_ecfp_01(
    rxn_str = "OCC1OC(O)C(O)C(O)C1O + OCC1OC(O)C(O)C(O)C1O = OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O + O"):

    reaction_object = featurizations.reaction(rxn_str)

    # we use our reaction fingerprinting function to automatically generate a reaction fingerprint
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type = "ecfp4",
        radius = 2,
        nBits = 2048,
        is_folded = True,
        dim = 2048,
        max_species = 4,
        cofactor_positioning = "add_subtract", # the add-subtract configuration is used here
        all_cofactors_wo_stereo = all_cofactors_wo_stereo)

    # then, we manually create a reaction fingerprint using the add-concat configuration
    substrate_1_smiles = 'OCC1OC(O)C(O)C(O)C1O'
    substrate_2_smiles = 'OCC1OC(O)C(O)C(O)C1O'
    product_smiles = 'OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O'
    water_smiles = 'O'

    substrate_1_mol = Chem.MolFromSmiles(substrate_1_smiles)
    substrate_2_mol = Chem.MolFromSmiles(substrate_2_smiles)

    product_mol = Chem.MolFromSmiles(product_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)

    substrate_1_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_1_mol, radius=2, nBits=2048)
    substrate_2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_2_mol, radius=2, nBits=2048)

    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)
    product_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=2, nBits=2048)

    # we initialize the lhs_fp with 2048 zeros and do the same for the rhs_fp
    lhs_fp = np.zeros(2048)
    rhs_fp = np.zeros(2048)

    # all manually fingerprinted compound fingerprints are summed up amongst reactants and products
    lhs_fp = lhs_fp + substrate_1_fp + substrate_2_fp
    rhs_fp = rhs_fp + water_fp + product_fp

    # finally, to create the reaction fingerprint, the sum of reactant fingerprints is subtracted from ...
    # ... the sum of product fingerprints
    manual_rxn_fp = rhs_fp - lhs_fp

    # Check the length of the manually assembled fingerprint
    # the add-subtract configuration only has 2048 elements since the sum of substrate fingerprints is ..
    # ... subtracted from the sum of product fingerprints
    assert len(manual_rxn_fp) == 2048

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)

def test_rxn_fp_rule_0007_rxn_with_add_subtract_configuration_ecfp_01(
    rxn_str = "CC(=O)NC1C(OP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)OC(CO)C(O)C1O + O = O=c1ccn(C2OC(CO)C(O)C2O)c(=O)[nH]1 + CC(=O)NC1C(OP(=O)(O)OP(=O)(O)O)OC(CO)C(O)C1O"):

    reaction_object = featurizations.reaction(rxn_str)

    # we use our reaction fingerprinting function to automatically generate a reaction fingerprint
    reaction_fp = reaction_object.rxn_2_fp_w_positioning(
        fp_type = "ecfp4",
        radius = 2,
        nBits = 2048,
        is_folded = True,
        dim = 2048,
        max_species = 4,
        cofactor_positioning = "add_subtract", # the add-subtract configuration is used here
        all_cofactors_wo_stereo = all_cofactors_wo_stereo)

    # then, we manually create a reaction fingerprint using the add-concat configuration
    substrate_smiles = 'CC(=O)NC1C(OP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)OC(CO)C(O)C1O'
    water_smiles = 'O'
    product_1_smiles = 'O=c1ccn(C2OC(CO)C(O)C2O)c(=O)[nH]1'
    product_2_smiles = 'CC(=O)NC1C(OP(=O)(O)OP(=O)(O)O)OC(CO)C(O)C1O'

    substrate_mol = Chem.MolFromSmiles(substrate_smiles)
    water_mol = Chem.MolFromSmiles(water_smiles)
    product_1_mol = Chem.MolFromSmiles(product_1_smiles)
    product_2_mol= Chem.MolFromSmiles(product_2_smiles)

    substrate_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(substrate_mol, radius=2, nBits=2048)
    water_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(water_mol, radius=2, nBits=2048)

    product_1_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_1_mol, radius=2, nBits=2048)
    product_2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(product_2_mol, radius=2, nBits=2048)

    # we initialize the lhs_fp with 2048 zeros and do the same for the rhs_fp
    lhs_fp = np.zeros(2048)
    rhs_fp = np.zeros(2048)

    # all manually fingerprinted compound fingerprints are summed up amongst reactants and products
    lhs_fp = lhs_fp + substrate_fp + water_fp
    rhs_fp = rhs_fp + product_1_fp + product_2_fp

    # finally, to create the reaction fingerprint, the sum of reactant fingerprints is subtracted from ...
    # ... the sum of product fingerprints
    manual_rxn_fp = rhs_fp - lhs_fp

    # Check the length of the manually assembled fingerprint
    assert len(manual_rxn_fp) == 2048

    # then check the length of the fingerprint generated by our custom featurizations module
    assert len(reaction_fp) == 2048

    # finally, check if both methods yield the same fingerprint in terms of length (admittedly redundant check)
    assert len(reaction_fp) == len(manual_rxn_fp)

    # finally, check - elementwise - if both methods actually yield the same fingerprint
    for i in range(0, len(reaction_fp)):
        element_frm_rxn_fp = int(reaction_fp[i])
        element_frm_manually_constructed_fp = int(manual_rxn_fp[i])
        assert element_frm_rxn_fp == element_frm_manually_constructed_fp

    # another way of checking that both arrays are truly equal could be to use the np.array_equal() function
    assert np.array_equal(reaction_fp, manual_rxn_fp)