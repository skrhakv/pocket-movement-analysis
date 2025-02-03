def get_protein_backbone(protein, chain_id):
    """Get the protein backbone atoms from a protein structure.
    
    Parameters
    ----------
    protein : AtomArray
        The protein structure.
    
    Returns
    -------
    protein_backbone : AtomArray
        The protein backbone atoms.
    """
    protein_backbone = protein[(protein.atom_name == "CA") 
                       & (protein.element == "C") 
                       & (protein.chain_id == chain_id) 
                       & (
                             (protein.res_name == 'ALA')
                           | (protein.res_name == 'ARG')
                           | (protein.res_name == 'ASN')
                           | (protein.res_name == 'ASP')
                           | (protein.res_name == 'CYS')
                           | (protein.res_name == 'GLN')
                           | (protein.res_name == 'GLU')
                           | (protein.res_name == 'GLY')
                           | (protein.res_name == 'HIS')
                           | (protein.res_name == 'ILE')
                           | (protein.res_name == 'LEU')
                           | (protein.res_name == 'LYS')
                           | (protein.res_name == 'MET')
                           | (protein.res_name == 'PHE')
                           | (protein.res_name == 'PRO')
                           | (protein.res_name == 'SER')
                           | (protein.res_name == 'THR')
                           | (protein.res_name == 'TRP')
                           | (protein.res_name == 'TYR')
                           | (protein.res_name == 'VAL'))]
    return protein_backbone
    