import sys
import statistics
from scipy.spatial import distance
from typing import List, Union
from Bio.PDB import PDBParser
from Bio.PDB import MMCIFParser
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
from Levenshtein import distance as lev_distance

class Chain3dStructure:
    def __init__(self, 
                 protein_id, chain_id,
                 structure_filename, 
                 load=False):
        
        self.__cache_key__ = f'{protein_id.lower()}{chain_id.upper()}'
        self.__cashable__ = [
            self.get_residue_count.__name__,
            self.get_nearest_residue_indexes.__name__,
            self.all_residues_has_alpha_carbon.__name__,
            self.compute_sequence.__name__,
            self.get_protrusion_vector.__name__,
            self.get_SASA_vector.__name__
        ]

        self.protein_id = protein_id
        self.chain_id = chain_id
        self.structure_file = structure_filename
        
        self._3d_structure: Union[Structure, None] = None
        self._3d_chain = None

        if load:
            self.load()

    def load(self, preferred_sequence=None):
        if (self.structure_file.endswith('cif')):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        self._3d_structure = parser.get_structure(id, self.structure_file)

        found_chains = []

        for model in self._3d_structure:
            for chain in model.get_chains():
                if chain.id.lower() == self.chain_id.lower():
                    found_chains.append(chain)                    

        if len(found_chains) == 0:
            raise Exception(f"Chain with ID {self.chain_id} not found in the structure. Protein {self.protein_id} - file {self.structure_file}")

        if len(found_chains) == 1 or preferred_sequence is None:
            self._3d_chain = found_chains[0]
            return
        
        best_match_chain = found_chains[0]
        best_lev_dist = sys.maxsize
        
        for chain in found_chains:
            chain_seq = Chain3dStructure.compute_chain_sequence(chain)
            dist = lev_distance(chain_seq, preferred_sequence)

            if dist < best_lev_dist:
                best_lev_dist = dist
                best_match_chain = chain

        self._3d_chain = best_match_chain

    def get_residue_count(self) -> int:
        return len(self.get_residue_list())

    def get_bfactors(self, bfactor_algorithm='average_bfactor') -> List[int]:
        residues = self.get_residue_list()

        if bfactor_algorithm in BFactorFunctionsCollection.names_to_functions:
            bfactor_algorithm = BFactorFunctionsCollection.names_to_functions[bfactor_algorithm] 

        bfactors = []
        for residue in residues:
            bfactors.append(bfactor_algorithm(residue))
        return bfactors


    def get_nearest_residue_indexes(
            self, 
            residue_index: int, 
            
            # do not change these default values as some values may be stored in cache (alternatively delete cache)
            n_nearest: int = 10,
            distance_func='alpha_atoms') -> List[int]:

        if distance_func in DistancesCollection.names_to_functions:
            distance_func = DistancesCollection.names_to_functions[distance_func] 

        residue_list = self.get_residue_list()
        center_residue = residue_list[residue_index]

        distances = []

        for i in range(len(residue_list)):
            distances.append((
                i, distance_func(center_residue, residue_list[i])
            ))

        distances.sort(key=lambda x: x[1])

        return list([residue_index for residue_index, dist in distances[:n_nearest]])

    def free_memory(self):
        self._3d_structure = None
        self._3d_chain = None

    def get_residue_list(self) -> List[Residue]: 
        return Chain3dStructure.get_chain_residue_list(self._3d_chain) 

    def compute_sequence(self) -> str:
        return Chain3dStructure.compute_chain_sequence(self._3d_chain)

    def all_residues_has_alpha_carbon(self):
        for res in self.get_residue_list():
            if 'CA' not in res:
                return False
            
        return True

    def get_origin_filename(self):
        return self.structure_file

    def get_protrusion_vector(
            self, radius, 
            # do not change this default value as some values may be stored in cache (alternatively delete cache)
            protrusion_algorithm='atom_count_from_center_of_mass'):
        
        if protrusion_algorithm in ProtrusionFunctionsCollection.names_to_functions:
            protrusion_algorithm = ProtrusionFunctionsCollection.names_to_functions[protrusion_algorithm] 

        all_atoms = list(self._3d_chain.get_atoms())
        ns = NeighborSearch(all_atoms)

        protrusion = []

        for residue in self.get_residue_list():
            value = protrusion_algorithm(ns, residue, radius)
            protrusion.append(value)

        return protrusion
    
    def get_SASA_vector(self):
        sr = ShrakeRupley()
        sr.compute(self._3d_chain, level="R")
        
        return [
            res.sasa for res in self.get_residue_list()
        ]

    @staticmethod
    def get_chain_residue_list(biopython_chain) -> List[Residue]:
        # Exclude hetero residues or non-standard residues with insertion code ' '
        return list([residue for residue in biopython_chain if residue.get_id()[0] == " "])

    @staticmethod
    def compute_chain_sequence(biopython_chain) -> str:
        return ''.join([AminoAcidMapper.to_one_letter_code(
                            residue.get_resname())
                        for residue in Chain3dStructure.get_chain_residue_list(biopython_chain)])

allowed_file_types = ['pdb', 'ent', 'cif']


class AminoAcidMapper:
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
        'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V'
    }

    @staticmethod
    def to_one_letter_code(long_code) -> str:
        def map_name(name):
            if name in AminoAcidMapper.three_to_one:
                return AminoAcidMapper.three_to_one[name]
            return f'#({name})'
        
        return map_name(long_code)


class ResidueDistances:
    @staticmethod
    def distance_between_alpha_atoms_sq(residue_1: Residue, residue_2: Residue) -> float:
        if ('CA' not in residue_1) or ('CA' not in residue_2):
            return sys.maxsize

        coord_1 = residue_1['CA'].get_coord()
        coord_2 = residue_2['CA'].get_coord()

        return sum((a - b) * (a - b) for a, b in zip(coord_1, coord_2))

class CenterAtomFunctions:
    @staticmethod
    def get_furthermost_atom(residue):
        atoms = list(residue.get_atoms())

        if len(atoms) == 1: center = atoms[0]
        elif len(atoms) < 5: center = atoms[1] # take c_alpha atom if side chain is not observed
        elif len(atoms) >= 5: 
            c_alpha_atom = atoms[1]
            max_distance = -1
            for idx, atom in enumerate(residue.get_atoms()):
                if idx == 1: continue
                atom_distance = distance.euclidean(c_alpha_atom.get_coord(), atom.get_coord())
                if max_distance < atom_distance: 
                    max_distance = atom_distance
                    center = atom
        
        return center

class ProtrusionFunctions:
    @staticmethod
    def compute_neighboring_atoms_from_center(ns, residue, radius):
        center = residue.center_of_mass()
        neighbors = ns.search(center=center, radius=radius, level='R') 
        return len(neighbors)
    
    @staticmethod
    def compute_neighboring_atoms_from_last_atom_coords(ns, residue, radius):
        center = list(residue.get_atoms())[-1].get_coord()
        neighbors = ns.search(center=center, radius=radius, level='R') 
        return len(neighbors)

    @staticmethod
    def compute_neighboring_atoms_from_furthermost_atom_coords(ns, residue, radius):
        center = CenterAtomFunctions.get_furthermost_atom(residue).get_coord()
        # compute protrusion
        neighbors = ns.search(center=center, radius=radius, level='R') 
        return len(neighbors)

class BFactorFunctions:
    @staticmethod
    def average_bfactor(residue) -> float:
        return statistics.mean([atom.get_bfactor() for atom in residue.get_atoms()])
    
    @staticmethod
    def furthest_atom_bfactor(residue) -> float:
        center = CenterAtomFunctions.get_furthermost_atom(residue)

        return center.get_bfactor()

class DistancesCollection:
    names_to_functions = {
        'alpha_atoms': ResidueDistances.distance_between_alpha_atoms_sq
    }

class ProtrusionFunctionsCollection:
    names_to_functions = {
        'atom_count_from_center_of_mass': ProtrusionFunctions.compute_neighboring_atoms_from_center,
        'atom_count_from_last_atom_coords': ProtrusionFunctions.compute_neighboring_atoms_from_last_atom_coords,
        'atom_count_from_furthermost_atom_coords': ProtrusionFunctions.compute_neighboring_atoms_from_furthermost_atom_coords
    }

class BFactorFunctionsCollection:
    names_to_functions = {
        'average_bfactor': BFactorFunctions.average_bfactor,
        'furthest_atom_bfactor': BFactorFunctions.furthest_atom_bfactor
    }

        