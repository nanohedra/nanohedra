# Copyright 2020 Joshua Laniado and Todd O. Yeates.
__author__ = "Joshua Laniado and Todd O. Yeates"
__copyright__ = "Copyright 2020, Nanohedra"
__version__ = "1.0"


class Atom:
    def __init__(self, number, type, alt_location, residue_type, chain, residue_number, code_for_insertion, x, y, z, occ, temp_fact, element_symbol, atom_charge):
        self.number = number
        self.type = type
        self.alt_location = alt_location
        self.residue_type = residue_type
        self.chain = chain
        self.residue_number = residue_number
        self.code_for_insertion = code_for_insertion
        self.x = x
        self.y = y
        self.z = z
        self.occ = occ
        self.temp_fact = temp_fact
        self.element_symbol = element_symbol
        self.atom_charge = atom_charge

    def __str__(self):
        # prints Atom in PDB format
        return "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format("ATOM", self.number, self.type, self.alt_location, self.residue_type, self.chain, self.residue_number, self.code_for_insertion, self.x, self.y, self.z, self.occ, self.temp_fact, self.element_symbol, self.atom_charge)

    def is_backbone(self):
        # returns True if atom is part of the proteins backbone and False otherwise
        backbone_specific_atom_type = ["N", "CA", "C", "O"]
        if self.type in backbone_specific_atom_type:
            return True
        else:
            return False

    def is_CB(self, InclGlyCA=False):
        if InclGlyCA:
            return self.type == "CB" or (self.type == "CA" and self.residue_type == "GLY")
        else:
            return self.type == "CB" or (self.type == "H" and self.residue_type == "GLY")

    def is_CA(self):
        return self.type == "CA"

    def coords(self):
        return [self.x, self.y, self.z]

    def __eq__(self, other):
        return (self.number == other.number and self.chain == other.chain and self.type == other.type and self.residue_type == other.residue_type)

    def get_number(self):
        return self.number

    def get_type(self):
        return self.type

    def get_alt_location(self):
        return self.alt_location

    def get_residue_type(self):
        return self.residue_type

    def get_chain(self):
        return self.chain

    def get_residue_number(self):
        return self.residue_number

    def get_code_for_insertion(self):
        return self.code_for_insertion

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def get_z(self):
        return self.z

    def get_occ(self):
        return self.occ

    def get_temp_fact(self):
        return self.temp_fact

    def get_element_symbol(self):
        return self.element_symbol

    def get_atom_charge(self):
        return self.atom_charge
