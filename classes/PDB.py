from Atom import Atom
import subprocess
import os
from shutil import copyfile
from shutil import move


# Copyright 2020 Joshua Laniado and Todd O. Yeates.
__author__ = "Joshua Laniado and Todd O. Yeates"
__copyright__ = "Copyright 2020, Nanohedra"
__version__ = "1.0"


class PDB:
    def __init__(self):
        self.all_atoms = []  # python list of Atoms
        self.filepath = None  # PDB filepath if instance is read from PDB file
        self.chain_id_list = []  # list of unique chain IDs in PDB
        self.cb_coords = []
        self.bb_coords = []

    def set_all_atoms(self, atom_list):
        self.all_atoms = atom_list

    def set_chain_id_list(self, chain_id_list):
        self.chain_id_list = chain_id_list

    def set_filepath(self, filepath):
        self.filepath = filepath

    def get_all_atoms(self):
        return self.all_atoms

    def get_chain_id_list(self):
        return self.chain_id_list

    def get_filepath(self):
        return self.filepath

    def readfile(self, filepath, remove_alt_location=False):
        # reads PDB file and feeds PDB instance
        self.filepath = filepath

        f = open(filepath, "r")
        pdb = f.readlines()
        f.close()

        available_chain_ids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R',
                               'S', 'T',
                               'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
                               'm', 'n',
                               'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4']

        chain_ids = []
        multimodel = False
        start_of_new_model = False
        model_chain_index = -1
        model_chain_id = None
        curr_chain_id = None
        for line in pdb:
            line = line.rstrip()
            if line[0:5] == "MODEL":
                start_of_new_model = True
                multimodel = True
                model_chain_index += 1
                model_chain_id = available_chain_ids[model_chain_index]
            elif line[0:4] == "ATOM":
                number = int(line[6:11].strip())
                type = line[12:16].strip()
                alt_location = line[16:17].strip()
                residue_type = line[17:20].strip()
                if multimodel:
                    if line[21:22].strip() != curr_chain_id:
                        curr_chain_id = line[21:22].strip()
                        if not start_of_new_model:
                            model_chain_index += 1
                            model_chain_id = available_chain_ids[model_chain_index]
                    start_of_new_model = False
                    chain = model_chain_id
                else:
                    chain = line[21:22].strip()
                residue_number = int(line[22:26].strip())
                code_for_insertion = line[26:27].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occ = float(line[54:60].strip())
                temp_fact = float(line[60:66].strip())
                element_symbol = line[76:78].strip()
                atom_charge = line[78:80].strip()
                atom = Atom(number, type, alt_location, residue_type, chain, residue_number, code_for_insertion, x, y, z, occ, temp_fact, element_symbol, atom_charge)
                if remove_alt_location:
                    if alt_location == "" or alt_location == "A":
                        if atom.chain not in chain_ids:
                            chain_ids.append(atom.chain)
                        self.all_atoms.append(atom)
                else:
                    if atom.chain not in chain_ids:
                        chain_ids.append(atom.chain)
                    self.all_atoms.append(atom)
        self.chain_id_list = chain_ids

    def read_atom_list(self, atom_list, store_cb_and_bb_coords=False):
        # reads a python list of Atoms and feeds PDB instance
        if store_cb_and_bb_coords:
            chain_ids = []
            for atom in atom_list:
                self.all_atoms.append(atom)
                if atom.is_backbone():
                    [x, y, z] = [atom.x, atom.y, atom.z]
                    self.bb_coords.append([x, y, z])
                if atom.is_CB(InclGlyCA=False):
                    [x, y, z] = [atom.x, atom.y, atom.z]
                    self.cb_coords.append([x, y, z])
                if atom.chain not in chain_ids:
                    chain_ids.append(atom.chain)
            self.chain_id_list = chain_ids
        else:
            chain_ids = []
            for atom in atom_list:
                self.all_atoms.append(atom)
                if atom.chain not in chain_ids:
                    chain_ids.append(atom.chain)
            self.chain_id_list = chain_ids

    def chain(self, chain_id):
        # returns a python list of Atoms containing the subset of Atoms in the PDB instance that belong to the selected chain ID
        selected_atoms = []
        for atom in self.all_atoms:
            if atom.chain == chain_id:
                selected_atoms.append(atom)
        return selected_atoms

    def extract_all_coords(self):
        coords = []
        for atom in self.all_atoms:
            [x, y, z] = [atom.x, atom.y, atom.z]
            coords.append([x, y, z])
        return coords

    def extract_backbone_coords(self):
        coords = []
        for atom in self.all_atoms:
            if atom.is_backbone():
                [x, y, z] = [atom.x, atom.y, atom.z]
                coords.append([x, y, z])
        return coords

    def get_CA_atoms(self):
        ca_atoms = []
        for atom in self.all_atoms:
            if atom.is_CA():
                ca_atoms.append(atom)
        return ca_atoms

    def get_backbone_atoms(self):
        bb_atoms = []
        for atom in self.all_atoms:
            if atom.is_backbone():
                bb_atoms.append(atom)
        return bb_atoms

    def get_CB_coords(self, ReturnWithCBIndices=False, InclGlyCA=False):
        coords = []
        cb_indices = []
        for i in range(len(self.all_atoms)):
            if self.all_atoms[i].is_CB(InclGlyCA=InclGlyCA):
                [x, y, z] = [self.all_atoms[i].x, self.all_atoms[i].y, self.all_atoms[i].z]
                coords.append([x, y, z])
                if ReturnWithCBIndices:
                    cb_indices.append(i)
        if ReturnWithCBIndices:
            return coords, cb_indices
        else:
            return coords

    def replace_coords(self, new_cords):
        for i in range(len(self.all_atoms)):
            self.all_atoms[i].x, self.all_atoms[i].y, self.all_atoms[i].z = new_cords[i][0], new_cords[i][1], new_cords[i][2]

    def mat_vec_mul3(self, a, b):
        c = [0. for i in range(3)]
        for i in range(3):
            c[i] = 0.
            for j in range(3):
                c[i] += a[i][j] * b[j]
        return c

    def rotate_translate(self, rot, tx):
        for atom in self.all_atoms:
            coord = [atom.x, atom.y, atom.z]
            coord_rot = self.mat_vec_mul3(rot, coord)
            newX = coord_rot[0] + tx[0]
            newY = coord_rot[1] + tx[1]
            newZ = coord_rot[2] + tx[2]
            atom.x, atom.y, atom.z = newX, newY, newZ

    def write(self, out_path, cryst1=None):
        outfile = open(out_path, "w")
        if cryst1 is not None and isinstance(cryst1, str) and cryst1.startswith("CRYST1"):
            outfile.write(str(cryst1))
        for atom in self.all_atoms:
            outfile.write(str(atom))
        outfile.close()

    def orient(self, symm, output_dir, orient_executable_dir):
        valid_subunit_number = {"C2": 2, "C3": 3, "C4": 4, "C5": 5, "C6": 6,
                                "D2": 4, "D3": 6, "D4": 8, "D5": 10, "D6": 12,
                                "I": 60, "O": 24, "T": 12}

        number_of_subunits = len(self.chain_id_list)

        pdb_file_name = os.path.basename(self.filepath)

        if number_of_subunits != valid_subunit_number[symm]:
            orient_log = open('%s/%s' % (output_dir, 'orient_oligomer_log.txt'), 'a+')
            orient_log.write("%s\n Oligomer could not be oriented: It has %s subunits while %s are expected "
                             "for %s symmetry\n\n" % (pdb_file_name, str(number_of_subunits),
                                                       str(valid_subunit_number[symm]), symm))
            orient_log.close()

            raise ValueError('orient_oligomer could not orient %s '
                             'check %s/orient_oligomer_log.txt for more information\n' % (pdb_file_name, output_dir))

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if os.path.exists(orient_executable_dir + "/input.pdb"):
            os.remove(orient_executable_dir + "/input.pdb")
        if os.path.exists(orient_executable_dir + "/output.pdb"):
            os.remove(orient_executable_dir + "/output.pdb")

        copyfile('%s' % self.filepath, '%s/input.pdb' % orient_executable_dir)

        process = subprocess.Popen(['%s/orient_oligomer' % orient_executable_dir],
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   cwd=orient_executable_dir)
        in_symm_file = '%s/symm_files/%s_symm.txt' % (orient_executable_dir, symm)
        stdout, stderr = process.communicate(input=in_symm_file)
        stdout = pdb_file_name + stdout[28:]

        orient_log = open('%s/%s' % (output_dir, 'orient_oligomer_log.txt'), 'a+')
        orient_log.write(stdout)
        if stderr != '':
            orient_log.write(stderr + "\n")
        else:
            orient_log.write('\n')
        orient_log.close()

        if os.path.exists(orient_executable_dir + "/output.pdb") and os.stat(orient_executable_dir + "/output.pdb").st_size != 0:
            move('%s/output.pdb' % orient_executable_dir, '%s/%s' % (output_dir, pdb_file_name))

        if os.path.exists(orient_executable_dir + "/input.pdb"):
            os.remove(orient_executable_dir + "/input.pdb")
        if os.path.exists(orient_executable_dir + "/output.pdb"):
            os.remove(orient_executable_dir + "/output.pdb")

        if not os.path.exists('%s/%s' % (output_dir, pdb_file_name)):
            raise RuntimeError('orient_oligomer could not orient %s '
                               'check %s/orient_oligomer_log.txt for more information\n' % (pdb_file_name, output_dir))

    def get_surface_resdiue_info(self, free_sasa_exe_path, probe_radius=2.2, sasa_thresh=0):
        # only works for monomers or homo-oligomers
        assert_error_message = "Could not locate freesasa executable here: %s" % free_sasa_exe_path
        assert os.path.exists(free_sasa_exe_path), assert_error_message

        proc = subprocess.Popen([free_sasa_exe_path, '--format=seq', '--probe-radius', str(probe_radius), self.filepath]
                                , stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = proc.communicate()
        out_lines = out.split("\n")
        sasa_out = []
        for line in out_lines:
            if line != "\n" and line != "" and not line.startswith("#"):
                chain_id = line[4:5]
                res_num = int(line[5:10])
                sasa = float(line[17:])
                if sasa > sasa_thresh:
                    sasa_out.append((chain_id, res_num))
        return sasa_out

