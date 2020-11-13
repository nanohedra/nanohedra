from classes.FragDock import dock
from classes.EulerLookup import EulerLookup
from classes.Fragment import *
from classes.SymEntry import *
from utils.SamplingUtils import get_degeneracy_matrices
from utils.CmdLineArgParseUtils import *
from utils.NanohedraManualUtils import *
from utils.ExpandAssemblyUtils import *
import sys
import os
import subprocess


# Copyright 2020 Joshua Laniado and Todd O. Yeates.
__author__ = "Joshua Laniado and Todd O. Yeates"
__copyright__ = "Copyright 2020, Nanohedra"
__version__ = "1.0"


def main():
    cmd_line_in_params = sys.argv

    if len(cmd_line_in_params) > 1 and cmd_line_in_params[1] == '-dock':

        # Parse Command Line Input
        sym_entry_number, pdb_dir1_path, pdb_dir2_path, rot_step_deg1, rot_step_deg2, master_outdir, output_exp_assembly, output_uc, output_surrounding_uc, min_matched, init_match_type = get_docking_parameters(cmd_line_in_params)

        # Master Log File
        master_log_filepath = master_outdir + "/nanohedra_master_logfile.txt"

        # Making Master Output Directory
        if not os.path.exists(master_outdir):
            os.makedirs(master_outdir)

        # Getting PDB1 File paths
        pdb1_filepaths = []
        for root1, dirs1, files1 in os.walk(pdb_dir1_path):
            for file1 in files1:
                if '.pdb' in file1:
                    pdb1_filepaths.append(pdb_dir1_path + "/" + file1)
        if len(pdb1_filepaths) == 0:
            master_log_file = open(master_log_filepath, "a+")
            master_log_file.write("\nCOULD NOT FIND PDB FILE(S) IN THE INPUT DIRECTORY SPECIFIED FOR OLIGOMER 1:\n")
            master_log_file.write("%s\n" % pdb_dir1_path)
            master_log_file.write("NANOHEDRA DOCKING RUN ENDED\n")
            master_log_file.close()
            sys.exit()

        # Getting PDB2 File paths
        pdb2_filepaths = []
        for root2, dirs2, files2 in os.walk(pdb_dir2_path):
            for file2 in files2:
                if '.pdb' in file2:
                    pdb2_filepaths.append(pdb_dir2_path + "/" + file2)
        if len(pdb2_filepaths) == 0:
            master_log_file = open(master_log_filepath, "a+")
            master_log_file.write("\nCOULD NOT FIND PDB FILE(S) IN THE INPUT DIRECTORY SPECIFIED FOR OLIGOMER 2:\n")
            master_log_file.write("%s\n" % pdb_dir2_path)
            master_log_file.write("NANOHEDRA DOCKING RUN ENDED\n")
            master_log_file.close()
            sys.exit()

        try:
            # Nanohedra.py Path
            main_script_dir = os.path.dirname(os.path.realpath(__file__))

            # Free SASA Executable Path
            free_sasa_exe_path = "/usr/local/bin/freesasa"
            sasa_assert_error_message = "Could not locate freesasa executable here: %s\n" \
                                        "FreeSASA might not be (correctly) installed.\n" \
                                        "Go to: https://freesasa.github.io for a quick " \
                                        "installation guide." % free_sasa_exe_path
            assert os.path.exists(free_sasa_exe_path), sasa_assert_error_message

            sasa_v_proc = subprocess.Popen(['%s' % free_sasa_exe_path, '--version'], stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
            (sasa_v_out, sasa_v_err) = sasa_v_proc.communicate()
            free_sasa_v = sasa_v_out.split("\n")[0]
            assert_free_sasa_v_message = free_sasa_v + " not supported.\nInstall supported version: " \
                                                       "FreeSASA 2.0.3 at https://freesasa.github.io"
            assert free_sasa_v == "FreeSASA 2.0.3", assert_free_sasa_v_message

            # Orient Oligomer Fortran Executable Path
            orient_executable_path = main_script_dir + "/orient/orient_oligomer"
            orient_assert_error_message = "Could not locate orient_oligomer executable here: %s\n" \
                                          "Check README file for instructions on how to compile " \
                                          "orient_oligomer.f" % orient_executable_path
            assert os.path.exists(orient_executable_path), orient_assert_error_message
            orient_executable_dir = os.path.dirname(orient_executable_path)

            # Fragment Database Directory Paths
            monofrag_cluster_rep_dirpath = main_script_dir + "/fragment_database/Top5MonoFragClustersRepresentativeCentered"
            ijk_intfrag_cluster_rep_dirpath = main_script_dir + "/fragment_database/Top75percent_IJK_ClusterRepresentatives_1A"
            intfrag_cluster_info_dirpath = main_script_dir + "/fragment_database/IJK_ClusteredInterfaceFragmentDBInfo_1A"

            # SymEntry Parameters
            sym_entry = SymEntry(sym_entry_number)

            oligomer_symmetry_1 = sym_entry.get_group1_sym()
            oligomer_symmetry_2 = sym_entry.get_group2_sym()
            design_symmetry = sym_entry.get_pt_grp_sym()

            rot_range_deg_pdb1 = sym_entry.get_rot_range_deg_1()
            rot_range_deg_pdb2 = sym_entry.get_rot_range_deg_2()

            set_mat1 = sym_entry.get_rot_set_mat_group1()
            set_mat2 = sym_entry.get_rot_set_mat_group2()

            is_internal_zshift1 = sym_entry.is_internal_tx1()
            is_internal_zshift2 = sym_entry.is_internal_tx2()

            is_internal_rot1 = sym_entry.is_internal_rot1()
            is_internal_rot2 = sym_entry.is_internal_rot2()

            design_dim = sym_entry.get_design_dim()

            ref_frame_tx_dof1 = sym_entry.get_ref_frame_tx_dof_group1()
            ref_frame_tx_dof2 = sym_entry.get_ref_frame_tx_dof_group2()

            result_design_sym = sym_entry.get_result_design_sym()
            uc_spec_string = sym_entry.get_uc_spec_string()

            # Default Fragment Guide Atom Overlap Z-Value Threshold For Initial Matches
            init_max_z_val = 1.0

            # Default Fragment Guide Atom Overlap Z-Value Threshold For All Subsequent Matches
            subseq_max_z_val = 2.0

            master_log_file = open(master_log_filepath, "a+")

            # Default Rotation Step
            if is_internal_rot1 and rot_step_deg1 is None:
                rot_step_deg1 = 3  # If rotation step not provided but required, set rotation step to default
            if is_internal_rot2 and rot_step_deg2 is None:
                rot_step_deg2 = 3  # If rotation step not provided but required, set rotation step to default

            if not is_internal_rot1 and rot_step_deg1 is not None:
                rot_step_deg1 = 1
                master_log_file.write(
                    "Warning: Rotation Step 1 Specified Was Ignored. Oligomer 1 Does Not Have Internal Rotational DOF\n\n")
            if not is_internal_rot2 and rot_step_deg2 is not None:
                rot_step_deg2 = 1
                master_log_file.write(
                    "Warning: Rotation Step 2 Specified Was Ignored. Oligomer 2 Does Not Have Internal Rotational DOF\n\n")

            if not is_internal_rot1 and rot_step_deg1 is None:
                rot_step_deg1 = 1
            if not is_internal_rot2 and rot_step_deg2 is None:
                rot_step_deg2 = 1

            # For SCMs where the two oligomeric components have the same point group symmetry
            # make sure that there is at least 2 PDB files in the input PDB directory
            if (oligomer_symmetry_1 == oligomer_symmetry_2) and len(pdb1_filepaths) < 2:
                master_log_file.write("\nAT LEAST 2 OLIGOMERS ARE REQUIRED TO BE IN THE SPECIFIED INPUT DIRECTORY,\n")
                master_log_file.write("WHEN THE 2 COMPONENTS OF A SCM OBEY THE SAME POINT GROUP SYMMETRY ")
                master_log_file.write("(IN THIS CASE %s)\n" % oligomer_symmetry_1)
                master_log_file.write("%s PDB FILE(S) FOUND IN: %s\n" % (str(len(pdb1_filepaths)), pdb_dir1_path))
                master_log_file.write("NANOHEDRA DOCKING RUN ENDED\n")
                master_log_file.close()
                sys.exit()

            master_log_file.write("NANOHEDRA PROJECT INFORMATION\n")
            master_log_file.write("Oligomer 1 Input Directory: %s\n" % pdb_dir1_path)
            master_log_file.write("Oligomer 2 Input Directory: %s\n" % pdb_dir2_path)
            master_log_file.write("Master Output Directory: %s\n\n" % master_outdir)

            master_log_file.write("SYMMETRY COMBINATION MATERIAL INFORMATION\n")
            master_log_file.write("Nanohedra Entry Number: %s\n" % str(sym_entry_number))
            master_log_file.write("Oligomer 1 Point Group Symmetry: %s\n" % oligomer_symmetry_1)
            master_log_file.write("Oligomer 2 Point Group Symmetry: %s\n" % oligomer_symmetry_2)
            master_log_file.write("SCM Point Group Symmetry: %s\n" % design_symmetry)
            master_log_file.write("Oligomer 1 Internal ROT DOF: %s\n" % str(sym_entry.get_internal_rot1()))
            master_log_file.write("Oligomer 2 Internal ROT DOF: %s\n" % str(sym_entry.get_internal_rot2()))
            master_log_file.write("Oligomer 1 Internal Tx DOF: %s\n" % str(sym_entry.get_internal_tx1()))
            master_log_file.write("Oligomer 2 Internal Tx DOF: %s\n" % str(sym_entry.get_internal_tx2()))
            master_log_file.write("Oligomer 1 Setting Matrix: %s\n" % str(set_mat1))
            master_log_file.write("Oligomer 2 Setting Matrix: %s\n" % str(set_mat2))
            master_log_file.write("Oligomer 1 Reference Frame Tx DOF: %s\n" % (str(ref_frame_tx_dof1) if sym_entry.is_ref_frame_tx_dof1() else str(None)))
            master_log_file.write("Oligomer 2 Reference Frame Tx DOF: %s\n" % (str(ref_frame_tx_dof2) if sym_entry.is_ref_frame_tx_dof2() else str(None)))
            master_log_file.write("Resulting SCM Symmetry: %s\n" % result_design_sym)
            master_log_file.write("SCM Dimension: %s\n" % str(design_dim))
            master_log_file.write("SCM Unit Cell Specification: %s\n\n" % uc_spec_string)

            master_log_file.write("ROTATIONAL SAMPLING INFORMATION\n")
            master_log_file.write("Oligomer 1 ROT Sampling Range: %s\n" % (str(rot_range_deg_pdb1) if is_internal_rot1 else "N/A"))
            master_log_file.write("Oligomer 2 ROT Sampling Range: %s\n" % (str(rot_range_deg_pdb2) if is_internal_rot2 else "N/A"))
            master_log_file.write("Oligomer 1 ROT Sampling Step: %s\n" % (str(rot_step_deg1) if is_internal_rot1 else "N/A"))
            master_log_file.write("Oligomer 2 ROT Sampling Step: %s\n\n" % (str(rot_step_deg2) if is_internal_rot2 else "N/A"))

            # Orient Input Oligomers to Canonical Orientation
            if oligomer_symmetry_1 == oligomer_symmetry_2:
                master_log_file.write("ORIENTING INPUT OLIGOMER PDB FILES\n")
                master_log_file.close()
                oriented_pdb1_outdir = master_outdir + "/" + oligomer_symmetry_1 + "_oriented"
                if not os.path.exists(oriented_pdb1_outdir):
                    os.makedirs(oriented_pdb1_outdir)
                pdb1_oriented_filepaths = []
                pdb2_oriented_filepaths = []
                for pdb1_path in pdb1_filepaths:
                    pdb1 = PDB()
                    pdb1.readfile(pdb1_path, remove_alt_location=True)
                    pdb1_filename = os.path.basename(pdb1_path)
                    try:
                        pdb1.orient(oligomer_symmetry_1, oriented_pdb1_outdir, orient_executable_dir)
                        pdb1_oriented_filepaths.append(oriented_pdb1_outdir + "/" + pdb1_filename)
                        master_log_file = open(master_log_filepath, 'a+')
                        master_log_file.write("oriented: %s\n" % pdb1_filename)
                        master_log_file.close()
                    except ValueError as val_err:
                        master_log_file = open(master_log_filepath, 'a+')
                        master_log_file.write(str(val_err))
                        master_log_file.close()
                    except RuntimeError as rt_err:
                        master_log_file = open(master_log_filepath, 'a+')
                        master_log_file.write(str(rt_err))
                        master_log_file.close()
                if len(pdb1_oriented_filepaths) == 0:
                    master_log_file = open(master_log_filepath, "a+")
                    master_log_file.write("\nCOULD NOT ORIENT OLIGOMER INPUT PDB FILES\n")
                    master_log_file.write(
                        "CHECK %s/orient_oligomer_log.txt FOR MORE INFORMATION\n" % oriented_pdb1_outdir)
                    master_log_file.write("NANOHEDRA DOCKING RUN ENDED\n")
                    master_log_file.close()
                    sys.exit()
                elif len(pdb1_oriented_filepaths) == 1:
                    master_log_file = open(master_log_filepath, "a+")
                    master_log_file.write("\nAT LEAST 2 OLIGOMERS ARE REQUIRED WHEN THE ")
                    master_log_file.write("2 OLIGOMERIC COMPONENTS OF A SCM OBEY THE SAME POINT GROUP SYMMETRY ")
                    master_log_file.write("(IN THIS CASE: %s)\n" % oligomer_symmetry_1)
                    master_log_file.write("HOWEVER ONLY 1 INPUT OLIGOMER PDB FILE COULD BE ORIENTED\n")
                    master_log_file.write(
                        "CHECK %s/orient_oligomer_log.txt FOR MORE INFORMATION\n" % oriented_pdb1_outdir)
                    master_log_file.write("NANOHEDRA DOCKING RUN ENDED\n")
                    master_log_file.close()
                    sys.exit()
                else:
                    master_log_file = open(master_log_filepath, "a+")
                    master_log_file.write(
                        "Successfully Oriented %s out of the %s Oligomer Input PDB Files\n==> %s\n\n"
                        % (str(len(pdb1_oriented_filepaths)), str(len(pdb1_filepaths)), oriented_pdb1_outdir))
                    master_log_file.close()
            else:
                master_log_file.write("ORIENTING OLIGOMER 1 INPUT PDB FILE(S)\n")
                master_log_file.close()
                oriented_pdb1_outdir = master_outdir + "/" + oligomer_symmetry_1 + "_oriented"
                if not os.path.exists(oriented_pdb1_outdir):
                    os.makedirs(oriented_pdb1_outdir)
                pdb1_oriented_filepaths = []
                for pdb1_path in pdb1_filepaths:
                    pdb1 = PDB()
                    pdb1.readfile(pdb1_path, remove_alt_location=True)
                    pdb1_filename = os.path.basename(pdb1_path)
                    try:
                        pdb1.orient(oligomer_symmetry_1, oriented_pdb1_outdir, orient_executable_dir)
                        pdb1_oriented_filepaths.append(oriented_pdb1_outdir + "/" + pdb1_filename)
                        master_log_file = open(master_log_filepath, 'a+')
                        master_log_file.write("oriented: %s\n" % pdb1_filename)
                        master_log_file.close()
                    except ValueError as val_err:
                        master_log_file = open(master_log_filepath, 'a+')
                        master_log_file.write(str(val_err))
                        master_log_file.close()
                    except RuntimeError as rt_err:
                        master_log_file = open(master_log_filepath, 'a+')
                        master_log_file.write(str(rt_err))
                        master_log_file.close()
                if len(pdb1_oriented_filepaths) == 0:
                    master_log_file = open(master_log_filepath, "a+")
                    master_log_file.write("\nCOULD NOT ORIENT OLIGOMER 1 INPUT PDB FILE(S)\n")
                    master_log_file.write(
                        "CHECK %s/orient_oligomer_log.txt FOR MORE INFORMATION\n" % oriented_pdb1_outdir)
                    master_log_file.write("NANOHEDRA DOCKING RUN ENDED\n")
                    master_log_file.close()
                    sys.exit()
                else:
                    master_log_file = open(master_log_filepath, "a+")
                    master_log_file.write(
                        "Successfully Oriented %s out of the %s Oligomer 1 Input PDB File(s)\n==> %s\n"
                        % (str(len(pdb1_oriented_filepaths)), str(len(pdb1_filepaths)), oriented_pdb1_outdir))
                    master_log_file.close()

                master_log_file = open(master_log_filepath, 'a+')
                master_log_file.write("\nORIENTING OLIGOMER 2 INPUT PDB FILE(S)\n")
                master_log_file.close()
                oriented_pdb2_outdir = master_outdir + "/" + oligomer_symmetry_2 + "_oriented"
                if not os.path.exists(oriented_pdb2_outdir):
                    os.makedirs(oriented_pdb2_outdir)
                pdb2_oriented_filepaths = []
                for pdb2_path in pdb2_filepaths:
                    pdb2 = PDB()
                    pdb2.readfile(pdb2_path, remove_alt_location=True)
                    pdb2_filename = os.path.basename(pdb2_path)
                    try:
                        pdb2.orient(oligomer_symmetry_2, oriented_pdb2_outdir, orient_executable_dir)
                        pdb2_oriented_filepaths.append(oriented_pdb2_outdir + "/" + pdb2_filename)
                        master_log_file = open(master_log_filepath, 'a+')
                        master_log_file.write("oriented: %s\n" % pdb2_filename)
                        master_log_file.close()
                    except ValueError as val_err:
                        master_log_file = open(master_log_filepath, 'a+')
                        master_log_file.write(str(val_err))
                        master_log_file.close()
                    except RuntimeError as rt_err:
                        master_log_file = open(master_log_filepath, 'a+')
                        master_log_file.write(str(rt_err))
                        master_log_file.close()
                if len(pdb2_oriented_filepaths) == 0:
                    master_log_file = open(master_log_filepath, "a+")
                    master_log_file.write("\nCOULD NOT ORIENT OLIGOMER 2 INPUT PDB FILE(S)\n")
                    master_log_file.write(
                        "CHECK %s/orient_oligomer_log.txt FOR MORE INFORMATION\n" % oriented_pdb2_outdir)
                    master_log_file.write("NANOHEDRA DOCKING RUN ENDED\n")
                    master_log_file.close()
                    sys.exit()
                else:
                    master_log_file = open(master_log_filepath, "a+")
                    master_log_file.write(
                        "Successfully Oriented %s out of the %s Oligomer 2 Input PDB File(s)\n==> %s\n\n"
                        % (str(len(pdb2_oriented_filepaths)), str(len(pdb2_filepaths)), oriented_pdb2_outdir))
                    master_log_file.close()

            # Get Degeneracy Matrices
            master_log_file = open(master_log_filepath, "a+")
            master_log_file.write("SEARCHING FOR POSSIBLE DEGENERACIES\n")
            degeneracy_matrices_1, degeneracy_matrices_2 = get_degeneracy_matrices(oligomer_symmetry_1,
                                                                                   oligomer_symmetry_2,
                                                                                   design_dim,
                                                                                   design_symmetry)
            if degeneracy_matrices_1 is None:
                master_log_file.write("No Degeneracies Found for Oligomer 1\n")
            elif len(degeneracy_matrices_1) == 1:
                master_log_file.write("1 Degeneracy Found for Oligomer 1\n")
                master_log_file.write("%s\n\n" % str(degeneracy_matrices_1[0]))
            else:
                master_log_file.write("%s Degeneracies Found for Oligomer 1\n" % str(len(degeneracy_matrices_1)))
                for degen_mat_1 in degeneracy_matrices_1:
                    master_log_file.write("%s\n" % str(degen_mat_1))
                master_log_file.write("\n")

            if degeneracy_matrices_2 is None:
                master_log_file.write("No Degeneracies Found for Oligomer 2\n\n")
            elif len(degeneracy_matrices_2) == 1:
                master_log_file.write("1 Degeneracy Found for Oligomer 2\n")
                master_log_file.write("%s\n\n" % str(degeneracy_matrices_2[0]))
            else:
                master_log_file.write("%s Degeneracies Found for Oligomer 2\n" % str(len(degeneracy_matrices_2)))
                for degen_mat_2 in degeneracy_matrices_2:
                    master_log_file.write("%s\n" % str(degen_mat_2))
                master_log_file.write("\n")

            master_log_file.write("LOADING COMPLETE INTERFACE FRAGMENT REPRESENTATIVES DATABASE\n")
            # Create fragment database for all ijk cluster representatives
            ijk_frag_db = FragmentDB(monofrag_cluster_rep_dirpath,
                                     ijk_intfrag_cluster_rep_dirpath,
                                     intfrag_cluster_info_dirpath)
            # Get complete IJK fragment representatives database dictionaries
            ijk_monofrag_cluster_rep_pdb_dict = ijk_frag_db.get_monofrag_cluster_rep_dict()
            ijk_intfrag_cluster_rep_dict = ijk_frag_db.get_intfrag_cluster_rep_dict()
            ijk_intfrag_cluster_info_dict = ijk_frag_db.get_intfrag_cluster_info_dict()

            if init_match_type == "1_2":
                master_log_file.write("RETRIEVING HELIX-STRAND INTERFACE FRAGMENT REPRESENTATIVES FOR INITIAL MATCH\n\n")
                # Get Helix-Strand fragment representatives database dictionaries for initial interface fragment matching
                init_monofrag_cluster_rep_pdb_dict_1 = {"1": ijk_monofrag_cluster_rep_pdb_dict["1"]}
                init_monofrag_cluster_rep_pdb_dict_2 = {"2": ijk_monofrag_cluster_rep_pdb_dict["2"]}
                init_intfrag_cluster_rep_dict = {"1": {"2": ijk_intfrag_cluster_rep_dict["1"]["2"]}}
                init_intfrag_cluster_info_dict = {"1": {"2": ijk_intfrag_cluster_info_dict["1"]["2"]}}
            elif init_match_type == "2_1":
                master_log_file.write("RETRIEVING STRAND-HELIX INTERFACE FRAGMENT REPRESENTATIVES FOR INITIAL MATCH\n\n")
                # Get Strand-Helix fragment representatives database dictionaries for initial interface fragment matching
                init_monofrag_cluster_rep_pdb_dict_1 = {"2": ijk_monofrag_cluster_rep_pdb_dict["2"]}
                init_monofrag_cluster_rep_pdb_dict_2 = {"1": ijk_monofrag_cluster_rep_pdb_dict["1"]}
                init_intfrag_cluster_rep_dict = {"2": {"1": ijk_intfrag_cluster_rep_dict["2"]["1"]}}
                init_intfrag_cluster_info_dict = {"2": {"1": ijk_intfrag_cluster_info_dict["2"]["1"]}}
            elif init_match_type == "2_2":
                master_log_file.write("RETRIEVING STRAND-STRAND INTERFACE FRAGMENT REPRESENTATIVES FOR INITIAL MATCH\n\n")
                # Get Strand-Strand fragment representatives database dictionaries for initial interface fragment matching
                init_monofrag_cluster_rep_pdb_dict_1 = {"2": ijk_monofrag_cluster_rep_pdb_dict["2"]}
                init_monofrag_cluster_rep_pdb_dict_2 = {"2": ijk_monofrag_cluster_rep_pdb_dict["2"]}
                init_intfrag_cluster_rep_dict = {"2": {"2": ijk_intfrag_cluster_rep_dict["2"]["2"]}}
                init_intfrag_cluster_info_dict = {"2": {"2": ijk_intfrag_cluster_info_dict["2"]["2"]}}
            else:
                master_log_file.write("RETRIEVING HELIX-HELIX INTERFACE FRAGMENT REPRESENTATIVES FOR INITIAL MATCH\n\n")
                # Get Helix-Helix fragment representatives database dictionaries for initial interface fragment matching
                init_monofrag_cluster_rep_pdb_dict_1 = {"1": ijk_monofrag_cluster_rep_pdb_dict["1"]}
                init_monofrag_cluster_rep_pdb_dict_2 = {"1": ijk_monofrag_cluster_rep_pdb_dict["1"]}
                init_intfrag_cluster_rep_dict = {"1": {"1": ijk_intfrag_cluster_rep_dict["1"]["1"]}}
                init_intfrag_cluster_info_dict = {"1": {"1": ijk_intfrag_cluster_info_dict["1"]["1"]}}
            master_log_file.close()

            # Initialize Euler Lookup Class
            eul_lookup = EulerLookup()

            # Get Expand Matrices
            if design_dim == 0:
                expand_matrices = get_ptgrp_sym_op(result_design_sym)
            elif design_dim == 2:
                expand_matrices = get_sg_sym_op(result_design_sym.upper())
            elif design_dim == 3:
                expand_matrices = get_sg_sym_op(result_design_sym)
            else:
                master_log_file = open(master_log_filepath, "a+")
                master_log_file.write(
                    "\n%s is an Invalid Design Dimension. The Only Valid Dimensions are: 0, 2, 3\n" % str(design_dim))
                master_log_file.close()
                sys.exit()

            master_log_file = open(master_log_filepath, "a+")
            master_log_file.write("\nSTARTING FRAGMENT-BASED SYMMETRY DOCKING PROTOCOL\n\n")
            master_log_file.close()

            if oligomer_symmetry_1 == oligomer_symmetry_2:
                n = len(pdb1_oriented_filepaths)
                for i in range(n - 1):
                    pdb1_oriented_path = pdb1_oriented_filepaths[i]
                    pdb1_oriented_filename = os.path.splitext(os.path.basename(pdb1_oriented_path))[0]

                    for j in range(i + 1, n):
                        pdb2_oriented_path = pdb1_oriented_filepaths[j]
                        pdb2_oriented_filename = os.path.splitext(os.path.basename(pdb2_oriented_path))[0]

                        master_log_file = open(master_log_filepath, "a+")
                        master_log_file.write("DOCKING %s | %s\n" % (pdb1_oriented_filename, pdb2_oriented_filename))
                        master_log_file.close()

                        dock(init_intfrag_cluster_rep_dict, ijk_intfrag_cluster_rep_dict,
                             init_monofrag_cluster_rep_pdb_dict_1, init_monofrag_cluster_rep_pdb_dict_2,
                             init_intfrag_cluster_info_dict, ijk_monofrag_cluster_rep_pdb_dict,
                             ijk_intfrag_cluster_info_dict, free_sasa_exe_path, master_outdir, pdb1_oriented_path,
                             pdb2_oriented_path, set_mat1, set_mat2, ref_frame_tx_dof1, ref_frame_tx_dof2,
                             is_internal_zshift1, is_internal_zshift2, result_design_sym, uc_spec_string,
                             design_dim, expand_matrices, eul_lookup, init_max_z_val, subseq_max_z_val,
                             degeneracy_matrices_1, degeneracy_matrices_2, rot_step_deg1,
                             rot_range_deg_pdb1, rot_step_deg2, rot_range_deg_pdb2, output_exp_assembly,
                             output_uc, output_surrounding_uc, min_matched)

                        master_log_file = open(master_log_filepath, "a+")
                        master_log_file.write(
                            "COMPLETE ==> %s\n\n" %
                            (master_outdir + "/" + pdb1_oriented_filename + "_" + pdb2_oriented_filename))
                        master_log_file.close()

            else:
                for pdb1_oriented_path in pdb1_oriented_filepaths:
                    pdb1_oriented_filename = os.path.splitext(os.path.basename(pdb1_oriented_path))[0]

                    for pdb2_oriented_path in pdb2_oriented_filepaths:
                        pdb2_oriented_filename = os.path.splitext(os.path.basename(pdb2_oriented_path))[0]

                        master_log_file = open(master_log_filepath, "a+")
                        master_log_file.write("DOCKING %s | %s\n" % (pdb1_oriented_filename, pdb2_oriented_filename))
                        master_log_file.close()

                        dock(init_intfrag_cluster_rep_dict, ijk_intfrag_cluster_rep_dict,
                             init_monofrag_cluster_rep_pdb_dict_1, init_monofrag_cluster_rep_pdb_dict_2,
                             init_intfrag_cluster_info_dict, ijk_monofrag_cluster_rep_pdb_dict,
                             ijk_intfrag_cluster_info_dict, free_sasa_exe_path, master_outdir, pdb1_oriented_path,
                             pdb2_oriented_path, set_mat1, set_mat2, ref_frame_tx_dof1, ref_frame_tx_dof2,
                             is_internal_zshift1, is_internal_zshift2, result_design_sym, uc_spec_string,
                             design_dim, expand_matrices, eul_lookup, init_max_z_val, subseq_max_z_val,
                             degeneracy_matrices_1, degeneracy_matrices_2, rot_step_deg1,
                             rot_range_deg_pdb1, rot_step_deg2, rot_range_deg_pdb2, output_exp_assembly,
                             output_uc, output_surrounding_uc, min_matched)

                        master_log_file = open(master_log_filepath, "a+")
                        master_log_file.write(
                            "COMPLETE ==> %s\n\n" %
                            (master_outdir + "/" + pdb1_oriented_filename + "_" + pdb2_oriented_filename))
                        master_log_file.close()

            master_log_file = open(master_log_filepath, "a+")
            master_log_file.write("\nCOMPLETED FRAGMENT-BASED SYMMETRY DOCKING PROTOCOL\n\n")
            master_log_file.write("DONE\n")
            master_log_file.close()
            return 0

        except KeyboardInterrupt:
            master_log_file = open(master_log_filepath, "a+")
            master_log_file.write("\nRun Ended By KeyboardInterrupt\n")
            master_log_file.close()
            sys.exit()

    elif len(cmd_line_in_params) > 1 and cmd_line_in_params[1] == '-query':
        query_mode(cmd_line_in_params)

    elif len(cmd_line_in_params) > 1 and cmd_line_in_params[1] == '-postprocess':
        postprocess_mode(cmd_line_in_params)

    else:
        print_usage()


if __name__ == "__main__":
    main()

