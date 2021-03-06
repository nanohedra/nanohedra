import numpy as np
import math


# Copyright 2020 Joshua Laniado and Todd O. Yeates.
__author__ = "Joshua Laniado and Todd O. Yeates"
__copyright__ = "Copyright 2020, Nanohedra"
__version__ = "1.0"


# ROTATION RANGE DEG
C2 = 180
C3 = 120
C4 = 90
C5 = 72
C6 = 60
RotRangeDict = {"C2": C2, "C3": C3, "C4": C4, "C5": C5, "C6": C6}


def get_degeneracy_matrices(oligomer_symmetry_1, oligomer_symmetry_2, design_dimension, design_symmetry):
    valid_pt_gp_symm_list = ["C2", "C3", "C4", "C5", "C6", "D2", "D3", "D4", "D6", "T", "O", "I"]

    if oligomer_symmetry_1 not in valid_pt_gp_symm_list:
        raise ValueError("Invalid Point Group Symmetry")

    if oligomer_symmetry_2 not in valid_pt_gp_symm_list:
        raise ValueError("Invalid Point Group Symmetry")

    if design_symmetry not in valid_pt_gp_symm_list:
        raise ValueError("Invalid Point Group Symmetry")

    if design_dimension not in [0, 2, 3]:
        raise ValueError("Invalid Design Dimension")

    degeneracies = [None, None]

    for i in range(2):

        degeneracy_matrices = None

        oligomer_symmetry = oligomer_symmetry_1 if i == 0 else oligomer_symmetry_2

        # For cages, only one of the two oligomers need to be flipped. By convention we flip oligomer 2.
        if design_dimension == 0 and i == 1:
            degeneracy_matrices = [[[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]]  # ROT180y

        # For layers that obey a cyclic point group symmetry
        # and that are constructed from two oligomers that both obey cyclic symmetry
        # only one of the two oligomers need to be flipped. By convention we flip oligomer 2.
        elif design_dimension == 2 and i == 1 and (oligomer_symmetry_1[0], oligomer_symmetry_2[0], design_symmetry[0]) == ("C", "C", "C"):
            degeneracy_matrices = [[[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]]  # ROT180y

        elif oligomer_symmetry in ["D3", "D4", "D6"] and design_symmetry in ["D3", "D4", "D6", "T", "O"]:
            if oligomer_symmetry == "D3":
                degeneracy_matrices = [[[0.5, -0.86603, 0.0], [0.86603, 0.5, 0.0], [0.0, 0.0, 1.0]]]  # ROT60z
            elif oligomer_symmetry == "D4":
                # 45 degrees about z; z unaffected; x goes to [1,-1,0] direction
                degeneracy_matrices = [[[0.707107, 0.707107, 0.0], [-0.707107, 0.707107, 0.0], [0.0, 0.0, 1.0]]]
            elif oligomer_symmetry == "D6":
                degeneracy_matrices = [[[0.86603, -0.5, 0.0], [0.5, 0.86603, 0.0], [0.0, 0.0, 1.0]]]  # ROT30z

        elif oligomer_symmetry == "D2" and design_symmetry != "O":
            if design_symmetry == "T":
                degeneracy_matrices = [[[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]]  # ROT90z

            elif design_symmetry == "D4":
                degeneracy_matrices = [[[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                                       [[0.0, 1.0, 0.0], [0.0,  0.0, 1.0], [1.0, 0.0, 0.0]]]  # z,x,y and y,z,x

            elif design_symmetry == "D2" or design_symmetry == "D6":
                degeneracy_matrices = [[[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                                       [[0.0, 1.0, 0.0], [0.0,  0.0, 1.0], [1.0, 0.0, 0.0]],
                                       [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
                                       [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]],
                                       [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]]

        elif oligomer_symmetry == "T" and design_symmetry == "T":
            degeneracy_matrices = [[[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]]  # ROT90z

        degeneracies[i] = degeneracy_matrices

    return degeneracies


def parse_ref_tx_dof_str_to_list(ref_frame_tx_dof_string):
    s1 = ref_frame_tx_dof_string.replace('<', '')
    s2 = s1.replace('>', '')
    l1 = s2.split(',')
    l2 = [x.replace(' ', '') for x in l1]
    return l2


def get_tx_dof_ref_frame_var_vec(string_vec, var):
    return_vec = [0.0, 0.0, 0.0]
    for i in range(3):
        if var in string_vec[i] and '*' in string_vec[i]:
            return_vec[i] = float(string_vec[i].split('*')[0])
        elif "-" + var in string_vec[i]:
            return_vec[i] = -1.0
        elif var == string_vec[i]:
            return_vec[i] = 1.0
    return return_vec


def get_ext_dof(ref_frame_tx_dof1, ref_frame_tx_dof2):

    ext_dof = []

    parsed_1 = parse_ref_tx_dof_str_to_list(ref_frame_tx_dof1)
    parsed_2 = parse_ref_tx_dof_str_to_list(ref_frame_tx_dof2)

    e1_var_vec = get_tx_dof_ref_frame_var_vec(parsed_1, 'e')
    f1_var_vec = get_tx_dof_ref_frame_var_vec(parsed_1, 'f')
    g1_var_vec = get_tx_dof_ref_frame_var_vec(parsed_1, 'g')

    e2_var_vec = get_tx_dof_ref_frame_var_vec(parsed_2, 'e')
    f2_var_vec = get_tx_dof_ref_frame_var_vec(parsed_2, 'f')
    g2_var_vec = get_tx_dof_ref_frame_var_vec(parsed_2, 'g')

    e2e1_diff = (np.array(e2_var_vec) - np.array(e1_var_vec)).tolist()
    f2f1_diff = (np.array(f2_var_vec) - np.array(f1_var_vec)).tolist()
    g2g1_diff = (np.array(g2_var_vec) - np.array(g1_var_vec)).tolist()

    if e2e1_diff != [0, 0, 0]:
        ext_dof.append(e2e1_diff)

    if f2f1_diff != [0, 0, 0]:
        ext_dof.append(f2f1_diff)

    if g2g1_diff != [0, 0, 0]:
        ext_dof.append(g2g1_diff)

    return ext_dof


def get_optimal_external_tx_vector(ref_frame_tx_dof, optimal_ext_dof_shifts):

    ext_dof_variables = ['e', 'f', 'g']

    parsed_ref_tx_vec = parse_ref_tx_dof_str_to_list(ref_frame_tx_dof)

    optimal_external_tx_vector = np.array([0.0, 0.0, 0.0])
    for dof_shift_index in range(len(optimal_ext_dof_shifts)):
        dof_shift = optimal_ext_dof_shifts[dof_shift_index]
        var_vec = get_tx_dof_ref_frame_var_vec(parsed_ref_tx_vec, ext_dof_variables[dof_shift_index])
        shifted_var_vec = np.array(var_vec) * dof_shift
        optimal_external_tx_vector += shifted_var_vec

    return optimal_external_tx_vector.tolist()


def get_rot_matrices(step_deg, axis, rot_range_deg):
    rot_matrices = []
    if axis == 'x':
        for angle_deg in range(0, rot_range_deg, step_deg):
            rad = math.radians(float(angle_deg))
            rotmatrix = [[1, 0, 0], [0, math.cos(rad), -1 * math.sin(rad)], [0, math.sin(rad), math.cos(rad)]]
            rot_matrices.append(rotmatrix)
        return rot_matrices

    elif axis == 'y':
        for angle_deg in range(0, rot_range_deg, step_deg):
            rad = math.radians(float(angle_deg))
            rotmatrix = [[math.cos(rad), 0, math.sin(rad)], [0, 1, 0], [-1 * math.sin(rad), 0, math.cos(rad)]]
            rot_matrices.append(rotmatrix)
        return rot_matrices

    elif axis == 'z':
        for angle_deg in range(0, rot_range_deg, step_deg):
            rad = math.radians(float(angle_deg))
            rotmatrix = [[math.cos(rad), -1 * math.sin(rad), 0], [math.sin(rad), math.cos(rad), 0], [0, 0, 1]]
            rot_matrices.append(rotmatrix)
        return rot_matrices

    else:
        print "AXIS SELECTED FOR SAMPLING IS NOT SUPPORTED"
        return None


def get_degen_rotmatrices(degeneracy_matrices, rotation_matrices):
    if rotation_matrices == list() and degeneracy_matrices is not None:
        identity_matrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        return [[identity_matrix]] + [[degen_mat] for degen_mat in degeneracy_matrices]

    elif rotation_matrices != list() and degeneracy_matrices is None:
        return [rotation_matrices]

    elif rotation_matrices != list() and degeneracy_matrices is not None:
        degen_rotmatrices = [rotation_matrices]
        for degen in degeneracy_matrices:
            degen_list = []
            for rot in rotation_matrices:
                combined = np.matmul(rot, degen)
                degen_list.append(combined.tolist())
            degen_rotmatrices.append(degen_list)
        return degen_rotmatrices
    else:
        identity_matrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        return [[identity_matrix]]
