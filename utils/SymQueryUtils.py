from classes.SymEntry import sym_comb_dict


# Copyright 2020 Joshua Laniado and Todd O. Yeates.
__author__ = "Joshua Laniado and Todd O. Yeates"
__copyright__ = "Copyright 2020, Nanohedra"
__version__ = "1.0"


query_out_format_str = "{:>5s} {:>6s} {:^9s} {:^9s} {:^20s} {:>6s} {:^9s} {:^9s} {:^20s} {:>6s} {:>3s} {:>4s} {:>4s}"


def print_query_header():
    header_format_str = "{:5s} {:6s} {:^9s} {:^9s} {:^20s} {:6s} {:^9s} {:^9s} {:^20s} {:6s} {:3s} {:4s} {:4s}"
    print header_format_str.format("ENTRY",
                                   "PtGrp1",
                                   "IntRtDOF1",
                                   "IntTxDOF1",
                                   "ReferenceFrameDOF1",
                                   "PtGrp2",
                                   "IntRtDOF2",
                                   "IntTxDOF2",
                                   "ReferenceFrameDOF2",
                                   "RESULT",
                                   "DIM",
                                   "#DOF",
                                   "RING")


def query_combination(combination_list):
    if type(combination_list) == list and len(combination_list) == 2:
        matching_entries = []
        for entry_number in sym_comb_dict:
            group1 = sym_comb_dict[entry_number][1]
            group2 = sym_comb_dict[entry_number][6]
            if combination_list == [group1, group2] or combination_list == [group2, group1]:
                int_rot1 = "none"
                int_tx1 = "none"
                int_rot2 = "none"
                int_tx2 = "none"
                int_dof_group1 = sym_comb_dict[entry_number][3]
                int_dof_group2 = sym_comb_dict[entry_number][8]
                for int_dof in int_dof_group1:
                    if int_dof.startswith('r'):
                        int_rot1 = int_dof[2:]
                    if int_dof.startswith('t'):
                        int_tx1 = int_dof[2:]
                for int_dof in int_dof_group2:
                    if int_dof.startswith('r'):
                        int_rot2 = int_dof[2:]
                    if int_dof.startswith('t'):
                        int_tx2 = int_dof[2:]
                ref_frame_tx_dof_group1 = sym_comb_dict[entry_number][5]
                ref_frame_tx_dof_group2 = sym_comb_dict[entry_number][10]
                if ref_frame_tx_dof_group1 == '<0,0,0>':
                    ref_frame_tx_dof_group1 = 'none'
                if ref_frame_tx_dof_group2 == '<0,0,0>':
                    ref_frame_tx_dof_group2 = 'none'
                result = sym_comb_dict[entry_number][12]
                dim = sym_comb_dict[entry_number][13]
                tot_num_dof = sym_comb_dict[entry_number][15]
                ring_size = sym_comb_dict[entry_number][16]
                matching_entries.append(query_out_format_str.format(str(entry_number),
                                                                    group1,
                                                                    int_rot1,
                                                                    int_tx1,
                                                                    ref_frame_tx_dof_group1,
                                                                    group2,
                                                                    int_rot2,
                                                                    int_tx2,
                                                                    ref_frame_tx_dof_group2,
                                                                    result,
                                                                    str(dim),
                                                                    str(tot_num_dof),
                                                                    str(ring_size)))
        if matching_entries == list():
            print '\033[1m' + "NO MATCHING ENTRY FOUND" + '\033[0m'
            print ''
        else:
            print '\033[1m' + "POSSIBLE COMBINATION(S) FOR: %s & %s" % (combination_list[0], combination_list[1]) + '\033[0m'
            print_query_header()
            for match in matching_entries:
                print match
    else:
        print "INVALID ENTRY"


def query_result(desired_result):
    if type(desired_result) == str:
        matching_entries = []
        for entry_number in sym_comb_dict:
            result = sym_comb_dict[entry_number][12]
            if desired_result == result:
                group1 = sym_comb_dict[entry_number][1]
                group2 = sym_comb_dict[entry_number][6]
                int_rot1 = "none"
                int_tx1 = "none"
                int_rot2 = "none"
                int_tx2 = "none"
                int_dof_group1 = sym_comb_dict[entry_number][3]
                int_dof_group2 = sym_comb_dict[entry_number][8]
                for int_dof in int_dof_group1:
                    if int_dof.startswith('r'):
                        int_rot1 = int_dof[2:]
                    if int_dof.startswith('t'):
                        int_tx1 = int_dof[2:]
                for int_dof in int_dof_group2:
                    if int_dof.startswith('r'):
                        int_rot2 = int_dof[2:]
                    if int_dof.startswith('t'):
                        int_tx2 = int_dof[2:]
                ref_frame_tx_dof_group1 = sym_comb_dict[entry_number][5]
                ref_frame_tx_dof_group2 = sym_comb_dict[entry_number][10]
                if ref_frame_tx_dof_group1 == '<0,0,0>':
                    ref_frame_tx_dof_group1 = 'none'
                if ref_frame_tx_dof_group2 == '<0,0,0>':
                    ref_frame_tx_dof_group2 = 'none'
                dim = sym_comb_dict[entry_number][13]
                tot_num_dof = sym_comb_dict[entry_number][15]
                ring_size = sym_comb_dict[entry_number][16]
                matching_entries.append(query_out_format_str.format(str(entry_number),
                                                                    group1,
                                                                    int_rot1,
                                                                    int_tx1,
                                                                    ref_frame_tx_dof_group1,
                                                                    group2,
                                                                    int_rot2,
                                                                    int_tx2,
                                                                    ref_frame_tx_dof_group2,
                                                                    result,
                                                                    str(dim),
                                                                    str(tot_num_dof),
                                                                    str(ring_size)))
        if matching_entries == list():
            print '\033[1m' + "NO MATCHING ENTRY FOUND" + '\033[0m'
            print ''
        else:
            print '\033[1m' + "POSSIBLE COMBINATION(S) FOR: %s" % desired_result + '\033[0m'
            print_query_header()
            for match in matching_entries:
                print match
    else:
        print "INVALID ENTRY"


def query_counterpart(query_group):
    if type(query_group) == str:
        matching_entries = []
        for entry_number in sym_comb_dict:
            group1 = sym_comb_dict[entry_number][1]
            group2 = sym_comb_dict[entry_number][6]
            if query_group in [group1, group2]:
                int_rot1 = "none"
                int_tx1 = "none"
                int_rot2 = "none"
                int_tx2 = "none"
                int_dof_group1 = sym_comb_dict[entry_number][3]
                int_dof_group2 = sym_comb_dict[entry_number][8]
                for int_dof in int_dof_group1:
                    if int_dof.startswith('r'):
                        int_rot1 = int_dof[2:]
                    if int_dof.startswith('t'):
                        int_tx1 = int_dof[2:]
                for int_dof in int_dof_group2:
                    if int_dof.startswith('r'):
                        int_rot2 = int_dof[2:]
                    if int_dof.startswith('t'):
                        int_tx2 = int_dof[2:]
                ref_frame_tx_dof_group1 = sym_comb_dict[entry_number][5]
                ref_frame_tx_dof_group2 = sym_comb_dict[entry_number][10]
                if ref_frame_tx_dof_group1 == '<0,0,0>':
                    ref_frame_tx_dof_group1 = 'none'
                if ref_frame_tx_dof_group2 == '<0,0,0>':
                    ref_frame_tx_dof_group2 = 'none'
                result = sym_comb_dict[entry_number][12]
                dim = sym_comb_dict[entry_number][13]
                tot_num_dof = sym_comb_dict[entry_number][15]
                ring_size = sym_comb_dict[entry_number][16]
                matching_entries.append(query_out_format_str.format(str(entry_number),
                                                                    group1,
                                                                    int_rot1,
                                                                    int_tx1,
                                                                    ref_frame_tx_dof_group1,
                                                                    group2,
                                                                    int_rot2,
                                                                    int_tx2,
                                                                    ref_frame_tx_dof_group2,
                                                                    result,
                                                                    str(dim),
                                                                    str(tot_num_dof),
                                                                    str(ring_size)))
        if matching_entries == list():
            print '\033[1m' + "NO MATCHING ENTRY FOUND" + '\033[0m'
            print ''
        else:
            print '\033[1m' + "POSSIBLE COMBINATION(S) FOR: %s" % query_group + '\033[0m'
            print_query_header()
            for match in matching_entries:
                print match
    else:
        print "INVALID ENTRY"


def all_entries():
    all_entries_list = []
    for entry_number in sym_comb_dict:
        group1 = sym_comb_dict[entry_number][1]
        group2 = sym_comb_dict[entry_number][6]
        int_rot1 = "none"
        int_tx1 = "none"
        int_rot2 = "none"
        int_tx2 = "none"
        int_dof_group1 = sym_comb_dict[entry_number][3]
        int_dof_group2 = sym_comb_dict[entry_number][8]
        for int_dof in int_dof_group1:
            if int_dof.startswith('r'):
                int_rot1 = int_dof[2:]
            if int_dof.startswith('t'):
                int_tx1 = int_dof[2:]
        for int_dof in int_dof_group2:
            if int_dof.startswith('r'):
                int_rot2 = int_dof[2:]
            if int_dof.startswith('t'):
                int_tx2 = int_dof[2:]
        ref_frame_tx_dof_group1 = sym_comb_dict[entry_number][5]
        ref_frame_tx_dof_group2 = sym_comb_dict[entry_number][10]
        if ref_frame_tx_dof_group1 == '<0,0,0>':
            ref_frame_tx_dof_group1 = 'none'
        if ref_frame_tx_dof_group2 == '<0,0,0>':
            ref_frame_tx_dof_group2 = 'none'
        result = sym_comb_dict[entry_number][12]
        dim = sym_comb_dict[entry_number][13]
        tot_num_dof = sym_comb_dict[entry_number][15]
        ring_size = sym_comb_dict[entry_number][16]
        all_entries_list.append(query_out_format_str.format(str(entry_number),
                                                            group1,
                                                            int_rot1,
                                                            int_tx1,
                                                            ref_frame_tx_dof_group1,
                                                            group2,
                                                            int_rot2,
                                                            int_tx2,
                                                            ref_frame_tx_dof_group2,
                                                            result,
                                                            str(dim),
                                                            str(tot_num_dof),
                                                            str(ring_size)))
    print '\033[1m' + "ALL ENTRIES" + '\033[0m'
    print_query_header()
    for entry in all_entries_list:
        print entry


def dimension(dim):
    if dim in [0, 2, 3]:
        matching_entries = []
        for entry_number in sym_comb_dict:
            if sym_comb_dict[entry_number][13] == dim:
                group1 = sym_comb_dict[entry_number][1]
                group2 = sym_comb_dict[entry_number][6]
                int_rot1 = "none"
                int_tx1 = "none"
                int_rot2 = "none"
                int_tx2 = "none"
                int_dof_group1 = sym_comb_dict[entry_number][3]
                int_dof_group2 = sym_comb_dict[entry_number][8]
                for int_dof in int_dof_group1:
                    if int_dof.startswith('r'):
                        int_rot1 = int_dof[2:]
                    if int_dof.startswith('t'):
                        int_tx1 = int_dof[2:]
                for int_dof in int_dof_group2:
                    if int_dof.startswith('r'):
                        int_rot2 = int_dof[2:]
                    if int_dof.startswith('t'):
                        int_tx2 = int_dof[2:]
                ref_frame_tx_dof_group1 = sym_comb_dict[entry_number][5]
                ref_frame_tx_dof_group2 = sym_comb_dict[entry_number][10]
                if ref_frame_tx_dof_group1 == '<0,0,0>':
                    ref_frame_tx_dof_group1 = 'none'
                if ref_frame_tx_dof_group2 == '<0,0,0>':
                    ref_frame_tx_dof_group2 = 'none'
                result = sym_comb_dict[entry_number][12]
                dim = sym_comb_dict[entry_number][13]
                tot_num_dof = sym_comb_dict[entry_number][15]
                ring_size = sym_comb_dict[entry_number][16]
                matching_entries.append(query_out_format_str.format(str(entry_number),
                                                                    group1,
                                                                    int_rot1,
                                                                    int_tx1,
                                                                    ref_frame_tx_dof_group1,
                                                                    group2,
                                                                    int_rot2,
                                                                    int_tx2,
                                                                    ref_frame_tx_dof_group2,
                                                                    result,
                                                                    str(dim),
                                                                    str(tot_num_dof),
                                                                    str(ring_size)))

        print '\033[1m' + "ALL ENTRIES FOUND WITH DIMENSION " + str(dim) + ": " + '\033[0m'
        print_query_header()
        for entry in matching_entries:
            print entry
    else:
        print "DIMENSION NOT SUPPORTED, VALID DIMENSIONS ARE: 0, 2 or 3 "

