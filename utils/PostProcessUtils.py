import os
import shutil


# Copyright 2020 Joshua Laniado and Todd O. Yeates.
__author__ = "Joshua Laniado and Todd O. Yeates"
__copyright__ = "Copyright 2020, Nanohedra"
__version__ = "1.0"


def frag_match_count_filter(master_design_dirpath, min_frag_match_count, master_design_outdir_path):
    for root1, dirs1, files1 in os.walk(master_design_dirpath):
        for file1 in files1:
            if "docked_pose_info_file.txt" in file1:
                info_file_filepath = root1 + "/" + file1

                tx_filepath = root1
                rot_filepath= os.path.dirname(tx_filepath)
                degen_filepath = os.path.dirname(rot_filepath)
                design_filepath = os.path.dirname(degen_filepath)

                tx_fname = tx_filepath.split("/")[-1]
                rot_fname = rot_filepath.split("/")[-1]
                degen_fname = degen_filepath.split("/")[-1]
                design_fname = design_filepath.split("/")[-1]

                outdir = master_design_outdir_path + "/" + design_fname + "/" + degen_fname + "/" + rot_fname + "/" + tx_fname

                info_file = open(info_file_filepath, 'r')
                for line in info_file.readlines():
                    if "Unique Mono Fragments Matched:" in line:
                        frag_match_count = int(line[30:])
                        if frag_match_count >= min_frag_match_count:
                            shutil.copytree(tx_filepath, outdir)
                info_file.close()


def score_filter(master_design_dirpath, min_score, master_design_outdir_path):
    for root1, dirs1, files1 in os.walk(master_design_dirpath):
        for file1 in files1:
            if "docked_pose_info_file.txt" in file1:
                info_file_filepath = root1 + "/" + file1

                tx_filepath = root1
                rot_filepath = os.path.dirname(tx_filepath)
                degen_filepath = os.path.dirname(rot_filepath)
                design_filepath = os.path.dirname(degen_filepath)

                tx_fname = tx_filepath.split("/")[-1]
                rot_fname = rot_filepath.split("/")[-1]
                degen_fname = degen_filepath.split("/")[-1]
                design_fname = design_filepath.split("/")[-1]

                outdir = master_design_outdir_path + "/" + design_fname + "/" + degen_fname + "/" + rot_fname + "/" + tx_fname

                info_file = open(info_file_filepath, 'r')
                for line in info_file.readlines():
                    if "Nanohedra Score:" in line:
                        score = float(line[17:])
                        if score >= min_score:
                            shutil.copytree(tx_filepath, outdir)
                info_file.close()


def score_and_frag_match_count_filter(master_design_dirpath, min_score, min_frag_match_count, master_design_outdir_path):
    for root1, dirs1, files1 in os.walk(master_design_dirpath):
        for file1 in files1:
            if "docked_pose_info_file.txt" in file1:
                info_file_filepath = root1 + "/" + file1

                tx_filepath = root1
                rot_filepath = os.path.dirname(tx_filepath)
                degen_filepath = os.path.dirname(rot_filepath)
                design_filepath = os.path.dirname(degen_filepath)

                tx_fname = tx_filepath.split("/")[-1]
                rot_fname = rot_filepath.split("/")[-1]
                degen_fname = degen_filepath.split("/")[-1]
                design_fname = design_filepath.split("/")[-1]

                outdir = master_design_outdir_path + "/" + design_fname + "/" + degen_fname + "/" + rot_fname + "/" + tx_fname

                score = None
                frag_match_count = None
                info_file = open(info_file_filepath, 'r')
                for line in info_file.readlines():
                    if "Nanohedra Score:" in line:
                        score = float(line[17:])
                    if "Unique Mono Fragments Matched:" in line:
                        frag_match_count = int(line[30:])
                info_file.close()

                if score is not None and frag_match_count is not None:
                    if score >= min_score and frag_match_count >= min_frag_match_count:
                        shutil.copytree(tx_filepath, outdir)


def rank(master_design_dirpath, metric, outdir):

    if metric == 'score':
        metric_str = "Nanohedra Score:"
    elif metric == 'matched':
        metric_str = "Unique Mono Fragments Matched:"
    else:
        raise ValueError('\n%s is not a recognized ranking metric. '
                         'Recognized ranking metrics are: score and matched.\n' %str(metric))

    designpath_metric_tup_list = []

    for root1, dirs1, files1 in os.walk(master_design_dirpath):
        for file1 in files1:
            if "docked_pose_info_file.txt" in file1:
                info_file_filepath = root1 + "/" + file1

                tx_filepath = root1
                rot_filepath = os.path.dirname(tx_filepath)
                degen_filepath = os.path.dirname(rot_filepath)
                design_filepath = os.path.dirname(degen_filepath)

                tx_fname = tx_filepath.split("/")[-1]
                rot_fname = rot_filepath.split("/")[-1]
                degen_fname = degen_filepath.split("/")[-1]
                design_fname = design_filepath.split("/")[-1]

                design_path = "/" + design_fname + "/" + degen_fname + "/" + rot_fname + "/" + tx_fname

                if metric == 'score':
                    info_file = open(info_file_filepath, 'r')
                    for line in info_file.readlines():
                        if metric_str in line:
                            score = float(line[17:])
                            designpath_metric_tup_list.append((design_path, score))
                    info_file.close()

                elif metric == 'matched':
                    info_file = open(info_file_filepath, 'r')
                    for line in info_file.readlines():
                        if metric_str in line:
                            frag_match_count = int(line[30:])
                            designpath_metric_tup_list.append((design_path, frag_match_count))
                    info_file.close()

    designpath_metric_tup_list_sorted = sorted(designpath_metric_tup_list, key=lambda tup: tup[1], reverse=True)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outfile = open(outdir + "/ranked_designs_%s.txt" % metric, 'w')
    for p, m in designpath_metric_tup_list_sorted:
        outfile.write("%s\t%s\n" % (str(p), str(m)))
    outfile.close()

