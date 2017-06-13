"""

"""
import argparse
import os
import sys


def main():
    args = get_args()
    pipeline(**args.__dict__)
    return 0


def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i', '--input-dir', default='.', help='path to the input directory')
    arg_parser.add_argument('-o', '--output-dir', default='.', help='path to the output directory')
    arg_parser.add_argument('-c', '--core-count', default=1, type=int, help='number of cores to use')
    arg_parser.add_argument('--forward-primer', default='ATTAGAWACCCVNGTAGTCC', help='forward primer to be clipped')
    arg_parser.add_argument('--reverse-primer', default='TTACCGCGGCKGCTGGCAC', help='reverse primer to be clipped')
    arg_parser.add_argument('--uchime-ref-db-fp', default='/cluster_16S/pr2/pr2_gb203_version_4.5.fasta',
                            help='database for vsearch --uchime_ref')
    args = arg_parser.parse_args()
    return args


def pipeline(input_dir, output_dir, core_count, forward_primer, reverse_primer, uchime_ref_db_fp):

    step_01_output_dir = step_01_split_ITS_16S(input_dir=input_dir)
    step_02_output_dir = step_02_adjust_files(input_dir=step_01_output_dir)
    step_03_output_dir = step_03_remove_primers_with_cutadapt(input_dir=step_02_output_dir)
    step_04_output_dir = step_04_merge_forward_reverse_reads_with_pear(input_dir=step_03_output_dir)
    step_05_output_dir = step_05_qc_reads_with_vsearch(input_dir=step_04_output_dir)
    step_06_output_dir = step_06_combine_runs(input_dir=step_05_output_dir)
    step_07_output_dir = step_07_dereplicate_sort_remove_low_abundance_reads(input_dir=step_06_output_dir)
    step_08_output_dir = step_08_cluster_97_percent(input_dir=step_07_output_dir)
    step_09_output_dir = step_09_reference_based_chimera_detection(input_dir=step_08_output_dir)
    step_10_output_dir = step_10_map_raw_reads_to_otus(input_dir=step_09_output_dir)
    step_11_output_dir = step_11_write_otu_table(input_dir=step_10_output_dir)

    output_dir_list = [
        step_01_output_dir, step_02_output_dir, step_03_output_dir, step_04_output_dir,
        step_05_output_dir, step_06_output_dir, step_07_output_dir, step_08_output_dir,
        step_09_output_dir, step_10_output_dir, step_11_output_dir
    ]

    return output_dir_list


def step_01_split_ITS_16S(input_dir):
    output_dir = create_output_dir(input_dir, sys._getframe().f_code.co_name)
    return output_dir


def step_02_adjust_files(input_dir):
    output_dir = create_output_dir(input_dir, sys._getframe().f_code.co_name)
    return output_dir


def step_03_remove_primers_with_cutadapt(input_dir):
    output_dir = create_output_dir(input_dir, sys._getframe().f_code.co_name)
    return output_dir


def step_04_merge_forward_reverse_reads_with_pear(input_dir):
    output_dir = create_output_dir(input_dir, sys._getframe().f_code.co_name)
    return output_dir


def step_05_qc_reads_with_vsearch(input_dir):
    output_dir = create_output_dir(input_dir, sys._getframe().f_code.co_name)
    return output_dir


def step_06_combine_runs(input_dir):
    output_dir = create_output_dir(input_dir, sys._getframe().f_code.co_name)
    return output_dir


def step_07_dereplicate_sort_remove_low_abundance_reads(input_dir):
    output_dir = create_output_dir(input_dir, sys._getframe().f_code.co_name)
    return output_dir


def step_08_cluster_97_percent(input_dir):
    output_dir = create_output_dir(input_dir, sys._getframe().f_code.co_name)
    return output_dir


def step_09_reference_based_chimera_detection(input_dir):
    output_dir = create_output_dir(input_dir, sys._getframe().f_code.co_name)
    return output_dir


def step_10_map_raw_reads_to_otus(input_dir):
    output_dir = create_output_dir(input_dir, sys._getframe().f_code.co_name)
    return output_dir


def step_11_write_otu_table(input_dir):
    output_dir = create_output_dir(input_dir, sys._getframe().f_code.co_name)
    return output_dir


def create_output_dir(input_dir, output_dir_name):
    input_parent_dir, _ = os.path.split(input_dir)
    output_dir = os.path.join(input_parent_dir, output_dir_name)
    os.mkdir(output_dir)
    return output_dir


if __name__ == '__main__':
    main()