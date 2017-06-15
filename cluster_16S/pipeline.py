"""

"""
import argparse
import glob
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import traceback

from Bio.Alphabet import generic_dna
from Bio import SeqIO


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

    output_dir_list = []
    output_dir_list.append(step_01_copy_and_compress(input_dir=input_dir, work_dir=output_dir))
    ##output_dir_list.append(step_02_adjust_headers(input_dir=output_dir_list[-1]))
    output_dir_list.append(step_03_remove_primers(input_dir=output_dir_list[-1]))
    output_dir_list.append(step_04_merge_forward_reverse_reads_with_pear(input_dir=output_dir_list[-1]))
    output_dir_list.append(step_05_qc_reads_with_vsearch(input_dir=output_dir_list[-1]))
    output_dir_list.append(step_06_combine_runs(input_dir=output_dir_list[-1]))
    output_dir_list.append(step_07_dereplicate_sort_remove_low_abundance_reads(input_dir=output_dir_list[-1]))
    output_dir_list.append(step_08_cluster_97_percent(input_dir=output_dir_list[-1]))
    output_dir_list.append(step_09_reference_based_chimera_detection(input_dir=output_dir_list[-1]))
    output_dir_list.append(step_10_map_raw_reads_to_otus(input_dir=output_dir_list[-1]))
    output_dir_list.append(step_11_write_otu_table(input_dir=output_dir_list[-1]))

    return output_dir_list


def run_cmd(cmd_line_list, **kwargs):
    log = logging.getLogger(name=__name__)
    try:
        log.debug('executing "%s"', ' '.join((str(x) for x in cmd_line_list)))
        output = subprocess.run(
            cmd_line_list,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            **kwargs)
        log.debug(output)
        return output
    except subprocess.CalledProcessError as c:
        logging.exception(c)
        print(c.message)
        print(c.cmd)
        print(c.output)
        raise c
    except Exception as e:
        logging.exception(e)
        print('blarg!')
        print(e)
        traceback.print_exc()
        raise e


def step_01_copy_and_compress(input_dir, work_dir):
    function_name = sys._getframe().f_code.co_name
    log = logging.getLogger(name=function_name)
    output_dir = create_output_dir(output_dir_name=function_name, parent_dir=work_dir)
    log.debug('output_dir: %s', output_dir)

    input_file_glob = os.path.join(input_dir, '*.fastq*')
    log.debug('input_file_glob: %s', input_file_glob)
    input_fp_list = glob.glob(input_file_glob)
    log.info('input files: %s', input_fp_list)

    for input_fp in input_fp_list:
        destination_fp = os.path.join(output_dir, os.path.basename(input_fp))
        log.debug('copying file "%s" to "%s"', input_fp, destination_fp)
        shutil.copyfile(src=input_fp, dst=destination_fp)
        if input_fp.endswith('.fastq'):
            gzipped_fp = destination_fp + '.gz'
            with open(destination_fp, 'rt') as f, gzip.open(gzipped_fp, 'wt') as g:
                shutil.copyfileobj(f, g)
            os.remove(destination_fp)

    return output_dir


def step_02_adjust_headers(input_dir):
    function_name = sys._getframe().f_code.co_name
    log = logging.getLogger(name=function_name)
    output_dir = create_output_dir(output_dir_name=function_name, input_dir=input_dir)

    input_file_glob = os.path.join(input_dir, '*.fastq.gz')
    input_file_list = glob.glob(input_file_glob)
    log.info('input files: %s', input_file_list)

    # the FASTQ headers look like this:
    #   R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:17273:1842 1:N:0
    # where 1:N:0 is the 'description', not part of the id
    fastq_header_pattern = re.compile(
        r'(?P<id>'
            r'^(?P<run>R\d+)'                     #  R3
            r'-(?P<locus>(16S|ITS))'              # -16S
            r'-(?P<sample_name>(mockE-\d+))'      # -mockE-1
            r'__(?P<unknown_1>[\w\d\-]+)'         # __HWI-M01380
            r':(?P<unknown_2>(\d+))'              # :86
            r':(?P<unknown_3>([\d\w\-]+))'        # :000000000-ALK4C
            r':(?P<unknown_4>(\d+))'              # :1
            r':(?P<unknown_5>(\d+))'              # :1101
            r':(?P<unknown_6>(\d+))'              # :17273
            r':(?P<unknown_7>(\d+))'              # :1842
        r')$')

    for input_fp in input_file_list:
        input_basename = os.path.basename(input_fp)
        output_16S_fp = os.path.join(output_dir, re.sub(string=input_basename, pattern='\.fastq\.gz$', repl='_16S.fastq.gz'))
        output_ITS_fp = os.path.join(output_dir, re.sub(string=input_basename, pattern='\.fastq\.gz$', repl='_ITS.fastq.gz'))
        log.info('16S output file: %s', output_16S_fp)
        log.info('ITS output file: %s', output_ITS_fp)

        with gzip.open(input_fp, 'rt') as input_file, \
            gzip.open(output_16S_fp, 'wt') as output_16S_file, gzip.open(output_ITS_fp, 'wt') as output_ITS_file:
                output_locus_file = { '16S': output_16S_file, 'ITS': output_ITS_file }
                for i, fastq_record in enumerate(SeqIO.parse(input_file, format='fastq', alphabet=generic_dna)):
                    m = fastq_header_pattern.match(fastq_record.id)
                    if m is None:
                        raise ValueError(
                            'in file "%s", id "%s" for record %d could not be parsed',
                            input_fp, fastq_record.id, i)
                    else:
                        log.debug('original ID : %s', fastq_record.id)
                        fastq_record.id = '@{sample_name}__{run}'.format(**m.groupdict())
                        log.debug('new ID      : %s', fastq_record.id)
                        SeqIO.write(fastq_record, output_locus_file[m.group('locus')], 'fastq')

    return output_dir


def step_03_remove_primers(input_dir, cutadapt_script_name='cutadapt'):
    function_name = sys._getframe().f_code.co_name
    log = logging.getLogger(name=function_name)
    output_dir = create_output_dir(output_dir_name=function_name, input_dir=input_dir)

    forward_fastq_files = glob.glob(os.path.join(input_dir, '*_01.fastq.gz'))
    if len(forward_fastq_files) == 0:
        raise Exception('found no forward reads in directory {}'.format(input_dir))

    minimum_length = 100

    forward_primer = 'ATTAGAWACCCVNGTAGTCC'
    reverse_primer = 'TTACCGCGGCKGCTGGCAC'
    for forward_fastq_fp in forward_fastq_files:
        log.info('removing forward primers from "%s"', forward_fastq_fp)
        reverse_fastq_fp = re.sub(string=forward_fastq_fp, pattern='_01\.fastq\.gz$', repl='_02.fastq.gz')
        log.info('removing reverse primers from "%s"', reverse_fastq_fp)

        forward_fastq_basename = os.path.basename(forward_fastq_fp)
        trimmed_forward_fastq_fp = os.path.join(
            output_dir,
            re.sub(string=forward_fastq_basename, pattern='_01\.fastq\.gz$', repl='_trimmed_01.fastq.gz'))
        trimmed_reverse_fastq_fp = os.path.join(
            output_dir,
            re.sub(string=forward_fastq_basename, pattern='_01\.fastq\.gz$', repl='_trimmed_02.fastq.gz'))

        run_cmd([
            cutadapt_script_name,
            '-a', forward_primer,
            '-A', reverse_primer,
            '-o', trimmed_forward_fastq_fp,
            '-p', trimmed_reverse_fastq_fp,
            '-m', str(minimum_length),
            forward_fastq_fp,
            reverse_fastq_fp
        ])

    return output_dir


def step_04_merge_forward_reverse_reads_with_pear(input_dir):
    function_name = sys._getframe().f_code.co_name
    log = logging.getLogger(name=function_name)
    output_dir = create_output_dir(output_dir_name=function_name, input_dir=input_dir)
    return output_dir


def step_05_qc_reads_with_vsearch(input_dir):
    function_name = sys._getframe().f_code.co_name
    log = logging.getLogger(name=function_name)
    output_dir = create_output_dir(output_dir_name=function_name, input_dir=input_dir)
    return output_dir


def step_06_combine_runs(input_dir):
    function_name = sys._getframe().f_code.co_name
    log = logging.getLogger(name=function_name)
    output_dir = create_output_dir(output_dir_name=function_name, input_dir=input_dir)
    return output_dir


def step_07_dereplicate_sort_remove_low_abundance_reads(input_dir):
    function_name = sys._getframe().f_code.co_name
    log = logging.getLogger(name=function_name)
    output_dir = create_output_dir(output_dir_name=function_name, input_dir=input_dir)
    return output_dir


def step_08_cluster_97_percent(input_dir):
    function_name = sys._getframe().f_code.co_name
    log = logging.getLogger(name=function_name)
    output_dir = create_output_dir(output_dir_name=function_name, input_dir=input_dir)
    return output_dir


def step_09_reference_based_chimera_detection(input_dir):
    function_name = sys._getframe().f_code.co_name
    log = logging.getLogger(name=function_name)
    output_dir = create_output_dir(output_dir_name=function_name, input_dir=input_dir)
    return output_dir


def step_10_map_raw_reads_to_otus(input_dir):
    function_name = sys._getframe().f_code.co_name
    log = logging.getLogger(name=function_name)
    output_dir = create_output_dir(output_dir_name=function_name, input_dir=input_dir)
    return output_dir


def step_11_write_otu_table(input_dir):
    function_name = sys._getframe().f_code.co_name
    log = logging.getLogger(name=function_name)
    output_dir = create_output_dir(output_dir_name=function_name, input_dir=input_dir)
    return output_dir


def create_output_dir(output_dir_name, parent_dir=None, input_dir=None):
    if parent_dir is not None and input_dir is None:
        pass
    elif input_dir is not None and parent_dir is None:
        parent_dir, _ = os.path.split(input_dir)
    else:
        raise ValueError('exactly one of parent_dir and input_dir must be None')

    output_dir = os.path.join(parent_dir, output_dir_name)
    os.mkdir(output_dir)
    return output_dir


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()