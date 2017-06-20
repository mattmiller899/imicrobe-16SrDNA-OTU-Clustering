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
    logging.basicConfig(level=logging.INFO)
    args = get_args()

    Pipeline().run(input_dir=args.input_dir)
    return 0


def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i', '--input-dir', default='.',
                            help='path to the input directory')
    arg_parser.add_argument('-o', '--output-dir', default='.',
                            help='path to the output directory')

    arg_parser.add_argument('-c', '--core-count', default=1, type=int,
                            help='number of cores to use')

    arg_parser.add_argument('--forward-primer', default='ATTAGAWACCCVNGTAGTCC',
                            help='forward primer to be clipped')
    arg_parser.add_argument('--reverse-primer', default='TTACCGCGGCKGCTGGCAC',
                            help='reverse primer to be clipped')

    arg_parser.add_argument('--uchime-ref-db-fp', default='/cluster_16S/pr2/pr2_gb203_version_4.5.fasta',
                            help='database for vsearch --uchime_ref')

    arg_parser.add_argument('--cutadapt-script-name', default='cutadapt3',
                            help='cutadapt or cutadapt3')
    arg_parser.add_argument('--cutadapt-min-length', required=True, type=int,
                            help='min_length for cutadapt')

    arg_parser.add_argument('--pear-min-overlap', required=True, type=int,
                            help='-v/--min-overlap for pear')
    arg_parser.add_argument('--pear-max-assembly-length', required=True, type=int,
                            help='-m/--max-assembly-length for pear')
    arg_parser.add_argument('--pear-min-assembly-length', required=True, type=int,
                            help='-m/--min-assembly-length for pear')

    args = arg_parser.parse_args()
    return args


class Pipeline:
    def __init__(self,
                 work_dir,
                 core_count,
                 cutadapt_script_name, cutadapt_min_length,
                 forward_primer, reverse_primer,
                 pear_min_overlap, pear_max_assembly_length, pear_min_assembly_length,
                 uchime_ref_db_fp):

        self.work_dir = work_dir
        self.core_count = core_count

        self.cutadapt_script_name = cutadapt_script_name
        self.cutadapt_min_length = cutadapt_min_length
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer

        self.pear_min_overlap = pear_min_overlap
        self.pear_max_assembly_length = pear_max_assembly_length
        self.pear_min_assembly_length = pear_min_assembly_length

        self.uchime_ref_db_fp = uchime_ref_db_fp

    def run(self, input_dir):
        output_dir_list = []
        output_dir_list.append(self.step_01_copy_and_compress(input_dir=input_dir))
        ##output_dir_list.append(self.step_02_adjust_headers(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_03_remove_primers(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_04_merge_forward_reverse_reads_with_pear(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_05_qc_reads_with_vsearch(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_06_combine_runs(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_07_dereplicate_sort_remove_low_abundance_reads(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_08_cluster_97_percent(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_09_reference_based_chimera_detection(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_10_map_raw_reads_to_otus(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_11_write_otu_table(input_dir=output_dir_list[-1]))

        return output_dir_list


    def step_01_copy_and_compress(self, input_dir):
        function_name = sys._getframe().f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)
        log.debug('output_dir: %s', output_dir)

        input_file_glob = os.path.join(input_dir, '*.fastq*')
        log.debug('input_file_glob: %s', input_file_glob)
        input_fp_list = sorted(glob.glob(input_file_glob))
        log.info('input files: %s', input_fp_list)

        if len(input_fp_list) == 0:
            raise Exception('found no fastq files in directory "{}"'.format(input_dir))

        for input_fp in input_fp_list:
            destination_fp = os.path.join(output_dir, os.path.basename(input_fp))
            if input_fp.endswith('.gz'):
                with open(input_fp, 'rb') as f, open(destination_fp, 'wb') as g:
                    shutil.copyfileobj(fsrc=f, fdst=g)
            else:
                destination_fp = destination_fp + '.gz'
                with open(input_fp, 'rt') as f, gzip.open(destination_fp, 'wt') as g:
                    shutil.copyfileobj(fsrc=f, fdst=g)

        return output_dir


    def step_02_adjust_headers(self, input_dir):
        function_name = sys._getframe().f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)

        input_file_glob = os.path.join(input_dir, '*.fastq.gz')
        log.debug('input file glob: "%s"', input_file_glob)
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
            output_16S_fp = os.path.join(output_dir, re.sub(
                string=input_basename, pattern='\.fastq\.gz$', repl='_16S.fastq.gz'))
            output_ITS_fp = os.path.join(output_dir, re.sub(
                string=input_basename, pattern='\.fastq\.gz$', repl='_ITS.fastq.gz'))
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


    def step_03_remove_primers(self, input_dir):
        function_name = sys._getframe().f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)

        for forward_fastq_fp in get_forward_fastq_files(input_dir=input_dir):
            log.info('removing forward primers from file "%s"', forward_fastq_fp)
            forward_fastq_basename = os.path.basename(forward_fastq_fp)

            reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=forward_fastq_fp)
            log.info('removing reverse primers from file "%s"', reverse_fastq_fp)

            trimmed_forward_fastq_fp = os.path.join(
                output_dir,
                re.sub(
                    string=forward_fastq_basename,
                    pattern='_([0R])1',
                    repl=lambda m: '_trimmed_{}1'.format(m.group(1))))
            trimmed_reverse_fastq_fp = os.path.join(
                output_dir,
                re.sub(
                    string=forward_fastq_basename,
                    pattern='_([0R])1',
                    repl=lambda m: '_trimmed_{}2'.format(m.group(1))))

            run_cmd([
                self.cutadapt_script_name,
                '-a', self.forward_primer,
                '-A', self.reverse_primer,
                '-o', trimmed_forward_fastq_fp,
                '-p', trimmed_reverse_fastq_fp,
                '-m', str(self.cutadapt_min_length),
                forward_fastq_fp,
                reverse_fastq_fp
            ])

        return output_dir


    def step_04_merge_forward_reverse_reads_with_pear(self, input_dir):
        function_name = sys._getframe().f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)

        for forward_fastq_fp in get_forward_fastq_files(input_dir=input_dir):
            reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=forward_fastq_fp)

            # decompress fastq files because PEAR can't cope
            forward_fastq_uncompressed_fp = os.path.join(
                output_dir,
                os.path.basename(forward_fastq_fp)[:-3]
            )
            with gzip.open(forward_fastq_fp, 'rt') as src, open(forward_fastq_uncompressed_fp, 'wt') as dst:
                shutil.copyfileobj(fsrc=src, fdst=dst)

            reverse_fastq_uncompressed_fp = os.path.join(
                output_dir,
                os.path.basename(reverse_fastq_fp)[:-3]
            )
            with gzip.open(reverse_fastq_fp, 'rt') as src, open(reverse_fastq_uncompressed_fp, 'wt') as dst:
                shutil.copyfileobj(fsrc=src, fdst=dst)

            joined_fastq_basename = re.sub(
                string=os.path.basename(forward_fastq_uncompressed_fp),
                pattern=r'_([0R]1)',
                repl=lambda m: '_joined'.format(m.group(1))
            )[:-6]
            joined_fastq_fp_prefix = os.path.join(output_dir, joined_fastq_basename)
            log.info('writing joined paired-end reads to "%s"', joined_fastq_fp_prefix)
            run_cmd([
                'pear',
                '-f', forward_fastq_uncompressed_fp,
                '-r', reverse_fastq_uncompressed_fp,
                '-o', joined_fastq_fp_prefix,
                '--min-overlap', str(self.pear_min_overlap),
                '--max-assembly-length', str(self.pear_max_assembly_length),
                '--min-assembly-length', str(self.pear_min_assembly_length),
                '-j', str(self.core_count)
            ])

            os.remove(forward_fastq_uncompressed_fp)
            os.remove(reverse_fastq_uncompressed_fp)

            pear_output_file_glob = joined_fastq_fp_prefix + '.*.fastq'
            pear_output_files = glob.glob(pear_output_file_glob)
            for pear_output_fp in pear_output_files:
                with open(pear_output_fp, 'rt') as src, gzip.open(pear_output_fp + '.gz', 'wt') as dst:
                    log.info('compressing PEAR output "%s"', pear_output_fp)
                    shutil.copyfileobj(fsrc=src, fdst=dst)
                os.remove(pear_output_fp)

        return output_dir


    def step_05_qc_reads_with_vsearch(self, input_dir):
        function_name = sys._getframe().f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)
        return output_dir


    def step_06_combine_runs(self, input_dir):
        function_name = sys._getframe().f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)
        return output_dir


    def step_07_dereplicate_sort_remove_low_abundance_reads(self, input_dir):
        function_name = sys._getframe().f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)
        return output_dir


    def step_08_cluster_97_percent(self, input_dir):
        function_name = sys._getframe().f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)
        return output_dir


    def step_09_reference_based_chimera_detection(self, input_dir):
        function_name = sys._getframe().f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)
        return output_dir


    def step_10_map_raw_reads_to_otus(self, input_dir):
        function_name = sys._getframe().f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)
        return output_dir


    def step_11_write_otu_table(self, input_dir):
        function_name = sys._getframe().f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)
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


def get_forward_fastq_files(input_dir):
    log = logging.getLogger(name=__name__)
    input_glob = os.path.join(input_dir, '*_?1*.fastq.gz')
    log.info('searcing for forward read files with glob "%s"', input_glob)
    forward_fastq_files = glob.glob(input_glob)
    if len(forward_fastq_files) == 0:
        raise Exception('found no forward reads from glob "{}"'.format(input_glob))
    return forward_fastq_files


def get_associated_reverse_fastq_fp(forward_fp):
    forward_input_dir, forward_basename = os.path.split(forward_fp)
    reverse_fastq_basename = re.sub(
        string=forward_basename,
        pattern=r'_([0R])1',
        repl=lambda m: '_{}2'.format(m.group(1)))
    reverse_fastq_fp = os.path.join(forward_input_dir, reverse_fastq_basename)
    return reverse_fastq_fp


def run_cmd(cmd_line_list, **kwargs):
    log = logging.getLogger(name=__name__)
    try:
        log.info('executing "%s"', ' '.join((str(x) for x in cmd_line_list)))
        output = subprocess.run(
            cmd_line_list,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            **kwargs)
        log.info(output)
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


if __name__ == '__main__':
    main()