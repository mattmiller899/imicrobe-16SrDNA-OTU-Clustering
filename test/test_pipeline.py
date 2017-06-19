import gzip
import logging
import os.path
import tempfile

import cluster_16S.pipeline as pipeline


logging.basicConfig(level=logging.DEBUG)


def get_pipeline(work_dir='/work_dir'):
    return pipeline.Pipeline(
        work_dir=work_dir, core_count=1,
        cutadapt_script_name='cutadapt', cutadapt_min_length=100,
        forward_primer='ATTAGAWACCCVNGTAGTCC', reverse_primer='TTACCGCGGCKGCTGGCAC',
        pear_min_overlap=1, pear_max_assembly_length=270, pear_min_assembly_length=0,
        uchime_ref_db_fp='')


def test_create_output_dir__input_dir(fs):
    input_dir = '/input_dir'
    fs.CreateDirectory(input_dir)
    output_dir = pipeline.create_output_dir(input_dir=input_dir, output_dir_name='output_dir')
    assert output_dir == '/output_dir'
    assert os.path.exists(output_dir)


def test_create_output_dir__parent_dir(fs):
    parent_dir = '/parent_dir'
    fs.CreateDirectory(parent_dir)
    output_dir = pipeline.create_output_dir(parent_dir=parent_dir, output_dir_name='output_dir')
    assert output_dir == '/parent_dir/output_dir'
    assert os.path.exists(output_dir)


def test_step_01__text_input():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:

        write_forward_reverse_read_files(input_dir=input_dir, compress=False)
        output_dir = get_pipeline(work_dir=work_dir).step_01_copy_and_compress(input_dir=input_dir)

        assert os.path.exists(output_dir)
        output_file_list = sorted(os.listdir(output_dir))
        assert len(output_file_list) == 2
        assert output_file_list[0] == 'input_file_01.fastq.gz'
        assert output_file_list[1] == 'input_file_02.fastq.gz'

        output_1_fp = os.path.join(output_dir, output_file_list[0])
        with gzip.open(output_1_fp, 'rt') as output_1:
            assert output_1.read() == forward_fastq_records

        output_2_fp = os.path.join(output_dir, output_file_list[1])
        with gzip.open(output_2_fp, 'rt') as output_2:
            assert output_2.read() == reverse_fastq_records


def test_step_01__compressed_input():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:

        write_forward_reverse_read_files(input_dir=input_dir, compress=True)
        output_dir = get_pipeline(work_dir=work_dir).step_01_copy_and_compress(input_dir=input_dir)

        assert os.path.exists(output_dir)
        output_file_list = sorted(os.listdir(output_dir))
        assert len(output_file_list) == 2
        assert output_file_list[0] == 'input_file_01.fastq.gz'
        assert output_file_list[1] == 'input_file_02.fastq.gz'

        output_1_fp = os.path.join(output_dir, output_file_list[0])
        with gzip.open(output_1_fp, 'rt') as output_1:
            assert output_1.read() == forward_fastq_records

        output_2_fp = os.path.join(output_dir, output_file_list[1])
        with gzip.open(output_2_fp, 'rt') as output_2:
            assert output_2.read() == reverse_fastq_records


def test_step_02():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:
        write_forward_reverse_read_files(input_dir=input_dir)
        output_dir = get_pipeline(work_dir=work_dir).step_02_adjust_headers(input_dir=input_dir)

        assert output_dir == os.path.join(work_dir, 'step_02_adjust_headers')
        assert os.path.exists(output_dir)

        output_files = sorted(os.listdir(output_dir))
        assert len(output_files) == 4
        assert output_files[0] == 'input_file_01_16S.fastq.gz'
        assert output_files[1] == 'input_file_01_ITS.fastq.gz'
        assert output_files[2] == 'input_file_02_16S.fastq.gz'
        assert output_files[3] == 'input_file_02_ITS.fastq.gz'

        with gzip.open(os.path.join(output_dir, output_files[0]), 'rt') as output_file:
            assert '-16S-' in output_file.readline()
        with gzip.open(os.path.join(output_dir, output_files[1]), 'rt') as output_file:
            assert '-ITS-' in output_file.readline()
        with gzip.open(os.path.join(output_dir, output_files[2]), 'rt') as output_file:
            assert '-16S-' in output_file.readline()
        with gzip.open(os.path.join(output_dir, output_files[3]), 'rt') as output_file:
            assert '-ITS-' in output_file.readline()


def test_step_03():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:
        write_forward_reverse_read_files(input_dir=input_dir)
        assert len(os.listdir(input_dir)) == 2

        output_dir = get_pipeline(work_dir=work_dir).step_03_remove_primers(input_dir=input_dir)

        assert output_dir == os.path.join(work_dir, 'step_03_remove_primers')
        assert os.path.exists(output_dir)

        output_file_list = sorted(os.listdir(output_dir))
        assert len(output_file_list) == 2
        assert output_file_list[0] == 'input_file_trimmed_01.fastq.gz'
        assert output_file_list[1] == 'input_file_trimmed_02.fastq.gz'


def test_step_04():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:
        write_forward_reverse_read_files(input_dir=input_dir)
        assert len(os.listdir(input_dir)) == 2

        output_dir = get_pipeline(work_dir=work_dir).step_04_merge_forward_reverse_reads_with_pear(
            input_dir=input_dir
        )

        assert output_dir == os.path.join(work_dir, 'step_04_merge_forward_reverse_reads_with_pear')
        assert os.path.exists(output_dir)

        output_file_list = sorted(os.listdir(output_dir))
        assert len(output_file_list) == 4
        assert output_file_list[0] == 'input_file_joined.assembled.fastq.gz'
        assert output_file_list[1] == 'input_file_joined.discarded.fastq.gz'
        assert output_file_list[2] == 'input_file_joined.unassembled.forward.fastq.gz'
        assert output_file_list[3] == 'input_file_joined.unassembled.reverse.fastq.gz'


def test_pipeline():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:
        write_forward_reverse_read_files(input_dir=input_dir)

        output_dir_list = get_pipeline(work_dir=work_dir).run(input_dir=input_dir)

        assert len(output_dir_list) == 10

        for work_dir in output_dir_list:
            assert os.path.exists(work_dir)
            assert os.path.isdir(work_dir)


forward_fastq_records = '''\
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:17273:1842 1:N:0
TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGGGCGAGGTAGAATTCCACGTTTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCTATGGCGAAG
+
GGGGFGG>EGGGGDGDEGE,,CF<7,8@F<FED88E@F77,9CCBFDCG7+:=FGGGDFFE,,CFA<++:,B,,,AAAB,B7+3@:D@9FDFGGGGGGGGGEGGGGCGFG
@R3-ITS-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:14076:1869 1:N:0
TACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAAC
+
FFGCGD@FGGFG?FFGGGGGGGGGGG7FGGGGAFGGGGEDCFGGGGGGGGGGGGGGGGGGGF8<FDG9,<?F<9,,?EE@,B>FFCGGGGGGGGGCGGGGGFGGGGFGGG
'''


reverse_fastq_records = '''\
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:17273:1842 2:N:0
CCTGTTTGCTACCCACGCTTTCGGGCATGAACGTCAGTGTTGTCCCAGGAGGCTGCCTTCGCCATCGGTATTCCTCCACATCTCTACGCATTTCACTGCTACACGTGGAA
+
E<E,CEFC<F9@@C,@66@BFDF@+6+C,,,C6C@,,C6FE,C@6C,:,,,,B@,,BEDB,,CCB=,,B,?@EFE,,,,,599C@F;F7+?=@E,AEED==DB>8@,++,
@R3-ITS-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:14076:1869 2:N:0
CCTGTTTGATCCCCACGCTTTCGCACATCAGCGTCAGTTACAGACCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCATATCTCTACGCATTTCACCGCTACACATGGAA
+
GGGGGGGGFGGFDFCFEGGGGGGGGG7FFCEFEFGE@FE9CFC9EFFFC<EFFGDGGGGGGGE:EF?FFGGGG9AFF8F;CFGGGGGCFDGGGGGGGGGGCBCG,E,,9>
'''


def write_forward_reverse_read_files(input_dir, compress=True):
    log = logging.getLogger(name=__name__)
    input_01_fp = os.path.join(input_dir, 'input_file_01.fastq')
    input_02_fp = os.path.join(input_dir, 'input_file_02.fastq')

    if compress:
        input_01_fp = input_01_fp + '.gz'
        with gzip.open(input_01_fp, 'wt') as input_01_file:
            input_01_file.write(forward_fastq_records)

        input_02_fp = input_02_fp + '.gz'
        with gzip.open(input_02_fp, 'wt') as input_02_file:
            input_02_file.write(reverse_fastq_records)
    else:
        with open(input_01_fp, 'wt') as input_01_file:
            input_01_file.write(forward_fastq_records)

        with open(input_02_fp, 'wt') as input_02_file:
            input_02_file.write(reverse_fastq_records)

    log.info('wrote forward read file "%s"', input_01_fp)
    log.info('wrote reverse read file "%s"', input_02_fp)
