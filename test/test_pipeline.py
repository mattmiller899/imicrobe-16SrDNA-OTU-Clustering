import gzip
import logging
import os.path
import tempfile

import cluster_16S.pipeline as pipeline


logging.basicConfig(level=logging.DEBUG)


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


def test_step_01(fs):
    input_dir = '/input_dir'
    fs.CreateDirectory(input_dir)
    fs.CreateFile(file_path=os.path.join(input_dir, 'input_file_01.fastq'))
    fs.CreateFile(file_path=os.path.join(input_dir, 'input_file_02.fastq'))
    work_dir = '/work_dir'
    fs.CreateDirectory(work_dir)

    pipeline.step_01_copy_and_compress(input_dir=input_dir, work_dir=work_dir)

    output_dir_list = os.listdir(work_dir)
    assert len(output_dir_list) == 1
    output_file_list = os.listdir(os.path.join(work_dir, output_dir_list[0]))
    assert len(output_file_list) == 2
    assert output_file_list[0] == 'input_file_01.fastq.gz'
    assert output_file_list[1] == 'input_file_02.fastq.gz'


def test_step_02(fs):
    input_dir = '/work_dir/step_01_dir'
    fs.CreateDirectory(input_dir)
    fs.CreateFile(
        file_path=os.path.join(input_dir, 'input_file_01.fastq.gz'),
        contents=gzip.compress(bytes(forward_fastq_records, 'utf-8')))

    fs.CreateFile(
        file_path=os.path.join(input_dir, 'input_file_02.fastq.gz'),
        contents=gzip.compress(bytes(reverse_fastq_records, 'utf-8')))

    output_dir = pipeline.step_02_adjust_headers(input_dir=input_dir)

    assert output_dir == '/work_dir/step_02_adjust_headers'
    assert os.path.exists(output_dir)

    output_files = os.listdir(output_dir)
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
    with tempfile.TemporaryDirectory() as work_dir:
        with tempfile.TemporaryDirectory(dir=work_dir) as input_dir:
            write_forward_reverse_read_files(input_dir=input_dir)
            print(os.listdir(input_dir))
            assert len(os.listdir(input_dir)) == 2

            output_dir = pipeline.step_03_remove_primers(input_dir=input_dir)

            assert output_dir == os.path.join(work_dir, 'step_03_remove_primers')
            assert os.path.exists(output_dir)

            output_file_list = sorted(os.listdir(output_dir))
            assert len(output_file_list) == 2
            assert output_file_list[0] == 'input_file_trimmed_01.fastq.gz'
            assert output_file_list[1] == 'input_file_trimmed_02.fastq.gz'



def test_pipeline():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as output_dir:
        write_forward_reverse_read_files(input_dir=input_dir)

        output_dir_list = pipeline.pipeline(
            input_dir=input_dir,
            output_dir=output_dir,
            core_count=1,
            forward_primer='ATTAGAWACCCVNGTAGTCC',
            reverse_primer='TTACCGCGGCKGCTGGCAC',
            uchime_ref_db_fp='',
            cutadapt_script_name='cutadapt')

        assert len(output_dir_list) == 10

        for output_dir in output_dir_list:
            assert os.path.exists(output_dir)
            assert os.path.isdir(output_dir)


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


def write_forward_reverse_read_files(input_dir):
    input_01_fp = os.path.join(input_dir, 'input_file_01.fastq.gz')
    input_02_fp = os.path.join(input_dir, 'input_file_02.fastq.gz')
    with gzip.open(input_01_fp, 'wt') as input_01_file:
        input_01_file.write(forward_fastq_records)

    with gzip.open(input_02_fp, 'wt') as input_02_file:
        input_02_file.write(reverse_fastq_records)
