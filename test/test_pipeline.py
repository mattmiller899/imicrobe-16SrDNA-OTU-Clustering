import gzip
import logging
import os.path
import tempfile

import cluster_16S.pipeline as pipeline
import cluster_16S.pipeline_util

logging.basicConfig(level=logging.DEBUG)


def get_pipeline(work_dir='/work_dir'):
    return pipeline.Pipeline(
        work_dir=work_dir, core_count=1,
        cutadapt_min_length=10,
        vsearch_filter_maxee=1, vsearch_filter_trunclen=200,
        forward_primer='ATTAGAWACCCVNGTAGTCC', reverse_primer='TTACCGCGGCKGCTGGCAC',
        pear_min_overlap=1, pear_max_assembly_length=270, pear_min_assembly_length=0,
        vsearch_derep_minuniquesize=3,
        uchime_ref_db_fp='')


def test_create_output_dir__input_dir(fs):
    input_dir = '/input_dir'
    fs.CreateDirectory(input_dir)
    output_dir = cluster_16S.pipeline_util.create_output_dir(input_dir=input_dir, output_dir_name='output_dir')
    assert output_dir == '/output_dir'
    assert os.path.exists(output_dir)


def test_create_output_dir__parent_dir(fs):
    parent_dir = '/parent_dir'
    fs.CreateDirectory(parent_dir)
    output_dir = cluster_16S.pipeline_util.create_output_dir(parent_dir=parent_dir, output_dir_name='output_dir')
    assert output_dir == '/parent_dir/output_dir'
    assert os.path.exists(output_dir)


def test_step_01__text_input():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:

        write_forward_reverse_read_files(input_dir=input_dir, suffix='.fastq', compress=False)
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

        write_forward_reverse_read_files(input_dir=input_dir, suffix='.fastq', compress=True)
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
        write_forward_reverse_read_files(input_dir=input_dir, suffix='.fastq')
        assert len(os.listdir(input_dir)) == 2

        output_dir = get_pipeline(work_dir=work_dir).step_02_remove_primers(input_dir=input_dir)

        assert output_dir == os.path.join(work_dir, 'step_02_remove_primers')
        assert os.path.exists(output_dir)

        output_file_list = sorted(os.listdir(output_dir))
        assert len(output_file_list) == 2
        assert output_file_list[0] == 'input_file_trimmed_01.fastq.gz'
        assert output_file_list[1] == 'input_file_trimmed_02.fastq.gz'


def test_step_03_pear():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:
        write_forward_reverse_read_files(input_dir=input_dir, suffix='.fastq')
        assert len(os.listdir(input_dir)) == 2

        output_dir = get_pipeline(work_dir=work_dir).step_03_merge_forward_reverse_reads_with_pear(
            input_dir=input_dir
        )

        assert output_dir == os.path.join(work_dir, 'step_03_merge_forward_reverse_reads_with_pear')
        assert os.path.exists(output_dir)

        output_file_list = sorted(os.listdir(output_dir))
        assert len(output_file_list) == 4
        assert output_file_list[0] == 'input_file_merged.assembled.fastq.gz'
        assert output_file_list[1] == 'input_file_merged.discarded.fastq.gz'
        assert output_file_list[2] == 'input_file_merged.unassembled.forward.fastq.gz'
        assert output_file_list[3] == 'input_file_merged.unassembled.reverse.fastq.gz'


def test_step_03_vsearch():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:
        write_forward_reverse_read_files(input_dir=input_dir, suffix='.fastq')
        assert len(os.listdir(input_dir)) == 2

        output_dir = get_pipeline(work_dir=work_dir).step_03_merge_forward_reverse_reads_with_vsearch(
            input_dir=input_dir
        )

        assert output_dir == os.path.join(work_dir, 'step_03_merge_forward_reverse_reads_with_vsearch')
        assert os.path.exists(output_dir)

        output_file_list = sorted(os.listdir(output_dir))
        assert len(output_file_list) == 3
        assert output_file_list[0] == 'input_file_merged.fastq.gz'
        assert output_file_list[1] == 'input_file_notmerged_fwd.fastq.gz'
        assert output_file_list[2] == 'input_file_notmerged_rev.fastq.gz'


def test_step_04():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:
        write_forward_reverse_read_files(input_dir=input_dir, suffix='.assembled.fastq')
        assert len(os.listdir(input_dir)) == 2

        output_dir = get_pipeline(work_dir=work_dir).step_04_qc_reads_with_vsearch(
            input_dir=input_dir
        )

        assert output_dir == os.path.join(work_dir, 'step_04_qc_reads_with_vsearch')
        assert os.path.exists(output_dir)

        output_file_list = sorted(os.listdir(output_dir))
        assert len(output_file_list) == 2
        assert output_file_list[0] == 'input_file_01.assembled.ee1trunc200.fastq.gz'
        assert output_file_list[1] == 'input_file_02.assembled.ee1trunc200.fastq.gz'


def test_step_05():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:
        write_forward_reverse_read_files(input_dir=input_dir, suffix='_trimmed_merged_V4.assembled.ee1trunc200.fastq')
        assert len(os.listdir(input_dir)) == 2

        output_dir = get_pipeline(work_dir=work_dir).step_05_combine_runs(input_dir=input_dir)

        assert output_dir == os.path.join(work_dir, 'step_05_combine_runs')
        assert os.path.exists(output_dir)

        output_file_list = os.listdir(output_dir)
        assert len(output_file_list) == 1
        assert output_file_list[0] == 'input_file_01_02_trimmed_merged_V4.assembled.ee1trunc200.fastq.gz'


def test_step_06():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:
        write_forward_reverse_read_files(input_dir=input_dir, suffix='_trimmed_merged_V4.assembled.ee1trunc200.fastq')
        assert len(os.listdir(input_dir)) == 2

        output_dir = get_pipeline(work_dir=work_dir).step_06_dereplicate_sort_remove_low_abundance_reads(input_dir=input_dir)

        assert output_dir == os.path.join(work_dir, 'step_06_dereplicate_sort_remove_low_abundance_reads')
        assert os.path.exists(output_dir)

        output_file_list = sorted(os.listdir(output_dir))
        assert len(output_file_list) == 4
        assert output_file_list[0] == 'input_file_01_trimmed_merged_V4.assembled.ee1trunc200.derepmin3.fasta.gz'
        assert output_file_list[1] == 'input_file_01_trimmed_merged_V4.assembled.ee1trunc200.derepmin3.txt'
        assert output_file_list[2] == 'input_file_02_trimmed_merged_V4.assembled.ee1trunc200.derepmin3.fasta.gz'
        assert output_file_list[3] == 'input_file_02_trimmed_merged_V4.assembled.ee1trunc200.derepmin3.txt'


def test_step_07():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:
        write_forward_reverse_read_files(input_dir=input_dir, suffix='.fasta')
        assert len(os.listdir(input_dir)) == 2

        output_dir = get_pipeline(work_dir=work_dir).step_07_cluster_97_percent(input_dir=input_dir)

        assert output_dir == os.path.join(work_dir, 'step_07_cluster_97_percent')
        assert os.path.exists(output_dir)

        output_file_list = sorted(os.listdir(output_dir))
        assert len(output_file_list) == 4
        assert output_file_list[0] == 'input_file_01.rad3.fasta'
        assert output_file_list[1] == 'input_file_01.rad3.txt'
        assert output_file_list[2] == 'input_file_02.rad3.fasta'
        assert output_file_list[3] == 'input_file_02.rad3.txt'


def test_pipeline():
    with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as work_dir:
        write_forward_reverse_read_files(input_dir=input_dir, suffix='.fastq')

        output_dir_list = get_pipeline(work_dir=work_dir).run(input_dir=input_dir)

        assert len(output_dir_list) == 10

        for work_dir in output_dir_list:
            assert os.path.exists(work_dir)
            assert os.path.isdir(work_dir)


forward_fastq_records = '''\
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:14076:1869 1:N:0
TACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGAGACTTGAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTTCGAAAGCGTGGGGATCAACCAGGATTCGCTACCCCTGTAGTCCGTAT
+
FFGCGD@FGGFG?FFGGGGGGGGGGG7FGGGGAFGGGGEDCFGGGGGGGGGGGGGGGGGGGF8<FDG9,<?F<9,,?EE@,B>FFCGGGGGGGGGCGGGGGFGGGGFGGGFGGGGGGGGGFFCFFG;CCC@@CFF7:DFF?CCGGGF;FFGFGFBCFDEE*?BFFGGGDG*BCFFG??CGGGGGC@CCCF9?E*8*:8CCCCC;CFG?CFFGGGE7FF7CDC3*/>F4F*>><DDCC3CCDFD:DFFF4**995@F#####################
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:10696:1889 1:N:0
TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGGGCGCAGACGGTTACTTAAGCAGGATGTGAAATCCCCGGGCTCAACCCGGGAACTGCGTTCTGAACTGGGTGACTCGAGTGTGTCAGAGGGAGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGACTACCGATGGCGAAGGCAGCCTCCTGGGACAACACTGACGTTCATGCCCGAAAGCGTGGGTAGCAAACAGGATTATATACCCTCGTAGTCCTTAT
+
F9FCGCGGDGCFC7@FGGG@CEFGC@7DFCEFEFGGGG@,,6EFECGGF:CFGGGGG7?FG,<<?9=++BAG99,<F?F<A+7:FCFFFGGGGGGFGGGGGGGGGGGGF<FGGGGGGGGGGC*@7DFGG,><9CC7>@FDF,3CCFF;C:B,@:;CC:FF@9+?CFGE*C*?+5;C??+=C**?7F77B8;=C58***:8EE8@EFF?DG6FGFFG=FF*9DDGFDGFCFGEG3<7*9D7C3<63?=F)9>C<>FB#####################
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:22552:1905 1:N:0
TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGCGCGCAGGCGGTCTTTTAAGTCTGATGTGAAAGCCCCCGGCTTAACCGGGGAGGGTCATTGGAAACTGGAAGACTGGAGTGCAGTAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGATATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTCACTGACGCTGAGGCTCGAAAGCGTGGGTAGCGAACAGGATTAGATACCCCTGTAGTCCGTACGAC
+
GGGGEGGGGGGGGGGGGGGGGFGFG@,BEGGGGGGGFGGGGGGGGGGGGGGGDFEFFBF9?DFFE9FG9FEFE9FA9DEDCEGGEGGGFFFFGGGEEEGFGGGGFCFFGGGCA,3@<CFEFGGGGGCG,>?;FGFDGGGGGGGGGGDFDFGGGGGGGGGGGF+BFFGGGG8FFCGGGGGGGGGGGGGCEGGGCE=G=ECG5=EEG?GGGFGGGGFFGGGCF:DGGGFGGGFGGG>DFGGDGGDGFFD3>@34DCGCFFFF496<C<36:AAFF#######
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:11319:1946 1:N:0
TACGTAGGTCCCGAGCGTTGTCCGGATTTACTGGGCGTAAAGGGAGCGTAGGTGGATATTTAAGTGGGATGTGAAATACTCGGGCTTAACCTGGGTGCTGCATTCCAAACTGGATATCTAGAGTGCAGGAGAGGAAAGTAGAATTCCTAGTGTAGCGGTGAAATGCGTAGAGATTAGGAAGAATACCAGTGGCGAAGGCGACTTTCTGGACTGTAACTGACACTGAGGCTCGAAAGCGTGGGGAGCAACCAGGATTAGATACCCGTGTATTCCTT
+
GCGGE,C@FFFGGCGGEFGEGFGGGEGGGGGGGGGGGGG@F<FGDGGGGGDGGGG7<FFGFGGCG99F?FEFAFGFFGFGGGGGGGGG9CFGGGFFEGGGGGGGFFGGGGGGGGGDFGF9?CFGGG<FGGGGGE8BFGGGGGGGGGGGGGGGGFGGGGG>EGGGGGGG58:D,DGFGGDFFCD>DE6CFGFC5EGCEGDGGEGCE?FGGFFGGGGC+;+:F4CG46CGD:<C5>*75)>BDDD>DB>F*9*95>:96*7:F?),7::?:4:AF4:
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:8349:1952 1:N:0
TACGGAGGGTGCAAGCGTTACTCGGAATCACTGGGCGTAAAGAGCGCGTAGGCGGGATAGTCAGTCAGGTGTGAAATCCTATGGCTTAACCATAGAACTGCATTTGAAACTACTATTCTAGAGTGTGGGAGAGGTAGGTGGAATTCTTGGTGTAGGGGTAAAATCCGTAGAGATCAAGAGGAATACTCATTGCGAAGGCGACCTGCTGGAACATTACTGACGCTGATTGCGCGAAAGCGTGGGGAGCAAACAGGATTAGCAACCCCAGTAGTCCGTAGA
+
FGCGEGGGGGGE8<FGGGGFFGGGG>7FFGFGGGGGGGGCFCFGGGGGGGFGCFEG7FCEGGGGGGGGFGGGGFFEDGGGGGDGGGGGGGGGGGGFFGAFDFGFGG=CFFF<FG9FFGGEGCFFGGGDGGGGEEE:CGCGC1,,DFGGGGGGGGCFGGCFGGGGGGGGFGGGFFGFFFGGFB+BFFGFFEC9C=8E55EE=CEECDGGGF6AFGGGFGCCFDGDGDGGCF:<5))177DDDGGF>)5>FF>FDG4<==?)9@>7C@F?>AF4A2?F###
'''


reverse_fastq_records = '''\
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:14076:1869 2:N:0
CCTGTTTGATCCCCACGCTTTCGCACATCAGCGTCAGTTACAGACCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCATATCTCTACGCATTTCACCGCTACACATGGAATTCCACTTTCCTCTTCTGCACTCAAGTCTCCCCGTTTCCAATGACCCTCCACGGTTGAGCCTTGGGCTTTCACATCAGCCTTAAAAAACCGCCTCCGCGCGCTTTACGCCCACTAATTCCGGTTAACGCTTCCCACCTACGTATTACCGCGTCTTCTTGCCCCCAG
+
GGGGGGGGFGGFDFCFEGGGGGGGGG7FFCEFEFGE@FE9CFC9EFFFC<EFFGDGGGGGGGE:EF?FFGGGG9AFF8F;CFGGGGGCFDGGGGGGGGGGCBCG,E,,9>EFGF,@FGGGGGGDFED8D+A?D8,2DFEGFG+=@DGF?+?E?DDDC=DD6A76@EGG+CF5A*0*;C?DD7+?10*3*0;A7A+1?>D5659@########################################################################
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:10696:1889 2:N:0
CCTGTTTGCTACCCACGCTTTCGGGCATGAACGTCAGTGTTGTCCCAGGAGGCTGCCTTCGCCATCGGTATTCCTCCACATCTCTACGCATTTCACTGCTACACGTGGAATTCTACCTCCCTCTGACACACTCGCGTCACCCCTTTCCGAACGCAGTTCCCGTGTTGAGCCCGTCGATTTCACATCCTGCTTAAGTAACCGTCTTCGCCCGCTTTCCCCCCAGTAATTCCGATTAACGCCTGCCCCTTACGTTTTCCCGCTGCTTCTGGCCCCCA
+
GGGFFGGGGGAAFC9CFEGGGGGCCC@F8E@FFGGC@E,C@,CFFFDFF7:C@BFGDFGGGDCECFEEGGGG?FFGE,9AEGFGGFFCF@CFFF9EDFDG9EEE@E,+84E?;=EFD9,AC,,@6,@+=,==D+++@EEEDF++6=+6+=C6C+07:22E,*,@FGGF###########################################################################################################
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:22552:1905 2:N:0
CCTGTTCGCTACCCACGCTTTCGAGCCTCAGCGTCAGTGACAGACCAGAGAGTCGCCTTCGCCACTGGTGTTCCCCCACATATCTACGCATTTCCCCGCTACACGTGGACTTCCACTCTCCTCTTATGCACTCCAGTCTTCCAGTTTCCAATGACCCTCCCCGGTTAAGCCGGGGGCTTTCACATCATACTTACAAGACCGCCTGCGCGCGCTTTACGCCCCATAAATCCGGACCACGCTTTCCACCTACGTATTACCGCGGCTGCTGTCACGCTGTA
+
EFF<FF<,:FCGFF97@:@FFF@F@EGD@,,@C@:@FFEF<<,6,C,,,C<FFF>:F:FF7?F=FDFFFGGGFF,,8+>+A,4CECFF@+@?BF,5B>C4B+@=+@8=E,,C=E9EFGGGG9FG=,6,,,2=@DDDFD@FGGADFAA=@++59AF?++:0==@*@9C,@CGDG*9):180+12;8A++/8+7++21=5AFD99:?#########################################################################
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:11319:1946 2:N:0
CCTGTTTGCTCCCCACGCTTTCGAGCCTCAGTGTCAGTTACAGTCCAGAAAGTCGCCTTCGCCACTGGTATTCTTCCTAATCTCTACGCATTTCACCGCTACACTAGGAATTCTACTTTCCTCTCCTGCACTCTAGATATCCAGTTTGGAATGCAGCACCCAGGTTAAGCCCGAGTATTTAACATCCCACTTAAATATCCACCTACGCTCCCTTTCCGCCCAGTAAATCCGAACAACGCCCGGGACCTACGTATTACCGCGCCGTGGTGCCCTCAGTC
+
GGGGGGGGGGGGCEE@+CFCFGGCGGDGFCFFFGCEFGGGG9DFCEFGCFFGFFGGGGGGDDFDCFGGGCFGFGF9AFF9A=;,CEGG@@FEEGGGCC@EG6BDFGFCF9FEECDF,=E;DD=8=DG88,6?DFFFGGGGGDFGGDD8+=+?8EGGGGGGDDC+@==,DC+@D+*<7<C5+++25*;*)003=68<?9;+?):>EE44)):?F#################################################################
@R3-16S-mockE-1__HWI-M01380:86:000000000-ALK4C:1:1101:8349:1952 2:N:0
CCTGTTTGCTCCCCACGCTTTCGCGCAATCAGCGTCAGTAATGTTCCAGCAGGTCGCCTTCGCAATGAGTATTCCTCTTGATCTCTACGGATTTTACCCCTACACCAAGAATTCCACCTACCTCTCCCACACTCTAGAATAGTAGTTTCAAATGCAGTTCTATGGTTAAGCCATAGGATTTCACACCTGACTGACTATCCCGCCTACGCGCTCTTTACGCCCAGTGATTCCGAGTACCGCTTTCACCCTCCGTATTACCGCGGCGGCTGTCACCCCTGT
+
GGGGGGGFGGAFGFGGGGGGGGGGDEFGGGCFFEGGGGG9AFGGGGGGGGGGGGGGGGGGGGGFCFGBAFGGGGGCFGGGGGGGGGCGGC+?FFGFGGGGGFFGFGCF@G@GGGG9DGGGGGGGFEE@,D+@?DFGGGGGGCFCFFGCF7=C9FGGGGGGGGGGGGFG7FG9CGGFFGGFG7<CFFEF6CFFB>AGF>8*;9:5>5:>A92/?>ACCE?D9>AA#######################################################
'''


def write_forward_reverse_read_files(input_dir, suffix='', compress=True):
    log = logging.getLogger(name=__name__)
    input_01_fp = os.path.join(input_dir, 'input_file_01{}'.format(suffix))
    input_02_fp = os.path.join(input_dir, 'input_file_02{}'.format(suffix))

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
