import os
import tempfile

import pytest

import cluster_16S.pipeline as pipeline
import cluster_16S.pipeline_util


def test_create_output_dir():
    with tempfile.TemporaryDirectory() as work_dir:
        output_dir = cluster_16S.pipeline_util.create_output_dir(output_dir_name='test', parent_dir=work_dir)

        assert output_dir == os.path.join(work_dir, 'test')
        assert os.path.exists(output_dir)


def test_get_forward_fastq_files():
    with tempfile.TemporaryDirectory() as test_dir:
        # specify prefix to make sorted order reproducible
        tempfile.mkstemp(dir=test_dir, prefix='_R1', suffix='_blahblah.fastq.gz')
        tempfile.mkstemp(dir=test_dir, prefix='_01', suffix='_hahahaha.fastq')
        tempfile.mkstemp(dir=test_dir, prefix='_02', suffix='_notthisone.fastq.gz')
        tempfile.mkstemp(dir=test_dir, prefix='_R2', suffix='_northisone.fastq')

        forward_fastq_fps = sorted(cluster_16S.pipeline_util.get_forward_fastq_files(test_dir))
        assert len(forward_fastq_fps) == 2
        assert forward_fastq_fps[1].endswith('_blahblah.fastq.gz')
        assert forward_fastq_fps[0].endswith('_hahahaha.fastq')


def test_get_forward_fastq_files__exception():
    with tempfile.TemporaryDirectory() as test_dir:
        # specify prefix to make sorted order reproducible
        tempfile.mkstemp(dir=test_dir, prefix='_R_', suffix='_blahblah.fastq.gz')
        tempfile.mkstemp(dir=test_dir, prefix='_0_', suffix='_hahahaha.fastq')
        tempfile.mkstemp(dir=test_dir, prefix='_2_', suffix='_notthisone.fastq.gz')
        tempfile.mkstemp(dir=test_dir, prefix='_R_', suffix='_northisone.fastq')

        with pytest.raises(cluster_16S.pipeline_util.PipelineException):
            cluster_16S.pipeline_util.get_forward_fastq_files(test_dir)


def test_get_combined_file_name__0():
    with pytest.raises(cluster_16S.pipeline_util.PipelineException):
        pipeline.get_combined_file_name([])


def test_get_combined_file_name__1():
    combined_file_name = pipeline.get_combined_file_name([
        '/input/data/Mock_Run1_V4.assembled.fastq.gz',
    ])

    assert combined_file_name == 'Mock_Run1_V4.assembled.fastq.gz'


def test_get_combined_file_name__2():
    combined_file_name = pipeline.get_combined_file_name([
        '/input/data/Mock_Run1_V4.assembled.fastq.gz',
        '/input/data/Mock_Run3_V4.assembled.fastq.gz'
    ])

    assert combined_file_name == 'Mock_Run1_Run3_V4.assembled.fastq.gz'
