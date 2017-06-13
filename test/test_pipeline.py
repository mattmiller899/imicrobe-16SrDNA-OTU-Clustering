import os.path

import cluster_16S.pipeline as pipeline


def test_pipeline(fs):
    input_dir = '/input'
    fs.CreateDirectory(input_dir)

    output_dir = '/output'
    fs.CreateDirectory(output_dir)

    output_dir_list = pipeline.pipeline(
        input_dir=input_dir,
        output_dir=output_dir,
        core_count=1,
        forward_primer='ATTAGAWACCCVNGTAGTCC',
        reverse_primer='TTACCGCGGCKGCTGGCAC',
        uchime_ref_db_fp='')

    assert len(output_dir_list) == 11

    for output_dir in output_dir_list:
        assert os.path.exists(output_dir)
        assert os.path.isdir(output_dir)
