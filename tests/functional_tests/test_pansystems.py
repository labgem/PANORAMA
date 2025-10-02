import pytest
from tests.utils.run_command import run_command
from tests.functional_tests.test_utils_cmd import utils_hmm_list, utils_model_list


@pytest.mark.requires_test_data
def test_pansystems_command(pangenome_list_file, utils_hmm_list, utils_model_list, num_cpus):

    outdir = pangenome_list_file.parent / "pansystems_outdir"
    annotation_command = ( f"panorama pansystems "
              f"--pangenomes {pangenome_list_file} "
              f"--source defensefinder_wf "
              f"--hmm {utils_hmm_list} "
              f"--threads {num_cpus} "
              f"--models {utils_model_list} "
              "--association all "
              f"-o {outdir}" )
    
    run_command(annotation_command)
    
    # Validate output structure and files
    assert outdir.exists(), f"Output directory {outdir} was not created"
    
    # Get list of pangenome names from the pangenome list file
    pangenome_names = []
    with open(pangenome_list_file, 'r') as f:
        for line in f:
            if line.strip():
                pangenome_name = line.split('\t')[0]
                pangenome_names.append(pangenome_name)
    
    # Expected files in each source directory
    expected_files = [
        "association.tsv",
        "correlation_modules.html",
        "correlation_RGPs.html", 
        "correlation_spots.html",
        "module_to_systems.tsv",
        "rgp_to_systems.tsv",
        "spot_to_systems.tsv"
    ]
    
    # Check each pangenome directory
    for pangenome_name in pangenome_names:
        pangenome_dir = outdir / pangenome_name
        assert pangenome_dir.exists(), f"Pangenome directory {pangenome_dir} was not created"
        
        # Check source directory
        source_dir = pangenome_dir / "defensefinder_wf"
        assert source_dir.exists(), f"Source directory {source_dir} was not created"
        
        # Check all expected files exist
        for expected_file in expected_files:
            file_path = source_dir / expected_file
            assert file_path.exists(), f"Expected file {file_path} was not created"
            assert file_path.stat().st_size > 0, f"Expected file {file_path} is empty"

