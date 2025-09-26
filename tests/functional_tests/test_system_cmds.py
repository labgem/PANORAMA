from distutils import command
import pytest
from tests.utils.run_command import run_command
from pathlib import Path
from tests.functional_tests.test_utils_cmd import utils_hmm_list, utils_model_list


@pytest.fixture(scope="session")
def annotation_and_systems_cmds(pangenome_list_file, utils_hmm_list, utils_model_list, num_cpus):


    annotation_command = ( f"panorama annotation "
              f"--pangenomes {pangenome_list_file} "
              f"--source defensefinder "
              f"--hmm {utils_hmm_list} " # "--mode sensitive " this produce an error 
              "--k_best_hit 3 "
              f"--threads {num_cpus}")
    
    run_command(annotation_command)


    systems_command = ("panorama systems "
                f"--pangenomes {pangenome_list_file} "
                f"--models {utils_model_list} "
                "--source defensefinder "
                "--annotation_sources defensefinder "
                f"--threads {num_cpus}"
                )

    run_command(systems_command)


@pytest.mark.requires_test_data
def test_write_systems(annotation_and_systems_cmds, pangenome_list_file, utils_model_list, num_cpus):

    outdir = pangenome_list_file.parent / "write_systems_outdir" # Create temp dir using factory

    command = (f"panorama write_systems  --pangenomes {pangenome_list_file} "
               f"--output {outdir} --models {utils_model_list} "
               "--sources defensefinder --projection "
               f"--threads {num_cpus} "
               "--association all"
               )

    run_command(command)