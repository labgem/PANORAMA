# from pathlib import Path

# import pytest

# from tests.functional_tests.test_utils_cmd import utils_clustering, utils_model_list, utils_hmm_list  # noqa: F401
# from tests.utils.run_command import run_command
# from tests.functional_tests.test_system_cmds import annotation_and_systems_cmds_pangenome_list  # noqa: F401


# @pytest.mark.requires_test_data
# def test_compare_spots_command(
#     annotation_and_systems_cmds_pangenome_list,
#     utils_model_list,  # noqa: F811
#     utils_clustering,  # noqa: F811
#     num_cpus,
#     update_golden,
# ):
#     # TODO Manage the problem with reproducibility for better testing

#     outdir = annotation_and_systems_cmds_pangenome_list.parent / "compare_spots_outdir"
#     compare_spots_command = (
#         "panorama compare_spots "
#         f"--pangenomes {annotation_and_systems_cmds_pangenome_list} "
#         f"-o {outdir} "
#         "--graph_formats gexf graphml "
#         f"--cpus {num_cpus} "
#         f"--cluster {utils_clustering} "
#         "--systems "
#         "-s defensefinder "
#         f"-m {utils_model_list} "
#     )
#     run_command(compare_spots_command)

#     # Validate output structure and files
#     assert outdir.exists(), f"Output directory {outdir} was not created"

#     # Expected files in each source directory
#     assert_or_update_files = [
#         "all_conserved_spots.tsv",
#         # "conserved_spots/conserved_spot_3.tsv",
#         # "conserved_spots/conserved_spot_24.tsv",
#         # "conserved_spots/conserved_spot_33.tsv",
#         # "conserved_spots/conserved_spot_60.tsv",
#     ]
#     expected_files = map(
#         Path,
#         assert_or_update_files
#         + [
#             "conserved_spots.gexf",
#             "conserved_spots.graphml",
#             "systems_link_with_conserved_spots_louvain.gexf",
#             "systems_link_with_conserved_spots_louvain.graphml",
#         ],
#     )
#     # Check all expected files exist
#     for expected_file in expected_files:
#         file_path = outdir / expected_file
#         assert file_path.exists(), f"Expected file {file_path} was not created"
#         assert file_path.stat().st_size > 0, f"Expected file {file_path} is empty"
#         # if expected_file in assert_or_update_files:
#         #     assert_or_update_file(expected_file, file_path, update_golden)
