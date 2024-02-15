#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
from typing import Union

# installed libraries
import ppanggolin.metadata

# local libraries
from panorama.annotate.hmm_search import profile_gfs
from panorama.format.read_binaries import check_pangenome_info, load_pangenomes
from panorama.format.write_proksee import write_proksee
from panorama.utils import check_tsv_sanity
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenomes
from panorama.format.conserved_spot import *
from panorama.alignment.align import all_against_all

# For Quentin
from panorama.format.figure import *
from multiprocessing import Manager
import tempfile

need_annotations = False
need_families = False
need_graph = False
need_partitions = False
need_spots = False
need_rgp = False
need_modules = False
need_gene_sequences = False
need_metadata = False
need_systems = False
bool_rgp = False
bool_modules = False
bool_spots = False


def check_flat_parameters(args):
    if args.hmm and (args.msa is None or args.msa_format is None):
        raise Exception("To write HMM you need to give msa files and format")
    if args.systems or args.systems_asso is not None or (args.proksee is not None and "systems" in args.proksee):
        if args.models is None or args.sources is None:
            raise Exception("To read system and write flat related, "
                            "it's necessary to give models and source.")
        if len(args.models) != len(args.sources):
            raise Exception("To read systems from different sources you need to give "
                            "one models directory by corresponding source.")


def write_annotations_to_families(pangenome: Pangenome, output: Path, sources: Set[str], disable_bar: bool = False):
    """ Write a tsv file with all annotations and sources present in pangenomes

    :param pangenome: Pangenome with all annotations
    :param output: Path to output directory
    :param disable_bar: allow to disable the progress bar
    """
    nb_source = len(sources)
    source_list = list(sources)
    column_name = np.array(
        f"Pangenome,families,{','.join([f'Annotation_{source},Accession_{source},Secondary_names_{source}' for source in source_list])}".split(
            ','))
    array_list = []
    for gf in tqdm(pangenome.gene_families, unit='gene families', disable=disable_bar):
        if any(source in source_list for source in gf.sources):
            annot_array = np.empty((gf.max_metadata_by_source()[1], 2 + nb_source * 3), dtype=object)
            if annot_array.shape[0] > 0:
                annot_array[:, 0] = pangenome.name
                annot_array[:, 1] = gf.name
                index_source = 2
                for source in source_list:
                    index_annot = 0
                    if source in gf.sources:
                        for annotation in gf.get_metadata_by_source(source):
                            annotation: ppanggolin.metadata.Metadata
                            annot_array[index_annot, index_source] = annotation.protein_name
                            annot_array[index_annot, index_source + 1] = annotation.Accession
                            if ('secondary_name' in annotation.__dict__.keys() and
                                    (annotation.secondary_name is not None or annotation.secondary_name != pd.NA)):
                                annot_array[index_annot, index_source + 2] = annotation.secondary_name
                            else:
                                annot_array[index_annot, index_source + 2] = '-'
                            index_annot += 1
                    index_source += 3
                array_list.append(annot_array)
    out_df = pd.DataFrame(np.concatenate(array_list), columns=column_name)
    out_df = out_df.sort_values(by=['Pangenome', 'families'] + list(column_name[range(3, len(column_name), 2)]))
    logging.getLogger("PANORAMA").info(",".join(out_df.columns[1:]))
    out_df.to_csv(f"{output}/{pangenome.name}/families_annotations.tsv", sep="\t",
                  columns=out_df.columns[1:], header=True, index=False)


def write_annotation_to_families_mp(pangenome_name: str, pangenome_info: dict, output: Path, sources: set = None,
                                    disable_bar: bool = False):
    pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
    pangenome.add_file(pangenome_info["path"])
    sources = sources if sources is not None else pangenome.status['metasources']["families"]
    check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                         need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_rgp,
                         need_spots=need_spots, need_gene_sequences=need_gene_sequences,
                         need_modules=need_modules, need_metadata=need_metadata, sources=sources,
                         metatypes={"families"}, disable_bar=disable_bar)
    write_annotations_to_families(pangenome, output, sources=sources, disable_bar=disable_bar)


def write_hmm(gf: GeneFamily, output: Path):
    """ Write an HMM profile for a gene family

    :param gf: Gene family with HMM
    :param output: Path to output directory
    """
    with open(f"{output.absolute().as_posix()}/{gf.name}.hmm", 'wb') as hmm_file:
        gf.HMM.write(hmm_file)


def write_hmm_profile(pangenome: Pangenome, output: Path, threads: int = 1,
                      disable_bar: bool = False):
    """Write an HMM profile for all gene families in pangenome

    :param pangenome: Pangenome with gene families
    :param output: Path to output directory
    :param threads: number of threads available
    :param disable_bar: allow to disable progress bar
    """
    with ThreadPoolExecutor(max_workers=threads) as executor:
        with tqdm(total=pangenome.number_of_gene_families(), unit='file',
                  desc='write gene families hmm/profile', disable=disable_bar) as progress:
            futures = []
            for gf in pangenome.gene_families:
                if gf.HMM is not None:
                    future = executor.submit(write_hmm, gf, output)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

                for future in futures:
                    future.result()


def write_flat_files(pan_to_path: Dict[str, Dict[str, Union[int, str]]], pangenomes: Pangenomes, output: Path,
                     annotation: bool = False, hmm: bool = False, systems: bool = False, systems_asso: List[str] = None,
                     conserved_spot: int = None, draw_spot: bool = False, proksee: List[str] = None,
                     proksee_template: Path = None,
                     organisms_list: List[str] = None, threads: int = 1, force: bool = False, disable_bar: bool = False,
                     **kwargs):
    """Launcher to write flat file from pangenomes

    :param pan_to_path: Pangenome name as key and path to hdf5 file as value
    :param pangenomes: Group of pangenomes
    :param output: Path to output directory
    :param annotation: Launch annotation write function
    :param hmm: Launch hmm write function
    :param systems: Launch systems write function
    :param systems_asso: Launch systems_asso write function
    :param conserved_spot: Launch conserved_spot write function
    :param draw_spot: Launch draw_spot write function
    :param proksee: Launch proksee write function
    :param proksee_template: Template used in proksee write function
    :param organisms_list: Organisms used in proksee write function
    :param threads: Number of available threads
    :param force: Boolean to allow force write in files
    :param disable_bar: disable progress bar
    """
    global need_annotations
    global need_families
    global need_graph
    global need_partitions
    global need_spots
    global need_rgp
    global need_modules
    global need_gene_sequences
    global need_metadata
    global need_systems
    global bool_rgp
    global bool_modules
    global bool_spots

    if annotation:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            logging.getLogger("PANORAMA").info('Write annotation')
            with tqdm(total=len(pan_to_path.keys()), unit='pangenome', disable=disable_bar) as progress:
                need_families = True
                need_metadata = True
                futures = []
                for pangenome_name, pangenome_info in pan_to_path.items():
                    future = executor.submit(write_annotation_to_families_mp, pangenome_name,
                                             pangenome_info, output, kwargs["sources"], disable_bar)
                    future.add_done_callback(lambda p: progress.update())
                    futures.append(future)
                for future in futures:
                    future.result()

            logging.getLogger("PANORAMA").info(
                f"Annotation has been written in {output}/{pangenome_name}/families_annotations.tsv")

    if hmm:
        need_families = True
        need_annotations = True
        for pangenome_name, pangenome_info in pan_to_path.items():
            pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
            pangenome.add_file(pangenome_info["path"])
            check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                                 need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_rgp,
                                 need_spots=need_spots, need_gene_sequences=need_gene_sequences,
                                 need_modules=need_modules,
                                 need_metadata=need_metadata, need_systems=need_systems,
                                 models=kwargs["models"], sources=kwargs["sources"], metatype="families",
                                 disable_bar=disable_bar)
            msa_path = kwargs.get('msa_path', None)
            msa_format = kwargs.get('msa_format', 'afa')
            threads = kwargs.get('threads', 1)
            profile_gfs(pangenome, msa_path, msa_format, threads, disable_bar)
            write_hmm_profile(pangenome, output, threads, disable_bar)

    if systems or systems_asso or conserved_spot or draw_spot:
        need_annotations = True
        need_families = True
        need_graph = True
        need_systems = True
        need_metadata = True
        need_partitions = True

        global_systems_proj = pd.DataFrame()
        global_systems_distribution = pd.DataFrame()
        global_id = pd.DataFrame()
        global_total = pd.DataFrame()

        df_spot_global = pd.DataFrame()
        df_module_global = pd.DataFrame()
        df_borders_global = pd.DataFrame()
        number_org_per_spot_global = pd.DataFrame()

        for pangenome_name, pangenome_info in pan_to_path.items():
            pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
            pangenome.add_file(pangenome_info["path"])
            if systems and not (systems_asso or conserved_spot or draw_spot):
                check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                                     need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_rgp,
                                     need_spots=need_spots, need_gene_sequences=need_gene_sequences,
                                     need_modules=need_modules,
                                     need_metadata=need_metadata, need_systems=need_systems, models=kwargs["models"],
                                     sources=kwargs["sources"], metatypes={"families"}, disable_bar=disable_bar)
            if systems_asso or conserved_spot or draw_spot:
                need_rgp = True
                need_modules = True
                need_spots = True
                check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                                     need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_rgp,
                                     need_spots=need_spots, need_gene_sequences=need_gene_sequences,
                                     need_modules=need_modules,
                                     need_metadata=need_metadata, need_systems=need_systems, models=kwargs["models"],
                                     sources=kwargs["sources"], metatypes={"families"}, disable_bar=disable_bar)
            logging.getLogger("PANORAMA").info(f"Begin writing systems projection for {pangenome_name}")
            for source in kwargs["sources"]:
                systems_proj = write_systems_projection(name=pangenome_name, pangenome=pangenome, output=output,
                                                        source=source,
                                                        threads=threads, force=force, disable_bar=disable_bar)
                logging.getLogger("PANORAMA").info(f"Finish writing systems projection for {pangenome_name}")
                # global_systems_proj = pd.concat([global_systems_proj, systems_proj])
                # per_pan_heatmap(name=pangenome_name, system_projection=systems_proj, output=output)
                # logging.getLogger("PANORAMA").info(f"Finish drawing heatmap figure for {pangenome_name}")
                #
                # systems_distribution = pan_distribution_system(name=pangenome_name, systems_projection=systems_proj)
                # global_systems_distribution = pd.concat([global_systems_distribution, systems_distribution])
                #
                # dataframe_id, dataframe_total = pan_number_system(name=pangenome_name, systems_projection=systems_proj)
                # global_id = pd.concat([global_id, dataframe_id])
                # global_total = pd.concat([global_total, dataframe_total])
                #
                # hbar_id_total(name=pangenome_name, dataframe_id=dataframe_id, dataframe_total=dataframe_total,
                #               output=output)
                # logging.getLogger("PANORAMA").info(f"Finish drawing histogram (ID and total count) for {pangenome_name}")
                #
                # if systems_asso or conserved_spot or draw_spot:
                #     if systems_asso in ["rgp", "rgp-modules", "rgp-spots", "modules-spots", "all"] or conserved_spot or draw_spot:
                #         bool_rgp = True
                #     if systems_asso in ["modules", "rgp-modules", "modules-spots", "all"] or conserved_spot or draw_spot:
                #         bool_modules = True
                #     if systems_asso in ["rgp-spots", "modules-spots", "all"] or conserved_spot or draw_spot:
                #         bool_spots = True
                #
                #     logging.getLogger("PANORAMA").info("Begin write systems with features projection")
                #     df_sys2feat = systems_to_features(name=pangenome_name, pangenome=pangenome, systems_projection=systems_proj,
                #                                       output=output, source=source, bool_rgp=bool_rgp, bool_modules=bool_modules,
                #                                       bool_spots=bool_spots, threads=threads, disable_bar=disable_bar)
                #
                #     df_borders, number_org_per_spot = write_borders_spot(name=pangenome_name, pangenome=pangenome)
                #     df_borders_global = pd.concat([df_borders_global, df_borders])
                #     number_org_per_spot_global = pd.concat([number_org_per_spot_global, number_org_per_spot])
                #
                #     if systems_asso in ["rgp", "rgp-spots", "modules-spots", "all"] or conserved_spot or draw_spot:
                #         df_spot, dict_spot_org = spot2sys(name=pangenome_name, pangenome=pangenome, system_to_feature=df_sys2feat,
                #                                           df_borders=df_borders, output=output)
                #         df_spot_global = pd.concat([df_spot_global, df_spot])
                #
                #     if systems_asso in ["modules", "rgp-modules", "modules-spots", "all"] or conserved_spot or draw_spot:
                #         df_module = mod2sys(name=pangenome_name, system_to_feature=df_sys2feat, output=output)
                #         df_module_global = pd.concat([df_module_global, df_module])
                #
                #     logging.getLogger("PANORAMA").info("Projection with features written")
                #     upsetplot(name=pangenome_name, systems_projection=systems_proj, system_to_feature=df_sys2feat, output=output)
                #     logging.getLogger("PANORAMA").info(f"Finish drawing upsetplot for {pangenome_name}")
                #
                #     if draw_spot:
                #         draw_spots(name=pangenome_name, pangenome=pangenome, output=output, df_spot=df_spot,
                #                    dict_spot_org=dict_spot_org, systems_projection=systems_proj)

        if conserved_spot:
            logging.getLogger("PANORAMA").info("Begin write conserved spot")
            manager = Manager()
            lock = manager.Lock()
            pans = load_pangenomes(pangenome_list=pangenomes, need_info={"need_families": True}, max_workers=1,
                                   lock=lock, disable_bar=False)
            tmpdir = tempfile.TemporaryDirectory(dir=output)
            df_align = all_against_all(pangenomes=pans, output=output, lock=lock, tmpdir=tmpdir)
            identical_spot(df_borders_global=df_borders_global, number_org_per_spot_global=number_org_per_spot_global,
                           df_spot_global=df_spot_global, df_align=df_align, output=output, threshold=conserved_spot)

        # heatmap(global_systems_proj=global_systems_proj, output=output)
        # logging.getLogger("PANORAMA").info(f"Global heatmap figure created")
        # global_systems_distribution.to_csv(f"{output}/systems_distribution.tsv", sep="\t", index=False)
        # global_systems_proj.to_csv(f"{output}/global_systems.tsv", sep="\t", index=False)
        # global_id.to_csv(f"{output}/global_systems_number_id.tsv", sep="\t", index=False)
        # global_total.to_csv(f"{output}/global_systems_number_total.tsv", sep="\t", index=False)
        #
        # if systems_asso in ["rgp", "rgp-spots", "modules-spots", "all"] or conserved_spot or draw_spot:
        #     df_spot_global.to_csv(f"{output}/global_spot_to_system.tsv", sep="\t", index=False)
        #
        # if systems_asso in ["modules", "rgp-modules", "modules-spots", "all"] or conserved_spot or draw_spot:
        #     df_module_global.to_csv(f"{output}/global_module_to_system.tsv", sep="\t", index=False)

    if proksee is not None:
        need_annotations = True
        need_gene_sequences = True
        need_families = True
        need_partitions = True
        need_metadata = True
        if "rgp" in proksee or "all" in proksee:
            need_rgp = True
        if "spots" in proksee or "all" in proksee:
            need_rgp = True
            need_spots = True
        if "modules" in proksee or "all" in proksee:
            need_modules = True
        if "annotations" in proksee or "all" in proksee:
            need_annotations = True
        if "systems" in proksee or "all" in proksee:
            need_systems = True
        for pangenome_name, pangenome_info in pan_to_path.items():
            pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
            pangenome.add_file(pangenome_info["path"])
            check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                                 need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_rgp,
                                 need_spots=need_spots, need_gene_sequences=need_gene_sequences,
                                 need_modules=need_modules,
                                 need_metadata=need_metadata, need_systems=need_systems,
                                 models=kwargs["models"], sources=kwargs["sources"], metatype="families",
                                 disable_bar=disable_bar)
            write_proksee(pangenome=pangenome, output=output, features=proksee, template=proksee_template,
                          organisms_list=organisms_list, threads=threads, disable_bar=disable_bar)


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    check_flat_parameters(args)
    pan_to_path = check_tsv_sanity(args.pangenomes)
    mkdir(args.output, force=args.force)
    write_flat_files(pan_to_path, args.pangenomes, output=args.output, annotation=args.annotations,
                     systems=args.systems,
                     systems_asso=args.systems_asso, conserved_spot=args.conserved_spot, draw_spot=args.draw_spot,
                     models=args.models, sources=args.sources, proksee=args.proksee,
                     proksee_template=args.proksee_template,
                     organisms_list=args.organisms, hmm=args.hmm, msa_path=args.msa, msa_format=args.msa_format,
                     threads=args.threads, force=args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("write", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_write(parser)
    return parser


def parser_write(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument("-o", "--output", required=True, type=Path, nargs='?',
                          help='Output directory')
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--annotations", required=False, action="store_true",
                          help="Write all the annotations from families")
    optional.add_argument("--systems", required=False, action="store_true",
                          help="Write all the systems in pangenomes and project on genomes")
    optional.add_argument("--systems_asso", required=False, type=str, default=None,
                          choices=["all", "rgp-modules", "rgp-spots", "modules-spots", "modules", "rgp"],
                          help="Write association between systems and others pangenomes elements")
    optional.add_argument("--conserved_spot", required=False, type=int, default=None,
                          help="Write all conserved spot in pangenomes having common families as borders (need a value between"
                               "1 and the total number of families defining a spot in PPanGGOLiN (6 by default)")
    optional.add_argument("--draw_spot", required=False, action="store_true",
                          help="Draw spots containing at least one system")
    optional.add_argument('--models', required=False, type=Path, default=None, nargs="+",
                          help="Path to model list file.")
    optional.add_argument("--sources", required=False, type=str, nargs="+", default=None,
                          help='Name of the annotation source where panorama as to select in pangenomes')
    optional.add_argument("--proksee", required=False, type=str, default=None, nargs='+',
                          choices=["all", "base", "modules", "rgp", "spots", "annotations", "systems"])
    optional.add_argument("--proksee_template", required=False, type=Path, default=None, nargs='?')
    optional.add_argument("--organisms", required=False, type=str, default=None, nargs='+')
    optional.add_argument("--hmm", required=False, action="store_true",
                          help="Write an hmm for each gene families in pangenomes")
    optional.add_argument("--msa", required=False, type=Path, default=None,
                          help="To create a HMM profile for families, you can give a msa of each gene in families."
                               "This msa could be get from ppanggolin (See ppanggolin msa). "
                               "If no msa provide Panorama will launch one.")
    optional.add_argument("--msa_format", required=False, type=str, default="afa",
                          choices=["stockholm", "pfam", "a2m", "psiblast", "selex", "afa",
                                   "clustal", "clustallike", "phylip", "phylips"],
                          help="Format of the input MSA.")
    optional.add_argument("--threads", required=False, type=int, default=1)
