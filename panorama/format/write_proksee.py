#!/usr/bin/env python3
# coding:utf-8
import logging
# default libraries
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import json
from pathlib import Path
from tqdm import tqdm
from typing import Dict, List
from uuid import uuid4

# installed libraries
from ppanggolin.genome import Organism, Contig, Gene
from ppanggolin.region import Spot

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenome


def read_settings(settings_data: dict):
    if "format" not in settings_data:
        settings_data["format"] = "circular"
    if "geneticCode" not in settings_data:
        # TODO Manage genetic code
        settings_data["geneticCode"] = "11"


def write_legend_items(legend_data: dict, features: List[str]):
    legend_data["items"] = [{"name": "CDS", "swatchColor": "rgba(0,0,153,0.5)", "decoration": "arrow"},
                            {"name": "persistent", "swatchColor": "rgba(229,156,4,1)", "decoration": "arc"},
                            {"name": "shell", "swatchColor": "rgba(60,254,91,1)", "decoration": "arc"},
                            {"name": "cloud", "swatchColor": "rgba(0,255,255,1)", "decoration": "arc"}]
    if "rgp" in features or "all" in features:
        legend_data["items"].append({"name": "RGP", "swatchColor": "rgba(255,0,165, 1)", "decoration": "arc"}),
    if "spots" in features or "all" in features:
        legend_data["items"].append({"name": "Spot", "swatchColor": "rgba(255,0,0, 1)", "decoration": "arc"})
    if "modules" in features or "all" in features:
        legend_data["items"].append({"name": "Module", "swatchColor": "rgba(0,255,0,1)", "decoration": "arc"})
    if "systems" in features or "all" in features:
        legend_data["items"].append({"name": "Systems", "swatchColor": "rgba(0,0,255,1)", "decoration": "arc"})


def write_tracks(features: List[str]):
    tracks = [{"name": "CDS", "separateFeaturesBy": "None", "position": "inside", "thicknessRatio": 1,
               "dataType": "feature", "dataMethod": "type", "dataKeys": "CDS"},
              {"name": "Partition", "separateFeaturesBy": "None", "position": "inside", "thicknessRatio": 1,
               "dataType": "feature", "dataMethod": "tag", "dataKeys": "partition"}]
    if "rgp" in features or "all" in features:
        tracks.append({"name": "RGP", "separateFeaturesBy": "None", "position": "outside", "thicknessRatio": 1,
                       "dataType": "feature", "dataMethod": "type", "dataKeys": "RGP"}),
    if "spots" in features or "all" in features:
        tracks.append({"name": "Spots", "separateFeaturesBy": "None", "position": "outside", "thicknessRatio": 1,
                       "dataType": "feature", "dataMethod": "type", "dataKeys": "Spot"})
    if "modules" in features or "all" in features:
        tracks.append({"name": "Module", "separateFeaturesBy": "None", "position": "outside", "thicknessRatio": 1,
                       "dataType": "feature", "dataMethod": "type", "dataKeys": "Module"})
    if "systems" in features or "all" in features:
        tracks.append({"name": "System", "separateFeaturesBy": "None", "position": "outside", "thicknessRatio": 1,
                       "dataType": "feature", "dataMethod": "type", "dataKeys": "System"})
    return tracks


def read_data(template: Path, features: List[str]) -> dict:
    with open(template, "r") as template_file:
        proksee_data = json.load(template_file)
    now = datetime.now()
    if "created" in proksee_data["cgview"]:
        proksee_data["cgview"]["updated"] = now.strftime("%Y-%m-%d %H:%M:%S")
        last_version = proksee_data["cgview"]["version"].split('.')
        proksee_data["cgview"]["version"] = ".".join(last_version[:-1] + last_version[-1] + 1)
    else:
        proksee_data["cgview"]["created"] = now.strftime("%Y-%m-%d %H:%M:%S")
        proksee_data["cgview"]["version"] = "1.0"
    if "name" not in proksee_data["cgview"]:
        proksee_data["cgview"]["name"] = "PANORAMA annotations at genome levels"
        proksee_data["cgview"]["id"] = uuid4().hex
    read_settings(proksee_data["cgview"]["settings"])
    if "items" not in proksee_data["cgview"]["legend"]:
        write_legend_items(proksee_data["cgview"]["legend"], features)
    if "tracks" not in proksee_data["cgview"]:
        proksee_data["cgview"]["tracks"] = write_tracks(features)
    return proksee_data


def write_contig(organism: Organism):
    contigs_data_list = []
    for contig in tqdm(organism.contigs, unit="contig", disable=True):
        contig: Contig
        contigs_data_list.append({"name": contig.name,
                                  "length": contig.length,
                                  "orientation": "+",
                                  "seq": "".join([gene.dna for gene in contig.genes])
                                  })
    return contigs_data_list


def write_genes(organism: Organism):
    genes_data_list = []
    gf2gene = {}
    for gene in tqdm(organism.genes, total=organism.number_of_genes(), unit="genes", disable=True):
        gene: Gene
        gf: GeneFamily
        gf = gene.family
        if gf.name in gf2gene:
            gf2gene[gf.name].append(gene)
        else:
            gf2gene[gf.name] = [gene]
        # annotations = {source: "|".join(str(annot) for annot in gf.get_source(source)) for source in gf.sources}
        genes_data_list.append({"name": gene.name,
                                "type": gene.type,
                                "contig": gene.contig.name,
                                "start": gene.start,
                                "stop": gene.stop,
                                "strand": 1 if gene.strand == "+" else -1,
                                "product": gene.product,
                                "legend": gene.type,
                                # "tags": list(annotations.keys()),
                                # "meta": annotations
                                })
    return genes_data_list, gf2gene


def write_partition(organism: Organism):
    partition_data_list = []
    for gene in tqdm(organism.genes, total=organism.number_of_genes(), unit="genes", disable=True):
        partition_data_list.append({"name": gene.family.name,
                                    "type": gene.family.named_partition,
                                    "contig": gene.contig.name,
                                    "start": gene.start,
                                    "stop": gene.stop,
                                    "legend": gene.family.named_partition,
                                    "tags": ["partition"]})
    return partition_data_list


def write_rgp(pangenome: Pangenome, organism: Organism):
    rgp_data_list = []
    for rgp in tqdm(pangenome.regions, unit="RGP", disable=True):
        if rgp.organism == organism:
            rgp_data_list.append({"name": rgp.name,
                                  "type": "RGP",
                                  "contig": rgp.contig.name,
                                  "start": rgp.start,
                                  "stop": rgp.stop,
                                  "legend": "RGP",
                                  "tags": []})
    return rgp_data_list


def write_spots(pangenome: Pangenome, organism: Organism, gf2genes: Dict[str, List[Gene]]):
    spots_data_list = []
    for spot in tqdm(pangenome.spots,  unit="Spot", disable=True):
        spot: Spot
        spot_orgs = set()
        for gf in spot.families:
            spot_orgs |= gf.organisms
        if organism in spot_orgs:
            gf_intersection = organism.families & spot.families
            completion = round(len(gf_intersection) / len(spot.families), 2)
            for gf in gf_intersection:
                for gene in gf2genes[gf.name]:
                    spots_data_list.append({"name": str(spot.ID),
                                            "type": "Spot",
                                            "start": gene.start,
                                            "stop": gene.stop,
                                            "contig": gene.contig.name,
                                            "legend": "Spot",
                                            "tags": [],
                                            "meta": {
                                                "completion": completion
                                            }})
    return spots_data_list


def write_modules(pangenome: Pangenome, organism: Organism, gf2genes: Dict[str, List[Gene]]):
    modules_data_list = []
    for module in tqdm(pangenome.modules, unit="Module", disable=True):
        mod_orgs = set()
        for gf in module.families:
            mod_orgs |= gf.organisms
        if organism in mod_orgs:
            gf_intersection = organism.families & module.families
            completion = round(len(gf_intersection) / len(module.families), 2)
            for gf in gf_intersection:
                for gene in gf2genes[gf.name]:
                    modules_data_list.append({"name": str(module.ID),
                                              "type": "Module",
                                              "start": gene.start,
                                              "stop": gene.stop,
                                              "contig": gene.contig.name,
                                              "legend": "Module",
                                              "tags": [],
                                              "meta": {
                                                  "completion": completion
                                              }})
    return modules_data_list


def write_systems(pangenome: Pangenome, organism: Organism, gf2genes: Dict[str, List[Gene]]):
    systems_data_list = []
    for sys in tqdm(pangenome.systems, total=pangenome.number_of_systems(), unit="System", disable=True):
        sys_orgs = set()
        for gf in sys.gene_families:
            sys_orgs |= gf.organisms
        if organism in sys_orgs:
            families_names = [family.name for family in sys.model.families]
            gf_intersection = organism.families & sys.gene_families
            completion = round(len(gf_intersection) / len(sys.gene_families), 2)
            for gf in gf_intersection:
                for gene in gf2genes[gf.name]:
                    annotations = []
                    if sys.source in gene.family.sources:
                        annotations = [annot.name for annot in gene.family.get_source(sys.source) if
                                       annot.name in families_names]

                    systems_data_list.append({"name": str(sys.ID),
                                              "type": "System",
                                              "start": gene.start,
                                              "stop": gene.stop,
                                              "contig": gene.contig.name,
                                              "legend": "Systems",
                                              "tags": [sys.source, sys.name],
                                              "meta": {
                                                  "completion": completion,
                                                  "annotation": "|".join(annotations)
                                              }})
    return systems_data_list


def write_proksee_organism(pangenome: Pangenome, organism: Organism, output: Path, template: Path,
                           features: List[str] = None):
    proksee_data = read_data(template=template, features=features)
    if "name" not in proksee_data["cgview"]["captions"]:
        proksee_data["cgview"]["captions"][0]["name"] = f"{organism.name} annotated with PANORAMA"
    proksee_data["cgview"]["sequence"]["contigs"] = write_contig(organism)
    if "features" not in proksee_data["cgview"]:
        proksee_data["cgview"]["features"] = []
    genes_features, gf2genes = write_genes(organism)
    proksee_data["cgview"]["features"] += genes_features
    proksee_data["cgview"]["features"] += write_partition(organism)
    if "rgp" in features or "all" in features:
        proksee_data["cgview"]["features"] += write_rgp(pangenome=pangenome, organism=organism)
    if "spots" in features or "all" in features:
        proksee_data["cgview"]["features"] += write_spots(pangenome=pangenome, organism=organism, gf2genes=gf2genes)
    if "modules" in features or "all" in features:
        proksee_data["cgview"]["features"] += write_modules(pangenome=pangenome, organism=organism, gf2genes=gf2genes)
    if "systems" in features or "all" in features:
        proksee_data["cgview"]["features"] += write_systems(pangenome=pangenome, organism=organism, gf2genes=gf2genes)
    logging.getLogger().debug(f"Write proksee for {organism.name}")
    with open(output.joinpath(organism.name).with_suffix(".json"), "w") as out_json:
        json.dump(proksee_data, out_json, indent=2)


def write_proksee(pangenome: Pangenome, output: Path, features: List[str] = None, template: Path = None,
                  organisms_list: List[str] = None, threads: int = 1, disable_bar: bool = False):
    assert features is not None
    if template is None:
        template = Path(__file__).parent.joinpath("cgview_template").with_suffix(".json")
    if organisms_list is not None:
        organisms = [organism for organism in pangenome.organisms if organism.name in organisms_list]
    else:
        organisms = pangenome.organisms
    with ThreadPoolExecutor(max_workers=threads) as executor:
        with tqdm(total=pangenome.number_of_organisms(), unit='organisms', disable=disable_bar) as progress:
            futures = []
            for organism in organisms:
                future = executor.submit(write_proksee_organism, pangenome, organism, output, template, features)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()
