#!/usr/bin/env python3
# coding:utf-8

# default libraries
from datetime import datetime
import json
from pathlib import Path
from tqdm import tqdm
from uuid import uuid4

# installed libraries
from ppanggolin.genome import Organism, Contig, Gene

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenome


def read_settings(organism: Organism, settings_data: dict):
    if "format" not in settings_data:
        settings_data["format"] = "circular"
    if "geneticCode" not in settings_data:
        # TODO Manage genetic code
        settings_data["geneticCode"] = "11"


def read_captions(organim: Organism, captions_data: dict):
    if "name" not in captions_data:
        captions_data[0]["name"] = f"{organim.name} annotated with PANORAMA"


def read_data(organism: Organism, template: Path) -> dict:
    with open(template, "r") as template_file:
        cgview_data = json.load(template_file)

    now = datetime.now()
    if "created" in cgview_data["cgview"]:
        cgview_data["cgview"]["updated"] = now.strftime("%Y-%m-%d %H:%M:%S")
        last_version = cgview_data["cgview"]["version"].split('.')
        cgview_data["cgview"]["version"] = ".".join(last_version[:-1] + last_version[-1] + 1)
    else:
        cgview_data["cgview"]["created"] = now.strftime("%Y-%m-%d %H:%M:%S")
        cgview_data["cgview"]["version"] = "1.0"
    if "name" not in cgview_data["cgview"]:
        cgview_data["cgview"]["name"] = "PANORAMA annotations at genome levels"
        cgview_data["cgview"]["id"] = uuid4().hex
    read_settings(organism, cgview_data["cgview"]["settings"])
    read_captions(organism, cgview_data["cgview"]["captions"])
    return cgview_data


def write_contig(organism: Organism):
    contigs_data_list = []
    for contig in tqdm(organism.contigs):
        contig: Contig
        contigs_data_list.append({"name": contig.name,
                                  "length": contig.length,
                                  "orientation": "+",
                                  "seq": "".join([gene.dna for gene in contig.genes])
                                  })
    return contigs_data_list


def write_genes(organism: Organism):
    genes_data_list = []
    for gene in tqdm(organism.genes):
        gene: Gene
        genes_data_list.append({"name": gene.name,
                                "type": gene.type,
                                "contig": gene.contig.name,
                                "start": gene.start,
                                "stop": gene.stop,
                                "strand": 1 if gene.strand == "+" else -1,
                                "product": gene.product,
                                "legend": gene.type,
                                "tags": []})
    return genes_data_list


def write_partition(organism: Organism):
    partition_data_list = []
    for gene in tqdm(organism.genes):
        partition_data_list.append({"name": gene.family.name,
                                    "type": gene.family.named_partition,
                                    "contig": gene.contig.name,
                                    "start": gene.start,
                                    "stop": gene.stop,
                                    "strand": 1 if gene.strand == "+" else -1,
                                    "legend": gene.family.named_partition,
                                    "tags": ["partition"]})
    return partition_data_list


def write_rgp(pangenome: Pangenome, organism: Organism):
    rgp_data_list = []
    for rgp in tqdm(pangenome.regions):
        if rgp.organism == organism:
            rgp_data_list.append({"name": rgp.name,
                                  "type": "RGP",
                                  "contig": rgp.contig.name,
                                  "start": rgp.start,
                                  "stop": rgp.stop,
                                  "legend": "RGP",
                                  "tags": []})
    return rgp_data_list


def write_spots(pangenome: Pangenome, organism: Organism):
    spots_data_list = []
    for spot in tqdm(pangenome.spots):
        if all(rgp.organism == organism for rgp in spot.regions):
            spots_data_list.append({"name": str(spot.ID),
                                    "type": "Spot",
                                    "start": min([rgp.start for rgp in spot.regions]),
                                    "stop": max([rgp.stop for rgp in spot.regions]),
                                    "legend": "Spot",
                                    "tags": []})
    return spots_data_list


def write_modules(pangenome: Pangenome, organism: Organism):
    def find_start_stop(gf_intersection, module):
        for gf in gf_intersection:
            for gene in gf.genes:
                if gene.organism == organism:
                    genes_module = {gene}
                    contig = gene.contig
                    neighboors = contig.genes
                    pos_l, pos_r = (gene.position, gene.position)
                    neighboor_l, neighboor_r = (neighboors[pos_l], neighboors[pos_r])
                    while neighboor_l.family in gf_intersection and pos_l > 0:
                        genes_module.add(neighboor_l)
                        pos_l -= 1
                        neighboor_l = neighboors[pos_l]
                    while neighboor_r.family in gf_intersection and pos_r < len(neighboors) - 1:
                        genes_module.add(neighboor_r)
                        pos_r += 1
                        try:
                            neighboor_r = neighboors[pos_r]
                        except IndexError:
                            print("pika")
                    if len(genes_module) == len(module.families):
                        return neighboors[pos_l].start, neighboors[pos_r].stop, contig.name
        return 0, 0, ""

    modules_data_list = []
    cmpt = 0
    for module in tqdm(pangenome.modules):
        mod_orgs = set()
        for gf in module.families:
            mod_orgs |= gf.organisms
        if organism in mod_orgs:
            gf_intersection = organism.families & module.families
            completion = round(len(gf_intersection) / len(module.families), 2)
            start, stop, contig = find_start_stop(gf_intersection, module)
            if start != 0 and stop != 0:
                modules_data_list.append({"name": str(module.ID),
                                          "type": "Module",
                                          "start": start,
                                          "stop": stop,
                                          "contig": contig,
                                          "legend": "Module",
                                          "tags": [],
                                          "meta": {
                                              "completion": completion
                                          }})
            else:
                cmpt += 1
    print(cmpt)
    return modules_data_list


# {"name":"trbG","type":"mobileOG-db","start":1636,"stop":2670,"strand":1,"source":"mobile_og_db-1.1","legend":"transfer","contig":"NZ_AJPX01000002.1","tags":[],"meta":{"mobileOG_id":"mobileOG_000347679","best_hit_accession":"WP_128929822.1","major_category":"transfer","minor_category":"conjugation","source_db":"Plasmid RefSeq","evidence":"Homology"}}
def write_cgview(organism: Organism, template: Path):
    cgview_data = read_data(organism, template)
    cgview_data["cgview"]["sequence"]["contigs"] += write_contig(organism)
    if "features" not in cgview_data["cgview"]:
        cgview_data["cgview"]["features"] = []
    cgview_data["cgview"]["features"] += write_genes(organism)
    cgview_data["cgview"]["features"] += write_partition(organism)
    cgview_data["cgview"]["features"] += write_rgp(pangenome=pangenome, organism=organism)
    # cgview_data["cgview"]["features"] += write_spots(pangenome=pangenome, organism=organism)
    cgview_data["cgview"]["features"] += write_modules(pangenome=pangenome, organism=organism)
    with open("test_cgview.json", "w") as out_json:
        json.dump(cgview_data, out_json, indent=2)


if __name__ == "__main__":
    from panorama.format.read_binaries import check_pangenome_info

    template = Path("./", "cgview_template").with_suffix(".json")
    pangenome = Pangenome(name="test_cgview", taxid=0)
    pangenome.add_file("/home/jerome/Projets/PANORAMA/testingDataset/annot/pangenome/pangenome.h5")
    check_pangenome_info(pangenome, need_annotations=True, need_gene_sequences=True, need_families=True,
                         need_partitions=True, need_rgp=True, need_spots=False, need_modules=True)
    organism = pangenome.get_organism("GCF_000261545.1_15354_genomic")
    write_cgview(organism, template)
