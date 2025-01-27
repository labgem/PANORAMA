
# default libraries
from __future__ import annotations

import logging
from typing import Dict, List, Optional
from pathlib import Path
# installed libraries
from tqdm import tqdm
import pandas as pd
import json

# installed libraries
from ppanggolin.genome import Organism
from ppanggolin.region import Module

# local libraries
from ppanggolin.formats.writeFlatGenomes import (
    manage_module_colors,
palette
)
from ppanggolin.formats.write_proksee import (
write_contig,
write_genes,
write_rgp,
write_modules,
initiate_proksee_data
)

# local libraries
from panorama.utils import mkdir

def write_pangenome_systems_in_proksee(pangenome, output: Path, pangenome_projection: pd.DataFrame, organisms_projection: pd.DataFrame,
                             organisms: Optional[List[str]] = None, force: bool = False, disable_bar: bool = False):
    """
    Write the projected systems to output files.

    Args:
        output (Path): Path to the output directory.
        pangenome_projection (pd.DataFrame): DataFrame containing the pangenome projection.
        organisms_projection (pd.DataFrame): DataFrame containing the organism projections.
        organisms (List[str], optional): List of organisms to project (default is all organisms).
        force (bool, optional): Force write to the output directory (default is False).

    Returns:
        None
    """

    
    module_to_colors = manage_module_colors(set(pangenome.modules))


    proj_dir = mkdir(output / "proksee", force=force)

    if organisms is not None:
        pangenome_projection = pangenome_projection[~pangenome_projection["organism"].isin(organisms)]
        organisms_projection = organisms_projection[~organisms_projection["organism"].isin(organisms)]

    organims_to_process = pangenome_projection["organism"].unique()
    for organism_name in tqdm(organims_to_process, total=len(organims_to_process), disable=disable_bar):
        
        organism = pangenome.get_organism(organism_name)

        org_df = organisms_projection.loc[organisms_projection["organism"] == organism_name]
        
        proksee_data = write_proksee_output(organism, org_df, module_to_colors) 


        with open(proj_dir / f"{organism_name}.json", 'w') as out_json:
            json.dump(proksee_data, out_json , indent=2, sort_keys=True)

        
def write_proksee_output(organism: Organism, org_df:pd.DataFrame, module_to_colors: Dict[Module, str],  features:List = ['all'], metadata_sep: str = ">@<" ):
        """
        """

        proksee_data = initiate_proksee_data(features, organism, {})
        
        systems_ids = set(org_df['system number'])

        colors = palette(len(systems_ids))
        for i, systems_id in enumerate(systems_ids):
            proksee_data['cgview']['legend']['items'].append({"name": f'system {systems_id}',
                                                        'swatchColor': colors[i],
                                                        'decoration': 'arc'
                                                        })
        
        tracks = proksee_data['cgview']['tracks']

        tracks.append(
                    {
                        "name": "System",
                        "separateFeaturesBy": "None",
                        "position": "inside",
                        "thicknessRatio": 1,
                        "dataType": "feature",
                        "dataMethod": "source",
                        "dataKeys": "System",
                    }
                )
        
        proksee_data["cgview"]["sequence"]["contigs"] = write_contig(
        organism, genome_sequences={}, metadata_sep=metadata_sep
        )


        genes_features, gf2genes = write_genes(
            organism, multigenics=[], metadata_sep=metadata_sep
        )

        proksee_data["cgview"]["features"] = genes_features


                
        if ("rgp" in features or "all" in features) and organism.regions is not None:
            proksee_data["cgview"]["features"] += write_rgp(
                organism=organism, metadata_sep=metadata_sep
            )

        if module_to_colors is not None and ("modules" in features or "all" in features):
            proksee_data["cgview"]["features"] += write_modules(
                organism=organism, gf2genes=gf2genes, metadata_sep=metadata_sep
            )
    
        system_data_list = []

        for _, row in org_df.iterrows():

            metadata_col = ['annotation',
                            'secondary_names',
                            'genomic organization',
                            'product']
            metadata_for_proksee = {c:row[c] for c in metadata_col if c in row}

            system_data_list.append(
                {
                    "name": row["system name"],
                    "contig": row["contig"],
                    "start": row["start"],
                    "stop": row["stop"],
                    "legend": f"system {row['system number']}",
                    "source": "System",
                    "tags": [],
                    "meta": metadata_for_proksee,
                }
            )

        proksee_data["cgview"]["features"] += system_data_list

        return proksee_data