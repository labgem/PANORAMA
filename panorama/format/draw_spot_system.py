#!/usr/bin/env python3
# coding:utf-8

# default libraries
from collections import defaultdict
import random

import pandas as pd
from math import pi
import sys
from pathlib import Path
from typing import Dict

# installed libraries
from pandas import DataFrame
from scipy.spatial.distance import pdist
from scipy.sparse import csc_matrix
from scipy.cluster.hierarchy import linkage, dendrogram

from tqdm import tqdm
from bokeh.plotting import ColumnDataSource, figure, save
from bokeh.io import output_file
from bokeh.layouts import column, row
from bokeh.models import WheelZoomTool, LabelSet, Slider, CustomJS, HoverTool, RadioGroup, Div, Column, GlyphRenderer

# local libraries
from ppanggolin.utils import jaccard_similarities
from ppanggolin.RGP.spot import comp_border
from panorama.pangenomes import Pangenome
from panorama.utils import mkdir

def draw_spots(name: str, pangenome: Pangenome, output: Path, df_spot: DataFrame, dict_spot_org: Dict):
    """ Draw spot containing systems

    :param name: Name of the pangenome
    :param pangenome: Pangenome containing spot
    :param output: Path to output directory
    :param df_spot: Systems in each spot of the pangenome
    :param dict_spot_org: Dictionary with spot as key and organisms having the spot as values
    """

    spot_list = df_spot['Spot'].tolist()
    inter_spot = []
    for s in pangenome.spots:
        if s.ID in spot_list:
            inter_spot.append(s)

    multigenics = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])

    fam2mod = {}
    for mod in pangenome.modules:
        for fam in mod.families:
            fam2mod[fam] = f"module_{mod.ID}"

    fam2sys = {}
    system = pangenome.systems
    for sys in system:
        for fam in sys.gene_families:
            existing_sys = fam2sys.get(fam)
            if existing_sys:
                fam2sys[fam] = f"{existing_sys} / {sys.name}"
            else:
                fam2sys[fam] = f"{sys.name}"

    output_2 = output / name / "spot_figure"
    mkdir(output_2, force=True)

    for spot in tqdm(inter_spot, total=len(inter_spot), unit="spot"):
        df_rep_identical = pd.DataFrame(columns=['representative_rgp','representative_rgp_organism','identical_rgp','identical_rgp_organism'])
        fname = output / name / "spot_figure" / f"spot_{spot.ID}"
        # write rgps representatives and the rgps they are identical to
        out_struc = open(output / name / "spot_figure" / f"spot_{spot.ID}_identical_rgps.tsv", 'w')
        out_struc.write('representative_rgp\trepresentative_rgp_organism\tidentical_rgp\tidentical_rgp_organism\n')
        for keyRGP, otherRGPs in spot.get_uniq_to_rgp().items():
            for rgp in otherRGPs:
                new_row = pd.DataFrame([{'representative_rgp':keyRGP.name, 'representative_rgp_organism':keyRGP.organism.name,
                                         'identical_rgp':rgp.name, 'identical_rgp_organism':rgp.organism.name}])
                df_rep_identical = pd.concat([df_rep_identical, new_row], ignore_index=True)
                out_struc.write(f"{keyRGP.name}\t{keyRGP.organism.name}\t{rgp.name}\t{rgp.organism.name}\n")
        out_struc.close()

        fams = set()
        gene_lists = []
        for rgp in spot.regions:
            borders = rgp.get_bordering_genes(pangenome.parameters["spots"]["set_size"], multigenics)
            minpos = min([gene.position for border in borders for gene in border])
            maxpos = max([gene.position for border in borders for gene in border])
            gene_list = rgp.contig.genes[minpos:maxpos + 1]
            minstart = min([gene.start for border in borders for gene in border])
            maxstop = max([gene.stop for border in borders for gene in border])
            rnas_toadd = set()
            for rna in rgp.contig.RNAs:
                if minstart < rna.start < maxstop:
                    rnas_toadd.add(rna)
            gene_list.extend(rnas_toadd)
            gene_list = sorted(gene_list, key=lambda x: x.start)

            fams |= {gene.family for gene in gene_list if gene.type == "CDS"}

            gene_lists.append([gene_list, borders, rgp])
        famcolors = make_colors_for_iterable(fams)
        # order all rgps the same way, and order them by similarity in gene content
        gene_lists = order_gene_lists(gene_lists, pangenome.parameters["spots"]["overlapping_match"],
                                      pangenome.parameters["spots"]["exact_match"], pangenome.parameters["spots"]["set_size"])

        count_uniq = spot.count_uniq_ordered_set()

        # keep only the representative rgps for the figure
        uniq_gene_lists = []
        ordered_counts = []
        for genelist in gene_lists:
            curr_genelist_count = count_uniq.get(genelist[2], None)
            if curr_genelist_count is not None:
                # Remove organisms with no systems
                for gene in genelist:
                    if not isinstance(gene, list):
                        if any(str(gene.organism) in value for value in dict_spot_org.get(spot.ID, [])):
                            uniq_gene_lists.append(genelist)
                            ordered_counts.append(curr_genelist_count)

        draw_curr_spot(uniq_gene_lists, ordered_counts, fam2mod, fam2sys, famcolors, fname, df_rep_identical)


def draw_curr_spot(gene_lists: list, ordered_counts: list, fam_to_mod: dict, fam_to_sys: dict, fam_col: dict, file_name: str, df_rep_identical: DataFrame):
    """
    :param gene_lists:
    :param ordered_counts:
    :param fam_to_mod: Dictionnary which link families and modules
    :param fam_to_mod: Dictionnary which link families and systems
    :param fam_col: Dictionnary with for each family the corresponding color
    :param file_name: Path and name of the file
    :return:
    """

    # Prepare the source data
    output_file(f"{file_name}.html")

    # generate the figure and add some tools to it
    wheel_zoom = WheelZoomTool()
    fig = figure(title="spot graphic", plot_width=1600, plot_height=600,
                 tools=["pan", "box_zoom", "reset", "save", wheel_zoom, "ywheel_zoom", "xwheel_zoom"])
    fig.axis.visible = True
    fig.toolbar.active_scroll = wheel_zoom

    # genome rectangles
    genome_source, genome_tooltip = mk_genomes(gene_lists, ordered_counts, df_rep_identical)
    genome_recs = fig.rect(x='x', y='y', fill_color="dimgray", width="width", height=0.5, source=genome_source)
    genome_recs_hover = HoverTool(renderers=[genome_recs], tooltips=genome_tooltip, mode="mouse",
                                  point_policy="follow_mouse")
    fig.add_tools(genome_recs_hover)

    # gene rectanges
    gene_source, gene_tooltips = mk_source_data(gene_lists, fam_col, fam_to_mod, fam_to_sys)
    recs = fig.rect(x='x', y='y', line_color='line_color', fill_color='fill_color', width='width', height=2,
                    line_width=5, source=gene_source)
    recs_hover = HoverTool(renderers=[recs], tooltips=gene_tooltips, mode="mouse", point_policy="follow_mouse")
    fig.add_tools(recs_hover)
    # gene modification tools
    gene_tools = add_gene_tools(recs, gene_source)

    # label modification tools
    labels_tools, labels = add_gene_labels(fig, gene_source)

    # genome tool
    genome_tools = add_genome_tools(fig, recs, genome_recs, gene_source, genome_source, len(gene_lists), labels)

    save(column(fig, row(labels_tools, gene_tools), row(genome_tools)))


def mk_source_data(genelists: list, fam_col: dict, fam_to_mod: dict, fam_to_sys: dict) -> (ColumnDataSource, list):
    """

    :param genelists:
    :param fam_col: Dictionary with for each family the corresponding color
    :param fam_to_mod: Dictionary with the correspondance modules families
    :param fam_to_sys: Dictionary with the correspondance systems families
    :return:
    """
    partition_colors = {"shell": "#00D860", "persistent": "#F7A507", "cloud": "#79DEFF"}

    df = {'name': [], 'ordered': [], 'strand': [], "start": [], "stop": [], "length": [], 'module': [],
          'module_color': [], 'system': [], 'system_color': [], 'x': [], 'y': [], 'width': [], 'family_color': [],
          'partition_color': [], 'partition': [], "family": [], "product": [], "x_label": [], "y_label": [],
          "label": [], "gene_type": [], 'gene_ID': [], "gene_local_ID": []}

    for index, GeneList in enumerate(genelists):
        genelist = GeneList[0]

        if genelist[0].start < genelist[1].start:
            # if the order has been inverted, positionning elements on the figure is different
            ordered = True
            start = genelist[0].start
        else:
            ordered = False
            start = genelist[0].stop

        for gene in genelist:
            df["ordered"].append(str(ordered))
            df["strand"].append(gene.strand)
            df["start"].append(gene.start)
            df["stop"].append(gene.stop)
            df["length"].append(max([gene.stop, gene.start])-min([gene.stop, gene.start]))
            df["gene_type"].append(gene.type)
            df["product"].append(gene.product)
            df["gene_local_ID"].append(gene.local_identifier)
            df['gene_ID'].append(gene.ID)

            if "RNA" in gene.type:  # dedicated values for RNA genes
                df["name"].append(gene.product)
                df["family"].append(gene.type)
                df["partition"].append("none")
                df["family_color"].append("#A109A7")
                df["partition_color"].append("#A109A7")
                df["module"].append("none")
                df["system"].append("")
            else:
                df["name"].append(gene.name)
                df["family"].append(gene.family.name)
                df["partition"].append(gene.family.named_partition)
                df["family_color"].append(fam_col[gene.family])
                df["partition_color"].append(partition_colors[gene.family.named_partition])
                df["module"].append(fam_to_mod.get(gene.family, "none"))
                df["system"].append(fam_to_sys.get(gene.family, ""))

            df["x"].append((abs(gene.start - start) + abs(gene.stop - start)) / 2)
            df["width"].append(gene.stop - gene.start)
            df["x_label"].append(str(int(df["x"][-1]) - int(int(df["width"][-1]) / 2)))
            if ordered:
                if gene.strand == "+":
                    df["y"].append((index * 10) + 1)
                else:
                    df["y"].append((index * 10) - 1)
            else:
                if gene.strand == "+":
                    df["y"].append((index * 10) - 1)
                else:
                    df["y"].append((index * 10) + 1)
            df["y_label"].append(float(df["y"][-1]) + 1.5)
    df["label"] = df["name"]
    df["line_color"] = df["partition_color"]
    df["fill_color"] = df["family_color"]

    # define colors for modules
    mod2col = make_colors_for_iterable(set(df["module"]))
    mod_colors = []
    for mod in df["module"]:
        mod_colors.append(mod2col[mod])
    df["module_color"] = mod_colors

    # define colors for systems
    sys2col = make_colors_for_iterable(set(df["system"]))
    sys_colors = []
    for sys in df["system"]:
        sys_colors.append(sys2col[sys])
    df["system_color"] = sys_colors

    # defining things that we will see when hovering over the graphical elements
    tooltips = [
        ("start", "@start"),
        ("stop", "@stop"),
        ("length", "@length"),
        ("name", "@name"),
        ("product", "@product"),
        ("family", "@family"),
        ("module", "@module"),
        ("system", "@system"),
        ("partition", "@partition"),
        ("local identifier", "@gene_local_ID"),
        ("gene ID", "@gene_ID"),
        ("ordered", "@ordered"),
        ("strand", "@strand"),
    ]

    return ColumnDataSource(data=df), tooltips


def add_gene_tools(recs: GlyphRenderer, source_data: ColumnDataSource) -> Column:
    """
    Define tools to change the outline and fill colors of genes

    :param recs:
    :param source_data:
    :return:
    """

    def color_str(color_element: str) -> str:
        """ Javascript code to switch between partition, family and module color for the given 'color_element'

        :param color_element:

        :return: Javascript code
        """
        return f"""
            if(this.active == 0){{
                source.data['{color_element}'] = source.data['partition_color'];
            }}else if(this.active == 1){{
                source.data['{color_element}'] = source.data['family_color'];
            }}else if(this.active == 2){{
                source.data['{color_element}'] = source.data['module_color'];
            }}else if(this.active == 3){{
                source.data['{color_element}'] = source.data['system_color'];
            }}
            recs.{color_element} = source.data['{color_element}'];
            source.change.emit();
        """

    radio_line_color = RadioGroup(labels=["partition", "family", "module", "system"], active=0)
    radio_fill_color = RadioGroup(labels=["partition", "family", "module", "system"], active=1)

    radio_line_color.js_on_click(CustomJS(args=dict(recs=recs, source=source_data),
                                          code=color_str("line_color")))

    radio_fill_color.js_on_click(CustomJS(args=dict(recs=recs, source=source_data),
                                          code=color_str("fill_color")))

    color_header = Div(text="<b>Genes:</b>")
    line_title = Div(text="""Color to use for gene outlines:""",
                     width=200, height=100)
    fill_title = Div(text="""Color to fill genes with:""",
                     width=200, height=100)

    gene_outline_size = Slider(start=0, end=10, value=5, step=0.1, title="Gene outline size:")
    gene_outline_size.js_on_change('value', CustomJS(args=dict(other=recs),
                                                     code="""
                other.glyph.line_width = this.value;
                """
                                                     ))

    return column(color_header, row(column(line_title, radio_line_color), column(fill_title, radio_fill_color)),
                  gene_outline_size)

def add_gene_labels(fig, source_data: ColumnDataSource) -> (Column, LabelSet):
    """

    :param fig:
    :param source_data:
    :return:
    """
    labels = LabelSet(x='x_label', y='y_label', text='label', source=source_data, render_mode='canvas',
                      text_font_size="18px")
    slider_font = Slider(start=0, end=64, value=16, step=1, title="Gene label font size in px")
    slider_angle = Slider(start=0, end=pi / 2, value=0, step=0.01, title="Gene label angle in radian")

    radio_label_type = RadioGroup(labels=["system", "name", "product", "family", "local identifier", "gene ID", "none"],
                                  active=1)

    slider_angle.js_link('value', labels, 'angle')

    slider_font.js_on_change('value',
                             CustomJS(args=dict(other=labels),
                                      code="other.text_font_size = this.value+'px';"
                                      )
                             )

    radio_label_type.js_on_click(CustomJS(args=dict(other=labels, source=source_data),
                                          code="""
                if(this.active == 6){
                    source.data['label'] = [];
                    for(var i=0;i<source.data['system'].length;i++){
                        source.data['label'].push('');
                    }
                }else if(this.active == 0){
                    source.data['label'] = source.data['system'];
                }else if(this.active == 1){
                    source.data['label'] = source.data['name'];
                }else if(this.active == 2){
                    source.data['label'] = source.data['product'];
                }else if(this.active == 3){
                    source.data['label'] = source.data['family'];
                }else if(this.active == 4){
                    source.data['label'] = source.data['gene_local_ID'];
                }else if(this.active == 5){
                    source.data['label'] = source.data['gene_ID'];
                }
                else{
                    source.data['label'] = source.data[this.labels[this.labels]];
                }
                other.source = source;
                source.change.emit();
                """
                                          ))

    label_header = Div(text="<b>Gene labels:</b>")
    radio_title = Div(text="""Gene labels to use:""",
                      width=200, height=100)
    labels_block = column(label_header, row(slider_font, slider_angle), column(radio_title, radio_label_type))

    fig.add_layout(labels)

    return labels_block, labels


def add_genome_tools(fig, gene_recs: GlyphRenderer, genome_recs: GlyphRenderer, gene_source: ColumnDataSource,
                     genome_source: ColumnDataSource, nb: int, gene_labels: LabelSet):
    """

    :param fig:
    :param gene_recs:
    :param genome_recs:
    :param gene_source:
    :param genome_source:
    :param nb:
    :param gene_labels:
    :return:
    """
    # add genome labels
    genome_labels = LabelSet(x='x_label', y='y', x_offset=-20, text='name', text_align="right", source=genome_source,
                             render_mode='canvas', text_font_size="16px")
    fig.add_layout(genome_labels)

    slider_font = Slider(start=0, end=64, value=16, step=1, title="Genome label font size in px")
    slider_font.js_on_change('value',
                             CustomJS(args=dict(other=genome_labels),
                                      code="other.text_font_size = this.value+'px';"
                                      )
                             )

    slider_offset = Slider(start=-400, end=0, value=-20, step=1, title="Genome label offset")
    slider_offset.js_link('value', genome_labels, 'x_offset')

    slider_spacing = Slider(start=1, end=40, value=10, step=1, title="Genomes spacing")
    slider_spacing.js_on_change('value', CustomJS(
        args=dict(gene_recs=gene_recs, gene_source=gene_source, genome_recs=genome_recs, genome_source=genome_source,
                  nb_elements=nb, genome_labels=genome_labels, gene_labels=gene_labels),
        code="""
            var current_val = genome_source.data['y'][genome_source.data['y'].length - 1] / (nb_elements-1);
            for (let i=0 ; i < genome_source.data['y'].length ; i++){
                genome_source.data['y'][i] =  (genome_source.data['y'][i] * this.value) / current_val;
            }
            for (let i=0 ; i < gene_source.data['y'].length ; i++){
                if((gene_source.data['ordered'][i] == 'True' && gene_source.data['strand'][i] == '+') || (gene_source.data['ordered'][i] == 'False' && gene_source.data['strand'][i] == '-') ){
                    gene_source.data['y'][i] = (((gene_source.data['y'][i]-1) * this.value) / current_val) +1;
                    gene_source.data['y_label'][i] = (((gene_source.data['y_label'][i]-1-1.5) * this.value) / current_val) + 1 + 1.5;
                }else{
                    gene_source.data['y'][i] = (((gene_source.data['y'][i]+1) * this.value) / current_val) -1;
                    gene_source.data['y_label'][i] = (((gene_source.data['y_label'][i]+1-1.5) * this.value) / current_val) -1 + 1.5;

                }
            }
            gene_recs.source = gene_source;
            genome_recs.source = genome_source;
            gene_labels.source = gene_source;
            genome_labels.source = genome_source;
            gene_source.change.emit();
            genome_source.change.emit();
        """))

    genome_header = Div(text="<b>Genomes:</b>")
    return column(genome_header, slider_spacing, slider_font, slider_offset)


def mk_genomes(gene_lists: list, ordered_counts: list, df_rep_identical: DataFrame) -> (ColumnDataSource, list):
    """

    :param gene_lists:
    :param ordered_counts:
    :return:
    """
    df = {"name": [], "width": [], "occurrences": [], 'x': [], 'y': [], "x_label": []}

    for index, GeneList in enumerate(gene_lists):
        genelist = GeneList[0]
        df["occurrences"].append(ordered_counts[index])
        df["y"].append(index * 10)
        if genelist[0].start < genelist[1].start:
            # if the order has been inverted, positionning elements on the figure is different
            df["width"].append(abs(genelist[-1].stop - genelist[0].start))
        else:
            # order has been inverted
            df["width"].append(abs(genelist[0].stop - genelist[-1].start))
        df["x"].append((df["width"][-1]) / 2)
        df["x_label"].append(0)

        # Replace representative organism name by number of organisms
        df_rep_identical_f = df_rep_identical[(df_rep_identical['representative_rgp_organism'] == genelist[0].organism.name)]
        number_identical_org = df_rep_identical_f.shape[0]
        if number_identical_org == 1:
            df["name"].append(f"{number_identical_org} organism")
        else:
            df["name"].append(f"{number_identical_org} organisms")

    tooltip = [
        ("name", "@name"),
        ("occurrences", "@occurrences"),
    ]

    return ColumnDataSource(data=df), tooltip


def make_colors_for_iterable(it: set) -> dict:
    """
    Randomly picks a color for all elements of a given iterable

    :param it: Iterable families or modules

    :return: Dictionary with for each element a random color associate
    """

    famcol = {}
    for element in it:
        col = list(random.choices(range(256), k=3))
        if element == "none" or element == "":
            famcol[element] = "#D3D3D3"
        else:
            famcol[element] = '#%02x%02x%02x' % (col[0], col[1], col[2])
    return famcol


def order_gene_lists(gene_lists: list, overlapping_match: int, exact_match: int, set_size: int):
    """
    Order all rgps the same way, and order them by similarity in gene content.

    :param gene_lists: List of genes in rgps
    :param overlapping_match: Allowed number of missing persistent genes when comparing flanking genes
    :param exact_match: Number of perfectly matching flanking single copy markers required to associate RGPs
    :param set_size: Number of single copy markers to use as flanking genes for RGP

    :return: List of ordered genes
    """
    line_order_gene_lists(gene_lists, overlapping_match, exact_match, set_size)
    return row_order_gene_lists(gene_lists)


def row_order_gene_lists(gene_lists: list) -> list:
    """
    Row ordering of all rgps

    :param gene_lists:

    :return : An ordered genes list
    """
    fam_dict = defaultdict(set)
    #if there is only one, ordering is useless
    if len(gene_lists) == 1:
        return gene_lists

    if len(gene_lists) > sys.getrecursionlimit():
        sys.setrecursionlimit(len(gene_lists))#we need the recursion limit to be higher than the number of regions.

    for index, genelist in enumerate([genelist[0] for genelist in gene_lists]):
        for gene in genelist:
            if hasattr(gene, "family"):
                fam_dict[gene.family].add(index)
    all_indexes = []
    all_columns = []
    data = []
    for famIndex, RGPindexes in enumerate(fam_dict.values()):
        all_indexes.extend([famIndex] * len(RGPindexes))
        all_columns.extend(RGPindexes)
        data.extend([1.0] * len(RGPindexes))

    mat_p_a = csc_matrix((data, (all_indexes, all_columns)), shape=(len(fam_dict), len(gene_lists)), dtype='float')
    dist = pdist(1 - jaccard_similarities(mat_p_a, 0).todense())
    hc = linkage(dist, 'single')

    dendro = dendrogram(hc, no_plot=True)

    new_gene_lists = [gene_lists[index] for index in dendro["leaves"]]

    return new_gene_lists


def line_order_gene_lists(gene_lists: list, overlapping_match: int, exact_match: int, set_size: int):
    """
    Line ordering of all rgps

    :param gene_lists: list
    :param overlapping_match: Allowed number of missing persistent genes when comparing flanking genes
    :param exact_match: Number of perfectly matching flanking single copy markers required to associate RGPs
    :param set_size: Number of single copy markers to use as flanking genes for RGP
    """
    classified = {0}  # first gene list has the right order
    new_classify = set()

    to_classify = set(range(1, len(gene_lists)))  # the others may (or may not) have it

    while len(to_classify) != 0:
        for classIndex in classified:
            base_border1 = [gene.family for gene in gene_lists[classIndex][1][0]]
            base_border2 = [gene.family for gene in gene_lists[classIndex][1][1]]
            for unclassIndex in list(to_classify):
                border1 = [gene.family for gene in gene_lists[unclassIndex][1][0]]
                border2 = [gene.family for gene in gene_lists[unclassIndex][1][1]]
                if comp_border(base_border1, border1, overlapping_match, set_size, exact_match) and \
                        comp_border(base_border2, border2, overlapping_match, set_size, exact_match):
                    to_classify.discard(unclassIndex)
                    new_classify.add(unclassIndex)
                elif comp_border(base_border2, border1, overlapping_match, set_size, exact_match) and \
                        comp_border(base_border1, border2, overlapping_match, set_size, exact_match):
                    # reverse the order of the genes to match the 'reference'
                    gene_lists[unclassIndex][0] = gene_lists[unclassIndex][0][::-1]
                    # inverse the borders
                    former_border_1 = gene_lists[unclassIndex][1][0]
                    former_border_2 = gene_lists[unclassIndex][1][1]
                    gene_lists[unclassIndex][1][0] = former_border_2
                    gene_lists[unclassIndex][1][1] = former_border_1

                    # specify the new 'classified' and remove from unclassified
                    to_classify.discard(unclassIndex)
                    new_classify.add(unclassIndex)
        classified |= new_classify  # the newly classified will help to check the unclassified,
        # the formerly classified are not useful for what remains (if something remains)
        new_classify = set()
