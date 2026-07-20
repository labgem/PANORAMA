# How to Cite

## Cite PANORAMA

If you use PANORAMA in your research, please cite:

> Arnoux J, Mainguy J, Bry L, Fernandez de Grado Q, Hoblos Y, Vallenet D, Calteau A. (2026)
> **Panorama: A robust pangenome-based method for predicting and comparing biological systems across species.**
> *PLOS Computational Biology* 22(7): e1013856; doi: [https://doi.org/10.1371/journal.pcbi.1013856](https://doi.org/10.1371/journal.pcbi.1013856)

```{button-link} ../_static/citations/panorama.bib
:color: primary
:outline:
:class: new-tab

Download BibTeX
```

:::{dropdown} Show BibTeX
```bibtex
@article{arnoux_panorama_2026,
	title = {Panorama: {A} robust pangenome-based method for predicting and comparing biological systems across species},
	volume = {22},
	copyright = {All rights reserved},
	issn = {1553-7358},
	shorttitle = {Panorama},
	url = {https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013856},
	doi = {10.1371/journal.pcbi.1013856},
	abstract = {Over the last decade, the expansion in the number of available genomes has profoundly transformed the study of genetic diversity, evolution, and ecological adaptation in prokaryotes. However, traditional bioinformatic approaches based on the analysis of individual genomes are showing their limitations when faced with the sheer scale of the data. To overcome these constraints, the concept of pangenome has emerged, offering a comprehensive framework to capture the full genetic repertoire of a species. In this study, we present PANORAMA, an innovative pangenomic tool designed to exploit pangenome graphs, enabling their annotation and comparison to explore the genomic diversity of several species. Based on the PPanGGOLiN pangenome graphs, PANORAMA integrates advanced methods for rule-based prediction of macromolecular systems and comparative analysis of conserved features between different pangenomes, such as spots of insertion. We illustrate the use of PANORAMA on a dataset of 941 Pseudomonas aeruginosa genomes, evaluating its performance against reference defense system prediction tools such as PADLOC and DefenseFinder. The analysis was then extended to a larger set, including four species of Enterobacteriaceae ({\textgreater}6,000 genomes), demonstrating PANORAMA’s ability to annotate, compare, and explore the diversity and distribution of biological systems across multiple species. This work provides new methods for the large-scale comparative study of microbial genomes and highlights the relevance of pangenome approaches in deciphering their evolutionary dynamics. PANORAMA is freely available and accessible at: https://github.com/labgem/PANORAMA},
	language = {en},
	number = {7},
	urldate = {2026-07-20},
	journal = {PLOS Computational Biology},
	publisher = {Public Library of Science},
	author = {Arnoux, Jérôme and Mainguy, Jean and Bry, Laura and Grado, Quentin Fernandez de and Hoblos, Yazid and Vallenet, David and Calteau, Alexandra},
	month = jul,
	year = {2026},
	keywords = {Bacterial genomics, Escherichia coli, Genome analysis, Genome annotation, Genomics, Graphs, Hidden Markov models, Salmonella enterica},
	pages = {e1013856},
	file = {Full Text PDF:D\:\\Bibliography\\storage\\PN66NKPD\\Arnoux et al. - 2026 - Panorama A robust pangenome-based method for predicting and comparing biological systems across spe.pdf:application/pdf},
}
```
:::

---

## Cite PPanGGOLiN

PANORAMA is built on top of [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN). Please also cite:

> Gautreau G, Bazin A, Gachet M, Planel R, Burlot L, Dubois M, Perrin A, Médigue C, Calteau A, Cruveiller S, Mateus C, Gaspin C, Vallenet D, Mariadassou M.
> **PPanGGOLiN: Depicting microbial diversity via a partitioned pangenome graph.**
> *PLOS Computational Biology* 16(3): e1007732 (2020); doi: [https://doi.org/10.1371/journal.pcbi.1007732](https://doi.org/10.1371/journal.pcbi.1007732)

```{button-link} ../_static/citations/ppanggolin.bib
:color: primary
:outline:
:class: new-tab

Download BibTeX
```

:::{dropdown} Show BibTeX
```bibtex
@article{gautreau_ppanggolin_2020,
	title = {{PPanGGOLiN}: {Depicting} microbial diversity via a partitioned pangenome graph},
	volume = {16},
	issn = {1553-7358},
	shorttitle = {{PPanGGOLiN}},
	url = {https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732},
	doi = {10.1371/journal.pcbi.1007732},
	abstract = {The use of comparative genomics for functional, evolutionary, and epidemiological studies requires methods to classify gene families in terms of occurrence in a given species. These methods usually lack multivariate statistical models to infer the partitions and the optimal number of classes and don’t account for genome organization. We introduce a graph structure to model pangenomes in which nodes represent gene families and edges represent genomic neighborhood. Our method, named PPanGGOLiN, partitions nodes using an Expectation-Maximization algorithm based on multivariate Bernoulli Mixture Model coupled with a Markov Random Field. This approach takes into account the topology of the graph and the presence/absence of genes in pangenomes to classify gene families into persistent, cloud, and one or several shell partitions. By analyzing the partitioned pangenome graphs of isolate genomes from 439 species and metagenome-assembled genomes from 78 species, we demonstrate that our method is effective in estimating the persistent genome. Interestingly, it shows that the shell genome is a key element to understand genome dynamics, presumably because it reflects how genes present at intermediate frequencies drive adaptation of species, and its proportion in genomes is independent of genome size. The graph-based approach proposed by PPanGGOLiN is useful to depict the overall genomic diversity of thousands of strains in a compact structure and provides an effective basis for very large scale comparative genomics. The software is freely available at https://github.com/labgem/PPanGGOLiN.},
	language = {en},
	number = {3},
	urldate = {2022-04-04},
	journal = {PLOS Computational Biology},
	publisher = {Public Library of Science},
	author = {Gautreau, Guillaume and Bazin, Adelme and Gachet, Mathieu and Planel, Rémi and Burlot, Laura and Dubois, Mathieu and Perrin, Amandine and Médigue, Claudine and Calteau, Alexandra and Cruveiller, Stéphane and Matias, Catherine and Ambroise, Christophe and Rocha, Eduardo P. C. and Vallenet, David},
	month = mar,
	year = {2020},
	keywords = {Genomics, Genome analysis, Bacterial genomics, Algorithms, Evolutionary genetics, Structural genomics, Synthetic genomics, Taxonomy},
	pages = {e1007732},
	file = {Full Text PDF:D\:\\Bibliography\\storage\\8MS5EU6Z\\Gautreau et al. - 2020 - PPanGGOLiN Depicting microbial diversity via a pa.pdf:application/pdf;Snapshot:D\:\\Bibliography\\storage\\66X8AY7R\\article.html:text/html},
}
```
:::

---

## Cite PPanGGOLiN tools used by PANORAMA

Depending on which features you use, additional citations may be appropriate.

### Genomic islands (RGPs / spots)

If you use PANORAMA to study genomic islands, please also cite:

> Bazin A, Gautreau G, Médigue C, Vallenet D, Calteau A.
> **panRGP: a pangenome-based method to predict genomic islands and explore their diversity.**
> *Bioinformatics*, Volume 36, Issue Supplement_2, December 2020, Pages i651–i658; doi: [https://doi.org/10.1093/bioinformatics/btaa792](https://doi.org/10.1093/bioinformatics/btaa792)

```{button-link} ../_static/citations/panrgp.bib
:color: primary
:outline:
:class: new-tab

Download BibTeX
```

:::{dropdown} Show BibTeX
```bibtex
@article{bazin_panrgp_2020,
	title = {{panRGP}: a pangenome-based method to predict genomic islands and explore their diversity},
	volume = {36},
	issn = {1367-4803},
	shorttitle = {{panRGP}},
	url = {https://doi.org/10.1093/bioinformatics/btaa792},
	doi = {10.1093/bioinformatics/btaa792},
	abstract = {Horizontal gene transfer (HGT) is a major source of variability in prokaryotic genomes. Regions of genome plasticity (RGPs) are clusters of genes located in highly variable genomic regions. Most of them arise from HGT and correspond to genomic islands (GIs). The study of those regions at the species level has become increasingly difficult with the data deluge of genomes. To date, no methods are available to identify GIs using hundreds of genomes to explore their diversity.We present here the panRGP method that predicts RGPs using pangenome graphs made of all available genomes for a given species. It allows the study of thousands of genomes in order to access the diversity of RGPs and to predict spots of insertions. It gave the best predictions when benchmarked along other GI detection tools against a reference dataset. In addition, we illustrated its use on metagenome assembled genomes by redefining the borders of the leuX tRNA hotspot, a well-studied spot of insertion in Escherichia coli. panRPG is a scalable and reliable tool to predict GIs and spots making it an ideal approach for large comparative studies.The methods presented in the current work are available through the following software: https://github.com/labgem/PPanGGOLiN. Detailed results and scripts to compute the benchmark metrics are available at https://github.com/axbazin/panrgp\_supdata.},
	number = {Supplement\_2},
	urldate = {2022-04-04},
	journal = {Bioinformatics},
	author = {Bazin, Adelme and Gautreau, Guillaume and Médigue, Claudine and Vallenet, David and Calteau, Alexandra},
	month = dec,
	year = {2020},
	pages = {i651--i658},
	file = {Full Text PDF:D\:\\Bibliography\\storage\\CBR7WCR7\\Bazin et al. - 2020 - panRGP a pangenome-based method to predict genomi.pdf:application/pdf;Snapshot:D\:\\Bibliography\\storage\\AHB42DC5\\6055938.html:text/html},
}
```
:::

### Modules

If you use PANORAMA to study modules, please also cite:

> Bazin A, Vallenet D, Calteau A.
> **panModule: detecting conserved modules in the variable regions of a pangenome graph.**
> *bioRxiv* 2021.12.06.471380; doi: [https://doi.org/10.1101/2021.12.06.471380](https://doi.org/10.1101/2021.12.06.471380)

```{button-link} ../_static/citations/panmodule.bib
:color: primary
:outline:
:class: new-tab

Download BibTeX
```

:::{dropdown} Show BibTeX
```bibtex
@misc{bazin_panmodule_2021,
	title = {{panModule}: detecting conserved modules in the variable regions of a pangenome graph},
	copyright = {© 2021, Posted by Cold Spring Harbor Laboratory. This pre-print is available under a Creative Commons License (Attribution 4.0 International), CC BY 4.0, as described at http://creativecommons.org/licenses/by/4.0/},
	shorttitle = {{panModule}},
	url = {https://www.biorxiv.org/content/10.1101/2021.12.06.471380v1},
	doi = {10.1101/2021.12.06.471380},
	abstract = {The recent years have seen the rise of pangenomes as comparative genomic tools to better understand the evolution of gene content among microbial genomes in close phylogenetic groups such as species. While the core or persistent genome is often well-known as it includes essential or ubiquitous genes, the variable genome is usually less characterized and includes many genes with unknown functions even among the most studied organisms. It gathers important genes for strain adaptation that are acquired by horizontal gene transfer. Here, we introduce panModule, an original method to identify conserved modules in pangenome graphs built from thousands of microbial genomes. These modules correspond to synteny blocks composed of consecutive genes that are conserved in a subset of the compared strains. Identifying conserved modules can provide insights on genes involved in the same functional processes, and as such is a very helpful tool to facilitate the understanding of genomic regions with complex evolutionary histories. The panModule method was benchmarked on a curated dataset of conserved modules in Escherichia coli genomes. Its use was illustrated through a study of a high pathogenicity island in Klebsiella pneumoniae that allowed a better understanding of this region. panModule is freely available and accessible through the PPanGGOLiN software suite (https://github.com/labgem/PPanGGOLiN).},
	language = {en},
	urldate = {2022-04-04},
	publisher = {bioRxiv},
	author = {Bazin, Adelme and Medigue, Claudine and Vallenet, David and Calteau, Alexandra},
	month = dec,
	year = {2021},
	note = {Section: New Results
Type: article},
	file = {Full Text PDF:D\:\\Bibliography\\storage\\USJBDV63\\Bazin et al. - 2021 - panModule detecting conserved modules in the vari.pdf:application/pdf;Snapshot:D\:\\Bibliography\\storage\\NTYBG7WJ\\2021.12.06.html:text/html},
}
```
:::