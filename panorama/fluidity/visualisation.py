"""
The main target of this script is to develop an interface to visualise the result of comparing fluidity between
pangenomes thanks to panorama and PPanGGOLiN
"""

# #!/usr/bin/env python3
# # coding:utf-8
#
# from bokeh.plotting import figure, show
# from bokeh.sampledata.penguins import data
# from bokeh.transform import factor_cmap, factor_mark
#
# SPECIES = sorted(data.species.unique())
# MARKERS = ['hex', 'circle_x', 'triangle']
#
# p = figure(title = "Penguin size", background_fill_color="#fafafa")
# p.xaxis.axis_label = 'Flipper Length (mm)'
# p.yaxis.axis_label = 'Body Mass (g)'
#
# p.scatter("flipper_length_mm", "body_mass_g", source=data,
#           legend_group="species", fill_alpha=0.4, size=12,
#           marker=factor_mark('species', MARKERS, SPECIES),
#           color=factor_cmap('species', 'Category10_3', SPECIES))
#
# p.legend.location = "top_left"
# p.legend.title = "Species"
#
# show(p)

# if __name__ == '__main__':
#     dic = {'Pan_A': [0.1, 5, 'Sp_A', 'Gen_A', 'Fam_A'],
#            'Pan_B': [0.15, 7, 'Sp_B', 'Gen_A', 'Fam_A'],
#            'Pan_C': [0.17, 4, 'Sp_C', 'Gen_A', 'Fam_A'],
#            'Pan_D': [0.3, 10, 'Sp_D', 'Gen_B', 'Fam_A'],
#            'Pan_E': [0.28, 8, 'Sp_E', 'Gen_B', 'Fam_A'],
#            'Pan_F': [0.9, 17, 'Sp_F', 'Gen_C', 'Fam_B'],
#            'Pan_G': [0.45, 9, 'Sp_G', 'Gen_D', 'Fam_C'],
#            'Pan_H': [0.47, 10, 'Sp_H', 'Gen_D', 'Fam_C'],
#            'Pan_I': [0.52, 12, 'Sp_I', 'Gen_E', 'Fam_C'],
#            'Pan_J': [0.62, 14, 'Sp_J', 'Gen_F', 'Fam_D'],
#            'Pan_K': [0.65, 13, 'Sp_K', 'Gen_F', 'Fam_D'],
#            }
#     df = pd.DataFrame.from_dict(dic, orient='index')
#     df.index.name = 'Pangenome'
#     df.columns = ['Genomes Fluidity', 'Nb of Organisms', 'Species', 'Genus', 'Family']
#     export_graph(df)