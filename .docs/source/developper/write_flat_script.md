# Write flat files code
<style>h1 {color: black}</style>
<style>h2 {color: #829356}</style>
<style>h3 {color: #1287A8}</style>

```{warning}
This part of the documentation is in progress. Need to add the rest of the function of the script.
```

## Write flat files
<div style="text-align: justify"> 
The function "write_flat_files" triggers specific functions based on the command used by the user.
</div>
<br>

### Annotation
<div style="text-align: justify">
If the user selects the "annotation" option to write the annotations from families, the function 
"write_annotation_to_families_mp" will be triggered.
</div>
<br>

### HMM
<div style="text-align: justify"> 
Alternatively, if the user chooses the "hmm" option to create an HMM for each gene family in pangenomes, the 
functions "profile_gfs" and "write_hmm_profile" will be triggered.
</div>
<br>

### Systems
<div style="text-align: justify"> 
If the user selects the "systems" option to write all systems in pangenomes and project them on organisms, it will 
trigger the functions: </div>

["write_systems_projection"](./system_asso_script.md#w_systems_proj) 
, ["per_pan_heatmap"](./figure_script.md#per_pan_heat) 
, ["pan_distribution_system"](./system_asso_script.md#pan_distrib) 
, ["pan_number_system"](./system_asso_script.md#pan_nb) 
, ["hbar_ID_total"](./figure_script.md#hbar_ID_total)
 and ["heatmap"](./figure_script.md#heatmap).

### Systems association
<div style="text-align: justify"> 
If the user selects the "systems_asso" option to link systems and features (RGPS, spots and modules), it will trigger the 
functions from the "systems" option, as well as the functions: </div>

["systems_to_features"](./system_asso_script.md#systems2feat) 
, ["write_borders_spot"](./conserved_spot_script.md#border) 
, ["spot2sys"](./system_asso_script.md#spot2sys) 
, ["mod2sys"](./system_asso_script.md#mod2sys) 
 and ["upsetplot"](./figure_script.md#upsetplot).

### Conserved spot
<div style="text-align: justify"> 
If the user selects the "conserved_spot" option to detect conserved spots across pangenomes, it will trigger the 
functions from the "systems" and "systems_asso" options, as well as the functions: 
</div>

"all_against_all"
 and ["identical_spot"](./conserved_spot_script.md#identical).

### Draw spot
<div style="text-align: justify">
If the user selects the "draw_spot" option to draw systems in spots, it will trigger the functions from the "systems" 
and "systems_asso" options, as well as: 
</div>

["draw_spot"](./draw_spot_script.md) script.

<br>