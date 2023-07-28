# System projection and association script
<style>h1 {color: black}</style>
<style>h2 {color: #829356}</style>
<style>h3 {color: #1287A8}</style>

```{warning}
This part of the documentation is in progress.
```

<p id="w_systems_proj"></p>

## Systems projection
### write_systems_projection
<div style="text-align: justify"> 
The function "write_systems_projection" is triggered by the "systems", "systems_asso", "conserved_spot" or "draw_spot" 
options of PANORAMA.
That function writes all systems in pangenomes and project them on organisms.
In the first part of the function, each system is projected onto all organisms of the pangenome, and for this purpose, 
the <a href="#w_sys_proj" title="Go to write_system_projection function">"write_system_projection"</a> function is 
called. <br>
After obtaining the results, they are processed to create a final dataframe, which will be saved and utilized 
in other functions.
</div>

<p id="w_sys_proj"></p>

### write_system_projection
<div style="text-align: justify"> 
The function "write_system_projection" projects a system onto all organisms of the pangenome. <br>
Firstly, the mandatory and accessory families of the system model are stored into two separate lists. <br>
Next, a list of organisms containing at least one family from the system families is generated. That list will be 
iterated through, and for each organism, the function <a href="#w_sys_org" title="Go to write_sys_to_organism function">
"write_sys_to_organism"</a> will be called. The <a href="#w_sys_org" title="Go to write_sys_to_organism function">
"write_sys_to_organism"</a> function will return a boolean value indicating whether the system is present in the 
organism. <br> 
If the boolean value is False, then another loop will iterate through the canonical form of the system. 
The mandatory and accessory families of the canonical system model will be stored in two separate lists. The 
<a href="#w_sys_org" title="Go to write_sys_to_organism function">"write_sys_to_organism"</a> function will then be 
called once again on the canonical form.
</div>

<p id="w_sys_org"></p>

### write_sys_to_organism
<div style="text-align: justify"> 
The function "write_sys_to_organism" projects a system on an organism. <br>
</div>

```{warning}
This function documentation is in progress.
```

<p id="pan_distrib"></p>

## Distribution of system types
<div style="text-align: justify">
The function "pan_distribution_system" is triggered by the "systems", "systems_asso", "conserved_spot" or "draw_spot" 
options of PANORAMA. <br>
The function returns a dataframe that represents the distribution of system types across organisms of the pangenome.
</div>

<p id="pan_nb"></p>

## ID and total count per system types
<div style="text-align: justify">
The function "pan_number_system" is triggered by the "systems", "systems_asso", "conserved_spot" or "draw_spot" 
options of PANORAMA. <br>
The function returns two dataframes that contain the count of ID and systems per system type.
</div>

<p id="systems2feat"></p>

## Systems and features association
### systems_to_features
<div style="text-align: justify">
The function "systems_to_features" is triggered by the "systems_asso", "conserved_spot" or "draw_spot" 
options of PANORAMA. <br>
The function links systems with features such as RGPs, spots and modules detected by 
<a href="https://github.com/labgem/PPanGGOLiN" target="_blank">PPanGGOLiN</a>. <br>
During the initial stage of the function, each system is projected onto features such as RGPs and modules within the 
pangenome. To accomplish this, the <a href="#sys2feat" title="Go to system_to_features function">"system_to_features"</a>
function is called. <br>
Once the results are obtained, they undergo processing to generate five dictionaries: <br>
- Two dictionaries with systems as keys and their corresponding RGPs/modules as values. <br>
- Two dictionaries with systems as keys and the organisms of the corresponding RGPs/modules as values. <br>
- One dictionary with RGPs as keys and their associated spot IDs as values. <br>
<br>
Following this, the function <a href="#w_sys2feat" title="Go to write_systems_to_features function">
"write_systems_to_features"</a> is called to write the results of the systems/features association.
</div>

<p id="sys2feat"></p>

### system_to_features
<div style="text-align: justify">
The function "system_to_features" associates a system to features (RGPs and modules). <br>
Firstly, the mandatory and accessory families of the system model are stored into two separate lists. <br>
<br>
Next, RGPs/modules of the pangenome are iterated through, and for each RGP/module, the function 
<a href="#w_sys_in_feat" title="Go to write_sys_in_feature function">"write_sys_in_feature"</a> 
will be called. The <a href="#w_sys_in_feat" title="Go to write_sys_in_feature function">"write_sys_in_feature"</a> 
function will return a boolean value indicating whether the system is present in the RGP/module. <br>
The "system_to_features" function returns the system ID and two lists containing the detected RGPs/modules.
</div>

<p id="w_sys_in_feat"></p>

### write_sys_in_feature
<div style="text-align: justify"> 
The function "write_sys_in_feature" projects a system on a feature (RGP or module). <br>
</div>

```{warning}
This function documentation is in progress.
```

<p id="w_sys2feat"></p>

### write_systems_to_features
<div style="text-align: justify"> 
The function "write_systems_to_features" write the results of the systems/features association. <br>
A list containing system IDs associated with RGPs or modules is generated. This list is iterated through, and for each 
system ID, additional information is provided, including the system name, modules' organisms, module IDs, RGPs' 
organisms, Spot IDs, and RGP IDs. <br>
The content of the output file can vary depending on the chosen "system_asso" option.
</div>

<p id="spot2sys"></p>

## System Composition within Spots
### spot2sys
<div style="text-align: justify"> 
The function "spot2sys" write spots associated to systems and extract spot content. To achieve this, 
there are three loops: <br>
- In the first loop, we iterate through spots that are linked to systems. For each spot, we filter a dataframe to retain
the system IDs associated with that spot. <br>
- In the second loop, we iterate through the rows/system_IDs of the filtered dataframe. Each system ID contains a list 
of spots, which can be different. <br>
- In the third loop, we iterate through that list and check if the spot of the third loop is similar to the spot of the
first loop. If this is the case, we store information in dictionaries. <br>
Upon completing the loops, we obtain three dictionaries with spot IDs as keys and system IDs, system names, and 
organism names as values, respectively. Additionally, there is one dictionary with spot IDs as keys and border 
information as values. Subsequently, all of these dictionaries are converted into a dataframe.
<br>
Additional information, such as the number of systems and the number of organisms, is incorporated into the dataframe.
<br>
Following this, the function <a href="#extract_spot" title="Go to extract_spot_content function">"extract_spot_content"</a>
is called to include details about the content of each spot.
</div>

<p id="extract_spot"></p>

### extract_spot_content
<div style="text-align: justify"> 
The function "extract_spot_content" supplements the details about the content of each spot, such as the mean number of 
genes, as well as the maximum and minimum number of genes contained within the spot.
</div>

<p id="mod2sys"></p>

## System Composition within Modules
### mod2sys
<div style="text-align: justify"> 
The function "mod2sys" write modules associated to systems. To achieve this, there are three loops: <br>
- In the first loop, we iterate through modules that are linked to systems. For each module, we filter a dataframe to 
retain the system IDs associated with that module. <br>
- In the second loop, we iterate through the rows/system_IDs of the filtered dataframe. Each system ID contains a list 
of modules, which can be different. <br>
- In the third loop, we iterate through that list and check if the module of the third loop is similar to the module of the
first loop. If this is the case, we store information in dictionaries. <br>
Upon completing the loops, we obtain two dictionaries with module IDs as keys and corresponding system names and 
organism names as values, respectively. Following this, all of these dictionaries are converted into a dataframe.
<br>
Additional information, such as the number of systems and the number of organisms, is incorporated into the dataframe.
</div>
<br>