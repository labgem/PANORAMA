# Output figures code
```{warning}
This part of the documentation is in progress.
```
## Figure script
<div style="text-align: justify"> 
PANORAMA utilizes Bokeh in various functions to provide representations of the results. 
</div>
<br><br/>

### Functions to draw heatmap pangenome-specific figures
<div style="text-align: justify"> 
The pangenome-specific heatmap figures display the frequency of each system type and the partitioning of system 
families for each organism.
</div>
<br><br/>

#### per_pan_heatmap
<div style="text-align: justify"> 
The function "per_pan_heatmap" processes the data and generates the following: <br>
- A list with system types. <br>
- A list with organism names. <br>
- A matrix with rows representing organisms, columns representing system types and the count of each system type
as the value. <br>
After extracting the data, the function calls "figure_partition_heatmap" and "figure_count_heatmap" to draw 
the corresponding figures.
</div>
<br><br/>

#### figure_partition_heatmap
<div style="text-align: justify"> 
The function "figure_partition_heatmap" generates a heatmap figure with the partitions of gene families of systems in
organisms. The gray color corresponds to systems whose families are present in multiple partitions.
</div>
<br><br/>

#### figure_count_heatmap
<div style="text-align: justify"> 
The function "figure_count_heatmap" generates a heatmap figure with the count of system types in organisms.
</div>
<br><br/>

### Functions to draw heatmap and histogram inter-pangenomes figures
<div style="text-align: justify"> 
The inter-pangenomes figures display the proportion/count of organisms with a system type compared to the total number 
of organisms.
</div>
<br><br/>

#### heatmap
<div style="text-align: justify">  
The function "heatmap" processes the data and generates the following: <br>
- A list with pangenome names. <br> 
- A list with system types. <br>
- Two Dataframes with rows representing pangenome names, columns representing system types and the count/ratio
of each system type as the value. <br>
<i> Note: To create both Dataframes, two lists of dictionaries are generated for each pangenome, with system types
as keys and count/ratio as values. </i> <br>
After processing the data, the function calls "figure_histogram" and "figure_heatmap" to draw the corresponding figures.
</div>
<br><br/>

#### figure_histogram
<div style="text-align: justify">  
The function "figure_histogram" generates a stacked histogram figure with the count of organisms with a system type.
</div>
<br><br/>

#### figure_heatmap
<div style="text-align: justify">  
The function "figure_heatmap" generates a figure with the proportion of organisms with a system type compared to 
the total number of organisms.
A filter is applied, removing system types that are present in less than 10% of all pangenomes.
</div>
<br><br/>

### Functions to draw upsetplot figure
<div style="text-align: justify"> 
Upsetplot is an interactive figure that consists of three sections. The top histogram displays the system count for 
each type. The left histogram represents the system count for each organism. The central figure illustrates the 
presence/absence of systems in spots or modules predicted by
<a href="https://github.com/labgem/PPanGGOLiN" target="_blank">PPanGGOLiN</a>.
</div>
<br><br/>

#### upsetplot
<div style="text-align: justify">  
The function "upsetplot" processes the data and generates the figure. <br>
First, the data is processed to generate the following: <br>
- A list with system types. <br>
- A list with organism names. <br>
- A dictionary with organism name as key and organism ID as value. <br>
<i> Note: To overcome the limitation of drawing a figure with both axes being categorical (to my knowledge), 
it is necessary to use an ID (ranging from 1 to the length of organism names) to each organism. </i> <br>
Then, the function "count" is utilized to obtain two lists: one with the number of systems per type of system, and the
other with the number of systems per organism.
<br><br/>
Second, the data is processed to generate lists and dictionaries.
There is two loops: <br>
- The first loop iterates through system types. <br>
- The second loop iterate through organism IDs. <br>
During each iteration of the first loop, i.e. for each system type, we check in the second loop if organisms have that 
particular system type. If they do, we store the system type and the organism ID in two separate lists. At the end of 
both loops, these two lists contain all the organisms that have the respective system types.
Additionally, we verify whether the system type exists in the dataframe connecting system types with spots/modules. 
If they do, we add this information to four dictionaries, i.e. spot_ID/module_ID as key and system_type/organism_ID 
as value. At the end of both loops, these dictionaries will contain organisms that possess systems in a feature 
(either spot/module or both).
<br><br/>
Next, the three sections of the upsetplot are generated : <br>
- The top histogram displays the system count for each type. <br>
- The left histogram represents the system count for each organism. <br>
- The central figure illustrates the presence/absence of systems in spots or modules. The horizontal axis contains system
types and the vertical axis contains organism names. White squares represent all systems, while colored squares 
indicate systems present in a spot. Colored circles represent systems within a module. The figure includes a clickable 
legend that enables users to hide/show features.
</div>
<br><br/>

#### count
<div style="text-align: justify"> 
The function "count" generates two lists : one with the number of systems per type of system, and the
other with the number of systems per organism. The lists are utilized in the "upsetplot" function. 
</div>

### Function to draw histogram pangenome-specific figure

#### hbar_ID_total
<div style="text-align: justify"> 
The function "hbar_ID_total" generates a figure with the number of ID per system type and the total number of system 
per type. Systems with identical model can have multiple IDs if the gene families defining each system are different.
</div>
<br><br/>

## Conserved spot script

## Draw spot script
