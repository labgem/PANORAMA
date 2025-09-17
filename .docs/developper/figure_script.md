# Figures script
<style>h1 {color: black}</style>
<style>h2 {color: #829356}</style>
<style>h3 {color: #1287A8}</style>

<div style="text-align: justify"> 
PANORAMA utilizes <a href="https://github.com/bokeh/bokeh" target="_blank">Bokeh</a> in various functions to provide 
representations of the results. 
</div>
<br>

## Pangenome-specific heatmaps
<div style="text-align: justify"> 
The pangenome-specific heatmap figures display the frequency of each system type and the partitioning of system 
families for each organism.
</div>

<p id="per_pan_heat"></p>

### per_pan_heatmap
<div style="text-align: justify"> 
The function "per_pan_heatmap" processes the data and generates the following: <br>
- A list with system types. <br>
- A list with organism names. <br>
- A matrix with rows representing organisms, columns representing system types and the count of each system type
as the value. <br>
<br>
After extracting the data, the function calls <a href="#partition_heatmap" title="Go to figure_partition_heatmap 
function">"figure_partition_heatmap"</a> and <a href="#count_heatmap" title="Go to figure_count_heatmap function">
"figure_count_heatmap"</a> to draw the corresponding figures.
</div>

<p id="partition_heatmap"></p>

### figure_partition_heatmap
<div style="text-align: justify"> 
The function "figure_partition_heatmap" generates a heatmap figure with the partitions of gene families of systems in
organisms. The gray color corresponds to systems whose families are present in multiple partitions.
</div>

<p id="count_heatmap"></p>

### figure_count_heatmap
<div style="text-align: justify"> 
The function "figure_count_heatmap" generates a heatmap figure with the count of system types in organisms.
</div>
<br>

## Inter-pangenomes heatmap and histogram
<div style="text-align: justify"> 
The inter-pangenomes figures display the proportion/count of organisms with a system type compared to the total number 
of organisms.
</div>

<p id="heatmap"></p>

### heatmap
<div style="text-align: justify">  
The function "heatmap" processes the data and generates the following: <br>
- A list with pangenome names. <br> 
- A list with system types. <br>
- Two Dataframes with rows representing pangenome names, columns representing system types and the count/ratio
of each system type as the value. <br>
<i> Note: To create both Dataframes, two lists of dictionaries are generated for each pangenome, with system types
as keys and count/ratio as values. </i> <br>
After processing the data, the function calls <a href="#fig_histogram" title="Go to figure_histogram function">
"figure_histogram"</a> and <a href="#fig_heatmap" title="Go to figure_heatmap function">"figure_heatmap"</a> to draw 
the corresponding figures.
</div>

<p id="fig_histogram"></p>

### figure_histogram
<div style="text-align: justify">  
The function "figure_histogram" generates a stacked histogram figure with the count of organisms with a system type.
</div>

<p id="fig_heatmap"></p>

### figure_heatmap
<div style="text-align: justify">  
The function "figure_heatmap" generates a figure with the proportion of organisms with a system type compared to 
the total number of organisms.
A filter is applied, removing system types that are present in less than 10% of all pangenomes.
</div>
<br>

## Upsetplot
<div style="text-align: justify"> 
Upsetplot is an interactive figure that consists of three sections. The top histogram displays the system count for 
each type. The left histogram represents the system count for each organism. The central figure illustrates the 
presence/absence of systems in spots or modules predicted by
<a href="https://github.com/labgem/PPanGGOLiN" target="_blank">PPanGGOLiN</a>.
</div>

<p id="upsetplot"></p>

### upsetplot
<div style="text-align: justify">  
The function "upsetplot" processes the data and generates the figure. <br>
<br>
First, the data is processed to generate the following: <br>
- A list with system types. <br>
- A list with organism names. <br>
- A dictionary with organism name as key and organism ID as value. <br>
<i> Note: To overcome the limitation of drawing a figure with both axes being categorical data (to my knowledge), 
it is necessary to use an ID (ranging from 1 to the length of organism names) to each organism (one axis being value data).
</i> <br>
Then, the function <a href="#count" title="Go to count function">"count"</a> is utilized to obtain two lists: one with 
the number of systems per system type, and the other with the number of systems per organism.
<br><br/>
Second, the data is processed to generate lists and dictionaries.
There is two loops: <br>
- The first loop iterates through system types. <br>
- The second loop iterate through organism IDs. <br>
During each iteration of the first loop, i.e. for each system type, we check in the second loop if organisms have that 
particular system type. If they do, we store the system type and the organism ID in two separate lists. At the end of 
both loops, these two lists contain all the organisms that have the respective system types. <br>
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

<p id="count"></p>

### count
<div style="text-align: justify"> 
The function "count" generates two lists : one with the number of systems per system type, and the
other with the number of systems per organism. The lists are utilized in the <a href="#upsetplot" title="Go to upsetplot
function">"upsetplot"</a> function. 
</div>

<p id="hbar_ID_total"></p>

## Pangenome-specific histogram
### hbar_ID_total
<div style="text-align: justify"> 
The function "hbar_ID_total" generates a figure with the number of ID per system type and the total number of system 
per type. Systems with identical model can have multiple IDs if the gene families defining each system are different.
</div>
<br>
