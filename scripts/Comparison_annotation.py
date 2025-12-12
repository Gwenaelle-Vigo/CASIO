#! /usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np
import itertools
import re
from collections import defaultdict
from plotly_upset.plotting import plot_upset


# Make graph with node (Tool_OG) and edge (find associate with an other Tool_OG)

## Algo in python 

### Create an adjacency list : dictionary with as key each Tool_OG and as value all Tool_OG associate + itself (to do in each direction i.e Scipio_OG000001:{Scipio_OG000001, Busco_OG0000005, Miniprot_OG0000008}; Busco_OG0000005:{Scipio_OG000001, Busco_OG0000005, Miniprot_OG0000008} ; Miniprot_OG0000008:{Busco_OG0000005, Scipio_OG000001, Miniprot_OG0000008}.
### To append information of edge weight and keep information of query that give this relationship it's possible to create a dictionary Tool1_OG1:{Tool1_OG1:{Query1S1,Query2S1}, Tool2_OG5:{Query1S1,Query2S1}, Tool3_OG8:{Query1S1,Query2S1}}

#### Read table sys.argv[1]
#### Create dictionary for adjacency list


# In graph identify all Connected component.

## Algo in python

### Create a set current_component : contains the component currently create in the loop
### Create a set components : contains all components find
### Create a set Todo : for the current Tool_OG contains all linked Tool_OG relate. Will save all node (Tool_OG) to treat until the set is empty to finish the current component and save it in components.
### Create a set to_visit : set that contains at the beginning all Tool_OG of the dictionary adjacency list. During the loop all visited node are removed from this list. When this set is empty all connected component are found. 

### Append in todo a node in to_visit
### while True
### Pick in Todo a node
### Search it neighbour (associate node)
### For n in neighbour
###     if neighbour not in to visit:
###         continue
###     if n in todo:
###         continue
###     todo.insert(n)
###     tovisit.delete(n)
### insert current in current_component
### if not todo (todo is empty)
###     add curent component in components
###     current_comp=set()
###     if tovisit is empty stop the while true if not to_visit
###         break
###     append in todo a node from tovisit


# Keep only connected component that doesn't have 2 or more OG for a given tool

###### Code

# Create dictionary
global dadjacency_list
dadjacency_list = defaultdict(lambda: defaultdict(list))

# Function

## Create dictionary for adjacency list : dictionary with as key each Tool_OG and as value all Tool_OG associate
def add_og_todict(key, query):
    global dadjacency_list
    key=list(filter(lambda x: x != "X", key))
    for a in key:
        for b in key:
            dadjacency_list[a][b].append(query)

## Define function to write result in a define ordered list based on the order of tools
def myKey(x):
    global tools
    # create the 'order list' for starting pattern
    patternsList = tools
    for i in range(len(patternsList)): # iterate patterns in order
        pattern = patternsList[i]
        if x.find(pattern) == 0: # check if x starts with pattern
            # return order value i and x without the pattern
            return (i, x.replace(pattern, ''))

## Read file to create adjacency list
with open(sys.argv[1],"r") as ASSOCIATION_OG_file:
    for line in ASSOCIATION_OG_file:
        line = line.strip()
        if re.search("^#",line):
            tools=line.split("\t")[1:]
        else:
            query_id=line.split("\t")[0]
            OG=line.split("\t")[1:]
            KEY_ADJACENCY=[]
            for i,og in enumerate(OG):
                if og != "X":
                    key_adjacency=str(tools[i]) + "_" + str(og)
                    KEY_ADJACENCY.append(key_adjacency)
                else:
                    KEY_ADJACENCY.append("X")
            add_og_todict(KEY_ADJACENCY,query_id)



## Define connected component from graph (adjacency list)

# Create var
c_component=set()      # current connected componnent in while loop
c_query=set()          # current connected componnent queries list
components=list()           # list of all connected components
components_queries=list()   # list of all connected components queries list
components_concordancy=list()   # list of all connected components concordancy ("Identical" if for all queries OG are identical between tools, "Matching" if for all queries OG are concordant/similar between tools, "Incorrect" if for all queries OG are bad between tools (i.e. a tool have 2 or more OGs))
to_do=set()             # list of OG to process before to go to another connected component
to_visit=set(dadjacency_list.keys())    # list of og not process (while will end when this set is empty)

# Initialize (randomly pick an OG in list to_visit)
to_do.add(to_visit.pop())   

# Loop to browse the adjacency list and define all connected components
while True:
    c_og=to_do.pop()
    neighbour=dadjacency_list[c_og].keys()
    for n in neighbour:
        if n not in to_visit:
            continue
        if n in to_do:
            continue
        to_do.add(n)
        to_visit.remove(n)
    c_component.add(c_og)
    if c_query != set(dadjacency_list[c_og][c_og]):
        if len(c_query) == 0:
            concordancy="Identical"
            c_query.update(dadjacency_list[c_og][c_og])
        else:
            concordancy="Matching"
            c_query.update(dadjacency_list[c_og][c_og])
    if not to_do:
        components.append(c_component)
        components_queries.append(c_query)
        components_concordancy.append(concordancy)
        c_component=set()
        c_query=set()
        if not to_visit:
            break
        to_do.add(to_visit.pop())


# Keep only connected component that doesn't have 2 or more OG for a given tool
good_components=list()      # List all components that pass the filter
bad_components=list()       # List all components that doesn't pass the filter
dog_group=defaultdict(list)     # Dictionary with associate tool in key and list of good components as value e.g. Busco_Scipio_Miniprot:[{Scipio_OG000001, Busco_OG0000005, Miniprot_OG0000008},{Scipio_OG000002, Busco_OG0000025, Miniprot_OG0000010}], Busco_Scipio:[{Scipio_OG000003, Busco_OG0000006}], Busco_Miniprot[{Busco_OG0000007, Miniprot_OG0000014},{Busco_OG0000008, Miniprot_OG0000015}], Busco:[{Busco_OG0000008},{Busco_OG0000009},{Busco_OG0000010}], ...

# Read all components
with open(sys.argv[2]+"/NewListGene_OG_GroupByAnnotationTools.tsv",'a') as fo:
    df_upset = pd.DataFrame(columns=tools)
    loc_index=0
    for i,c in enumerate(components):
        n=list()
        og_group=list()
        all_og_group=list()
        upset=list()
        for motif in tools: # Define in component the encountered tools and occurence of each tools
            count = sum(len(re.findall(motif, s)) for s in c)    
            n.append(count)
            if count>0:
                og_group.append(motif)
                all_og_group.append(motif)
                upset.append(1)
            else:
                upset.append(0)
                all_og_group.append("X")
        if max(n)>=2:
            components_concordancy[i]="Incorrect"
            bad_components.append(c)
        else:
            good_components.append(c)
            dog_group["_".join(og_group)].append(c)
            df_upset.loc[loc_index] = upset
            loc_index+=1 
        lc=list(c)
        lc.sort(key = myKey)
        fo.write("OG"+str(i+1).zfill(7)+ "_" + "".join(all_og_group)  + "\t" + "_".join(og_group) + "\t" + "|".join(lc) + "\t" + "|".join(components_queries[i]) + "\t" + str(len(components_queries[i])) + "\t" + components_concordancy[i] + "\n")

fo.close()


fig = plot_upset(
    dataframes=[df_upset],
    legendgroups=["OGs"],
    marker_size=16,
    exclude_zeros=True,
    sorted_x="d",
    sorted_y="a",
)

fig.write_image(sys.argv[2]+ "/UpsetPlot_OG_AnnotationTools.png")











#########################################################################################
#########################################################################################

# OLD VERSION 
### Test with pandas

# df = pd.read_csv(sys.argv[1], sep="\t")

# # Dynamically determine the columns to group by (all except the first column). Each column is an annotation tool.
# group_cols = df.columns[1:]

# # Query name
# id_col = df.columns[0]

# # Create all combinations of the tool list
# for L in range(len(group_cols) + 1):
    # for subset in itertools.combinations(group_cols, L):
        # list(subset)

# # Perform the aggregation
# result = (
    # df.groupby(list(group_cols), sort=False).agg(ID_List=(id_col, lambda x: ';'.join(x)), Count=(id_col, 'count')).reset_index()
# )

# # Save or print the result
# result.to_csv(sys.argv[2], sep="\t", index=False, encoding='utf-8')


# # Identify tool columns (exclude ID_List and Count)
# tool_cols = [col for col in result.columns if col not in ["ID_List", "Count"]]



# # For each tool column, find values that occur only once
# map_df = pd.DataFrame(index=result.index)

# # For each tool column, compute the uniqueness map
# for col in tool_cols:
    # counts = result[col].value_counts() #Count each value in each column 
    # col_map = result[col].map(counts) == 1  # Set as False if count is more than 1 (not unique)
    # col_map = col_map.where(result[col] != "X", np.nan)  # Replace False where the value was 'X' (=NaN) 
    # map_df[col] = col_map

# # Nombre d'orthogroup qui sont True partout (sans compter NaN)
# map_df_true=map_df[~map_df.eq(False).any(axis=1)]
# # Display initial row of this row
# df_true=result[~map_df.eq(False).any(axis=1)]


# # Orthogroup qui ont un False
# map_df_false=map_df[map_df.eq(False).any(axis=1)]
# # Display initial row of this row
# df_false=result[map_df.eq(False).any(axis=1)]



