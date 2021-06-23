# shastaHapPath
Creates haploid assemblies from [shasta](https://github.com/chanzuckerberg/shasta) assembly graphs. Shasta assembly graphs are bidirectional graphs where nodes represent DNA sequence and the edges reprensent links between different DNA segments. Bubbles represent variation in the diploid DNA assembly graph.   
The ( [vg toolkit](https://github.com/vgteam) ) provides tools for working with graphs that have blunt edges. Shasta graphs are produced with overlapping edges and need to be bluntified before it can be converted into a vg graph.   
This program uses the boundry nodes, identified using vg snarls, to travers the differnt structures within the graph and identify haploid walks through the graph. Those walks are output as haploid fasta files.  
 
If parental kmer alignment files are provided haploid paths can be resolved. 

## Algorithm
![image](https://user-images.githubusercontent.com/28329271/122816002-0882a980-d28b-11eb-9b01-d5823c6a45f2.png)

## Python Dependencies
  python bindings: [bdsg](https://github.com/vgteam/libbdsg#from-pip-python-bindings-only)  
  ``` pip install bdsg ```  
  And:  
  ```argparse, numpy, json```  

## Input files
1) a bluntified vg PackedGraph file   
  bluntifiers: [getBlunted](https://github.com/vgteam/GetBlunted), gimbricate/seqwish, stark)  
  TODO: Add more about going from shasta overlappy gfa to chunk_1.pg
2) vg snarls json ( [vg toolkit](https://github.com/vgteam/vg#command-line-interface) )  
  TODO: Add snarls command  
  
~ not necessary because this function is not turned on yet ~  
3) Alignment files for parental Kmers  

## Running the Command:
```python3 parentPath.1.05.noParentAlns.py -g graph.pg -j graph.pg.snarls.json > graph.pg.log.txt```  

## Output Files
1)  ``` graph.path.csv ``` A csv file that colors a gfa graph when viewed in [bandage](https://rrwick.github.io/Bandage/) graph viewer
2)  ``` graph.pg.log.txt ``` A log file for development, graph stats (homozygous node lengths, mean heterozygous node lengths, etc.) are at the final lines of this file for now.   
3)  Two haploid fasta files - not turned on right now

## Toy Examples
The test directory has some toy graph examples that illustrate traversals of the graph that avoid tips in the graph, navigate through snarls in the graph, and traverse one side of bi-allelic bubbles in the graph.    
  
Tests are graphs that contain different snarl types:
- directed acyclic graphs (dag)  
- parent snarls : A top level snarl that has a snarl within it's boundires and which has no parents itself
- bi-allelic bubbles : bubbles with just two nodes
  
The test graphs are a set of 4 files:  
1) Two Input Files ```test_graph.pg``` & ```test_graph.pg.snarls.json ```    
2) One Output ```test_graph.path.csv ```    
3) A ```test_graph.pg.gfa``` file for viewing the graph in bandage. Here the haploid path is red ( currently the code doesn't attempt to traverse the top level Parent snarl so it is grey ).  
  
This example graph has a parent snarl (not covered in grey) and dag. The tip is avoided and the graph is covered from end to end. 
```python3 parentPath.1.05.noParentAlns.py -g test_dag.2.pg -j test_dag.2.pg.snarls.json > test_dag.2.pg.log.txt```  
![test_dag 2 path](https://user-images.githubusercontent.com/28329271/122820486-84cbbb80-d290-11eb-8747-44c2c6348148.png)
  
This graph has a dag that is an inversion just after the parent snarl. The inversion is included in the path.  
```python3 parentPath.1.05.noParentAlns.py -g test_dag.inv.pg -j test_dag.inv.pg.snarls.json > test_dag.inv.pg.log.txt```  
![test_dag inv path](https://user-images.githubusercontent.com/28329271/122843356-1ac50d80-d2b4-11eb-8eac-44b18a4bfbdf.png)

  
The nodes in these bidirectional graph represent DNA sequence and the edges reprensent links between different DNA segments. Bubbles represent variation in the diploid DNA assembly graph. The red path is an arbitrary walk through the graph representing a haploid path through the graph.   
  
Phasing the biallelic bubbles is the next step once traversal is solid.
