# shastaHapPath
Creates haploid assemblies from shasta assembly graphs.

## Algorithm
![image](https://user-images.githubusercontent.com/28329271/122816002-0882a980-d28b-11eb-9b01-d5823c6a45f2.png)

## Python Dependencies
Required Python libries:  
  python bindings: [bdsg](https://github.com/vgteam/libbdsg#from-pip-python-bindings-only)  
  ``` pip install bdsg ```  
  And:  
  ```argparse, numpy, json```  

## Input files
The python program requries:
1) a bluntified vg PackedGraph file 
  bluntifiers: [getBlunted](https://github.com/vgteam/GetBlunted), gimbricate/seqwish, stark)
  TODO: Add more about going from shasta overlappy gfa to chunk_1.pg
2) vg snarls json ( [vg toolkit](https://github.com/vgteam/vg#command-line-interface) )
  TODO: Add snarls command  
  
~ not necessary because this function is not turned on yet ~  
3) Alignment files for parental Kmers  

## Running the Command:
```python3 parentPath.1.05.noParentAlns.py -g graph.pg -j graph.pg.snarls.json >> graph.pg.log.txt```  
The test directory has some toy graph examples that illustrate traversals of the graph that avoid tips in the graph, navigate through snarls in the graph, and traverse one side of bi-allelic bubbles in the graph:  
Directed Acyclic Graph (dag) and Parent Snarls (A snarl that has a snarl within it's boundires )
The test graphs have a set of 4 files:  
1) Input ```graph.pg & graph.pg.snarls.json ```    
2) Output ``` graph.path.csv ```    
3) A gfa file for viewing the graph in bandage. Here the haploid path is red ( currently the code doesn't attempt to traverse the top level Parent snarl so it is grey )
![test_dag 2 path](https://user-images.githubusercontent.com/28329271/122820486-84cbbb80-d290-11eb-8747-44c2c6348148.png)

## Output Files
1)  ``` graph.path.csv ``` A csv file that colors a gfa file when viewed in [bandage](https://rrwick.github.io/Bandage/) graph viewer

2)  ``` graph.pg.log.txt ``` A log file for development

3)  Two haploid fasta files - not turned on right now
