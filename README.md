# shastaHapPath
Creates haploid assemblies from shasta assembly graphs.

## Algorithm
![image](https://user-images.githubusercontent.com/28329271/122816002-0882a980-d28b-11eb-9b01-d5823c6a45f2.png)

## Python Dependencies
Required Python libries:
[bdsg](https://github.com/vgteam/libbdsg#from-pip-python-bindings-only)
``` pip install bdsg ```
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

## Output Files
A csv file that colors a gfa file when viewed in [bandage](https://rrwick.github.io/Bandage/) graph viewer
``` graph.path.csv ```
A log file for development
``` graph.pg.log.txt ```
Two haploid fasta files - not turned on right now
