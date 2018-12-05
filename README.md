# Blaster

Blaster uses Blast(tblastn) to perform homologus proteins search in a set of target genomes. 

## Require

python modules:
Tkinter, biopython, pillow, seaborn 

other binary programs:
1. ncbi blast
2. mview
3. muscle

Important Note: Orthoprok.gif file has to be at the same path that Blaster.py.

## Usage

<br />

`python Blaster.py`

<br />

## Input

Proteins and genomes sequences have to be in separated fasta (.fa, .fna, .fasta) files (multifasta is not supported).


## Output

Blaster analysis generates two different folders. /Result folder that contains all the analysis result and /tmp folder in which several intermediate-tables analysis are stored.

/result/ file description:

- Analysis_presence_summary.csv = table containning number of copies found by blast for query proteins in every genome.
- Analysis_pseudo_summary.csv = table containning number of pseudogenes found by blast for query proteins in every genome.
- Protein_ident.csv = table containning % identity value of proteins found by blast in every genome.
- Protein_size.csv = table containning % coverage value of proteins found by blast in every genome.
- summarystats.csv = number of total queries proteins found in all genomes, this data are graphical representaed in the grpahic_stats.pdf file.

In /result/blast/ you can find all the blast hits found for every query in the target genomes according to coverage and identity parameters used in the analysis.


Alignment analysis:

all the protein sequences found according to coverage and identity values set are going to be alignned with MUSCLE and reformated MVIEW.


## Graphic analysis


In this section different graphic analysis can be performed and width, height and font scale of the output  figure could be controlled.

### Protein plot analysis:
In these graphics the coverage of all homologous proteins found in the target genomes for a determinated query protein will be plotted. Resulting graphic will be pop up but it will also be saved at /result/proteins/


### Heatmap analysis:

next files present in /tmp:
- groupbyserovar_1normalized_absence.csv
- groupbyserovar_1normalized_merge.csv
- groupbyserovar_1normalized_pseudo.csv

are used to plot presence/absence of query proteins in all genomes. In the case collapse genomes option is set, all genome from the same species will be collapsed in the heatmap.
This analysis will generate three differente graphic files:

- clustermap_presence.pdf
- clustermap_pseudogenes.pdf
- clustermap_presence_pseudogenes.pdf



### Blast plot:

This function plot the coverage and the identity of all query protein making easy the visual inspection of the whole Blast analysis.

This figure is stored in blast_plot.tif file




