#
# Plotter.py
# is part of Orthoprok
#MIT License
#
#Copyright (c) 2018 Joaquin Giner Lamia
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

import os
import pandas as pd
import numpy as np
import pylab as plt
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt

import Analyzer as analyzer

PROTEIN_FIGURES = "/protein_figures/"
PROTEIN_SIZE = "protein_size_"
BLAST_PLOT_TIF = 'blast_plot.tif'
CLUSTER_MAP_TIF = 'clustermap_presence.tif'
CLUSTER_MAP_PDF = 'clustermap_presence.pdf'
CLUSTER_MAP_PRESENCE = 'clustermap_presence_pseudogenes.pdf'
CLUSTER_MAP_PSEUDO = 'clustermap_pseudogenes.pdf'

def get_protein_figures(savepath):
	return analyzer.get_resultspath(savepath)+PROTEIN_FIGURES

def get_protein_size_pdf(savepath, prot_name):
	return "".join([analyzer.get_resultspath(savepath), '/', PROTEIN_SIZE, '_', str(prot_name), ".pdf"])

def get_protein_size_tif(savepath, prot_name):
	return "".join([analyzer.get_tmpfiles(savepath), '/', PROTEIN_SIZE, '_', str(prot_name), ".tif"])

def get_blast_plot_tif(savepath):
	return analyzer.get_resultspath(savepath)+"/"+BLAST_PLOT_TIF

def get_cluster_map_tif(savepath):
	return analyzer.get_tmpfiles(savepath)+"/"+CLUSTER_MAP_TIF

def get_cluster_map_pdf(savepath):
	return analyzer.get_resultspath(savepath)+"/"+CLUSTER_MAP_PDF

def get_cluster_map_presence(savepath):
	return analyzer.get_resultspath(savepath)+"/"+CLUSTER_MAP_PRESENCE

def get_cluster_map_pseudo(savepath):
	return analyzer.get_resultspath(savepath)+"/"+CLUSTER_MAP_PSEUDO

############

def protein_length_graphic(prot_name, savepath, height, width, font_scale, err):

	""" 
	creates a seaborn graphic using the coverage of proteins stored at
	Proteins_size.csv file 
	"""
	
	resultpath = get_protein_figures(savepath)
	
	if not os.path.exists(resultpath):
		os.makedirs(resultpath)
	
	inputFile = analyzer.get_file5(savepath)#get_proteins_size_csv(savepath)
	tmpfiles = analyzer.get_tmpfiles(savepath)
	df = pd.read_csv(inputFile)
	
	sns.set(font_scale=font_scale)
	
	# index = protein_listbox1.curselection()[0]
	# n =  protein_listbox1.get(index)
	
	dfname = df['Name']
	dfprot = df[prot_name]
	
	#df2 = pd.concat([dfname, dfprot])
	sns.set_style("whitegrid")
	xticks = list(dfname)
	f, ax = plt.subplots(figsize=(width,height))
	plt.tight_layout()
	
	pal = sns.color_palette('Set2') # con esto elegimos la paleta de color, hay varias para elegir
	#ax = sns.barplot(data=df2, y=dfname, x=dfprot, errwidth=0.5, palette = pal, capsize=0.5) # Error bars represent standard deviations
	ax = sns.barplot(data=df, y="Name", x=prot_name, errwidth=0.5, palette = pal, capsize=0.5) # Error bars represent standard deviations
	
	plt.xticks(np.arange(0, 100+1, 25))
	ax.set_yticklabels(ax.get_yticklabels(), rotation = 0)
	plt.tight_layout()
	
	plt.savefig(get_protein_size_pdf(savepath, prot_name))
	plt.savefig(get_protein_size_tif(savepath, prot_name))
	plt.close()
	
	return

def stripplot (savepath, height, width, err):
	
	sns.set(style="whitegrid", color_codes=True)
	
	df = pd.read_csv(analyzer.get_file5(savepath))
	df2 = pd.read_csv(analyzer.get_file6(savepath))

	df= df.drop('Name', 1).drop('Assembly', 1)  
	df2= df2.drop('Name', 1).drop('Assembly', 1)  

	variable_count = len(df.columns)

	df = df.stack().reset_index().rename(columns={'level_1' : 'Query sequences', 0 : 'value', 'Coverage':'Coverage'})
	df2 = df2.stack().reset_index().rename(columns={'level_1' : 'Query sequences', 0 : 'value', 'Identity':'Identity'})

	df["Type"] = str("cov")
	df2["Type"] = str("iden")

	frames = [df, df2]
	result = pd.concat(frames) # concatenate both dataframe in only one df, result.
	#result = result[result.value != 0]

	err.write(str(result)+"\n")
	flatui = ["#3498db", "#e74c3c"]
	colors = sns.color_palette(flatui)

	#colors = ["windows blue",  "amber"]
	#palete = sns.xkcd_palette(colors)
	
	fig, g = plt.subplots(figsize=(width,height))
	g= sns.stripplot(x="Query sequences", y="value", hue="Type",edgecolor="black", data= result, jitter=True, dodge=True, palette=colors)
	#pl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=1,
           ncol=2, borderaxespad=0.)
	
	g.set_xticklabels(g.get_xticklabels(), rotation = 90)
	file_saved = get_blast_plot_tif(savepath)
	plt.tight_layout()
	plt.savefig(file_saved)
	
	plt.close()
	
	return

def heatmap(savepath, height, width, font_scale, right_scale, bottom_scale, err):
	
	absfile = analyzer.get_absfile(savepath)
	df = pd.read_csv(absfile)
	mergefile = analyzer.get_mergefile(savepath)
	df_merge = pd.read_csv(mergefile)
	pseudofile = analyzer.get_pseudofile(savepath)
	df_pseudo = pd.read_csv(pseudofile)
	
	err.write(str(df)+"\n")
	
	df = df.set_index('Serovar')
	df_merge = df_merge.set_index('Serovar')
	df_pseudo = df_pseudo.set_index('Serovar')

	listindex = list(df.index)
	listindex_merge = list(df_merge.index)
	listindex_pseudo = list(df_pseudo.index)
	
	#**df_absence
	xticks = list(df.columns)
	yticks = listindex

	sns.set(font_scale=font_scale)
	#--- Construct clustermap out of the data, no estoy seguro si hacerlo por method average o no!!.#
	
	g = sns.clustermap(df, method='average',
					   xticklabels=xticks, yticklabels=yticks,
					   cmap="Blues", robust=True, fmt="d", vmin=0, vmax=1,
					   linewidths=.05, linecolor='grey', figsize=(width, height))
	
	plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90) 
	plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0) 
	g.fig.suptitle('genes')
	plt.subplots_adjust(right=float(right_scale),bottom=float(bottom_scale))
	ax = g.ax_heatmap
	ax.set_ylabel("")
	
	clustermap_tif = get_cluster_map_tif(savepath)
	plt.savefig(clustermap_tif)
	clustermap_pdf = get_cluster_map_pdf(savepath)
	plt.savefig(clustermap_pdf)
	
	plt.close()
	
	#**df_merge
	xticks = list(df.columns)
	yticks = listindex_merge

	#--- Construct clustermap out of the data, no estoy seguro si hacerlo por method average o no!!.#
	g = sns.clustermap(df_merge, method='average',
					   xticklabels=xticks,yticklabels=yticks,
					   cmap="Blues",robust=True, fmt="d",
					   vmin=0, vmax=1, linewidths=.05,
					   linecolor='grey', figsize=(width, height))
	plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90) 
	plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0) 
	g.fig.suptitle('genes and pseudogenes')
	plt.subplots_adjust(right=float(right_scale),bottom=float(bottom_scale))
	ax = g.ax_heatmap
	ax.set_ylabel("")
	
	clustermap_presence = get_cluster_map_presence(savepath)
	#pl.savefig(resultspath +'clustermap_presence_pseudogenes.pdf')
	plt.savefig(clustermap_presence)
	
	plt.close()
	
	#**df_pseudo
	xticks = list(df.columns)
	yticks = listindex_pseudo

	#--- Construct clustermap out of the data, no estoy seguro si hacerlo por method average o no!!.#
	g = sns.clustermap(df_pseudo,method='average', xticklabels=xticks,yticklabels=yticks,cmap="Blues",robust=True, fmt="d", vmin=0, vmax=1, linewidths=.05, linecolor='grey', figsize=(width, height))

	plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90) 
	plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0) 
	g.fig.suptitle('Pseudogenes')
	plt.subplots_adjust(right=float(right_scale),bottom=float(bottom_scale))
	ax = g.ax_heatmap
	ax.set_ylabel("")

	clustermap_pseudo = get_cluster_map_pseudo(savepath)
	#pl.savefig(resultspath+'clustermap_pseudogenes.pdf')
	plt.savefig(clustermap_pseudo)
	plt.close()
	
	return

## END