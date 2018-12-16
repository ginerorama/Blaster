#!/home/gargamelle/anaconda2/bin/python
#/usr/bin/python
#Blaster v0.1
#
# Blaster.py
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

##########################################################################
#
#                   Import  Modules
#
##########################################################################

import Tkinter, tkFileDialog, Tkconstants 
from tkFileDialog import *
from tkFileDialog import askopenfilename
from Tkinter import *
from ttk import *
import tkMessageBox
import os,sys
import numpy as np
import pandas as pd
import re
from Bio import SeqIO
from PIL import Image, ImageTk

# CPC2018
import gzip

import subprocess

##########################################################################
#
#                   Initzialize main Tkinter window
#
##########################################################################

window=Tkinter.Tk()
window.geometry("600x800")
window.title("Blaster")
dir_path = os.path.dirname(os.path.realpath(__file__))
imagen1=PhotoImage(file=dir_path+"/orthoprok.gif") 
label1 = Label(window, image=imagen1).pack() 

#s = Style()
#s.theme_names()
#s.theme_use('aqua')

##########################################################################
#
#                   Tkinter components
#
##########################################################################

OpenName = StringVar() # variable that stores target sequences path to show in GUI
SaveName = StringVar()  # variable that stores project directory to show in GUI
QueryName = StringVar() # variable that stores query sequences to show in GUI
displayedText = StringVar()
height = IntVar()
width = IntVar()
font_scale = DoubleVar()
right_scale = DoubleVar()
bottom_scale = DoubleVar()
cov = IntVar() 
ident = IntVar() 
function = StringVar()
redo = IntVar()
collapsed =IntVar()
absence = IntVar()

##########################################################################
#
#                   Commands
#
##########################################################################

import src.Facade as facade
import src.GenomeParser as genomeParser
import src.Analyzer as analyzer
import src.utils as utils

paths_dict = {}
known_query_list = []
known_target_list = []

OPENPATH = "openpath"
DEFOPENPATH = "/home/gargamelle/Blaster_data"
SAVEPATH = "savepath"
DEFSAVEPATH = "/home/gargamelle/Blaster_results"
QUERYPATH = "querypath"
DEFQUERYPATH = "/home/gargamelle/Blaster_data/proteins"

def _ask_directory(FieldPath, paths_dict, FIELDPATH):
	dirpath = askdirectory()
	#dirpath = "/".join(dirpath.split("/")[2::])
	FieldPath.set(dirpath)
	paths_dict[FIELDPATH] = dirpath
	return dirpath

def _display(text):
	displayedText.set(text)
	output_text.update_idletasks()
	window.update()
	return

def iterateDirNs(paths_dict, UPDATE):
	"""
	Generate the log file input for blaster v2.0
	automatic detection of .gz compressed genomes from genbank an uncompress
	log file also contain in the 3th column from fasta header file
	Arguments: /fasta_files_folder log_basename_file

	"""
	
	subjectDir = paths_dict[OPENPATH]+"/"
	outputDir = paths_dict[SAVEPATH]+"/"
	
	try:
		dict_file = facade.updatelog(outputDir, subjectDir, UPDATE, sys.stderr) # CPC2018
		
		if os.path.exists(dict_file):
			_display("Target sequences database done.")
			sys.stderr.write("Target sequences database done.\n")
		else:
			raise Exception("Error in target sequences database!!")
		
	except Exception as e:
		sys.stderr.write(str(e)+"\n")
		_display(str(e))
	
	return

def Blast (paths_dict, REDO):
	
	sys.stderr.write("Starting Blast command...\n")
	
	querypath = paths_dict[QUERYPATH]
	savepath = paths_dict[SAVEPATH]
	openpath = paths_dict[OPENPATH]
	
	utils.change_extension_fasta(openpath)
	utils.change_extension_fasta(querypath)

	queryDir = querypath+"/"
	subjectDir = openpath+"/"
	outputDir = savepath+"/output/"
	
	#creates output directory in project directory
	newpath = savepath+"/output"
	if not os.path.exists(newpath):
		os.makedirs(newpath)
	
	facade.Blast(queryDir, subjectDir, outputDir, 
				known_query_list, known_target_list,
				REDO , _display, sys.stderr) # CPC2018
	
	tkMessageBox.showinfo("ORTHOPROK","Blast is finished!!")				
	displayedText.set('Blast is done!!')
	
	sys.stderr.write("Finished Blast command\n")
	
	return

def analysis(paths_dict, cov, ident, collapsed, absence, height, width): #statsA + statisticsA
	cov = int(cov.get())
	ident = int(ident.get())
	collapsed = int(collapsed.get())
	absence = int(absence.get())
	height =int(height.get())
	width = int(width.get())
	savepath = paths_dict[SAVEPATH]
	#openpath = paths_dict[OPENPATH]
	
	try:
		facade.analysis(savepath, cov, ident, collapsed, absence, sys.stderr)
	except Exception as e:
		_display(str(e))
		sys.stderr.write(str(e)+"\n")

	tkMessageBox.showinfo("ORTHOPROK","analysis is done!!")

	# load the protein listbox with data
	df_protein = pd.read_csv(analyzer.get_file5(savepath))

	protein_list= df_protein.ix[:,2::].columns.tolist()

	for item in protein_list:
		protein_listbox1.insert(END,item)

	novi = Toplevel()
	image = Image.open(analyzer.get_tmpfiles(savepath) +'graphic_stats.tif') 
	photo = ImageTk.PhotoImage(image)
	#image not visual
	
	canvas = Canvas(novi, width = (width*100), height = (height*100),
					scrollregion = (0,0,(width*100),(height*100)))
	canvas.create_image(0, 0,image = photo, anchor="nw")
	canvas.photo = photo
		
	hbar=Scrollbar(novi,orient=HORIZONTAL)
	hbar.pack(side=BOTTOM,fill=X)
	hbar.config(command=canvas.xview)
	vbar=Scrollbar(novi,orient=VERTICAL)
	vbar.pack(side=RIGHT,fill=Y)
	vbar.config(command=canvas.yview)
		
	canvas.config(width=900,height=600)
	canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
	canvas.pack(side=LEFT,expand=True,fill=BOTH)
	canvas.pack(expand = YES, fill = BOTH)
	
	return

def protein_length_graphic(savepath,height,width,font_scale):

	""" 
	creates a seaborn graphic using the coverage of proteins stored at
	Proteins_size.csv file 
	"""
	import pylab as plt
	import seaborn as sns; sns.set()


	height = int(height.get())
	width = int(width.get())
	font_scale = float(font_scale.get())


	outputDir = savepath+"/results/"
	resultpath= outputDir+"/protein_figures/"
	
	if not os.path.exists(resultpath):
		os.makedirs(resultpath)

	inputFile = outputDir+"Proteins_size.csv"
	tmpfiles = savepath+"/tmp/"
	df = pd.read_csv(inputFile)

	sns.set(font_scale=font_scale)


	try:
		index = protein_listbox1.curselection()[0]
		n =  protein_listbox1.get(index)
		
		df2 = pd.concat([df['Name'], df[n]])
		sns.set_style("whitegrid")
		xticks = list(df['Name'])
		f, ax = plt.subplots(figsize=(width,height))
		plt.tight_layout()	

		pal = sns.color_palette('Set2') # con esto elegimos la paleta de color, hay varias para elegir
		ax = sns.barplot(data=df2, y=df['Name'], x=df[n], errwidth=0.5, palette = pal, capsize=0.5) # Error bars represent standard deviations
		plt.xticks(np.arange(0, 100+1, 25))
		ax.set_yticklabels(ax.get_yticklabels(), rotation = 0)
		plt.tight_layout()

		#fig = plt.gcf()
		#size = fig.get_size_inches()*fig.dpi # con esto obtenemos los dpi de la figura
		#print size

		plt.savefig(resultpath+'/protein_size_'+str(n)+".pdf")
		plt.savefig(tmpfiles+'/protein_size_'+str(n)+".tif")
		plt.close()
	
		
		novi = Toplevel()
		image = Image.open(tmpfiles +'/protein_size_'+str(n)+".tif")
		photo = ImageTk.PhotoImage(image)

		
		
		canvas = Canvas(novi, width = (width*100), height = (height*100),scrollregion = (0,0,(width*100),(height*100)))
		canvas.create_image(0, 0,image = photo, anchor="nw")
		canvas.photo = photo
		
		hbar=Scrollbar(novi,orient=HORIZONTAL)
		hbar.pack(side=BOTTOM,fill=X)
		hbar.config(command=canvas.xview)
		vbar=Scrollbar(novi,orient=VERTICAL)
		vbar.pack(side=RIGHT,fill=Y)
		vbar.config(command=canvas.yview)
	

		canvas.config(width=800,height=1200)
		canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
		canvas.pack(side=LEFT,expand=True,fill=BOTH)
		canvas.pack(expand = YES, fill = BOTH)





	except IndexError:
		displayedText.set('all proteins will be plotted!!')
		column_names = list(df)
		proteins =column_names[2::]
	
		for n in proteins:

			df2 = pd.concat([df['Name'], df[n]])
			sns.set_style("whitegrid")
			xticks = list(df['Name'])
			f, ax = plt.subplots(figsize=(width,height))
			plt.tight_layout()	

			pal = sns.color_palette('Set2') # con esto elegimos la paleta de color, hay varias para elegir
			ax = sns.barplot(data=df2, y=df['Name'], x=df[n], errwidth=0.5, palette = pal, capsize=0.5) # Error bars represent standard deviations
			plt.xticks(np.arange(0, 100+1, 25))
			ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize = 4.5)
		
			plt.savefig(resultpath+'/protein_size_'+str(n)+".pdf")
			plt.close()

	displayedText.set('Protein plot done!!')
	
	return

def stripplot (savepath,height,width):

	import matplotlib.pyplot as plt
	import seaborn as sns

	height = int(height.get())
	width = int(width.get())

	resultspath = savepath+"/results/"


	sns.set(style="whitegrid", color_codes=True)


	df = pd.read_csv(resultspath+"Proteins_size.csv")
	df2 = pd.read_csv(resultspath+"Proteins_ident.csv")

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

	print result
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
	file_saved = resultspath+'blast_plot.tif'
	plt.tight_layout()
	plt.savefig(file_saved)
	
	plt.close()



	novi = Toplevel()
	image = Image.open(file_saved)
	photo = ImageTk.PhotoImage(image)

	canvas = Canvas(novi, width = (width*100), height = (height*100),scrollregion = (0,0,(width*100),(height*100)))
	canvas.create_image(0, 0,image = photo, anchor="nw")
	canvas.photo = photo
	hbar=Scrollbar(novi,orient=HORIZONTAL)
	hbar.pack(side=BOTTOM,fill=X)
	hbar.config(command=canvas.xview)
	vbar=Scrollbar(novi,orient=VERTICAL)
	vbar.pack(side=RIGHT,fill=Y)
	vbar.config(command=canvas.yview)
	
	
	canvas.config(width=600,height=1200)
	canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
	canvas.pack(side=LEFT,expand=True,fill=BOTH)
	canvas.pack(expand = YES, fill = BOTH)
	
	return

def heatmap(savepath,height,width,font_scale,right_scale, bottom_scale):
	


	height = int(height.get())
	width = int(width.get())
	font_scale = float(font_scale.get())
	bottom_scale = float(bottom_scale.get())
	right_scale = float(right_scale.get())


	import pylab as pl
	import seaborn as sns; sns.set()

	tmpfiles = savepath+"/tmp/"
	resultspath = savepath+"/results/"

	df = pd.read_csv(tmpfiles +'groupbyserovar_1normalized_absence.csv')
	df_merge = pd.read_csv(tmpfiles +'groupbyserovar_1normalized_merge.csv')
	df_pseudo= pd.read_csv(tmpfiles+'groupbyserovar_1normalized_pseudo.csv')

	
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
	


	g = sns.clustermap(df,method='average', xticklabels=xticks,yticklabels=yticks,cmap="Blues",robust=True, fmt="d", vmin=0, vmax=1, linewidths=.05, linecolor='grey', figsize=(width, height))


	pl.setp(g.ax_heatmap.get_xticklabels(), rotation=90) 
	pl.setp(g.ax_heatmap.get_yticklabels(), rotation=0) 
	g.fig.suptitle('genes')
	pl.subplots_adjust(right=float(right_scale),bottom=float(bottom_scale))
	ax = g.ax_heatmap
	ax.set_ylabel("")

	pl.savefig(tmpfiles+'clustermap_presence.tif')
	pl.savefig(resultspath+'clustermap_presence.pdf')

	
	pl.close()


	#**df_merge
	xticks = list(df.columns)
	yticks = listindex_merge

	#--- Construct clustermap out of the data, no estoy seguro si hacerlo por method average o no!!.#
	g = sns.clustermap(df_merge,method='average', xticklabels=xticks,yticklabels=yticks,cmap="Blues",robust=True, fmt="d", vmin=0, vmax=1, linewidths=.05, linecolor='grey', figsize=(width, height))
	pl.setp(g.ax_heatmap.get_xticklabels(), rotation=90) 
	pl.setp(g.ax_heatmap.get_yticklabels(), rotation=0) 
	g.fig.suptitle('genes and pseudogenes')
	pl.subplots_adjust(right=float(right_scale),bottom=float(bottom_scale))
	ax = g.ax_heatmap
	ax.set_ylabel("")

	pl.savefig(resultspath +'clustermap_presence_pseudogenes.pdf')
	
	pl.close()

	#**df_pseudo
	xticks = list(df.columns)
	yticks = listindex_pseudo

	#--- Construct clustermap out of the data, no estoy seguro si hacerlo por method average o no!!.#
	g = sns.clustermap(df_pseudo,method='average', xticklabels=xticks,yticklabels=yticks,cmap="Blues",robust=True, fmt="d", vmin=0, vmax=1, linewidths=.05, linecolor='grey', figsize=(width, height))

	pl.setp(g.ax_heatmap.get_xticklabels(), rotation=90) 
	pl.setp(g.ax_heatmap.get_yticklabels(), rotation=0) 
	g.fig.suptitle('Pseudogenes')
	pl.subplots_adjust(right=float(right_scale),bottom=float(bottom_scale))
	ax = g.ax_heatmap
	ax.set_ylabel("")


	pl.savefig(resultspath+'clustermap_pseudogenes.pdf')
	pl.close()

	

	novi = Toplevel()
	image = Image.open(tmpfiles +'clustermap_presence.tif')
	photo = ImageTk.PhotoImage(image)

		
		
	canvas = Canvas(novi, width = (width*100), height = (height*100),scrollregion = (0,0,(width*100),(height*100)))
	canvas.create_image(0, 0,image = photo, anchor="nw")
	canvas.photo = photo
		
	hbar=Scrollbar(novi,orient=HORIZONTAL)
	hbar.pack(side=BOTTOM,fill=X)
	hbar.config(command=canvas.xview)
	vbar=Scrollbar(novi,orient=VERTICAL)
	vbar.pack(side=RIGHT,fill=Y)
	vbar.config(command=canvas.yview)
	

	canvas.config(width=800,height=1200)
	canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
	canvas.pack(side=LEFT,expand=True,fill=BOTH)
	canvas.pack(expand = YES, fill = BOTH)
	
	return

def sequence_alignment(savepath,openpath,querypath, cov, ident):
	from Bio.Align.Applications import MuscleCommandline
	import subprocess

	cov = int(cov.get())
	ident = int(ident.get())	

	seq_resultspath = savepath+"/alignment/"
	if not os.path.exists(seq_resultspath):
		os.makedirs(seq_resultspath)


	resultsDir = savepath+"/output/"
	tmpfiles = savepath+"/tmp/"
	
	#dictionaries
	presence_dict = {}
	genome_file = {}
	sequences_dict={}

	#parameters for blast protein analysis retrieving, identity and coverage

	for i in os.listdir(resultsDir):
		if i.endswith(".txt"):
			fileName = resultsDir+i
			f = open(fileName,'r')	
			print f,ident,cov
			sequences = analyzer.parseResults(f,ident,cov,"retrieve_sequences")					
			f.close()
			tmp = os.path.splitext(i)[0]
			query = tmp.split('##')[0]
			subject = tmp.split('##')[1]
			genome_file[subject] = subject+".fasta"
			if subject not in presence_dict:
				presence_dict[subject] = {}
			presence_dict[subject][query] = sequences #[0]
			




	outputString = ""
	count = 0
	genomes =[]		
	# Generate the Sequence_dict containing the sequence for align clasified by protein instead of genome sequence.							
	BD_dict = np.load(savepath+'/BD_.dict.npy').item()
	for genome_seq in sorted(presence_dict,key=utils.natural_keys):
		if genome_seq in BD_dict:	
			genome = [genome_seq][0]
			value_dict = BD_dict[genome_seq].split("##")
			Strain =value_dict[0].strip("\n").replace(" ","_")
			assembly= "_".join(genome.split("_")[0:2])
			line = [genome_seq, Strain]		
			for protein in sorted(presence_dict[genome_seq],key=utils.natural_keys):		
				if protein not in sequences_dict:
					sequences_dict[protein]={}		
				if presence_dict[genome_seq][protein][0] > 0:						
					seq = []	
					head = Strain+"_"+assembly
					for n in range(len(presence_dict[genome_seq][protein][1])):
						sequence = presence_dict[genome_seq][protein][1][n]				
						seq.append(sequence.rstrip("\n"))
					sequences_dict[protein][head]=seq
			
	#for each protein, all genome sequence are extracted and write in a fasta file     
	for protein in 	sequences_dict:		
		filename = seq_resultspath+protein+".fasta"


		handle = querypath+"/"+protein+".fasta"
		Seq = SeqIO.read(handle,"fasta")
		query_sequence = Seq.seq
		query_name = Seq.name
		f = open(filename,'w')
		f.write(">"+query_name+"\n")	
		f.write(str(query_sequence)+"\n\n")		
		for genome_seq in sequences_dict[protein]:
			genome = genome_seq.replace("_genomic","")	
			for n in range(len(sequences_dict[protein][genome_seq])):
				count += 1
				seq = sequences_dict[protein][genome_seq][n]
				f.write(">"+genome+"_"+str(count)+"\n")	
				if seq =="NN":					#podemos modificar esto para controlar que salgan o no los nombresde los genomas que no contienen las proteinas query
					seq=""
				f.write(str(seq)+"\n\n")			
			count = 0	
		f.close()




	outputDir  = seq_resultspath+"muscle_results/"
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)

	outputDir2  = seq_resultspath+"mview_results/"
	if not os.path.exists(outputDir2):
		os.makedirs(outputDir2)

	outputDir3  = seq_resultspath+"ete_results/"
	if not os.path.exists(outputDir3):
		os.makedirs(outputDir3)
	




	""" align sequence from blast analysis in fasta 
	format(.fas), clustal format (.txt and .html) and generated .dat file
	for mview program """	
	displayedText.set('Alignment with muscle running!!')
	output_text.update_idletasks()
	window.update()	
	for file in os.listdir(seq_resultspath):
		if file.endswith(".fasta"):
			filename = os.path.splitext(file)[0]			
			res_fasta= outputDir+filename+"_muscle_align.fas"
			res_txt = outputDir+filename+"_muscle_align.txt"
			res_html = outputDir+filename+"_muscle_align.html"
			res_dat = tmpfiles+filename+"_muscle_align.dat"
			file2 = seq_resultspath+file
			print file2
			muscle_cline = MuscleCommandline(input=file2,fastaout=res_fasta, clwout=res_txt, htmlout=res_html)
			os.system(str(muscle_cline) +">/dev/null 2>&1") 				
			muscle_cline2 = MuscleCommandline(input=file2,fastaout=res_dat)		
			os.system(str(muscle_cline2)+">/dev/null 2>&1")	
			
	

		

	""" align sequence alignment dat files from muscle with mview """
	displayedText.set('mview alignment running!!')
	output_text.update_idletasks()
	window.update()	
	for file in os.listdir(tmpfiles):   
		if file.endswith(".dat"):
			filename = os.path.splitext(file)[0]			
			res_html = outputDir2+filename+"_mview_align.html"
			with open(res_html,"wb") as out, open("stderr.txt","wb") as err:
	  			subprocess.Popen(["mview","-in","fasta","-html","head","-css","on","-coloring","any","-threshold", "90",tmpfiles+file],stdout=out,stderr=err)					
			
			
	
	#try:
	#	subprocess.call(["ete3"])
	#except OSError as e:
	#	if e.errno == os.errno.ENOENT:
        # handle file not found error.
	#		tkMessageBox.showinfo("ORTHOPROK can not access to ete3!!")
	#		displayedText.set('ORTHOPROK can not access to ete3!!')	



	#for file in os.listdir(seq_resultspath):		
	#		if file.endswith(".fasta"):
	#			fileName = resultsDir+"/Seq_for_alignments/"+file
	#			query = os.path.splitext(file)[0]
	#			ResultFolder= outputDir3+query	
	#			subprocess.Popen(["ete3", "build", "-w", "standard_fasttree","--rename-dup-seqnames", "-a", fileName,"-o", ResultFolder+"/", "--clearall"])

			
	

	tkMessageBox.showinfo("ORTHOPROK","Alignment is done!!")


	displayedText.set('Sequence alignment is done!!')
	output_text.update_idletasks()
	window.update()	

##########################################################################
#
#                   Tkinter  Widgets and Labels
#
##########################################################################

#Openbutton = Button(window, text="Browse", command=Openfunc)
Openbutton = Button(window, text="Browse", command=lambda: _ask_directory(OpenName, paths_dict, OPENPATH))
Openbutton.pack()
Openbutton.place(y= 125, x=25)

grpahic_label = Label(window,text="/Genomes folder")
grpahic_label.pack()
grpahic_label.place(y=150,x=50)

pathName = Entry(window, textvariable=OpenName)
pathName.update()
pathName.focus_set()
pathName.pack( anchor='e')
pathName.place(y = 125, x = 110, width = 400, height = 25)
OpenName.set(DEFOPENPATH)
paths_dict[OPENPATH] = DEFOPENPATH

#Savebutton = Button(window, text="Browse", command=Savefunc)
Savebutton = Button(window, text="Browse", command=lambda: _ask_directory(SaveName, paths_dict, SAVEPATH))
Savebutton.pack()
Savebutton.place(y= 175, x=25) 

grpahic_label = Label(window,text="/Save folder")
grpahic_label.pack()
grpahic_label.place(y=200,x=50)

saveName = Entry(window, textvariable=SaveName)
saveName.update()
saveName.focus_set()
saveName.pack( anchor='e')
saveName.place(y = 175, x = 110, width = 400, height = 25)
SaveName.set(DEFSAVEPATH)
paths_dict[SAVEPATH] = DEFSAVEPATH

#querybutton = Button(window, text="Browse", command=Queryfunc)
querybutton = Button(window, text="Browse", command=lambda: _ask_directory(QueryName, paths_dict, QUERYPATH))
querybutton.pack()
querybutton.place(y= 225, x=25) 

grpahic_label = Label(window,text="/Queries folder")
grpahic_label.pack()
grpahic_label.place(y=250,x=50)


queryName = Entry(window, textvariable=QueryName)
queryName.update()
queryName.focus_set()
queryName.pack( anchor='e')
queryName.place(y = 225, x = 110, width = 400, height = 25)
QueryName.set(DEFQUERYPATH)
paths_dict[QUERYPATH] = DEFQUERYPATH

Blastbutton = Button(window, text="Blast", command= lambda:  Blast(paths_dict, int(redo.get())))	
Blastbutton.pack()
Blastbutton.place(y= 325, x=25)


redoblastbutton = Checkbutton(window, text="redo Blast", variable=redo)	
redoblastbutton.pack()
redoblastbutton.place(y= 325, x=145)
redo.set(0)



logbutton = Button(window, text="log", command= lambda:  iterateDirNs(paths_dict,
																	  genomeParser.UPDATELOGNO))
logbutton.pack()
logbutton.place(y= 280, x=25)

log_update = Button(window, text="log update", command= lambda:  iterateDirNs(paths_dict,
																			  genomeParser.UPDATELOGYES))
log_update.pack()
log_update.place(y= 280, x=150)




Analysisbutton = Button(window, text="Analysis", command= lambda:
	analysis(paths_dict,cov,ident,collapsed,absence,height,width))
Analysisbutton.pack()
Analysisbutton.place(y= 370, x=25)


coverage_description= Label(window, text="Coverage")
coverage_description.pack()
coverage_description.place(y=370,x=125)
coverage_entry = Entry(window, textvariable=cov)
coverage_entry.pack()
coverage_entry.place(y=370,x=195, width=50)
cov.set(51)

identity_description= Label(window, text="Identity")
identity_description.pack()
identity_description.place(y=370,x=250)
identity_entry = Entry(window, textvariable=ident)
identity_entry.pack()
identity_entry.place(y=370,x=310, width=50)
ident.set(90)


collapsebutton = Checkbutton(window, text="collapse genomes", variable=collapsed)	
collapsebutton.pack()
collapsebutton.place(y= 370, x=380)
collapsed.set(0)

absencebutton = Checkbutton(window, text="absence analysis", variable=absence)
absencebutton.pack()
absencebutton.place(y= 415, x=380)
absence.set(0)



grpahic_label = Label(window,text="Graphic analysis")
grpahic_label.pack()
grpahic_label.place(y=480,x=30)


#Figures parameters

width_description= Label(window, text="Width")
width_description.pack()
width_description.place(y=515,x=135)
width_entry = Entry(window, textvariable=width)
width_entry.pack()
width_entry.place(y=515,x=180, width=50)
width.set(6)


height_description= Label(window, text="Height")
height_description.pack()
height_description.place(y=515,x=240)
height_entry = Entry(window, textvariable=height)
height_entry.pack()
height_entry.place(y=515,x=290, width=50)
height.set(12)


Font_scale= Label(window, text="font scale")
Font_scale.pack()
Font_scale.place(y=515,x=350)
Font_scale_entry = Entry(window, textvariable=font_scale)
Font_scale_entry.pack()
Font_scale_entry.place(y=515,x=420, width=50)
font_scale.set(0.8)


Right_scale= Label(window, text="x_adj")
Right_scale.pack()
Right_scale.place(y=550,x=135)
Right_scale_entry = Entry(window, textvariable=right_scale)
Right_scale_entry.pack()
Right_scale_entry.place(y=550, x=180, width=50)
right_scale.set(0.5)


Bottom_scale= Label(window, text="y_adj")
Bottom_scale.pack()
Bottom_scale.place(y=550,x=240)
Bottom_scale_entry = Entry(window, textvariable=bottom_scale)
Bottom_scale_entry.pack()
Bottom_scale_entry.place(y=550, x=290, width=50)
bottom_scale.set(0.28)


heatmapbutton = Button(window, text="Heatmap", command= lambda:  heatmap(savepath,height,width,font_scale,right_scale,bottom_scale))
heatmapbutton.pack()
heatmapbutton.place(y= 515, x=25)

sizebutton = Button(window, text="protein plot", command= lambda:  protein_length_graphic(savepath,height,width,font_scale))
sizebutton.pack()
sizebutton.place(y= 630, x=25)


blast_graphic_button = Button(window, text="blast plot", command= lambda:  stripplot (savepath,height,width))
blast_graphic_button.pack()
blast_graphic_button.place(y= 555, x=25)



Alignbutton = Button(window, text="Alignment", command= lambda:  sequence_alignment(savepath,openpath,querypath,cov,ident))
Alignbutton.pack()
Alignbutton.place(y= 415, x=25)



# create the listbox for protein plot
protein_listbox1= Listbox(window, width = 25, height = 6)
protein_listbox1.pack()
protein_listbox1.place(y=600, x= 150)


# parece que no funciona, da error, cambiar por otro tipo de widget
# create a vertical scrollbar to the right of the protein listbox
#yscroll = Scrollbar(command=protein_listbox1.yview, orient=VERTICAL)
#yscroll.grid(row=0, column=1, sticky=N+S)
#protein_listbox1.configure(yscrollcommand=yscroll.set)

output_text = Tkinter.Label(window, textvariable=displayedText) # console information
output_text.pack()
output_text.configure(background="white",fg="black")
output_text.place(y= 750, x=125, width = 400, height = 35)

QuitButton = Button(window, text = "Quit", command=lambda: window.destroy())#lambda: QuitButtonCommand(window))
QuitButton.pack(anchor='center')
QuitButton.place(x = 25, y = 750, width = 80, height = 25)

window.configure(background='snow2')
window.mainloop()
