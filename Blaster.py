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

import Tkinter, tkFileDialog, Tkconstants 
from tkFileDialog import *
from tkFileDialog import askopenfilename
from Tkinter import *
from ttk import *
import tkMessageBox
import os, sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from PIL import Image, ImageTk

import subprocess

import src.Facade as facade
import src.GenomeParser as genomeParser
import src.Analyzer as analyzer
import src.Plotter as plotter
import src.utils as utils

##################
# Global variables

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

def _show_image(imagepath, height, width, hcanvas, wcanvas):
	
	novi = Toplevel()
	image = Image.open(imagepath)
	photo = ImageTk.PhotoImage(image)
	
	canvas = Canvas(novi,
					width = (width*100), height = (height*100),
					scrollregion = (0,0,(width*100),(height*100)))
	canvas.create_image(0, 0,image = photo, anchor="nw")
	canvas.photo = photo
	
	hbar=Scrollbar(novi,orient=HORIZONTAL)
	hbar.pack(side=BOTTOM,fill=X)
	hbar.config(command=canvas.xview)
	vbar=Scrollbar(novi,orient=VERTICAL)
	vbar.pack(side=RIGHT,fill=Y)
	vbar.config(command=canvas.yview)
	
	canvas.config(width=wcanvas,height=hcanvas)
	canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
	canvas.pack(side=LEFT,expand=True,fill=BOTH)
	canvas.pack(expand = YES, fill = BOTH)
	
	return

##########################################################################
#
#                   Commands
#
##########################################################################

def iterateDirNs(paths_dict, UPDATE):
	"""
	Generate the log file input for blaster v2.0
	automatic detection of .gz compressed genomes from genbank an uncompress
	log file also contain in the 3th column from fasta header file
	Arguments: /fasta_files_folder log_basename_file

	"""
	
	savepath = paths_dict[SAVEPATH]
	openpath = paths_dict[OPENPATH]
	
	try:
		dict_file = facade.updatelog(savepath, openpath, UPDATE, sys.stderr) # CPC2018
		
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
	outputDir = analyzer.get_outputfolder(savepath)
	
	#creates output directory in project directory
	# newpath = savepath+"/output"
	# if not os.path.exists(newpath):
	# 	os.makedirs(newpath)
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
	
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

def sequence_alignment(paths_dict, cov, ident):

	cov = int(cov.get())
	ident = int(ident.get())
	
	savepath = paths_dict[SAVEPATH]
	querypath = paths_dict[QUERYPATH]
	
	facade.sequence_alignment(savepath, querypath, cov, ident, _display, sys.stderr)
	
	tkMessageBox.showinfo("ORTHOPROK","Alignment is done!!")
	
	_display('Sequence alignment is done')
	
	return

def protein_length_graphic(paths_dict, height, width, font_scale):
	
	height = int(height.get())
	width = int(width.get())
	font_scale = float(font_scale.get())
	
	savepath = paths_dict[SAVEPATH]
	
	prot_index = 0
	if len(protein_listbox1.curselection()) > 0:
		prot_index = protein_listbox1.curselection()[0]
	else:
		raise Exception("At least one protein must be selected")
	
	prot_name = protein_listbox1.get(prot_index)
	
	facade.protein_length_graphic(prot_name, savepath, height, width, font_scale, sys.stderr)
	
	protein_size_tif = plotter.get_protein_size_tif(savepath, prot_name)
	
	_show_image(protein_size_tif, height, width, 1200, 800)
	
	displayedText.set('Protein plot done!!')
	
	return

def stripplot (paths_dict, height, width):
	
	height = int(height.get())
	width = int(width.get())
	
	savepath = paths_dict[SAVEPATH]
	
	facade.stripplot(savepath, height, width, sys.stderr)
	
	blast_plot_tif = plotter.get_blast_plot_tif(savepath)
	
	_show_image(blast_plot_tif, height, width, 1200, 600)
	
	return

def heatmap(paths_dict, height, width, font_scale, right_scale, bottom_scale):
	
	height = int(height.get())
	width = int(width.get())
	font_scale = float(font_scale.get())
	bottom_scale = float(bottom_scale.get())
	right_scale = float(right_scale.get())
	
	savepath = paths_dict[SAVEPATH]
	
	facade.heatmap(savepath, height, width, font_scale, right_scale, bottom_scale, sys.stderr)
	
	clustermap = plotter.get_cluster_map_tif(savepath)
	
	_show_image(clustermap, height, width, 1200, 800)
	
	return

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

heatmapbutton = Button(window, text="Heatmap", command= lambda:
	heatmap(paths_dict, height, width, font_scale, right_scale, bottom_scale))

heatmapbutton.pack()
heatmapbutton.place(y= 515, x=25)

sizebutton = Button(window, text="protein plot", command= lambda:
	protein_length_graphic(paths_dict ,height ,width, font_scale))

sizebutton.pack()
sizebutton.place(y= 630, x=25)


blast_graphic_button = Button(window, text="blast plot", command= lambda:
	stripplot(paths_dict, height, width))

blast_graphic_button.pack()
blast_graphic_button.place(y= 555, x=25)

Alignbutton = Button(window, text="Alignment", command= lambda:
	sequence_alignment(paths_dict, cov, ident))

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

## END