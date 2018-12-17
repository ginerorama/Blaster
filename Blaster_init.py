# -*- coding: utf-8 *-*

"""initialize the Blaster project"""	



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


window=Tkinter.Tk()
window.geometry("300x250")
window.title("Blaster")

OpenName = StringVar() # variable that stores target sequences path to show in GUI
SaveName = StringVar()  # variable that stores project directory to show in GUI
QueryName = StringVar() 
Pro_NAME = StringVar()	


global old_project # la uso para saber si tengo que crear un boton de load previous project en new_projectcommand
old_project = "no"




def create_project(project_Name,OpenName,queryName,saveName,fasta_genomes_count,fasta_query_count,log_file):	
	#escribe el archivo del proyecto con todos los parametros
	print project_Name,OpenName,queryName,saveName,fasta_genomes_count,fasta_query_count,log_file
	create_log_file(OpenName,saveName,log_file)
	project_file = open(project_Name,"w")
	genome_folder = "genome_folder: "+OpenName
	query_folder = "query_folder: "+queryName
	save_folder = "save_folder: "+saveName
	genomes_count = "genomes: "+str(fasta_genomes_count)
	query_count = "queries: "+str(fasta_query_count)
	log_file = "log_file: "+str(log_file)
	tab="\n"
	project_file.write(genome_folder+tab+query_folder+tab+save_folder+tab+genomes_count+tab+query_count+tab+log_file)
	
	return tkMessageBox.showinfo("INFO","Project has been created!!")




def check_if_project_exist(saveName):
	#chequea si existe un proyecto previo en la carpeta Save con datos e importa los datos del proyecto. Sino new project
	exist ="no"

	for file in os.listdir(saveName):
		if file.endswith(".pro"):
			project = open(saveName+"/"+file,"r")
			project_name = file
			for line in project:
				line = line.split(": ")
				if line[0] =="genome_folder":
					OpenName = line[1].rstrip("\n")
		

				if line[0] =="query_folder":
					queryName = line[1].rstrip("\n")
		

				if line[0] =="genomes":
					fasta_genomes_count = line[1].rstrip("\n")
	

				if line[0] =="queries":
					fasta_query_count = line[1].rstrip("\n")
		

			exist= "yes"		
	

	if exist == "yes":		
		tkMessageBox.showinfo("INFO","A previous project file exist in that folder!!")	
		try:
			return [exist, project_name,saveName, OpenName, queryName]
		except:
			return "no"

	else:
		return exist
	





def open_project():
	global old_project
	count = 0
	#open an load a previous project
	open_project_file= askopenfilename()
	if open_project_file.endswith(".pro"):
			project = open(open_project_file,"r")
			project_name = open_project_file
			for line in project:
				line = line.split(": ")
				if line[0] =="genome_folder":
					openpath = line[1].rstrip("\n")
				if line[0] =="query_folder":
					querypath = line[1].rstrip("\n")
				if line[0] =="save_folder":
					savepath = line[1].rstrip("\n")
				if line[0] =="genomes":
					fasta_genomes_count = line[1].rstrip("\n")
				if line[0] =="queries":
					fasta_query_count = line[1].rstrip("\n")
			
			#checking for log and dictionary log existance
			for file in os.listdir(savepath):		
				if file.endswith(".log"):
					create_log = "no"
				else:
					create_log = "yes"


			old_project = "yes"
			return old_project, New_projectCommand(project_name,querypath,savepath,openpath,"no_window",create_log) 
			

	else:
		tkMessageBox.showerror("ERROR","Please load a project.pro file!!")
	
	



def check_fasta(path_file):
	#check if folder contain fasta files
	test_value = ""
	for file in os.listdir(path_file):
		if file.endswith(".fasta"):
			test_value = "pass"
			break
		else:
			test_value = "not_pass"	

	return test_value	




def count_fasta_files(path_file):
	#count the number of fasta files in a folder
	count = 0

	for file in os.listdir(path_file):
		if file.endswith(".fasta"):
			count += 1
	return count		




def change_extension_fasta(path_file):
	#Change file extension to .fasta
	for file in os.listdir(path_file):
		if file.endswith(".fa") or file.endswith(".fna"):
			file_base = os.path.splitext(file)[0]
			os.rename(path_file+"/"+file, str(path_file+"/"+file_base) + ".fasta")




def Openfunc(): 
	# Open directory where target sequence for blast are stored.
	global openpath
	openpath = askdirectory() # variable that stores target sequences path
	print openpath
	fileopen = "/".join(openpath.split("/")[2::])
	OpenName.set(fileopen)




def Savefunc():
	# Directory in which the project will be stored.
	global savepath
	global old_project
	savepath = askdirectory() # variable that stores project directory
	filesave = "/".join(savepath.split("/")[2::])
	SaveName.set(filesave)
	previous_project = check_if_project_exist(savepath)
	if previous_project[0] == "yes":
		project_name = previous_project[1]
		savepath = previous_project[2]
		openpath = previous_project[3]
		querypath = previous_project[4]
		old_project = "yes"
		#checking for log and dictionary log existance
		for file in os.listdir(savepath):		
			if file.endswith(".log"):
				create_log = "no"
				break
		print project_name	
		New_projectCommand(project_name,querypath,savepath,openpath,novi,create_log)




def Queryfunc():
	# Directory in which the query sequences are stored.
	global querypath
	querypath = askdirectory() # variable that stores query sequences
	Querysave = "/".join(querypath.split("/")[2::])
	QueryName.set(Querysave)




def Projectfunc():
	#ask for project name and generate a file with .pro extension
	global project_Name
	project_Name = asksaveasfile(mode='w', defaultextension=".pro").name
	Pro_NAME.set(project_Name)




def project_window():
	""" ventana principal para seleccionar las variables del proyecto"""
	global novi
	#primero chekea que no exista un proyecto previo




	novi = Toplevel()
	canvas = Canvas(novi, width = 550, height = 400,background='snow2') # modify this value depending on number of variables
	canvas.pack(expand = YES, fill = BOTH)
	novi.title('OrthoProk: Create a new project')


	
	Open_label = Label(novi, text="Target genomes sequences folder:")
	Open_label.pack()
	Open_label.place(y= 25, x=25)
	
	Openbutton = Button(novi, text="Browse", command=Openfunc)	
	Openbutton.pack()
	Openbutton.place(y= 50, x=25)

	pathName = Entry(novi, textvariable=OpenName)
	pathName.update()
	pathName.focus_set()
	pathName.pack( anchor='e')
	pathName.place(y = 50, x = 110, width = 400, height = 25)




	Query_label = Label(novi, text="Query sequences folder:")
	Query_label.pack()
	Query_label.place(y= 75, x=25)

	querybutton = Button(novi, text="Browse", command=Queryfunc)	
	querybutton.pack()
	querybutton.place(y= 100, x=25) 

	queryName = Entry(novi, textvariable=QueryName)
	queryName.update()
	queryName.focus_set()
	queryName.pack( anchor='e')
	queryName.place(y = 100, x = 110, width = 400, height = 25)
	



	Output_label = Label(novi, text="save folder:")
	Output_label.pack()
	Output_label.place(y= 125, x=25)

	Savebutton = Button(novi, text="Browse", command=Savefunc)	
	Savebutton.pack()
	Savebutton.place(y= 150, x=25) 

	saveName = Entry(novi, textvariable=SaveName)
	saveName.update()
	saveName.focus_set()
	saveName.pack( anchor='e')
	saveName.place(y = 150, x = 110, width = 400, height = 25)


	

	projectName_label = Label(novi, text="Project Name:")
	projectName_label.pack()
	projectName_label.place(y= 250, x=25)

	projectNamebutton = Button(novi, text="Browse", command=Projectfunc)	
	projectNamebutton.pack()
	projectNamebutton.place(y= 275, x=25) 


	projectName = Entry(novi, textvariable=Pro_NAME)
	projectName.update()
	projectName.focus_set()
	projectName.pack( anchor='e')
	projectName.place(y = 275, x = 110, width = 400, height = 25)



	next_button = Button(novi, text="Next >", command= lambda: New_projectCommand(project_Name,querypath,savepath,openpath,novi,"yes"))
	next_button.pack()
	next_button.place(x=450,y=350) 


	cancel_button = Button(novi, text="Cancel", command = lambda: QuitButtonCommand(novi))
	cancel_button.pack()
	cancel_button.place(x=350,y=350)





def New_projectCommand(ProjectName,QueryName,SaveName,Open_Name,window,create_log):
	global old_project

	#close nove window	
	if window =="no_window":
		pass
	else:
		window.destroy()
	

	#getting variables values
	ProjectName
	queryName = QueryName
	saveName = SaveName
	OpenName= Open_Name


	#create a new window for summary project view
	summary_window = Toplevel()
	canvas = Canvas(summary_window, width = 550, height = 400,background='snow2') # modify this value depending on number of variables
	canvas.pack(expand = YES, fill = BOTH)
	summary_window.title('Project summary')


	#change files extension to .fasta in either genome and query folders
	change_extension_fasta(OpenName)
	change_extension_fasta(queryName)

	

	#Count fasta files in either genome and query folders
	fasta_genomes_count = count_fasta_files(OpenName)
	fasta_query_count = count_fasta_files(queryName)


	if fasta_genomes_count == "0":
		fasta_genomes_count_color = "red"
	else:
		fasta_genomes_count_color = "blue"	

	if fasta_query_count == "0":
		fasta_query_count_color = "red"
	else:
		fasta_query_count_color = "blue"

	#check if query and open folder exist fasta files
	if check_fasta(queryName) == "not_pass":
		queryName = "Error: this folder does not contain Fasta files"
		query_color ="red"
	else:
		query_color ="black"

	if check_fasta(OpenName) == "not_pass":
		OpenName = "Error: this folder does not contain Fasta files"
		genome_color="red"
	else:
		genome_color="black"

	
	#Create log file
	if create_log =="yes":
		log_file = ProjectName.replace(".pro","")+"_target_sequences.log"
		create_log_file(OpenName,saveName,log_file)
	else:
		pass




	#######
	# Tkinter GUI paramaters
	
	project_name = Label(summary_window, text="Project:" )
	project_name.pack()
	project_name.place(y= 55, x=25)
	project_name = Label(summary_window, text=ProjectName )
	project_name.pack()
	project_name.place(y= 55, x=145)

	genome_folder = Label(summary_window, text="Loaded genomes:" )
	genome_folder.pack()
	genome_folder.place(y= 75, x=25)
	genome_folder = Label(summary_window, text=OpenName, foreground = genome_color)
	genome_folder.pack()
	genome_folder.place(y= 75, x=145)

	protein_folder = Label(summary_window, text="Loaded proteins:" )
	protein_folder.pack()
	protein_folder.place(y= 95, x=25)
	protein_folder = Label(summary_window, text=queryName , foreground = query_color)
	protein_folder.pack()
	protein_folder.place(y= 95, x=145)

	save_folder = Label(summary_window, text="Save folder:" )
	save_folder.pack()
	save_folder.place(y= 115, x=25)
	save_folder = Label(summary_window, text=saveName )
	save_folder.pack()
	save_folder.place(y= 115, x=145)



	genome_count = Label(summary_window, text="Number of Genomes: " )
	genome_count.pack()
	genome_count.place(y= 155, x=25)
	genome_count = Label(summary_window, text=fasta_genomes_count , foreground = fasta_genomes_count_color)
	genome_count.pack()
	genome_count.place(y= 155, x=165)

	query_count = Label(summary_window, text="Number of proteins: " )
	query_count.pack()
	query_count.place(y= 175, x=25)
	query_count = Label(summary_window, text=fasta_query_count , foreground = fasta_query_count_color)
	query_count.pack()
	query_count.place(y= 175, x=165)


	
	if check_fasta(QueryName)== "pass" and check_fasta(Open_Name) == "pass" and old_project == "no":

		next_button = Button(summary_window, text="Create project >", command= lambda: create_project(project_Name,OpenName,queryName,saveName,fasta_genomes_count,fasta_query_count,log_file))
		next_button.pack()
		next_button.place(x=395,y=350) 
	else:
		load_button = Button(summary_window, text="load project >", command= lambda: load_project())
		load_button.pack()
		load_button.place(x=395,y=350) 


		
	back_button = Button(summary_window, text="< back", command= lambda: project_window_destroy(summary_window))
	back_button.pack()
	back_button.place(x=300,y=350) 





def create_log_file(openpath,savepath,log_file):
	"""
	Generate the log file input for blaster v2.0
	automatic detection of .gz compressed genomes from genbank an uncompress
	log file also contain in the 3th column from fasta header file
	Arguments: /fasta_files_folder log_basename_file

	"""
	outputDir = savepath+"/"
	subjectDir = openpath+"/"
	
	



	output = open(log_file,"wr")
	for i in os.listdir(subjectDir):
		if i !=".DS_Store":
			path = os.path.abspath(subjectDir+i)
			name = os.path.splitext(i)[0]
			extension = os.path.splitext(i)[1]
		
			#check for gz files from ncbi
			if extension == ".gz":
				import gzip
				print "comprised gz file found\nuncompressing...(This could take some time!!"
				file = gzip.open(path,"rb")
				file = file.read()
				uncompressed_file = open(subjectDir+name,"w")
				uncompressed_file.write(str(file))
			else:	
				file = (open(path,"r")).readline()

			#read the files and check for fasta files
			if ">" not in file[0]:
				displayedText.set("ERROR: " + name+ " is not a fasta file")
				output_text.update_idletasks()
				window.update()



			else:
				fasta_header = file.split("\n")[0][1::]
				name_fasta= (" ".join(fasta_header.split(" ")[1::]).split(",")[0].split(".")[0].split("_")[0])
				name_fasta= name_fasta.replace("circular","").replace("chromosome","").replace("complete","").replace("genome","").replace("sequence","").replace("assembly","")
			#	output.write(str(name)+"\t"+str(path)+"\t"+str(fasta_header)+"\n")		
				output.write(str(name_fasta)+"\t"+str(path)+"\t"+str(name)+"\n")		
			

	output.close()




	
def create_log_dict(savepath):
	"""
	Code to create target sequences dictionary
	"""	

	outputDir = savepath+"/"



	selected_genomes_dir = outputDir+"target_sequences.log" 
	selected_genomes= open(selected_genomes_dir, 'r')

	BD_dict= {}
	for genome in selected_genomes:
		ID = genome.split("\t")[2].strip("\n").replace(",","") # nombre del genoma sin la extension, que sera usado en el diccionario
		name= genome.split("\t")[0]	
		#name= name.replace("circular","").replace("chromosome","").replace("complete","").replace("genome","").replace("sequence","").replace("assembly","")
		genome = genome.split("\t")[1].strip("\n").replace(".fna","").replace(".fasta","").replace(".fa","").replace("faa","").replace("frn","").replace("ffn","")
		BD_dict[ID]=name+"##"+genome



	np.save(outputDir+'BD_.dict.npy', BD_dict) 

	#if os.path.exists(outputDir+'BD_.dict.npy'):
	#	tkMessageBox.showinfo("INFO","target sequences database done!!")
	#	output_text.update_idletasks()
	#	window.update()		
	#else:
	#	displayedText.set("Error in target sequences database!!")
	#	tkMessageBox.showinfo("ERROR","Error in target sequences database!!")
	#	output_text.update_idletasks()
	#	window.update()	





def project_window_destroy(window):
   
    window.destroy()
    project_window()




def QuitButtonCommand(window):
	window.destroy()












project_buttom = Button(window, text="New project", command = lambda: project_window())
project_buttom.pack()
project_buttom.place(x=25,y=100)

load_project_buttom = Button(window, text="Load project", command = lambda: open_project())
load_project_buttom.pack()
load_project_buttom.place(x=25,y=140)

QuitButton = Button(window, text = "Quit", command=lambda: QuitButtonCommand(window))
QuitButton.pack(anchor='center')
QuitButton.place(x = 25, y = 180, width = 80, height = 25)



window.configure(background='snow2')
window.mainloop()

