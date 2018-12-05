

#!/usr/bin/python
# primergenerator v0.2 
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
#                   Global variables
#
##########################################################################


current_directory = os.path.dirname(os.path.abspath(__file__))

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
query_list = []
target_list = []



##########################################################################
#
#                   Utility functions
#
##########################################################################



def Openfunc(): 
	# Open directory where target sequence for blast are stored.
	global openpath
	openpath = askdirectory() # variable that stores target sequences path
	fileopen = "/".join(openpath.split("/")[2::])
	OpenName.set(fileopen)

def Savefunc():
	# Directory in which the project will be stored.
	global savepath
	savepath = askdirectory() # variable that stores project directory
	filesave = "/".join(savepath.split("/")[2::])
	SaveName.set(filesave)

def Queryfunc():
	# Directory in which the query sequences are stored.
	global querypath
	querypath = askdirectory() # variable that stores query sequences
	Querysave = "/".join(querypath.split("/")[2::])
	QueryName.set(Querysave)


def change_extension_fasta(path):

	for file in os.listdir(path):
		if file.endswith(".fa") or file.endswith(".fna") :
			file_base = os.path.splitext(file)[0]
			os.rename(path_file+"/"+file, str(path_file+"/"+file_base) + ".fasta")





def atoi(text):
	return int(text) if text.isdigit() else text

def natural_keys(text):
	'''
	alist.sort(key=natural_keys) sorts in human order
	http://nedbatchelder.com/blog/200712/human_sorting.html
	(See Toothy's implementation in the comments)
	'''
	return [ atoi(c) for c in re.split('(\d+)', text) ]



def parseResults(file,cov,ident,function): #statsA
	function = function
	cov = int(cov)
	ident = int(ident)
	numPass = 0
	line_number=0
	nCount=0
	protein_size = 0
	iden = 0
	if function == "retrieve_sequences":  	
		seq = []
		for line in file:
			if line_number >= 1:  
				fields = line.split("\t")
				qcovhsp = int(fields[13])
				pident = int(round(float(fields[2])))
				sequence = fields[15]
				if int(qcovhsp) > cov and int(pident) > ident: # ajustar los stats al menos 85% #
					numPass += 1			
					if "*" in sequence: #take only the sequence before codon stop	
						sequence = sequence.split("*")[0]
						seq.append(sequence)
					else:
						seq.append(sequence)
			line_number = line_number +1		
		if numPass == 0:
			numPass =1
			seq.append("NN")	
		return [numPass,seq] 


	if function == "Analysis":  	
		seq = []
		coverage = []
		identity= []			
		for line in file:
			if line_number >= 1:        #this is to remove columns legends
				fields = line.split("\t")
				qcovhsp = int(fields[13])
				pident = int(round(float(fields[2])))
				sequence = fields[15]
				if int(qcovhsp) > cov and int(pident) > ident: 
					
					if	'*' in sequence:     # si el stop codon no se encuentra en el 20% final de la proteina es un pseudogen
						position = int(sequence.find('*'))
						
						if int(round(position/len(sequence)))*100 < 90: 
							nCount = sequence.count("*") 

						protein_size = int(round(position/len(sequence)))*100 
					
					else:				
						protein_size = qcovhsp	
					
					iden = pident
					numPass += 1			
					seq.append(sequence)
					coverage.append(qcovhsp)
					identity.append(pident)
			line_number = line_number +1		
		
		return [numPass,coverage,identity,seq,nCount,protein_size,iden]	





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





def analysis(savepath,openpath, cov, ident, collapsed, absence,height,width): #statsA + statisticsA
	cov = int(cov.get())
	ident = int(ident.get())
	collapsed = int(collapsed.get())
	absence = int(absence.get())
	height =int(height.get())
	width = int(width.get())


	resultspath = savepath+"/results/"
	if not os.path.exists(resultspath):
		os.makedirs(resultspath)

	blast_folder = resultspath+"/blast/"
	if not os.path.exists(blast_folder):
		os.makedirs(blast_folder)	

	tmpfiles = savepath+"/tmp/"
	if not os.path.exists(tmpfiles):
		os.makedirs(tmpfiles)


	outputDir = savepath+"/output/"
	outputFile = resultspath+"/Analysis_presence_summary.csv"
	outputFile2 = resultspath+"/Analysis_pseudo_summary.csv"
	outputFile3 = tmpfiles+"/HitsPresence.csv"
	outputFile4 = tmpfiles+"/HitsPseudo.csv"
	outputFile5 = resultspath+"/Proteins_size.csv"
	outputFile6 = resultspath+"/Proteins_ident.csv"

	resHash = {}

	for i in os.listdir(outputDir):
		if i.endswith(".txt"):
			fileName = outputDir+i
			f = open(fileName,'r')
			numCov = parseResults(f,cov,ident,"Analysis")					
			f.close()
			tmp = os.path.splitext(i)[0]
			query = tmp.split('##')[0]
			subject = tmp.split('##')[1]
			if subject not in resHash:
				resHash[subject] = {}
			resHash[subject][query] = numCov #[0]
			


	outputString_presence = ""
	outputString_pseudo = ""
	outputString_presence_tmp = ""
	outputString_pseudo_tmp = ""
	outputString_psize= ""
	outputString_ident= ""
	genomes =[]		
	que_list=[]

	#print resHash
	BD_dict = np.load(savepath+'/BD_.dict.npy').item()

	first = True
	head = ["Assembly,Name"]
	head_tmp =[]
	count=1

	for assembly in sorted(resHash,key=natural_keys):
		if assembly in BD_dict:	
			genome = [assembly][0]
			value_dict = BD_dict[assembly].split("##")
			Strain = value_dict[0].strip("\n")
			line_presence = [assembly, Strain]
			line_pseudo = [assembly, Strain]
			line_presence2 = []
			line_pseudo2 = []
			line_psize = [assembly, Strain]
			line_ident = [assembly, Strain]

			for query in sorted(resHash[assembly],key=natural_keys):		
				if query not in que_list:
					que_list.append(query)
				if first:
					head.append(query)
					head_tmp.append(query)
				line_presence.append(str(resHash[assembly][query][0]))
				line_pseudo.append(str(resHash[assembly][query][4])) 				
				line_presence2.append(str(resHash[assembly][query][0]))
				line_pseudo2.append(str(resHash[assembly][query][4]))	
				line_psize.append(str(resHash[assembly][query][5]))	
				line_ident.append(str(resHash[assembly][query][6]))	

			if first:
				first = False
				outputString_presence+=','.join(head)+"\n"
				outputString_pseudo+=','.join(head)+"\n"
				outputString_presence_tmp+=','.join(head_tmp)+","+"Serovar"+","+"Genome_file"+","+"Count"+"\n"
				outputString_pseudo_tmp+=','.join(head_tmp)+","+"Serovar"+","+"Genome_file"+","+"Count"+"\n"
				outputString_psize +=','.join(head)+"\n"
				outputString_ident +=','.join(head)+"\n"

			value_dict = BD_dict[assembly].split("##")
			Strain = value_dict[0].strip("\n")
			Assembly = value_dict[1].strip("\n")
			outputString_presence+=','.join(line_presence)+"\n"
			outputString_pseudo+=','.join(line_pseudo)+"\n"
			outputString_presence_tmp+=','.join(line_presence2)+","+Strain+","+Assembly+","+str(count)+"\n"
			outputString_pseudo_tmp+=','.join(line_pseudo2)+","+Strain+","+Assembly+","+str(count)+"\n"
			outputString_psize+=','.join(line_psize)+"\n"
			outputString_ident+=','.join(line_ident)+"\n"


	##### output files ######

	f = open(outputFile,'w')
	f.write(outputString_presence)
	f.close()			
	f2 = open(outputFile2,'w')
	f2.write(outputString_pseudo)
	f2.close()		
	f3 = open(outputFile3,'w')
	f3.write(outputString_presence_tmp)
	f3.close()			
	f4 = open(outputFile4,'w')
	f4.write(outputString_pseudo_tmp)
	f4.close()	
	f5 = open(outputFile5,'w')
	f5.write(outputString_psize)
	f5.close()
	f6 = open(outputFile6,'w')
	f6.write(outputString_ident)
	f6.close()

	
	for que in que_list:
		outputFile2 = blast_folder+"/"+que+".csv"
		f2 = open(outputFile2,'w')

		f2.write(
			str(que.upper())+
			","+"Coverage"+
			","+"Identity"+
			","+"Sequence"+"\n"
			)
		
		for sub in sorted(resHash,key=natural_keys):
			if sub in BD_dict:	
				genome = [sub][0]
				value_dict = BD_dict[sub].split("##")
				Strain = value_dict[0].strip("\n")
				if resHash[sub][que][0] > 0:						
					valor = 0							
					f2.write(str(Strain)+"\n")
					for n in range(len(resHash[sub][que][1])):
						valor +=1
						seq = resHash[sub][que][3][n]
						cov = resHash[sub][que][1][n]
						iden= resHash[sub][que][2][n]
						f2.write(
							str(int(valor))+
							","+str(cov)+
							","+str(iden)+
							","+ str(seq)
							)
						f2.write("\n")			
					f2.write("\n\n")
		f2.close()			



	#####################################
	#####statisctisA

	df1 = pd.read_csv(outputFile3) #'HitsPresence.csv
	df2 = pd.read_csv(outputFile4) #'HitsPseudo.csv

	df1_query_number = len(list(df1.columns.values))-3
	df2_query_number = len(list(df2.columns.values))-3

	DBlength = len(df1.index) #calculate number of rows in the dataframe

	absence_protein = df1.columns.values[0:df1_query_number] # to retrieve column names of absence proteins
	stop_codon_pseudo = df2.columns.values[0:df2_query_number] # to retrieve column names of stop_codons gene proteins
	columns_names = ['absence','STOP_codon']

	# creamos df_result donde guardamos los resultados de analizar las tablas donde guardar los datos
	df_result = pd.DataFrame(np.random.randint(0,100,size=(2,len(absence_protein))),columns=absence_protein , index=columns_names)


	# modify df1 to conver values greater than 1 in 1, and 0 in 0 for absence_protein
	for protein in absence_protein:
		count = 0
		for n in df1[protein]:
			if n == 0:
				df1.ix[count,protein]=0
			if n > 0:
				df1.ix[count,protein]=1
			df1 = df1
			count += 1	

	for protein in absence_protein:
		total_presence = df1[protein].sum()
		df_result.ix[0,protein]=total_presence

	### esta tabla se usara en csv_groupby para generar tabla de ausencia collapse
	df1.to_csv(tmpfiles+'absence_for_collapse.csv')


	# modify df1 to conver values greater than 1 in 1, and 0 in 0 for stop_codon_pseudo
	for protein in stop_codon_pseudo:
		count = 0
		for n in df2[protein]:
			if n == 0:
				df2.ix[count,protein]=0
			if n > 0:
				df2.ix[count,protein]=1
			df2 = df2
			count += 1	

	for protein in stop_codon_pseudo:
		total_pseudo = df2[protein].sum()
		df_result.ix[1,protein]=total_pseudo
	

	### esta tabla se usara en csv_groupby para generar tabla de pseudogenes for collapse
	df2.to_csv(tmpfiles+'pseudo_for_collapse.csv')

	## save the df_result
	df_result.to_csv(resultspath+'summarystats.csv')	



	## Graphication phase
	import matplotlib.pyplot as plt

	#volvemos a cargar el archivo para graficar porque sino no se ajustan las coordenadas a la hora de graficar
	df_result = pd.read_csv(resultspath+'summarystats.csv')


	fig1 = plt.figure(figsize = (9,7))
	plt.subplots_adjust(hspace=1.4)
	p1 = plt.subplot(2,1,1)
	l1 = plt.bar(range(len(df_result.ix[0,1::])),df_result.ix[0,1::])
	ly = plt.ylabel('Numbers')
	lx = plt.xticks(range(len(df_result.ix[0,1::])), absence_protein, size='small', rotation='vertical')
	if absence == 1:
		ttl = plt.title("Absence genes/proteins")
	else:
		ttl = plt.title("Proteins founds")

	p4 = plt.subplot(2,1,2)
	l4 = plt.bar(range(len(df_result.ix[1,1::])),df_result.ix[1,1::], color='m')
	ly = plt.ylabel('Numbers')
	lx = plt.xticks(range(len(df_result.ix[1,1::]+1)), absence_protein, size='small', rotation='vertical')
	ttl = plt.title("Pseudogenes(STOP)")
	if absence == 1:
		ttl = plt.suptitle("Absent proteins and pseudogens distribution in DataBase")
	else:
		ttl = plt.suptitle("Proteins presents and pseudogens distribution in DataBase")


	
	figure = resultspath+"graphic_stats.pdf"
	figure_2 = tmpfiles+"graphic_stats.tif"
	if not os.path.exists(figure):
		fig1.savefig(resultspath+'graphic_stats.pdf')
		fig1.savefig(tmpfiles+"graphic_stats.tif")
		
		plt.close()
	else:
		os.remove(figure)
		fig1.savefig(resultspath+'graphic_stats.pdf')
		fig1.savefig(tmpfiles+"graphic_stats.tif")
		
		plt.close()
	


	#####################
	#csvgroupby

	df = pd.read_csv(tmpfiles+'absence_for_collapse.csv', index_col = 0)
	df_pseudo = pd.read_csv(tmpfiles+'pseudo_for_collapse.csv', index_col = 0)




	#### suma los serovares apilando la tabla		
	if collapsed ==1:
		df = df.groupby(['Serovar']).sum()
		df_pseudo = df_pseudo.groupby(['Serovar']).sum()

# esto es para calcular las ausencias en vez de las presencias
	if absence == 1:
		for n in list(df.index):
			for x in (df.columns):
				if x != 'Count':
					value = df.ix[n,x]
					inverse_value = df.ix[n,'Count']- value
					df.ix[n,x] = inverse_value
					df= df

	### ahora sumamos ambos dataframes el de ausencia y el de pseudogenes

	df_merge = df.add(df_pseudo, fill_value=0)
	count = 0
	for n in df_merge['Count']:
		new_n= int(n/2)
		df_merge.ix[count,'Count']= new_n
		df_merge = df_merge
		count +=1
		



	#df_merge = df_merge.drop('Count', 1) # elimina la columna count


	df_merge = df_merge.reset_index()
	df_merge['Count']= df_merge['Count'].astype(str)
	df_merge['Serovar'] = df_merge[['Serovar', 'Count']].apply(lambda x: '_'.join(x), axis=1)
	df_merge = df_merge.set_index('Serovar')
	df_merge= df_merge.drop('Count', 1) 
	df_merge.to_csv(tmpfiles+'groupbyserovar_1normalized_merge.csv')

	

	#normaliza los valores de 1-0 
	if collapsed == 1:
		for protein in list(df.columns)[1:-2]:
			if protein != "Count":
				count = 0
				for n in df[protein]:
					if n == 0:
						df.ix[count,protein]=0
						df = df
						count += 1	
					else:
						df.ix[count,protein]=(float(n)/float(df.ix[count,'Count']))# +1
						df = df
						count += 1	

	
	df = df.reset_index()
	if collapsed == 1:
		df['Count']= df['Count'].astype(str)
		df['Serovar'] = df[['Serovar', 'Count']].apply(lambda x: '_'.join(x), axis=1)
	
	df = df.set_index('Serovar')
	df= df.drop('Count', 1) 
	df.to_csv(tmpfiles+ 'groupbyserovar_1normalized_absence.csv')



	df_pseudo = df_pseudo.reset_index()
	df_pseudo['Count']= df_pseudo['Count'].astype(str)
	df_pseudo['Serovar'] = df_pseudo[['Serovar', 'Count']].apply(lambda x: '_'.join(x), axis=1)
	df_pseudo = df_pseudo.set_index('Serovar')
	df_pseudo= df_pseudo.drop('Count', 1) 
	df_pseudo.to_csv(tmpfiles+'groupbyserovar_1normalized_pseudo.csv')

	tkMessageBox.showinfo("ORTHOPROK","analysis is done!!")

	if os.path.exists(outputFile) and os.path.exists(outputFile2) :
		displayedText.set("analysis is done!!")
		output_text.update_idletasks()
		window.update()		
	else:
		displayedText.set("Error: presence analysis!!")
		output_text.update_idletasks()
		window.update()	




# load the protein listbox with data
	df_protein = pd.read_csv(outputFile5)

	protein_list= df_protein.ix[:,2::].columns.tolist()

	for item in protein_list:
		protein_listbox1.insert(END,item)

	novi = Toplevel()
	image = Image.open(tmpfiles +'graphic_stats.tif') 
	photo = ImageTk.PhotoImage(image)
	#image not visual
	

	canvas = Canvas(novi, width = (width*100), height = (height*100),scrollregion = (0,0,(width*100),(height*100)))
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






def Blast (openpath,savepath,querypath, redo):
	from Bio.Blast.Applications import NcbiblastnCommandline
	import subprocess
	global query_list
	global target_list
	redo = int(redo.get())


	change_extension_fasta(openpath)
	change_extension_fasta(querypath)

	queryDir = querypath+"/"
	subjectDir = openpath+"/"
	outputDir = savepath+"/output/"




	#creates output directory in project directory
	newpath = savepath+"/output"
	if not os.path.exists(newpath):
		os.makedirs(newpath)
	
	count =0
	totalfiles = (len([name for name in os.listdir(queryDir) if os.path.isfile(os.path.join(queryDir, name))]))-1
	
	for i in  os.listdir(queryDir):
		if i != ".DS_Store":
			count += 1	
			displayedText.set("Running blast: "+str(count)+" of "+str(totalfiles)+" sequences")
			output_text.update_idletasks()
			window.update()
			queryFile = queryDir+i
			try: 		
				handle = open(queryFile,"rU")
				record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
				handle.close()
				queryFile = '"'+queryDir+i+'"'
				queryName = os.path.splitext(i)[0]# cuidado si interaccionas con los ficheros en finder y aparece DS_store aqui fallara hayq eliminar los DS_store files antes
				for j in os.listdir(subjectDir):
					if j.endswith(".fna") or j.endswith(".fasta") or j.endswith(".fa"):	
						if redo == 1:
							if i not in query_list and j not in target_list:
								subjectFile = '"'+subjectDir+j+'"'
								blastn_cline = NcbiblastnCommandline(cmd='tblastn',query=queryFile,subject=subjectFile,outfmt='"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qseq sseq"',evalue='0.00001',num_threads='4')
								stdout, stderr = blastn_cline()
								outFileName = outputDir+queryName+"##"+os.path.splitext(j)[0]+".txt"
								f = open(outFileName,'w')
								f.write("qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qcovs	qcovhsp	qseq	sseq\n"+stdout)
								f.close()
						if redo == 0:
							subjectFile = '"'+subjectDir+j+'"'
							blastn_cline = NcbiblastnCommandline(cmd='tblastn',query=queryFile,subject=subjectFile,outfmt='"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qseq sseq"',evalue='0.00001',num_threads='4')
							stdout, stderr = blastn_cline()
							outFileName = outputDir+queryName+"##"+os.path.splitext(j)[0]+".txt"
							f = open(outFileName,'w')
							f.write("qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qcovs	qcovhsp	qseq	sseq\n"+stdout)
							f.close()
			except:
				count = count-1
				continue

	for seq in os.listdir(queryDir):
		if i != ".DS_Store":
			queryName = os.path.splitext(i)[0]
			query_list.append(queryName)
	for target in os.listdir(subjectDir):
		if i != ".DS_Store":
			targetName = os.path.splitext(i)[0]
			target_list.append(targetName)


	tkMessageBox.showinfo("ORTHOPROK","Blast is finished!!")				
	displayedText.set('Blast is done!!')
	return query_list,target_list




def iterateDirNs(openpath,savepath,update):
	"""
	Generate the log file input for blaster v2.0
	automatic detection of .gz compressed genomes from genbank an uncompress
	log file also contain in the 3th column from fasta header file
	Arguments: /fasta_files_folder log_basename_file

	"""
	outputDir = savepath+"/"
	subjectDir = openpath+"/"
	
	

	if update =="No":

		output = open(outputDir+"target_sequences.log","wr")
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

	

	"""
		Code to create target sequences dictionary
	"""	


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

	if os.path.exists(outputDir+'BD_.dict.npy'):
		displayedText.set("target sequences database done!!")
		output_text.update_idletasks()
		window.update()		
	else:
		displayedText.set("Error in target sequences database!!")
		output_text.update_idletasks()
		window.update()	





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
			sequences = parseResults(f,ident,cov,"retrieve_sequences")					
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
	for genome_seq in sorted(presence_dict,key=natural_keys):
		if genome_seq in BD_dict:	
			genome = [genome_seq][0]
			value_dict = BD_dict[genome_seq].split("##")
			Strain =value_dict[0].strip("\n").replace(" ","_")
			assembly= "_".join(genome.split("_")[0:2])
			line = [genome_seq, Strain]		
			for protein in sorted(presence_dict[genome_seq],key=natural_keys):		
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




def QuitButtonCommand(window):
    window.destroy()





##########################################################################
#
#                   Tkinter  Widgets and Labels
#
##########################################################################



Openbutton = Button(window, text="Browse", command=Openfunc)	
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


Savebutton = Button(window, text="Browse", command=Savefunc)	
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


querybutton = Button(window, text="Browse", command=Queryfunc)	
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



Blastbutton = Button(window, text="Blast", command= lambda:  Blast(openpath,savepath,querypath,redo))	
Blastbutton.pack()
Blastbutton.place(y= 325, x=25)


redoblastbutton = Checkbutton(window, text="redo Blast", variable=redo)	
redoblastbutton.pack()
redoblastbutton.place(y= 325, x=145)
redo.set(0)



logbutton = Button(window, text="log", command= lambda:  iterateDirNs(openpath,savepath,"No"))
logbutton.pack()
logbutton.place(y= 280, x=25)

log_update = Button(window, text="log update", command= lambda:  iterateDirNs(openpath,savepath,"Yes"))
log_update.pack()
log_update.place(y= 280, x=150)




Analysisbutton = Button(window, text="Analysis", command= lambda:  analysis(savepath,openpath,cov,ident,collapsed,absence,height,width))
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

QuitButton = Button(window, text = "Quit", command=lambda: QuitButtonCommand(window))
QuitButton.pack(anchor='center')
QuitButton.place(x = 25, y = 750, width = 80, height = 25)





window.configure(background='snow2')
window.mainloop()