#
# Analyzer.py
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

import os, sys
import numpy as np
import pandas as pd

import GenomeParser as genomeParser
import utils as utils

OUTPUT_FOLDER = "/output/"
RESULTS_PATH = "/results/"
BLAST_FOLDER = "/blast/"
TMP_FOLDER = "/tmp/"
OUTPUT_FILE = "/Analysis_presence_summary.csv"
FILE2 = "/Analysis_pseudo_summary.csv"
FILE3 = "/HitsPresence.csv"
FILE4 = "/HitsPseudo.csv"
FILE5 = "/Proteins_size.csv"
FILE6 = "/Proteins_ident.csv"

def get_resultspath(savepath):
	return savepath+RESULTS_PATH

def get_blastfolder(savepath):
	return savepath+BLAST_FOLDER

def get_tmpfiles(savepath):
	return savepath+TMP_FOLDER

def get_outputfolder(savepath):
	return savepath+OUTPUT_FOLDER

def get_outputfile(savepath):
	return get_resultspath(savepath)+OUTPUT_FILE

def get_file2(savepath):
	return get_resultspath(savepath)+FILE2

def get_file3(savepath):
	return get_tmpfiles(savepath)+FILE3

def get_file4(savepath):
	return get_tmpfiles(savepath)+FILE4

def get_file5(savepath):
	return get_resultspath(savepath)+FILE5

def get_file6(savepath):
	return get_resultspath(savepath)+FILE6

def get_dictfile(savepath):
	return savepath+"/"+genomeParser.BD_DICT_FILE

def analysis(savepath, cov, ident, collapsed, absence, err):
	
	err.write("Analyzer: analysis\n")
	
	resultspath = get_resultspath(savepath)
	if not os.path.exists(resultspath):
		os.makedirs(resultspath)
		
	blast_folder = get_blastfolder(savepath)
	if not os.path.exists(blast_folder):
		os.makedirs(blast_folder)	

	tmpfiles = get_tmpfiles(savepath)
	if not os.path.exists(tmpfiles):
		os.makedirs(tmpfiles)
	
	outputDir = get_outputfolder(savepath)
	outputFile = get_outputfile(savepath)
	outputFile2 = get_file2(savepath)
	outputFile3 = get_file3(savepath)
	outputFile4 = get_file4(savepath)
	outputFile5 = get_file5(savepath)
	outputFile6 = get_file6(savepath)
	
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
	BD_dict = np.load(get_dictfile(savepath)).item()

	first = True
	head = ["Assembly,Name"]
	head_tmp =[]
	count=1

	for assembly in sorted(resHash,key=utils.natural_keys):
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

			for query in sorted(resHash[assembly],key=utils.natural_keys):		
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
		
		for sub in sorted(resHash,key=utils.natural_keys):
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
	
	if os.path.exists(outputFile) and os.path.exists(outputFile2) :
		pass
	else:
		raise Exception("Error: presence analysis!!")
	
	err.write("Analyzer: analysis finished\n")

	return

def parseResults(file, cov, ident,function): #statsA
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
	
## END