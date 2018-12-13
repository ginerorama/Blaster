#
# GenomeParser.py
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

UPDATELOGYES = "Yes"
UPDATELOGNO = "No"

BD_DICT_FILE = 'BD_.dict.npy'
LOG_FILE = "target_sequences.log"

########## GENOMES LOG

# TODO rename "parse_genomes"
def updatelog(outputDir, subjectDir, UPDATE, err = sys.stderr):
	
	err.write("Orthoprok_Facade updatelog\n")
	
	if UPDATE == UPDATELOGNO:
		_update_log(outputDir, subjectDir)
		err.write("Finished parsing target genomes.\n")
		
	dict_file = _create_seqs_dir(outputDir, err)
	
	return dict_file

# TODO: refactor
def _update_log(outputDir, subjectDir):
	
	sys.stderr.write("\tupdate --> NO\n")
	
	output = open(outputDir+"target_sequences.log","wr")
	
	for i in os.listdir(subjectDir):
		
		if i == ".DS_Store": continue # CPC2018
		
		sys.stderr.write("\tdir --> "+str(i)+"\n")
		
		path = os.path.abspath(subjectDir+i)
		name = os.path.splitext(i)[0]
		extension = os.path.splitext(i)[1]

		#check for gz files from ncbi
		# TODO: rename "file" vars
		if extension == ".gz":
			
			print "comprised gz file found\nuncompressing...(This could take some time!!"
			file = gzip.open(path,"rb")
			file = file.read()
			uncompressed_file = open(subjectDir+name,"w")
			uncompressed_file.write(str(file))
		else:	
			file = (open(path,"r")).readline()

		#read the files and check for fasta files
		if ">" not in file[0]:
			# TODO: raise exception
			displayedText.set("ERROR: " + name+ " is not a fasta file")
			output_text.update_idletasks()
			window.update()
		else:
			fasta_header = file.split("\n")[0][1::]
			name_fasta= (" ".join(fasta_header.split(" ")[1::]).split(",")[0].split(".")[0].split("_")[0])
			name_fasta= name_fasta.replace("circular","").replace("chromosome","").replace("complete","").replace("genome","").replace("sequence","").replace("assembly","")
		#	output.write(str(name)+"\t"+str(path)+"\t"+str(fasta_header)+"\n")		
			output.write(str(name_fasta)+"\t"+str(path)+"\t"+str(name)+"\n")
			sys.stderr.write(str(name_fasta)+"\t"+str(path)+"\t"+str(name)+"\n")
		
	output.close()
	
	return

def _replace_fasta(text):
	reptext = ""	
	mapping = [
		('.fna', ''),
		('.fasta', ''),
		('.fa', ''),
		('faa', ''),
		('frn', ''),
		('ffn', '')
		]
	
	for k, v in mapping:
		reptext = text.replace(k, v)
		
	return reptext
	
def _create_seqs_dir(outputDir, err):
	"""
		Code to create target sequences dictionary
	"""
	
	err.write("Creating target sequences dictionary...\n")
	
	dict_file = outputDir+BD_DICT_FILE
	log_file = outputDir+LOG_FILE
	
	try:
		with open(log_file, 'r') as selected_genomes:
			BD_dict = {}
			for genome in selected_genomes:
				
				err.write("\t"+str(genome)+"\n")
				
				genome_data = genome.strip().split("\t")
				genome_name = genome_data[0]
				genome_path = _replace_fasta(genome_data[1])
				genome_ID = genome_data[2].replace(",","") # nombre del genoma sin la extension, que sera usado en el diccionario
				BD_dict[genome_ID] = genome_name+"##"+genome_path
			
			np.save(dict_file, BD_dict)
			
	except IOError as e:
		err.write("Could not open "+log_file+"\n")
		raise e
	
	err.write("Created target sequences dictionary.\n")
	# TODO makeblastdb of all the targets
	
	return dict_file

## END