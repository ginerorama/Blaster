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

import utils

UPDATELOGYES = "Yes"
UPDATELOGNO = "No"

BD_DICT_FILE = 'BD_.dict.npy'
LOG_FILE = "target_sequences.log"
DS_STORE = ".DS_Store"

def get_bd_dict_file(savepath):
	return savepath+"/"+BD_DICT_FILE

def get_log_file(savepath):
	return savepath+"/"+LOG_FILE

def get_ds_store(savepath):
	return savepath+"/"+DS_STORE

########## GENOMES LOG

# TODO rename "parse_genomes"
def updatelog(savepath, openpath, UPDATE, err = sys.stderr):
	
	if UPDATE == UPDATELOGNO:
		err.write("Start parsing target genomes...\n")
		_update_log(savepath, openpath, err)
		err.write("Finished parsing target genomes.\n")
		
	dict_file = _create_seqs_dir(savepath, err)
	
	return dict_file

def _update_log(outputDir, subjectDir, err):
	
	log_file = get_log_file(outputDir)
	with open(log_file,"wr") as output:
		
		for i in os.listdir(subjectDir):
			
			if i == DS_STORE: continue # CPC2018
			if not utils.is_fasta(i): continue
			
			fullpath = subjectDir+"/"+i
			
			if not os.path.isfile(fullpath): continue
			
			fafile = os.path.abspath(fullpath)
			name = os.path.splitext(i)[0]
			extension = os.path.splitext(i)[1]
	
			#check for gz files from ncbi
			# TODO: rename "file" vars
			if extension == ".gz":
				fafile = utils.uncompress_file(fafile, subjectDir+"/"+name)
			
			fafileheader = ""
			with open(fafile,"r") as fafileopen:
				fafileheader = fafileopen.readline().strip()
	
			#read the files and check for fasta files
			if ">" not in fafileheader[0]:
				raise Exception("ERROR: " + name + " is not a fasta file")
			else:
				target_name = _get_target_name(fafileheader)
				
				fasta_record = "\t".join([target_name, fafile, name])
				#if VERBOSE: err.write(fasta_record+"\n")
				
				output.write(fasta_record+"\n")
	
	return

def _get_target_name(fafileheader):
	fasta_header = fafileheader[1:]
	name_fasta = (" ".join(fasta_header.split(" ")[1:]).split(",")[0].split(".")[0].split("_")[0])
	
	replist = [("circular",""),
		("chromosome",""),
		("complete",""),
		("genome",""),
		("sequence",""),
		("assembly","")]
	
	name_fasta = utils.replace_list(name_fasta, replist)
	
	return name_fasta
	
def _create_seqs_dir(outputDir, err):
	"""
		Code to create target sequences dictionary
	"""
	
	err.write("Creating target sequences dictionary...\n")
	
	dict_file = get_bd_dict_file(outputDir)
	log_file = get_log_file(outputDir)
	
	try:
		with open(log_file, 'r') as selected_genomes:
			BD_dict = {}
			for genome in selected_genomes:
				
				genome_data = genome.strip().split("\t")
				genome_name = genome_data[0]
				genome_path = utils.replace_fasta(genome_data[1])
				genome_ID = genome_data[2] #.replace(",","") CPC2018 
				# nombre del genoma sin la extension, que sera usado en el diccionario
				#if VERBOSE: err.write("\t"+genome_ID+"\n")
				genome_value = "##".join([genome_name, genome_path])
				BD_dict[genome_ID] = genome_value
			
			np.save(dict_file, BD_dict)
			
	except IOError as e:
		err.write("Could not open "+log_file+"\n")
		raise e
	
	err.write("Created target sequences dictionary.\n")
	# TODO makeblastdb of all the targets
	
	return dict_file
	
## END