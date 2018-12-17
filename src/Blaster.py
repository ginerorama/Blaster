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

import os, sys
from Bio.Blast.Applications import NcbiblastnCommandline

import GenomeParser as genomeParser
import utils

# CPC2018
BLAST_OUTFMT='"6 qseqid sseqid pident length \
mismatch gapopen qstart qend sstart send \
evalue bitscore qcovs qcovhsp qseq sseq"'
BLAST_EVALUE='0.00001'
BLAST_NUMTHREADS='4'
BLAST_HEADER="qseqid	sseqid	pident	length	mismatch	\
gapopen	qstart	qend	sstart	send	\
evalue	bitscore	qcovs	qcovhsp	qseq	sseq"

######## BLAST

def Blast(queryDir, subjectDir, outputDir, 
		  known_query_list, known_target_list,
		  redo, displayf, err = sys.stderr):
	
	err.write("\tReading list of queries...\n")
	# CPC2018
	queries_list = [name for name in os.listdir(queryDir)
					if name != genomeParser.DS_STORE and
					os.path.isfile(os.path.join(queryDir, name)) and
					utils.is_fasta(name)]
	
	totalfiles = len(queries_list) # -1 # CPC2018
	
	err.write("\tTotal queries found: "+str(totalfiles)+"\n")
	
	err.write("\tReading list of targets...\n")
	targets_list = os.listdir(subjectDir) # CPC2018
	err.write("\tTotal targets found: "+str(len(targets_list))+"\n")
	
	target_tmp_list = []
	
	#for i in os.listdir(queryDir):
	count = 0
	for i in queries_list: # CPC2018
		#if i == ".DS_Store": continue # CPC2018
		
		count += 1
		displayf("Running blast: "+str(count)+" of "+str(totalfiles)+" sequences")
		#queryFile = queryDir+i
		try: 		
			queryFile = '"'+queryDir+"/"+i+'"'
			queryName = os.path.splitext(i)[0]# cuidado si interaccionas con los ficheros en finder y aparece DS_store aqui fallara hayq eliminar los DS_store files antes
			#for j in os.listdir(subjectDir):
			for j in targets_list:
				if utils.is_fasta(j) and _redo(redo, i, known_query_list, j, known_target_list): # CPC2018
					targetName = os.path.splitext(j)[0] # CPC2018
					subjectFile = '"'+subjectDir+"/"+j+'"'
					
					err.write("\tRunning blast: "+str(queryName)+" --> "+str(targetName)+"\n")
					
					_tblastn(queryName, queryFile, targetName, subjectFile, outputDir)
					
					err.write("\tBlast finished: "+str(queryName)+" --> "+str(targetName)+"\n")
					
					# CPC2018					
					if targetName not in set(target_tmp_list):
						target_tmp_list.append(targetName)
				
			# CPC2018
			if queryName not in set(known_query_list):
				known_query_list.append(queryName)
						
		except Exception as e: # CPC2018
			print e # CPC2018
	
	# CPC2018
	for target in target_tmp_list:
		if target not in set(known_target_list):
			known_target_list.append(target)
	
	err.write("Finished blast of all queries and targets.\n")
	
	return

def _redo(redo, query, known_query_list, target, known_target_list):
	return (redo == 0 or \
			(redo == 1 and
			 query not in set(known_query_list) and
			 target not in set(known_target_list)))
	
def _tblastn(queryName, queryFile, targetName, subjectFile, outputDir):
	blastn_cline = NcbiblastnCommandline(cmd='tblastn',
										query=queryFile,
										subject=subjectFile,
										outfmt=BLAST_OUTFMT,
										evalue=BLAST_EVALUE,
										num_threads=BLAST_NUMTHREADS)
	stdout, stderr = blastn_cline()
	
	outFileName = outputDir+"/"+queryName+"##"+targetName+".txt"
	with open(outFileName,'w') as f:
		f.write(BLAST_HEADER+"\n"+stdout)
	
	return

## END