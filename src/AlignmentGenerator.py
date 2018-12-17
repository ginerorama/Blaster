#
# AlignmentGenerator.py
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
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
import subprocess

import Analyzer as analyzer
import GenomeParser as genomeParser
import utils as utils

ALIGNMENT_FOLDER = "/alignment/"
MUSCLE_RESULTS = "muscle_results/"
MVIEW_RESULTS = "mview_results/"
ETE_RESULTS = "ete_results/"
MUSCLE_ALIGN_FASTA = "_muscle_align.fas"
MUSCLE_ALIGN_TEXT = "_muscle_align.txt"
MUSCLE_ALIGN_HTML = "_muscle_align.html"
MUSCLE_ALIGN_DAT = "_muscle_align.dat"

def get_alignment_folder(savepath):
	return savepath+ALIGNMENT_FOLDER

def get_muscle_results(savepath):
	return get_alignment_folder(savepath)+MUSCLE_RESULTS

def get_mview_results(savepath):
	return get_alignment_folder(savepath)+MVIEW_RESULTS

def get_ete_results(savepath):
	return get_alignment_folder(savepath)+ETE_RESULTS

def get_muscle_align_fasta(savepath, filename):
	return get_muscle_results(savepath)+filename+MUSCLE_ALIGN_FASTA

def get_muscle_align_text(savepath, filename):
	return get_muscle_results(savepath)+filename+MUSCLE_ALIGN_TEXT

def get_muscle_align_html(savepath, filename):
	return get_muscle_results(savepath)+filename+MUSCLE_ALIGN_HTML

def get_muscle_align_dat(savepath, filename):
	return analyzer.get_tmpfiles(savepath)+filename+MUSCLE_ALIGN_DAT

def sequence_alignment(savepath, querypath, cov, ident, display, err = sys.stderr):
	
	seq_resultspath = get_alignment_folder(savepath)
	if not os.path.exists(seq_resultspath):
		os.makedirs(seq_resultspath)
	
	resultsDir = analyzer.get_outputfolder(savepath)
	tmpfiles = analyzer.get_tmpfiles(savepath)
	
	#dictionaries
	presence_dict = {}
	genome_file = {}
	sequences_dict={}
	
	#parameters for blast protein analysis retrieving, identity and coverage
	# for each file of blast results
	for i in os.listdir(resultsDir):
		if i.endswith(".txt"):
			fullpath = resultsDir+i
			with open(fullpath,'r')	as f:
				err.write(str(f)+str(ident)+str(cov)+"\n")
				sequences = analyzer.parseResults(f, ident, cov, "retrieve_sequences")					
				
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
	BD_dict = np.load(savepath+genomeParse.BD_DICT_FILE).item()
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
	
	outputDir = get_muscle_results(savepath)
	#outputDir  = seq_resultspath+"muscle_results/"
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
	
	outputDir2 = get_mview_results(savepath)
	#outputDir2  = seq_resultspath+"mview_results/"
	if not os.path.exists(outputDir2):
		os.makedirs(outputDir2)
	
	outputDir3 = get_ete_results(savepath)
	#outputDir3  = seq_resultspath+"ete_results/"
	if not os.path.exists(outputDir3):
		os.makedirs(outputDir3)
	
	""" align sequence from blast analysis in fasta 
	format(.fas), clustal format (.txt and .html) and generated .dat file
	for mview program """
	
	display('Alignment with muscle running!!')
	# displayedText.set('Alignment with muscle running!!')
	# output_text.update_idletasks()
	# window.update()
	
	for file in os.listdir(seq_resultspath):
		if file.endswith(".fasta"):
			filename = os.path.splitext(file)[0]
			res_fasta = get_muscle_align_fasta(savepath, filename)			
			res_txt = get_muscle_align_text(savepath, filename)
			res_html = get_muscle_align_html(savepath, filename)
			res_dat = get_muscle_align_dat(savepath, filename)
			
			file2 = seq_resultspath+file
			print file2
			muscle_cline = MuscleCommandline(input=file2,fastaout=res_fasta, clwout=res_txt, htmlout=res_html)
			os.system(str(muscle_cline) +">/dev/null 2>&1") 				
			muscle_cline2 = MuscleCommandline(input=file2,fastaout=res_dat)		
			os.system(str(muscle_cline2)+">/dev/null 2>&1")	
	
	""" align sequence alignment dat files from muscle with mview """
	
	display('mview alignment running!!')
	# displayedText.set('mview alignment running!!')
	# output_text.update_idletasks()
	# window.update()	
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
	
	return

## END