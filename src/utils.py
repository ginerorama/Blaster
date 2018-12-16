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

import os, sys, re
import gzip

def replace_list(text, replist):
	reptext = text
	
	for k, v in replist:
		reptext = reptext.replace(k, v)
		
	return reptext

def replace_fasta(text):
	reptext = ""
	mapping = [
		('.fna', ''),
		('.fasta', ''),
		('.fa', ''),
		('faa', ''),
		('frn', ''),
		('ffn', '')
		]
	
	reptext = replace_list(text, mapping)
		
	return reptext

def uncompress_file(gzfile, ungzfile):
	with gzip.open(path, 'rb') as gzfile:
		with open(ungzfile, 'w') as uncompressed_file:
			uncompressed_file.write(str(gzfile.read()))
			
	return ungzfile

def is_short_fasta(filen):
	return filen.endswith(".fna") or filen.endswith(".fa")

def is_fasta(filen):
	return filen.endswith(".fna") or filen.endswith(".fasta") or filen.endswith(".fa")

def change_extension_fasta(path):

	for filen in os.listdir(path):
		#if file.endswith(".fa") or file.endswith(".fna") :
		if is_short_fasta(filen):
			file_base = os.path.splitext(filen)[0]
			#os.rename(path+"/"+filen, str(path+"/"+file_base) + ".fasta")
			os.rename(path+"/"+filen, path+"/"+file_base+".fasta")
			
	return

def atoi(text):
	return int(text) if text.isdigit() else text

def natural_keys(text):
	'''
	alist.sort(key=natural_keys) sorts in human order
	http://nedbatchelder.com/blog/200712/human_sorting.html
	(See Toothy's implementation in the comments)
	'''
	return [ atoi(c) for c in re.split('(\d+)', text) ]
