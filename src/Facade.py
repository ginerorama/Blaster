#
# Facade.py
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

import sys

import GenomeParser as genomeParser
import Blaster as blaster
import Analyzer as analyzer
import AlignmentGenerator as alignmentGenerator
import Plotter as plotter

########## GENOMES LOG

# TODO Rename "parse_genomes"
def updatelog(savepath, openpath, UPDATE, err = sys.stderr):
	# returns dict_file
	return genomeParser.updatelog(savepath, openpath, UPDATE, err)

######## BLAST

def Blast(queryDir, subjectDir, outputDir, 
		  known_query_list, known_target_list,
		  redo, displayf, err = sys.stderr):
	
	return blaster.Blast(queryDir, subjectDir, outputDir, 
		  known_query_list, known_target_list,
		  redo, displayf, err)

####### Analysis

def analysis(savepath, cov, ident, collapsed, absence, err = sys.stderr):
	
	return analyzer.analysis(savepath, cov, ident, collapsed, absence, err)

###### Sequence alignment

def sequence_alignment(savepath, querypath, cov, ident, display, err = sys.stderr):
	
	return alignmentGenerator.sequence_alignment(savepath, querypath, cov, ident, display, err)

#######

def protein_length_graphic(prot_name, savepath, height, width, font_scale, err = sys.stderr):
	
	return plotter.protein_length_graphic(prot_name, savepath, height, width, font_scale, err)

#######

def stripplot(savepath, height, width, err = sys.stderr):
	return plotter.stripplot(savepath, height, width, err)

#######

def heatmap(savepath, height, width, font_scale, right_scale, bottom_scale, err = sys.stderr):
	
	return plotter.heatmap(savepath, height, width,
						   font_scale, right_scale, bottom_scale, err)

## END