import os
import re
import pandas as pd
import numpy as np
import datefinder
import datetime
import hgvs.parser


def open_dated_hgvs_file(gene, type):
	'''

	Parameters:
		gene (str)

	Returns:

	'''	
	gene_pattern = re.compile(gene)

	if type == 'gdc':
		folder = "./data_gdc/"
	elif type == 'clinvar':
		folder = "./data_clinvar/"
	
	datestamps = []
	data = pd.DataFrame()
	files = os.listdir(folder)
	
	for file in files:
		try:
			gene_file = gene_pattern.search(file).group()
		except:
			gene_file = gene_pattern.search(file)
		if gene_file == gene:
			date_match = str(list(datefinder.find_dates(file))[0].date())
			datestamps.append(date_match)

	if gene_file != None:
		most_recent = sorted(datestamps, reverse = True)[0]
	
		if type == 'gdc':
			file = "./data_gdc/{}_{}_GDC.tsv".format(most_recent,gene)
		if type == 'clinvar':
			file = "./data_clinvar/{}_{}_ClinVar.tsv".format(most_recent,gene)

		data = pd.read_csv(file,sep='\t')

	else:
		print ("No file found for {} in {}".format(gene, type))

	return data

def extract_var_g(data):
	'''

	Parameters:
		gene (str)

	Returns:

	'''	
	hp = hgvs.parser.Parser()
	all_var_g = []
	data['var_g'] = ''
	
	for hgvs_g in data['hgvs']:
		var_g = hp.parse_hgvs_variant(hgvs_g)
		all_var_g.append(str(var_g.posedit))

	data['var_g'] = all_var_g
	return data

def obtain_ncid(hgvs):
	'''
	
	
	Parameters:
		hgvs (str)

	Returns:

	'''
	hp = hgvs.parser.Parser()
	chr2grch38 = pd.read_csv('chr2grch38.txt', sep = "\t")
	hgvs = hp.parse_hgvs_variant(output_gdc['hgvs']).ac
	ncid = hgvs.split('.')[0]
	ncid_final = chr2grch38['grch38'][chr2grch38['grch38'].str.contains(ncid)].to_string(index = False)
	return ncid_final


def intersection_var_g(data1, data2):
	'''
	
	
	Parameters:
		gene (str)

	Returns:

	'''	
	intersection = np.intersect1d(data1,data2)
	return intersection


def union_var_g(data1, data2):
	'''
	
	
	Parameters:
		gene (str)

	Returns:

	'''	
	union = np.union1d(data1,data2)
	return union

def difference_var_g(data1,data2):
	'''
	
	
	Parameters:
		gene (str)

	Returns:

	'''	
	difference = np.setdiff1d(data1,data2)
	return difference


def venn2_groups(data1, data2, label1, label2):
	'''
	
	
	Parameters:
		gene (str)

	Returns:

	'''	
	venn2_df = pd.DataFrame('variant','source')
	varg1 = np.array(set(data1['var_g'].values))
	varg2 = np.array(set(data2['var_g'].values))
	varg1_only = difference_var_g(varg1, varg2)
	label1_only = np.repeat(label1,varg1_only.size)
	varg2_only = difference_var_g(varg2, varg1)
	label2_only = np.repeat(label2,varg2_only.size)
	varg12_intersect = intersection_var_g(var_g1, var_g2)
	label12_intersect = np.repeat("{}_{}".format(label1,label2),varg12_intersect.size)
	venn2_df['variant'] = pd.Series[varg1_only, varg2_only, varg12_intersect]
	venn2_df['source'] = pd.Series[label1_only, label2_only, label12_intersect]
	return venn2_df

def venn3_groups(data1, data2, data3):
	'''
	
	
	Parameters:
		gene (str)

	Returns:

	'''	
	var_g1 = np.array(set(data1['var_g'].values))
	var_g2 = np.array(set(data2['var_g'].values))
	var_g3 = np.array(set(data3['var_g'].values))
	var_g12 = union(var_g1, var_g2)
	var_g13 = union(var_g1, var_g3)
	var_g23 = union(var_g2, var_g3)
	var_g1_only = difference_var_g(var_g1, var_g23)
	var_g2_only = difference_var_g(var_g1, var_g13)
	var_g3_only = difference_var_g(var_g1, var_g12)
	var_g123_intersect = intersection_var_g(var_g1, var_g23)
	var_g12_intersect_only = difference_var_g(intersection_var_g(var_g1, var_g2), var_g123_intersect)
	var_g13_intersect_only = difference_var_g(intersection_var_g(var_g1, var_g3), var_g123_intersect)
	var_g23_intersect_only = difference_var_g(intersection_var_g(var_g2, var_g3), var_g123_intersect)
	venn3 = [var_g1_only, var_g2_only, var_g3_only, var_g123_intersect, var_g12_intersect_only, var_g13_intersect_only, var_g23_intersect_only]
	return venn3



def main():
	'''
	Run as:

	import sys
	sys.path.append('/path/to/mutation-gatherer/mutation-gatherer/')
	import miscellaneous as misc

	'''
	pass    


if __name__ == "__main__":
	main()