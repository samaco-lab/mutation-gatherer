import os
import re
import pandas as pd
import numpy as np
import datefinder
import datetime
import hgvs.parser
from datetime import datetime

def open_dated_hgvs_file(gene, file_type):
	'''

	Parameters:
		gene (str)

	Returns:

	'''	
	gene_pattern = re.compile(gene)

	if file_type == 'gdc':
		folder = "./data_gdc/"
	elif file_type == 'clinvar':
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

	if len(datestamps) > 0:
		most_recent = sorted(datestamps, reverse = True)[0]
	
		if file_type == 'gdc':
			file = "./data_gdc/{}_{}_GDC.tsv".format(most_recent,gene)
		if file_type == 'clinvar':
			file = "./data_clinvar/{}_{}_ClinVar.tsv".format(most_recent,gene)

		print ("Retrieved {}".format(file))
		data = pd.read_csv(file,sep='\t')

	else:
		print ("No file found for {} in {}".format(gene, file_type))

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
		try:
			var_g = hp.parse_hgvs_variant(hgvs_g)
			all_var_g.append(str(var_g.posedit))
		except:
			all_var_g.append('')

	data['var_g'] = all_var_g
	return data

def obtain_ncid(hgvs, parsee):
	'''
	
	
	Parameters:
		hgvs (str)

	Returns:

	'''
	chr2grch38 = pd.read_csv('chr2grch38.txt', sep = "\t")
	hgvs = parsee.parse_hgvs_variant(hgvs).ac
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


def venn2_groups(data1, data2, label1, label2, ncid):
	'''
	
	
	Parameters:
		gene (str)

	Returns:

	'''	
	venn2_df = pd.DataFrame(columns=['hgvs','source'])
	
	varg1 = np.array(list(set(data1['var_g'].values)))
	varg2 = np.array(list(set(data2['var_g'].values)))
	
	varg1_only = difference_var_g(varg1, varg2)
	label1_only = np.repeat(label1,varg1_only.size)
	varg2_only = difference_var_g(varg2, varg1)
	label2_only = np.repeat(label2,varg2_only.size)
	varg12_intersect = intersection_var_g(varg1, varg2)
	label12_intersect = np.repeat("{}_{}".format(label1,label2),varg12_intersect.size)
	
	print (ncid)
	venn2_df['hgvs'] = list(varg1_only) + list(varg2_only) + list(varg12_intersect)
	venn2_df['hgvs'] = ncid + venn2_df['hgvs'].astype(str)
	venn2_df['source'] = list(label1_only) + list(label2_only) + list(label12_intersect)
	return venn2_df

def venn3_groups(data1, data2, data3, label1, label2, label3, ncid):
	'''
	
	
	Parameters:
		gene (str)

	Returns:

	'''	
	venn3_df = pd.DataFrame(columns=['hgvs','source'])
	
	varg1 = np.array(list(set(data1['var_g'].values)))
	varg2 = np.array(list(set(data2['var_g'].values)))
	varg3 = np.array(list(set(data3['var_g'].values)))
	varg12 = union(varg1, varg2)
	varg13 = union(varg1, varg3)
	varg23 = union(varg2, varg3)
	
	varg1_only = difference_var_g(varg1, varg23)
	label1_only = np.repeat(label1,varg1_only.size)
	varg2_only = difference_var_g(varg2, varg13)
	label2_only = np.repeat(label2,varg2_only.size)
	varg3_only = difference_var_g(varg3, varg12)
	label3_only = np.repeat(label3,varg3_only.size)
	varg123_intersect = intersection_var_g(var_g1, var_g23)
	label123_intersect = np.repeat("{}_{}_{}".format(label1,label2,label3),varg122_intersect.size)
	varg12_intersect_only = difference_var_g(intersection_var_g(varg1, varg2), varg123_intersect)
	label12_intersect_only = np.repeat("{}_{}".format(label1,label2),varg12_intersect_only.size)
	varg13_intersect_only = difference_var_g(intersection_var_g(varg1, varg3), varg123_intersect)
	label13_intersect_only = np.repeat("{}_{}".format(label1,label3),varg13_intersect_only.size)
	varg23_intersect_only = difference_var_g(intersection_var_g(varg2, varg3), varg123_intersect)
	label23_intersect_only = np.repeat("{}_{}".format(label2,label3),varg23_intersect_only.size)
	
	venn3_df['hgvs'] = list(varg1_only) + list(varg2_only) + list(varg3_only) + list(varg123_intersect) + list(varg12_intersect_only) + list(varg13_intersect_only) + list(varg23_intersect_only)
	venn3_df['hgvs'] = ncid + venn3_df['hgvs'].astype(str)
	venn3_df['source'] = list(label1_only) + list(label2_only) + list(label3_only) + list(label123_intersect) + list(label12_intersect_only) + list(label13_intersect_only) + list(label23_intersect_only)
	
	return venn3


def print_hgvs_venn_file(gene, venn):
	'''
	
	
	Parameters:
		gene (str)

	Returns:

	'''	
	hgvs_venn_output_file = "./data_hgvs_venn/{}_{}_hgvs_venn.tsv".format(str(datetime.date(datetime.now())),gene)
	hgvs_venn_output = open(hgvs_venn_output_file,'w+')
	hgvs_venn_output.write("hgvs\tsource\n")
	for index,variant in venn.iterrows():
		output_tsv = "{}\t{}\n".format(variant['hgvs'].strip(),variant['source'].strip())
		hgvs_venn_output.write(output_tsv)
	hgvs_venn_output.close()
	print ("{} file written".format(hgvs_venn_output_file))


def create_hgvs_venn_list(gene, venn_groups = 2):
	'''
	
	
	Parameters:
		gene (str)

	Returns:

	'''	
	output_gdc = open_dated_hgvs_file(gene = gene, file_type = 'gdc')
	output_clinvar = open_dated_hgvs_file(gene = gene, file_type = 'clinvar')
	output_gdc = extract_var_g(data = output_gdc)
	output_clinvar = extract_var_g(data = output_clinvar)

	hp = hgvs.parser.Parser()

	try:
		ncid = obtain_ncid(hgvs = output_clinvar['hgvs'][0], parsee = hp)
	except:
		ncid = ''

	if ncid != '':
		print ("Doing overlap of hgvs data")
		if venn_groups == 2:
			venn = venn2_groups(data1 = output_gdc, data2 = output_clinvar, label1 = 'gdc', label2='clinvar', ncid = ncid)
		elif venn_groups ==3:
			venn = venn3_groups(data1 = output_gdc, data2 = output_clinvar, data3 = output_sfari, label1 = 'gdc', label2='clinvar', label3 = 'sfari', ncid = ncid)
		else:
			venn = pd.DataFrame()
		print_hgvs_venn_file(gene,venn)

	else:
		print ("Error on NCID {}".format(ncid))


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