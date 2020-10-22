import re
import json
import requests
import pandas as pd
import hgvs.parser
from datetime import datetime

def query_mutalyzer_for_hgvs_g(refseq_transcript, variant):
	'''
	Parameters:
	
	Returns:
		response (json): JSON object with query response.
	'''
	
	hgvs_c = "{}:{}".format(refseq_transcript, variant)
	mutalyzer_url = "https://mutalyzer.nl/json/numberConversion?build=hg38;variant={}".format(hgvs_c)
	hgvs_g = requests.get(mutalyzer_url).json()[0]
	
	return hgvs_g


def find_valid_transcript_subversion(data,refseq_transcript):
	'''
	Parameters:
	
	Returns:
		response (json): JSON object with query response.
	'''
	for index,variant in enumerate(data.iloc[:,0]):
		if re.search ("^c.",variant):
			c_dot = variant
			break

	non_version_transcript = refseq_transcript.split('.')[0]
	current_subversion = int(refseq_transcript.split('.')[1])

	for subversion in range(current_subversion,-1,-1):
		adjusted_transcript = '.'.join([non_version_transcript,str(subversion)])
		hgvs_g = query_mutalyzer_for_hgvs_g(adjusted_transcript, c_dot)
		if hgvs_g != None:
			break

	return adjusted_transcript


def process_cdna_file(data, refseq_transcript):
	'''
	Parameters:
	
	Returns:
		response (json): JSON object with query response.
	'''
	variants = []
	for index,variant in enumerate(data.iloc[:,0]):
		if re.search ("^c.",variant):
			hgvs_g = query_mutalyzer_for_hgvs_g(refseq_transcript, variant)
			variants.append(hgvs_g)
	return variants

def extract_var_g(variants, parsee):
	'''

	Parameters:
		gene (str)

	Returns:

	'''	
	all_var_g = []
	data = pd.DataFrame(columns = ['hgvs', 'var_g'])

	for hgvs_g in variants:
		try:
			var_g = parsee.parse_hgvs_variant(hgvs_g)
			all_var_g.append(str(var_g.posedit))
		except:
			all_var_g.append('')

	data['hgvs'] = variants
	data['var_g'] = all_var_g
	
	return data

def write_hgvs_sfari_file(gene, data):
	'''

	Parameters:
		gene (str)

	Returns:

	'''	
	sfari_output_file = "./data_sfari/{}_{}_SFARI.tsv".format(str(datetime.date(datetime.now())),gene)
	sfari_output = open(sfari_output_file,'w+')
	sfari_output.write("hgvs\tsource\tvar_g\n")

	for index,mutation in data.iterrows():
		try:
			output_tsv = "{}\tsfari\t{}\n".format(mutation['hgvs'].strip(),mutation['var_g'].strip())
			sfari_output.write (output_tsv)
		except:
			pass
			
	sfari_output.close()
	print ("Written {}".format(sfari_output_file))



def process_sfari_cdna_per_gene(gene, refseq_transcript,parsee):
	'''
	Parameters:
	
	Returns:

	'''
	input_sfari_file = "./input_sfari/{}.tsv".format(gene)
	input_sfari = pd.read_csv(input_sfari_file, sep="\t")
	adjusted_transcript = find_valid_transcript_subversion(data = input_sfari,refseq_transcript = refseq_transcript)
	print ("Adjusting transcript to {} for Mutalyzer".format(adjusted_transcript))
	variants = process_cdna_file(data = input_sfari, refseq_transcript = adjusted_transcript)
	output = extract_var_g(variants = variants, parsee = parsee)
	write_hgvs_sfari_file(gene = gene, data = output)


def main():
	'''
	Run as:

	import sys
	sys.path.append('/path/to/mutation-gatherer/mutation-gatherer/')
	import process_sfari_cdna as sfari

	'''
	pass    


if __name__ == "__main__":
	main()