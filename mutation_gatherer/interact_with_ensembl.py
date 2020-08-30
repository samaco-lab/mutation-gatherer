import requests
import json
import re
import pandas as pd

def query_ensembl_xrefs_symbol(gene):
	'''
	Interact with ENSEMBL REST "id" based APIs (https://rest.ensembl.org/; the ones with [endpoint]/id/[id]).
	'''
	url = 'https://rest.ensembl.org/xrefs/symbol/homo_sapiens'

	gene_url = "{}/{}".format(url, gene)

	parameters = {
	'expand' : 1
	}

	response = requests.get(gene_url, params = parameters,  headers = { "Content-Type" : "application/json"})

	return response


def query_ensembl_by_id(ensembl_id, endpoint):
	'''
	Interact with ENSEMBL REST "id" based APIs (https://rest.ensembl.org/; the ones with [endpoint]/id/[id]).
	
	Parameters:
		ensembl_id (str):
		endpoint (str): 
	
	Returns:
		response (json): JSON object with query response.
	'''
	url = 'https://rest.ensembl.org'
	
	id_url = "{}/{}/id/{}".format(url, endpoint, ensembl_id)

	parameters = {
	'expand' : 1
	}

	response = requests.get(id_url, params = parameters, headers={ "Content-Type" : "application/json"})

	return response

def query_ensembl_for_transcripts(gene):
	'''
	Using ENSEMBL API endpoints obtain the transcripts that have valid protein translations and RefSeq ids
	
	Parameters:
		gene (str):
	
	Returns:
		ensembl_valid_transcripts (list): #should it be dataframe?
	'''

	response = query_ensembl_xrefs_symbol(gene)
	decoded = response.json()

	ensembl_gene_ids = []
	for entry in decoded:
		if re.search("^ENS",entry['id']):
			ensembl_gene_ids.append(entry['id'])

	ensembl_valid_transcripts = []

	for ensembl_gene_id in ensembl_gene_ids:
		gene_output = query_ensembl_by_id(ensembl_gene_id,'lookup')
		gene_decoded = gene_output.json()

		for transcript_by_gene_data in gene_decoded['Transcript']:
			ensembl_transcript_id = transcript_by_gene_data['id'].strip()
			ensembl_transcript_start = transcript_by_gene_data['start']
			ensembl_transcript_end = transcript_by_gene_data['end']
			ensembl_transcript_strand = transcript_by_gene_data['strand']
			orientation = None
			if ensembl_transcript_strand == 1:
				orientation = 'forward'
			elif ensembl_transcript_strand == -1:
				orientation = 'reverse'
			
			transcript_id_response = query_ensembl_by_id(ensembl_transcript_id,'lookup')
			transcript_id_decoded = transcript_id_response.json()
			try:
				ensembl_protein_id = transcript_id_decoded['Translation']['id'].strip()
			except:
				ensembl_protein_id = 'NA'
			try:
				ensembl_protein_start = transcript_id_decoded['Translation']['start']
			except:
				ensembl_protein_start = -1
			try:
				ensembl_protein_end = transcript_id_decoded['Translation']['end']
			except:
				ensembl_protein_end = -1
			
			if ensembl_protein_id != 'NA':
				transcript_databases_response = query_ensembl_by_id(ensembl_transcript_id,'xrefs')
				transcript_databases_decoded = transcript_databases_response.json()

				for transcript_data in transcript_databases_decoded:          
					if transcript_data['dbname'] == 'RefSeq_mRNA':

						try:
							xref_id_score = transcript_data['xref_identity']
						except:
							xref_id_score = -1
						try:
							ens_id_score = transcript_data['ensembl_identity']
						except:
							ens_id_score = -1
						try:
							nm_id = transcript_data['display_id'].strip()
						except:
							nm_id = 'NA'
						
						ensembl_valid_transcripts.append([ensembl_gene_id, ensembl_transcript_id, nm_id, orientation, ensembl_transcript_start, ensembl_transcript_end, ensembl_protein_id, ensembl_protein_start, ensembl_protein_end])

	output = pd.DataFrame(ensembl_valid_transcripts, columns = ('gene_id', 'transcript_id', 'refseq_id', 'orientation', 'transcript_start', 'transcript_end', 'protein_id', 'protein_start', 'protein_end'))

	return output


def obtain_exon_coordinates(transcript):
	'''
	Using a valid ENSEMBL transcript id, obtain exon boundaries
	
	Parameters:
		transcript (str):
	
	Returns:
		output (dataframe): 
	'''
	transcript_id_response = query_ensembl_by_id(transcript,'lookup')
	transcript_id_decoded = transcript_id_response.json()
	exon_list = transcript_id_decoded['Exon']
	exon_coordinates = []
	for exon in exon_list:
			exon_coordinates.append([transcript,exon['id'], exon['start'], exon['end']])

	output = pd.DataFrame(exon_coordinates, columns = ('transcript_id', 'exon_id', 'exon_start', 'exon_end'))

	return output


def main():
	'''
	Run as:

	import sys
	sys.path.append('/path/to/mutation-gatherer/mutation-gatherer/')
	import interact_with_ensembl as ensembl

	'''
	#usage 
	pass    


if __name__ == "__main__":
	main()