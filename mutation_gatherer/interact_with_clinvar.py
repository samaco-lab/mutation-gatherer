import requests
import xmltodict
import pandas as pd
import numpy as np
from datetime import datetime

def query_clinvar_api_gene_for_ids(term):
	'''
	Interact with ClinVar API (https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/, https://www.ncbi.nlm.nih.gov/clinvar/docs/linking/).

	Parameters:
		gene (str)

	Returns:

	'''	

	url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
	query = "{}db=clinvar&term={}[gene]+AND+single_gene[prop]&retmax=500".format(url,term)
	response = requests.get(query)
	try:
		ids = xmltodict.parse(response.content)['eSearchResult']['IdList']['Id']
	except:
		ids = ''
		print ("Error with ClinVar API\n")
	return ids
	

def parse_clinvar_response_old(gene):
	'''
	DEPRECATED

	Parameters:
		gene (str)

	Returns:

	'''	

	clinvar_output_file = "{}_ClinVar.tsv".format(gene)
	clinvar_output = open(clinvar_output_file,'w+')
	time.sleep(10)
	for i in range (0,len(ids)):
		clinvar_id = ids[i].strip()
		query = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={}".format(clinvar_id)
		clinvar_by_id = requests.get(query)
		title = xmltodict.parse(clinvar_by_id.content)['eSummaryResult']['DocumentSummarySet']["DocumentSummary"]['title'].strip()
		try:
			clinical_significance = xmltodict.parse(clinvar_by_id.content)['eSummaryResult']['DocumentSummarySet']["DocumentSummary"]["clinical_significance"]["description"]
		except AttributeError:
			clinical_significance = ''
		try:
			canonical_spdi = xmltodict.parse(clinvar_by_id.content)['eSummaryResult']['DocumentSummarySet']["DocumentSummary"]['variation_set']['variation']['canonical_spdi'].strip()
		except AttributeError:
			canonical_spdi = ''
		output = "{}\t{}\t{}\t{}".format(i,clinvar_id,title,canonical_spdi)
		if clinical_significance == "Pathogenic" or clinical_significance == "Likely Pathogenic":
			print (output)
			clinvar_output.write(output)
			clinvar_output.write ("\n")
		time.sleep(10)
	return


def query_ids_for_data(ids):
	'''
	Interact with ClinVar API (https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/, https://www.ncbi.nlm.nih.gov/clinvar/docs/linking/).

	Parameters:
		gene (str)

	Returns:

	'''		
	ids_query = ','.join(ids)
	url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id="
	query_ids_url= "{}{}&from_esearch=true".format(url,ids_query)
	output = requests.get(query_ids_url)
	output_df = pd.DataFrame(columns=['clinvar_id','accession','interpreted_condition','clinical_significance','method_type','origin','hgvs','dataset'])
	for i in range(0,len(ids)):
		output_ordered_dict = xmltodict.parse(output.content)
		try:
			hgvs = output_ordered_dict['ClinVarResult-Set']['VariationArchive'][i]['InterpretedRecord']['SimpleAllele']['HGVSlist']['HGVS'][0]['NucleotideExpression']['Expression'].strip()
		except:
			hgvs = ''
		try:
			method_type = output_ordered_dict['ClinVarResult-Set']['VariationArchive'][i]['InterpretedRecord']['ClinicalAssertionList']['ClinicalAssertion']['ObservedInList']['ObservedIn']['Method']['MethodType'].strip()
		except:
			method_type = ''
		try:
			origin = output_ordered_dict['ClinVarResult-Set']['VariationArchive'][i]['InterpretedRecord']['ClinicalAssertionList']['ClinicalAssertion']['ObservedInList']['ObservedIn']['Sample']['Origin'].strip()
		except:
			origin = ''
		try:
			clinical_significance = output_ordered_dict['ClinVarResult-Set']['VariationArchive'][i]['InterpretedRecord']['RCVList']['RCVAccession']['@Interpretation']
		except:
			clinical_significance = ''
		try:
			accession = output_ordered_dict['ClinVarResult-Set']['VariationArchive'][i]['@Accession']
		except:
			accession = ''
		try:
			interpreted_condition = output_ordered_dict['ClinVarResult-Set']['VariationArchive'][i]['InterpretedRecord']['RCVList']['RCVAccession']['InterpretedConditionList']['InterpretedCondition']['#text']
		except:
			interpreted_condition = ''

		id_output = [ids[i], accession, interpreted_condition, clinical_significance, method_type, origin, hgvs, 'clinvar']
		output_df.loc[i] = id_output
		
	return output_df


def filter_clinvar_variants(data, filter_clinical='original', filter_origin = 'yes', filter_method = 'yes'):
	'''
	

	Parameters:
		

	Returns:
		

	'''	
	filtered_data = pd.DataFrame()

	pathogenic = data['clinical_significance'] == 'Pathogenic'
	likely_pathogenic = data['clinical_significance'] == 'Likely pathogenic'
	vus = data['clinical_significance'] == 'Uncertain significance'
	
	if filter_origin == 'yes':
		germline = data['origin'] == 'germline'
	else:
		germline = pd.Series(np.repeat(True,data.shape[0]))
	
	if filter_method == 'yes':
		clinical_testing = data['method_type'] == 'clinical testing'
	else:
		clinical_testing = pd.Series(np.repeat(True,data.shape[0]))	
	
	if filter_clinical == 'original':
		filtered_data = data[(pathogenic | likely_pathogenic) & clinical_testing & germline]
	
	elif filter_clinical == 'vus':
		filtered_data = data[(pathogenic | likely_pathogenic | vus) & clinical_testing & germline]
	
	return filtered_data


def write_hgvs_clinvar_file(gene, data):
	'''
	

	Parameters:
		

	Returns:
		

	'''

	clinvar_output_file = "./data_clinvar/{}_{}_ClinVar.tsv".format(str(datetime.date(datetime.now())),gene)
	clinvar_output = open(clinvar_output_file,'w+')
	clinvar_output.write("variant_id\tdisease\thgvs\tdataset\n")
	for index,mutation in data.iterrows():
		output_tsv = "{}\t{}\t{}\tclinvar\n".format(mutation['accession'].strip(),mutation['interpreted_condition'].strip(),mutation['hgvs'].strip())
		clinvar_output.write(output_tsv)
	clinvar_output.close()
	print ("{} file written".format(clinvar_output_file))


def process_clinvar_per_gene(gene, eda = 'no'):
	'''
	

	Parameters:
		

	Returns:
		

	'''

	#query ClinVar for the variant IDs per gene
	ids = query_clinvar_api_gene_for_ids(term = gene)

	if len(ids) > 0:
		print ("{} has ClinVar response".format(gene))
		#query ClinVar for the information per variant ids, and filter by germline, clinical testing, and pathogenic/likelypathogenic
		output = query_ids_for_data(ids=ids)

		#filter ClinVar entries by variants' clinical significance, method of obtention type, and variant origin
		filtered_output = filter_clinvar_variants(data = output, filter_clinical='original', filter_origin = 'yes', filter_method = 'yes')

		#write clinvar file
		write_hgvs_clinvar_file(gene= gene, data = filtered_output)

		
def main():
	'''
	Run as:

	import sys
	sys.path.append('/path/to/mutation-gatherer/mutation-gatherer/')
	import interact_with_clinvar as clinvar


	'''
	pass    


if __name__ == "__main__":
	main()