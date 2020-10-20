import requests
import json
import pandas as pd
import hgvs.parser
from datetime import datetime


def query_gdc_api(endpoint,field,value,fields='NA'):
	'''
	Interact with GDC API (https://gdc.cancer.gov/developers/gdc-application-programming-interface-api).

	Parameters:
		endpoint (str): The GDC API endpoint to interact with (https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/#api-endpoints).
		field (str): The field to use to filter the GDC API call (https://docs.gdc.cancer.gov/API/Users_Guide/Python_Examples/#a-filtered-query).
		value (str):The used to use to filter the GDC API call.
		fields (str): (OPTIONAL) Specific fields you'd want to retrieve, that can be found via the _mapping endpoint (https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#endpoints). WARNING: It only works with one field.

	Returns:
		response (json): JSON object with query response.
	'''

	url =  'https://api.gdc.cancer.gov'
	
	endpoint_url = "{}/{}".format(url,endpoint)

	filters = {
	'op':'in',
	'content' : {
	'field' : field,
	'value' : value
	}
	}

	if fields != 'NA':
		fields = [fields]
		
		parameters = {
		'filters' : json.dumps(filters),
		'fields' : fields,
		'format' : 'JSON',
		'size' : 10000 #how to get all?
		}

	else:
		
		parameters = {
		'filters' : json.dumps(filters),
		'format' : 'JSON',
		'size' : 10000 #how to get all?
		}

	response = requests.get(endpoint_url, params = parameters)

	return response


def gdc_response_to_dataframe(response = 'NA' , list_to_strings = False):
	'''
	Convert GDC JSON response hits to dataframe using the response.json()['data']['hits'] format.

	Parameters:
		response (json): JSON object with query response.
		list_to_strings (boolean): NOT WORKING NOW. Will be used to convert lists inside columns to strings separated by commas.

	Returns:
		output (dataframe): The response object as a pandas dataframe.
	'''
	
	output = pd.DataFrame(response.json()['data']['hits'])

	return output

	
def merge_field_to_standard_response(standard_response, field_response, field_root, field, field_subroot = 'NA', standard_format = 'json', field_format = 'json'):
	'''
	Return the output dataframe with an extra column with case_id information.

	Parameters:
		standard_response (json): JSON object with query response. Obtained from: query_gdc_api(endpoint,field,value,fields='NA').
		field_response (json): JSON object with query response with one field. Obtained from: query_gdc_api(endpoint,field,value,fields='fields').
		field_root (str): Main root of field_response field. Use _mapping endpoint to identify (https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#endpoints). WARNING: This program can only do one or two levels deep (root and subroot), separated by dots from field.
		field (str): Field name without root or subroot information. Use _mapping endpoint to identify. 
		field_subroot (str): Subroot of field_response field. Use _mapping endpoint to identify. 
		standard_format (str): (OPTIONAL) Indicate whether the standard response is a json object or has already been transformed to a dataframe.
		field_format (str): Indicate whether the field response is a json object or has already been transformed to a dataframe.

	Returns:
		output (dataframe): Standard_response dataframe with field_response extra column.

	'''

	if standard_format == 'json':
		standard_dataframe = gdc_response_to_dataframe(standard_response)
	elif standard_format == 'dataframe':
		standard_dataframe = standard_response
	if field_format == 'json':
		field_tmp_dataframe = gdc_response_to_dataframe(field_response)
	elif field_format == 'dataframe':
		field_tmp_dataframe == field_response

	field_values = []
	for i in range(0, len(field_tmp_dataframe)):
		standard_id = field_tmp_dataframe['id'][i].strip()
		if field_subroot == 'NA':
			field_value = field_tmp_dataframe[field_root][i][0][field].strip()
		else:
			field_value = field_tmp_dataframe[field_root][i][0][field_subroot][field].strip()
		field_values.append((standard_id,field_value))

	field_dataframe =  pd.DataFrame(field_values, columns = ('id',field))

	output = pd.merge(left = standard_dataframe, right = field_dataframe, left_on = 'id', right_on = 'id')

	return output


def standardize_hgvs(data, parsee):
	'''
	

	Parameters:
		data (dataframe): 

	Returns:
		data (dataframe): 

	'''
	
	chr2grch38 = pd.read_csv('chr2grch38.txt', sep = "\t")
	data['hgvs'] = ''
	all_hgvs = []

	for index,value in enumerate(data['genomic_dna_change']):
		hgvs_data =  parsee.parse_hgvs_variant(value)
		chrom = hgvs_data.ac.strip()
		ncid = chr2grch38['grch38'][chr2grch38['chromosome'] == chrom].to_string(index=False)
		hgvs = value.replace(chrom,ncid)
		all_hgvs.append(hgvs)
	
	data['hgvs'] = all_hgvs
	
	return data


def query_based_on_results(data, field, fields, endpoint):
	'''
	

	Parameters:
		data (dataframe): 

	Returns:
		data (dataframe): 

	'''

	data[fields] = ''
	all_fields = []
	for index, case in enumerate(data[field]):
		case_id = case.strip()
		response_case = query_gdc_api(endpoint = endpoint, field = field, value = case_id, fields = fields)
		decoded_case = response_case.json()['data']['hits'][0][fields].strip()
		all_fields.append(decoded_case)
		
	data[fields] = all_fields
	return data


def write_hgvs_gdc_file(gene, data):
	'''
	

	Parameters:
		

	Returns:
		

	'''
	gdc_output_file = "./data_gdc/{}_{}_GDC.tsv".format(str(datetime.date(datetime.now())),gene)
	gdc_output = open(gdc_output_file,'w+')
	gdc_output.write("variant_id\tdisease\thgvs\tdataset\n")
	for index,mutation in data.iterrows():
		output_tsv = "{}\t{}\t{}\tgdc\n".format(mutation['ssm_id'].strip(),mutation['disease_type'].strip(),mutation['hgvs'].strip())
		gdc_output.write (output_tsv)
	gdc_output.close()
	print ("Written {}".format(gdc_output_file))

def process_gdc_per_gene(gene, eda = 'no'):
	'''
	

	Parameters:
		

	Returns:
		

	'''

	#get all the mutations
	response = query_gdc_api(endpoint = 'ssms', field = 'consequence.transcript.gene.symbol', value = gene)

	if response.json()['data']['pagination']['total'] > 0:
		print ("{} has GDC response".format(gene))
		#get all the case identifiers to mutations
		response_case_ids = query_gdc_api(endpoint = 'ssms', field = 'consequence.transcript.gene.symbol', value = gene, fields = 'occurrence.case.case_id')

		#merge into one table
		gene_data = merge_field_to_standard_response(standard_response = response, field_response = response_case_ids, field_root = 'occurrence', field_subroot = 'case', field = 'case_id')

		#change hgvs to standard expression
		hp = hgvs.parser.Parser()
		gene_data = standardize_hgvs(data = gene_data, parsee = hp)

		#add disease_type data
		gene_data = query_based_on_results(data = gene_data, field = 'case_id', fields = 'disease_type', endpoint = 'cases')

		#write GDC file
		write_hgvs_gdc_file(gene= gene, data = gene_data)
		if eda == 'yes':
			return gene_data

	else:
		print ("{} doesn't have GDC response".format(gene))
		if eda== 'yes':
			gene_data = pd.DataFrame()
			return gene_data


			
def main():
	'''
	Run as:

	import sys
	sys.path.append('/path/to/mutation-gatherer/mutation-gatherer/')
	import interact_with_gdc as gdc

	'''
	pass
	


if __name__ == "__main__":
	main()