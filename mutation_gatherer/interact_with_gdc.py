import requests
import json
import pandas as pd


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