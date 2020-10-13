import requests
import xmltodict
import hgvs.parser as hgvs

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
	all_hgvs = []
	ids_query = ','.join(ids)
	url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id="
	query_ids_url= "{}{}&from_esearch=true".format(url,ids_query)
	data = requests.get(query_ids_url)
	for i in range(0,len(ids)):
		try:
			hgvs = xmltodict.parse(data.content)['ClinVarResult-Set']['VariationArchive'][i]['InterpretedRecord']['SimpleAllele']['HGVSlist']['HGVS'][0]['NucleotideExpression']['Expression'].strip()
		except:
			hgvs = ''
		try:
			method_type = xmltodict.parse(data.content)['ClinVarResult-Set']['VariationArchive'][i]['InterpretedRecord']['ClinicalAssertionList']['ClinicalAssertion']['ObservedInList']['ObservedIn']['Method']['MethodType'].strip()
		except:
			method_type = ''
		try:
			origin = xmltodict.parse(data.content)['ClinVarResult-Set']['VariationArchive'][i]['InterpretedRecord']['ClinicalAssertionList']['ClinicalAssertion']['ObservedInList']['ObservedIn']['Sample']['Origin'].strip()
		except:
			origin = ''
		try:
			clinical_significance = xmltodict.parse(data.content)['ClinVarResult-Set']['VariationArchive'][i]['InterpretedRecord']['RCVList']['RCVAccession']['@Interpretation']
		except:
			clinical_significance = ''
		
		if hgvs.startswith('NC',0) and method_type == 'clinical testing' and origin == 'germline' and (clinical_significance == 'Pathogenic' or clinical_significance == 'Likely Pathogenic'):
			all_hgvs.append(hgvs)
	return all_hgvs


def write_hgvs_clinvar_file(gene, data):
	'''
	

	Parameters:
		

	Returns:
		

	'''

	clinvar_output_file = "{}_ClinVar.tsv".format(gene)
	clinvar_output = open(clinvar_output_file,'w+')
	clinvar_output.write("ssm_id\tdisease_type\thgvs\tdataset\n")

	for hgvs in enumerate(data):
		output = "NA\tNA\t{}\tclinvar\n".format(hgvs.strip())
		clinvar_output.write(output)

	clinvar_output.close()


def process_gdc_per_gene(gene, eda = 'no'):
	'''
	

	Parameters:
		

	Returns:
		

	'''

	#get all the mutations
	response = query_gdc_api(endpoint = 'ssms', field = 'consequence.transcript.gene.symbol', value = gene)

	#get all the case identifiers to mutations
	response_case_ids = query_gdc_api(endpoint = 'ssms', field = 'consequence.transcript.gene.symbol', value = gene, fields = 'occurrence.case.case_id')

	#merge into one table
	gene_data = merge_field_to_standard_response(standard_response = response, field_response = response_case_ids, field_root = 'occurrence', field_subroot = 'case', field = 'case_id')

	#change hgvs to standard expression
	hp  = hgvs.Parser()
	gene_data = standardize_hgvs(data = gene_data, parsee = hp)

	#add disease_type data
	gene_data = query_based_on_results(data = gene_data, field = 'case_id', fields = 'disease_type', endpoint = 'cases')

	#write GDC file
	write_hgvs_gdc_file(gene= gene, data = gene_data)
	if eda == 'yes':
		return gene_data


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