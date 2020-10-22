import re
import pandas as pd
import miscellaneous as misc
from xlsxwriter.workbook import Workbook
from datetime import datetime

def consume_input_vai(gene):
	'''

	Parameters:
		gene (str)

	Returns:

	'''

	vai_input_file = "./input_vai/{}_VAI.txt".format(gene)
	output_vai =  pd.DataFrame()

	with open (vai_input_file , 'r') as vai_input:
		for line in vai_input:
			if line.startswith("#"):
				pass
			elif line.startswith("Uploaded"):
				output_vai = pd.DataFrame(columns = line.strip().split('\t'))
			else:
				row = pd.Series(line.strip().split('\t'), index = output_vai.columns)
				output_vai = output_vai.append(row, ignore_index=True)

	return output_vai

def extract_info_from_vai_extra(output_vai, refseq_transcript):
	'''

	Parameters:
		gene (str)

	Returns:

	'''

	non_versioned_transcript = refseq_transcript.split('.')[0]	
	output_vai_filtered = output_vai[(output_vai.Feature.str.match(non_versioned_transcript)) & ((output_vai.Consequence == 'missense_variant') | (output_vai.Consequence == 'frameshift_variant') | (output_vai.Consequence == 'stop_gained') | (output_vai.Consequence == 'inframe_insertion'))]

	all_transcripts = []
	all_cdots = []
	all_cdot_numbers = []
	all_proteins = []
	all_pdots = []
	all_pdot_numbers = []
	all_exons = []
	all_interpros = []

	output_vai_filtered['transcript'] = ''
	output_vai_filtered['cdot'] = ''
	output_vai_filtered['cdot_number'] = ''
	output_vai_filtered['protein'] = ''
	output_vai_filtered['pdot'] = ''
	output_vai_filtered['pdot_number'] = ''
	output_vai_filtered['exon'] = ''


	for line in output_vai_filtered.Extra.items():
		hgvs_c_start = line[1].find('HGVSCN') + 7
		hgvs_c_end = hgvs_c_start + line[1][hgvs_c_start:].find(';')
		hgvs_p_start = line[1].find('HGVSP') + 6
		hgvs_p_end = hgvs_p_start + line[1][hgvs_p_start:].find(';')
		exon_start = line[1].find('EXON') + 5
		exon_end = exon_start + line[1][exon_start:].find('/')
		
		hgvs_c = line[1][hgvs_c_start:hgvs_c_end]
		hgvs_p = line[1][hgvs_p_start:hgvs_p_end]
		exon = int(line[1][exon_start:exon_end])
		transcript = hgvs_c[:hgvs_c.find(':')]
		cdot = hgvs_c[hgvs_c.find(':') + 1:]
		cdot_number = int(''.join(re.findall(r'\d+', cdot)))
		protein = hgvs_p[:hgvs_p.find(':')]
		pdot = hgvs_p[hgvs_p.find(':') + 1:].replace('(','').replace(')','')
		pdot_number = int(''.join(re.findall(r'\d+', pdot)))
		all_transcripts.append(transcript)
		all_cdots.append(cdot)
		all_cdot_numbers.append(cdot_number)
		all_proteins.append(protein)
		all_pdots.append(pdot)
		all_pdot_numbers.append(pdot_number)
		all_exons.append(exon)

		interpro_start = line[1].find('INTERPRO') + 9
		if interpro_start > -1:
			interpro_end = interpro_start + line[1][interpro_start:].find(';')
			interpro = line[1][interpro_start:interpro_end]
		else:
			interpro = ''

		all_interpros.append(interpro)

	output_vai_filtered['transcript'] = all_transcripts
	output_vai_filtered['cdot'] = all_cdots
	output_vai_filtered['cdot_number'] = all_cdot_numbers
	output_vai_filtered['protein'] = all_proteins
	output_vai_filtered['pdot'] = all_pdots
	output_vai_filtered['pdot_number'] = all_pdot_numbers
	output_vai_filtered['exon'] = all_exons
	output_vai_filtered['interpro'] = all_interpros

	return output_vai_filtered

def merge_data_frames(gene, output_vai_filtered):
	'''

	Parameters:
		gene (str)

	Returns:

	'''
	output_gdc = misc.open_dated_hgvs_file(gene = gene, file_type = 'gdc')
	output_gdc = output_gdc[['hgvs','disease']]
	output_clinvar = misc.open_dated_hgvs_file(gene = gene, file_type = 'clinvar')
	output_clinvar = output_clinvar[['hgvs','disease']]
	output_combined = pd.concat([output_gdc, output_clinvar])
	output_venn = misc.open_dated_hgvs_file(gene = gene, file_type = 'venn')
	output_vai_final = pd.merge(left = output_vai_filtered, right = output_combined, left_on = 'Uploaded Variation', right_on = 'hgvs')
	output_vai_final = pd.merge(left = output_vai_final, right = output_venn, left_on = 'Uploaded Variation', right_on = 'hgvs')

	return output_vai_final

def write_vai_final_file(gene, data):
	'''
	

	Parameters:
		

	Returns:
		

	'''
	vai_final_output_file = "./data_vai_final/{}_{}_VAI_final.tsv".format(str(datetime.date(datetime.now())),gene)
	vai_final_output = open(vai_final_output_file,'w+')
	header_tsv = "hgvs\tsource\tgene\tfeature\tconsequence\ttranscript\tcdot\tcdot_number\tprotein\tpdot\tpdot_number\texon\tinterpro\tdisease_type\n"
	vai_final_output.write(header_tsv)
	
	xlsx_file = "./data_vai_final/{}_overlap_analysis.xlsx".format(gene)
	workbook = Workbook(xlsx_file)
	worksheet = workbook.add_worksheet()
	header_xlsx = ('hgvs','source','gene','feature','consequence','transcript','cdot','cdot_number','protein','pdot','pdot_number','exon','interpro','disease_type')
	worksheet.write_row(0, 0, header_xlsx)

	for index,mutation in data.iterrows():
		if isinstance(mutation['disease'],float):
			mutation['disease']=''
		output_tsv = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(mutation['Uploaded Variation'].strip(),mutation['source'].strip(),mutation['Gene'].strip(),mutation['Feature type'].strip(),mutation['Consequence'].strip(),mutation['Feature'].strip(),mutation['cdot'].strip(),mutation['cdot_number'],mutation['protein'].strip(),mutation['pdot'].strip(),mutation['pdot_number'],mutation['exon'],mutation['interpro'].strip(),mutation['disease'].strip())
		output_xlsx = (mutation['Uploaded Variation'].strip(),mutation['source'].strip(),mutation['Gene'].strip(),mutation['Feature type'].strip(),mutation['Consequence'].strip(),mutation['Feature'].strip(),mutation['cdot'].strip(),mutation['cdot_number'],mutation['protein'].strip(),mutation['pdot'].strip(),mutation['pdot_number'],mutation['exon'],mutation['interpro'].strip(),mutation['disease'].strip())
		vai_final_output.write(output_tsv)
		worksheet.write_row(index+1, 0, output_xlsx)
	
	vai_final_output.close()
	workbook.close()

	print ("Written {}".format(vai_final_output_file))
	print ("Written {}\n\nPLEASE PROCEED TO UPLOAD\n\n".format(xlsx_file))

def process_vai_to_data_final(gene,refseq_transcript):
	'''

	Parameters:
		gene (str)

	Returns:

	'''
	output_vai = consume_input_vai(gene = gene)
	output_vai_filtered = extract_info_from_vai_extra(output_vai = output_vai, refseq_transcript = refseq_transcript)
	output_vai_final = merge_data_frames(gene = gene, output_vai_filtered = output_vai_filtered)
	write_vai_final_file(gene = gene, data = output_vai_final)

	return output_vai_final

def main():
	'''
	Run as:

	import sys
	sys.path.append('/path/to/mutation-gatherer/mutation-gatherer/')
	import process_vai as vai


	'''
	pass    


if __name__ == "__main__":
	main()