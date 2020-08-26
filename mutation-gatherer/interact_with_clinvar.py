import requests
import xmltodict

def query_clinvar_api_gene(term =):
    '''
    Interact with ClinVar API (https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/, https://www.ncbi.nlm.nih.gov/clinvar/docs/linking/).

    Parameters:
    	gene (str)

    Returns:

    '''	

    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'
    query = "{}db=clinvar&term={}[gene]+AND+single_gene[prop]&retmax=500".format(url,term)
    response = requests.get(response)
    pass

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