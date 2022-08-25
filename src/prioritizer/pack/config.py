from pathlib import Path
import os

class Config:
    #Config file of Ilyome
    mainpack =  os.path.join(Path.home(), '.ilyome')
    #ORPHA
    url_orpha_GARD = 'http://www.orphadata.org/data/xml/en_product6.xml'
    url_orpha_PARD = 'http://www.orphadata.org/data/xml/en_product4.xml'
    out_orpha_GARD = os.path.join(mainpack, 'Orpha_GARD')
    out_orpha_PARD = os.path.join(mainpack, 'Orpha_PARD')
    json_orpha_GARD = os.path.join(mainpack, 'Orpha_GARD.json')
    json_orpha_PARD = os.path.join(mainpack, 'Orpha_PARD.json')

    #HPO url and outputs
    url_HPO_GP = 'http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt'
    url_HPO_PG = 'http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt'
    out_hpo_PG = os.path.join(mainpack, 'phenotype_to_genes.txt')
    out_hpo_GP = os.path.join(mainpack, 'genes_to_phenotype.txt')
    json_hpo_GP = os.path.join(mainpack, 'genes_to_phenotype.json')
    json_hpo_PG = os.path.join(mainpack, 'phenotype_to_genes.json')
    json_hpo_inheritance = os.path.join(mainpack, 'hpo_inheritance.json')
    json_string_gene_associations = os.path.join(mainpack, 'GeneAssociations.json')

    #clinvar
    url_clinvar = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    url_clinvar_tbi = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
    clinvar_file = os.path.join(mainpack, 'clinvar.vcf.gz')
    clinvar_tbi_file = os.path.join(mainpack, 'clinvar.vcf.gz.tbi')

    #Omim
    Json_omim = os.path.join(mainpack, 'omimdict.json')


    def __init__(self):
        if not os.path.isdir(self.mainpack):
            os.mkdir(self.mainpack)