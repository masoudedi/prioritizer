import json
import os
from .clinic import Scores
from .config import Config
from .log import *
import pandas as pd

class Read_input_vep(Config):
    def __init__(self, TabularFile, ListOfClinics, DisorderInheritancePedigree, ROHnumber):
        super().__init__()
        self.ListOfClinics = ListOfClinics
        self.input = TabularFile
        self.DisorderInheritancePedigree = DisorderInheritancePedigree
        self.GeneP = {}
        self.hpo = {}
        self.omim = {}
        self.inheritance = {}
        self.associations = {}
        self.ROHnumber = ROHnumber
        self.allele_freq = 0.02
        self.HomoCount = 4
        self.header = {}
        self.variants = []
        self.setup()
    
    def setup(self):
        if self.isValid():
            self.Read_database()
            self.Read_header()
            self.read_variants()
            self.FilterVariants()
    
    def isValid(self):
        """
        The input file is ckecked to make sure if it's valid or not
        """
        if os.path.isfile(self.input):
            with open(self.input, 'r') as T:
                if T.readline().startswith('## ENSEMBL VARIANT EFFECT PREDICTOR'):
                    return True
                else:
                    logging.error(f'Vep input file is not valid! The input file must be tab deliminated file with header.')
                    return False
        logging.error(f"The VEP file not found: {self.input}")
        return False

    def Read_database(self):
        """
        This function reads in local databases from ~/.ilyome as dictionary.
        If the provided input file is valid, databases will be loaded.
        """
        try:
            with open(self.json_hpo_GP) as T:
                self.GeneP = json.load(T)

        except Exception as e:
            logging.error(f'Failed in reading HPO database: {e}')
        try:
            with open(self.Json_omim) as T:
                self.omim = json.load(T)
        except Exception as e:
            logging.error(f'Failed in reading Omim database: {e}')
        try:
            with open(self.json_hpo_PG) as T:
                self.hpo = json.load(T)

        except Exception as e:
            logging.error(f'Failed in reading HPO database: {e}')
        try:
            with open(self.json_hpo_inheritance) as T:
                self.inheritance = json.load(T)
        except Exception as e:
            logging.error(f'Failed in reading HPO database: {e}')
        try:
            with open(self.json_string_gene_associations) as F:
                self.associations = json.load(F)
        except Exception as e:
            logging.error(f'Failed in reading HPO database: {e}')

    def Read_header(self):
        #Uploaded_variation
        with open(self.input, 'r', encoding='latin1') as T:
            for line in T:
                if line.startswith('#Uploaded_variation'):
                    for index, value in enumerate(line.strip('\n').split('\t')):
                        self.header[value] =  index
                    break
    
    def read_variants(self):
        with open(self.input, 'r', encoding='latin1') as T:
            for line in T:
                if not line.startswith('#'):
                    line = line.strip('\n').split('\t')
                    if len(self.header) == len(line):
                        #remove high allele frequency
                        max_allele = line[self.header["MAX_AF"]]
                        max_allele = float(max_allele) if max_allele != "-" else 0
                        if max_allele <= self.allele_freq:
                            self.variants.append(line)
    
    def dividecomma(self, colm):
        if ',' not in colm:
            if 'E' in colm or 'e' in colm:
                return float(colm)
            if colm == '.':
                return float(colm.replace('.', '0'))
            return float(colm.replace('-', '0'))
        return max([float(el) if el != '.' else 0 for el in colm.split(',')])

    def FilterVariants(self):
        df = pd.DataFrame(self.variants, columns=self.header)
        gnomAD_AF = (df['nomadExome_G_Exome_AF'].apply(self.dividecomma) < self.allele_freq)
        gnomAD_Homo = (df['nomadExome_nhomalt'].apply(self.dividecomma) < self.HomoCount)
        gnomad_Genome_AF = (df['nomadGenome_AF'].apply(self.dividecomma) < self.allele_freq)
        gnomad_Genome = (df['nomadGenome_nhomalt'].apply(self.dividecomma) < self.HomoCount)
        filteredDB = df.loc[ gnomAD_AF & gnomAD_Homo & gnomad_Genome & gnomad_Genome_AF ]
        self.variants = filteredDB.values.tolist()

    def flatten(self, x):
        result = []
        for el in x:
            if hasattr(el, "__iter__") and not isinstance(el, str):
                result.extend(self.flatten(el))
            else:
                result.append(el)
        return result

    def Annotate_variant(self):
        result = []
        for variant in self.variants:
            variant = {key:variant[indexx] for key,indexx in self.header.items()}
            gene = variant["SYMBOL"]
            ##
            ## Need to change main input for hpo to string not list
            variant["HPO_features"] = ", ".join(self.hpo[gene]) if gene in self.hpo.keys() else "-"
            variant["Inheritance"] = ", ".join(self.inheritance[gene]) if gene in self.inheritance.keys() else "-"
            variant["Omim"] = self.omim[gene][0] if gene in self.omim else "-"
            
            """
            check which genes have assosication and examine their phenotypes.
            gene = [list of genes with associations]
            [list of genes with associations] = [list of phenotypes from hpo]
            assosication score = like match score but for all phenotypes of associated genes
            Need to correct approved gene
            """
            if gene in self.associations.keys():
                associated_genes = self.associations[gene]["association"]
                associated_phenotypes = []
                for A_gene in associated_genes:
                    if A_gene in self.hpo.keys():
                        associated_phenotypes.append(self.hpo[A_gene])
                flat_associated_phenotypes = self.flatten(associated_phenotypes)
                variant["Associations_genes"] = ", ".join(associated_genes)
                variant["Associations_phenotypes"] = ", ".join(flat_associated_phenotypes)
            else:
                variant["Associations_genes"] = "-"
                variant["Associations_phenotypes"] = "-"
            score = Scores(variant, self.ListOfClinics, self.DisorderInheritancePedigree)
            variant["Total_allele_count"] = score.AlleleCount
            variant["Total_homos"] = score.Nhomalt
            variant["Ilyome_score"] = score.run()
            variant["Roh_number"] = self.ROHnumber
            result.append(variant)
        return result

if __name__ == "__main__":
    inn = "/home/edi/ngs/09-reports/734-G-MeErg_S23_hg38/734-G-MeErg_S23_hg38.vep"
    listOfClinic = ['Abnormality of retinal pigmentation', 'Macular dystrophy']
    a = Read_input_vep(inn, listOfClinic)
    print(a.Annotate_variant())