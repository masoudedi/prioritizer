import json
import xmltodict
import requests
import os
from pathlib import Path
from config import Config

class database:
    def __init__(self, config):
        self.conf = Config
    def download(self):
        # self.down(self.conf.url_orpha_GARD, self.conf.out_orpha_GARD)
        # self.down(self.conf.url_orpha_PARD, self.conf.out_orpha_PARD)
        self.down(self.conf.url_HPO_GP, self.conf.out_hpo_PG)
        self.down(self.conf.url_HPO_PG, self.conf.out_hpo_GP)
        self.down(self.conf.url_clinvar, self.conf.clinvar_file)
        self.down(self.conf.url_clinvar_tbi, self.conf.clinvar_tbi_file)
         
    def update(self):
        print("Updating the database ...")
        self.updat(self.conf.out_hpo_GP, self.conf.json_hpo_GP)
        self.updat(self.conf.out_hpo_PG, self.conf.json_hpo_PG)
        #self.xmlToJson(self.conf.out_orpha_GARD, self.conf.json_orpha_GARD)
        #self.xmlToJson(self.conf.out_orpha_PARD, self.conf.json_orpha_PARD)
    
    def down(self, url, outputname):
        
        file = requests.get(url)
        print("Downloading file: ", url.split("/")[-1])
        print("Server response: ",file.status_code)
        with open(os.path.join(self.conf.mainpack, outputname), 'wb') as T:
            T.write(file.content)
    
    def updat(self, input, output):
        result = {}
        with open(input) as T:
            for line in T:
                line = line.strip('\n').split('\t')
                if len(line) > 1:
                    if not line[1] in result:
                        result[line[1]] = [line[3]]
                    else:
                        result[line[1]].append(line[3])
        with open(output, 'w') as T:
            json.dump(result, T, indent=1)
        return result
    
    def Create_inheritance(self):
        with open(self.conf.json_hpo_PG) as T:
            hpo = json.load(T)
        hpo_inheritance = {}
        for gene, value in hpo.items():
            inheritances = []
            for pheno in value:
                if "inheritance" in pheno:
                    inheritances.append(pheno)
            inheritance = []
            for inh in inheritances:
                if "Autosomal recessive inheritance" == inh:
                    inheritance.append("AR")
                if "Autosomal dominant inheritance" == inh:
                    inheritance.append("AD")
                if "X-linked recessive inheritance" == inh:
                    inheritance.append("XLR")
                if "X-linked dominant inheritance" == inh:
                    inheritance.append("XLD")
                if "X-linked inheritance" == inh:
                    inheritance.append("XL")
                if "Y-linked inheritance" == inh:
                    inheritance.append("YL")
                if "Oligogenic inheritance" == inh:
                    inheritance.append("Oligogenic")
                if "Multifactorial inheritance" == inh:
                    inheritance.append("MF")
                if "Mitochondrial inheritance" == inh:
                    inheritance.append("MT")
                if "Digenic inheritance" == inh:
                    inheritance.append("Digenic")
            
            inheritance = list(set(inheritance))
            hpo_inheritance[gene] = inheritance
        with open(self.conf.json_hpo_inheritance, 'w') as T:
            json.dump(hpo_inheritance, T, indent=1) 
    
    def xmlToJson(self, file, output):
        with open(file, encoding="ISO-8859-1") as T:
            data_dict = xmltodict.parse(T.read())
            data_dict = data_dict['JDBOR']['DisorderList']['Disorder']
        result = {}
        for key in data_dict:
            OrphaCode = key['OrphaCode']
            DisorderName = key['Name']['#text']
            DisorderType = key['DisorderType']['Name']['#text']
            DisorderGroup = key['DisorderGroup']['Name']['#text']
            DisAsso = key['DisorderGeneAssociationList']['DisorderGeneAssociation']
            GeneName = DisAsso['Gene']['Name']['#text']
            GeneSymbol = DisAsso['Gene']['Symbol']
            DisorderGAT = DisAsso['DisorderGeneAssociationType']['Name']['#text']

        with open(output, 'w') as T:
            json.dump(data_dict, T, indent=1)


## READ XML files
def xmlToJson(file, output):
    with open(file, encoding="ISO-8859-1") as T:
        data_dict = xmltodict.parse(T.read())
    with open(output, 'w') as T:
        json.dump(data_dict, T, indent=1)
    # disorders = {}
    # with open(output) as T:
    #     data = json.load(T)
    # for key in data['JDBOR']['DisorderList']['Disorder']:
    #     disorders[key['OrphaCode']] = {}
    #     disorders[key['OrphaCode']] = []
    #     print(key['Name']['#text'])
    #     break
if __name__ == "__main__":
    conf = Config()
    x = database(conf)
    x.download()
    x.update()
    x.Create_inheritance()
    
