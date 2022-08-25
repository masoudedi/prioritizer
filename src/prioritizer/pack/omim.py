import re
import json

seprator = '/'
cardinal_seprator = '\t'
gap_filer = '.'
omimfile = "/home/edi/ngs/04.III/III/core/omim/genemap2.txt"
strin = "Incontinentia pigmenti, 308300 (3), X-linked dominant; Ectodermal dysplasia and immunodeficiency 1, 300291 (3), X-linked recessive; Immunodeficiency 33, 300636 (3), X-linked recessive"

class Omim:
    def __init__(self, omim):
        self.omim = omim
        self.bedFile = []
        self.run()
        self.writeBedFile()
    
    def run(self):
        genemap = self.openOmim()
        IlyomeWebSearchDatabase = []
        for values in genemap.values():
            values = values[0].split(";")
            for val in values:
                IlyomeWebSearchDatabase.append(val.strip())
        IlyomeWebSearchDatabase = list(set(IlyomeWebSearchDatabase))
        print("Total Omim genes ", len(genemap))
        with open("/home/edi/ngs/10-prioritizer/prioritizer/src/prioritizer/pack/data/omimdict.json", "w") as T:
            json.dump(genemap, T, indent=1)
        with open("/home/edi/ngs/10-prioritizer/prioritizer/src/prioritizer/pack/data/Webomimdict.json", "w") as T:
            json.dump(IlyomeWebSearchDatabase, T, indent=1)

    
    def Search_feature(self, search_pattern, SearchTaregt):
        pattern = re.compile(search_pattern)
        return pattern.findall(SearchTaregt)
     
    def Omim_divider(self, st):
        st = st.replace("?", "").replace("{","").replace("}","").replace("[","").replace("]","")
        phenos = st.split(";")
        res = []
        pattern = r"^\d\d\d*"
        respheno = []
        for pheno in phenos:
            res = []
            for el in pheno.split(","):
                if self.Search_feature(pattern, el.strip()):
                    pass
                else:
                    listt = {"X-linked recessive": "XLR",
                            "X-linked dominant":"XLD",
                            "Pseudoautosomal dominant": "PAD",
                            "Pseudoautosomal recessive": "PAR",
                            "X-linked": "XL",
                            "Somatic mosaicism": "SOM",
                            "Autosomal recessive": "AR",
                            "Autosomal dominant":"AD",
                            "susceptibility to": "SuT"
                            }
                    if el.strip() in listt:
                        el = listt[el.strip()]
                        res.append(el)
                    else:
                        res.append(el)
            respheno.append(" ".join(res))
        return ";".join(respheno)

    def openOmim(self):
        genemap = {}
        
        with open(self.omim) as T:
            omimdb = T.read()
        for line in omimdb.split('\n'):
            line = line.split('\t')
            if len(line) > 1:
                chro = line[0]
                start = line[1]
                end = line[2]
                genename = line[6]
                approvedgene = line[8]
                genecode = line[9]
                ensemble = line[10]
                clinic = line[12]
                # Select only approved genes for bed file and genes with phenotypes
                if approvedgene != '' and start != '0' and clinic != '':
                    clinic = self.Omim_divider(clinic)
                    phenotyps = clinic.split(';')
                    inheritence = []
                    for phenotype in phenotyps:
                        if 'AD' in phenotype:
                            inheritence.append('AD')
                        if 'AR' in phenotype:
                            inheritence.append('AR')
                        elif 'XD' in phenotype:
                            inheritence.append('XD')
                        elif 'XR' in phenotype:
                            inheritence.append('XR')
                        elif 'SuT' in phenotype:
                            inheritence.append('SuT')
                        elif 'YL' in phenotype:
                            inheritence.append('YL')
                        elif 'PAR' in phenotype:
                            inheritence.append('PAR')
                        elif 'PAD' in phenotype:
                            inheritence.append('PAD')
                        elif 'SOM' in phenotype:
                            inheritence.append('SOM')
                        else:
                            inheritence.append('.')
                    if len(inheritence) > 1:
                        inheritence = set(inheritence)
                        # for n in inheritence:
                        if '.' in inheritence:
                            inheritence.remove('.')
                        if 'XL' in inheritence:
                            inheritence.remove('XL')

                    inheritence = seprator.join(inheritence)

                    # omimbedfile.write(f'{chro}\t{start}\t{end}\t{approvedgene}\n')
                    # , ensemble, genename, chro, start, end,
                    genemap[approvedgene] = [clinic, inheritence]
                    if not chro.startswith("#"):
                        self.bedFile.append(f"{chro}\t{start}\t{end}\t{clinic}&{inheritence}\n")


        return genemap
    
    def writeBedFile(self):
        with open("OmimBedFile.bed", "w") as T:
            for el in self.bedFile:
                T.write(el)
Omim(omimfile)