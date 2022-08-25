import json
file = "/home/edi/.ilyome/phenotype_to_genes.json"


with open(file) as T:
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
            inheritance.append("Multifactorial")
        if "Mitochondrial inheritance" == inh:
            inheritance.append("MT")
        if "Digenic inheritance" == inh:
            inheritance.append("Digenic")
    
    inheritance = list(set(inheritance))
    hpo_inheritance[gene] = inheritance
