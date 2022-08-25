import re 

class Scores:
    def __init__(self, variant, listofClinics, DisorderInheritancePedigree):
        self.variant = variant
        self.diseaseInheritance = DisorderInheritancePedigree
        ## Need to check if we have two genes or not
        ## Geneselect function from Ilyome class
        self.gene = variant["SYMBOL"]
        self.AlleleCount = self.AlleleCount(variant)
        self.Nhomalt = self.Nhomalt(variant)
        self.HPO_features = variant["HPO_features"].split(", ")
        self.Associations_phenotypes = variant["Associations_phenotypes"].split(", ")
        self.gene_inheritances = variant["Inheritance"].split(", ")
        self.phenotype = listofClinics

    def run(self):
        MS = self.match_score() 
        AS = self.Assosication_score() 
        IS = self.impact_score() 
        IHS = self.inheritance_score() 
        CS = self.clinvar_score() 
        ## Need to check if IHS is suitable then consider clinvar
        ##If MS is zero give two times wieght to AS.
        AS = AS * 2 if MS == 0 else AS
        # if self.gene ==  "PCARE" or self.gene == "PHKA2":
        #     print(MS, CS, self.HPO_features, self.phenotype)
        return MS + AS + IS + IHS + CS 
    
    def Search_feature(self, search_item, HPO_list):
        pattern = re.compile(search_item.strip().lower())
        for feature in HPO_list:
            matches = pattern.findall(feature.lower())
            if matches:
                return 1
        return 0 

    def match_score(self):
        sc = sum([self.Search_feature(feature, self.HPO_features) for feature in self.phenotype])
        return (sc/len(self.phenotype)) * 3 if sc else 0

    def Assosication_score(self):
        sc = sum([self.Search_feature(feature, self.Associations_phenotypes) for feature in self.phenotype ])
        return (sc/len(self.phenotype)) * 1.5 if sc else 0

    def impact_score(self):
        _category = {'HIGH': 3, 'MODERATE': 1.5, 'MODIFIER': 0, 'LOW': -3}
        try:
            return _category[self.variant["IMPACT"]]
        except Exception:
            return 0

    def inheritance_score(self):
        variant_zygosity = self.variant["ZYG"] 
        zygosity_scores = [-3,]
        for gene_inheritance in self.gene_inheritances:
            if self.Nhomalt <= 10:
                if variant_zygosity == "HOM" and gene_inheritance == "AR" and self.diseaseInheritance == "AR":
                    zygosity_scores.append(3)
                if variant_zygosity == "HOM" and gene_inheritance == "XLR" and self.diseaseInheritance == "XLR":
                    zygosity_scores.append(3)

            if self.AlleleCount == 0:
                if variant_zygosity == "HET" and gene_inheritance == "AD" and self.diseaseInheritance == "AD":
                    zygosity_scores.append(3)
                if variant_zygosity == "HET" and gene_inheritance == "XLD" and self.diseaseInheritance == "XLD":
                    zygosity_scores.append(3)
            
            elif self.AlleleCount <= 20:
                if variant_zygosity == "HET" and gene_inheritance == "AD" and self.diseaseInheritance == "AD":
                    zygosity_scores.append(0)
                if variant_zygosity == "HET" and gene_inheritance == "XLD" and self.diseaseInheritance == "XLD":
                    zygosity_scores.append(0)

        return max(zygosity_scores)

    def clinvar_score(self):

        clinvar1 = self.variant["CLIN_SIG"].split(",")
        clinvar2 = self.variant["ClinVar_CLNSIG"].split(",")
        clinvar = clinvar1 + clinvar2
        pathogenics  = ["pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"]
        benigns = ["benign", "likely_benign", "benign/likely_benign"]
        for feature in clinvar:
            if feature.lower() in pathogenics:
                return 3
        for feature in clinvar:
            if feature.lower() in benigns:
                return -3
        return 0
    
    def dividecomma(self, colm):
        if ',' not in colm:
            if 'E' in colm or 'e' in colm:
                return float(colm)
            if colm == '.':
                return float(colm.replace('.', '0'))
            return float(colm.replace('-', '0'))
        return max([float(el) if el != '.' else 0 for el in colm.split(',')])

    def AlleleCount (self, variant):
        nomadExome = variant["nomadExome_AC"]
        nomadGenomeAC = variant["nomadGenome_AC"]
        return self.dividecomma(nomadExome) + self.dividecomma(nomadGenomeAC)
    
    def Nhomalt (self, variant):
        nomadExome = variant["nomadExome_nhomalt"]
        nomadGenomeAC= variant["nomadGenome_nhomalt"]
        return self.dividecomma(nomadExome) + self.dividecomma(nomadGenomeAC)
        

