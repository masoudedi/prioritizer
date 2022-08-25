import pandas as pd
class Filter_write:
    def __init__(self, Dataframe, output):
        self.df = pd.DataFrame(Dataframe)
        self.output = output
        self.run()

    def getConsequence(self, inn):
        consequence = ['transcript_ablation', 'splice_acceptor_variant', 
                        'splice_donor_variant', 'stop_gained', 'frameshift_variant',
                        'stop_lost', 'start_lost', 'transcript_amplification',
                        'inframe_insertion','inframe_deletion','missense_variant',
                        'protein_altering_variant','splice_region_variant','incomplete_terminal_codon_variant',
                        'start_retained_variant','stop_retained_variant']
        for el in inn.split(','):
            if el in consequence:
                return True
        return False
    
    def clinvar(self, x):
        for el in x.split(","):
            if el == "pathogenic" or el == "likely_pathogenic" or el == "pathogenic/likely_pathogenic":
                return True
        return False

    def run(self):
        #recolumn
        firstColumns = ["Omim","Inheritance", "HPO_features", "Ilyome_score", "ZYG", "#Uploaded_variation", "SYMBOL", "HGVSc", "HGVSp", "ClinVar_CLNSIG", "Total_allele_count", "Total_homos", "Roh_number", "Associations_genes", "Associations_phenotypes"]
        columns = firstColumns + [cl for cl in self.df.columns.values if cl not in firstColumns]
        df = self.df.reindex(columns=columns)
        df.sort_values(by=['Ilyome_score'], inplace=True, ascending=False)
        pathogenic_df = df.loc[(df['CLIN_SIG'].apply(self.clinvar) | df['ClinVar_CLNSIG'].apply(self.clinvar))]
        consequence = df['Consequence'].apply(self.getConsequence)
        df = df.loc[consequence]
        df1 = df.loc[consequence]
        Most_sever_df = df1.head(20)
        pheno = (df['HPO_features'].apply(lambda value: value != '-'))
        Genes_without_pheno = df.loc[~pheno]
        Homos = df.loc[pheno & df['ZYG'].str.contains('HOM')]
        Heteros_all = df.loc[pheno & df['ZYG'].str.contains('HET')]
        genes = Heteros_all['SYMBOL'].values.tolist()
        def duplicate_genes(hetro_genes):
            return list(set([gene for gene in hetro_genes if hetro_genes.count(gene) > 1]))
        duplicate_list = duplicate_genes(genes)

        compund_df = Heteros_all.loc[Heteros_all['SYMBOL'].isin(duplicate_list)]
        other_HET_all = Heteros_all.loc[~Heteros_all['SYMBOL'].isin(duplicate_list)]
        Het_AD_df = other_HET_all.loc[other_HET_all['Inheritance'].str.contains('AD', na=False) | other_HET_all['Inheritance'].str.contains('XD', na=False)]
        Het_AR_df = other_HET_all.loc[~ other_HET_all['Inheritance'].str.contains('AD', na=False) | other_HET_all['Inheritance'].str.contains('XD', na=False)]

        with pd.ExcelWriter(f'{self.output}_Ilyome.xlsx') as writer:
            Most_sever_df.to_excel(writer, sheet_name="Candidates", index=False, header=True)
            Homos.to_excel(writer, sheet_name='HOMO', index=False, header=True)
            compund_df.to_excel(writer, sheet_name='Compound_Heterozygous', index=False, header=True)
            Het_AD_df.to_excel(writer, sheet_name='HET_AD', index=False, header=True)
            Het_AR_df.to_excel(writer, sheet_name='HET_AR', index=False, header=True)
            pathogenic_df.to_excel(writer, sheet_name='Pathogenic', index=False, header=True)
            Genes_without_pheno.to_excel(writer, sheet_name='Unknown_Genes', index=False, header=True)