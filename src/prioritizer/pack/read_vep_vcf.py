import gzip

class vep_vcf_input:
    # class variables:
    seprator = '\t'
    gap_filer = ""

    def __init__(self, vcf):
        self.vcf = vcf
        self.variants = []
        self.vcfColumns = None
        self.ColumnsLength = 0
        self.info = []
        self.CSQ = None
        self.Format = []
        self.readHeader()
        self.readVariants()

    def readHeader(self):
        with gzip.open(self.vcf, 'rb') as T:
            for line in T:
                line = line.decode().strip('\n')
                if line.startswith('#'):
                    if line.startswith("##INFO="):
                        self.info.append(line.split(",")[0].split("=<ID=")[-1])
                    if line.startswith("##INFO=<ID=CSQ"):
                        self.CSQ = line.split("Format:")[-1].strip().replace('"', "").replace(">","").split("|")
                    if line.startswith("#CHROM"):
                        self.vcfColumns = line.split(self.seprator)
                        self.ColumnsLength += len(self.vcfColumns)
                    if line.startswith("##FORMAT"):
                        self.Format.append(line.split(",")[0].split("=<ID=")[-1])
                else:
                    break
    
    def readVariants(self):
        with gzip.open(self.vcf, 'rb') as T:
            for line in T:
                line = line.decode().strip('\n')
                if not line.startswith('#'):
                    line = line.split(self.seprator)
                    resDict = {self.vcfColumns[0]:line[0],
                                self.vcfColumns[1]:line[1],
                                self.vcfColumns[2]:line[2],
                                self.vcfColumns[3]:line[3],
                                self.vcfColumns[4]:line[4],
                                self.vcfColumns[5]:line[5],
                                self.vcfColumns[6]:line[6],
                                }
                    infoList = [i.split("=") for i in line[7].split(";")]
                    infoDict = {i[0]:i[1] for i in infoList if len(i) == 2}

                    for infoHeader in self.info:
                        if infoHeader == "CSQ":
                            csq = infoDict["CSQ"].split("|")
                            for index,csqHeader in enumerate(self.CSQ):
                                resDict[csqHeader] = csq[index]
                        else:
                            try:
                                resDict[infoHeader] = infoDict[infoHeader]
                            except Exception:
                                 resDict[infoHeader] = self.gap_filer

                    formatDict = {i[0]:i[1] for i in zip(line[8].split(":"), line[9].split(":"))}
                    for formatHeader in self.Format:
                        if formatHeader == "AD":
                            try:
                               fractions = [int(i) for i in formatDict[formatHeader].split(",")]
                               resDict["AD"] = fractions
                               resDict["alleleFraction"] = fractions[1]/sum(fractions) * 100
                            except Exception:
                                resDict["alleleFraction"] = self.gap_filer
                                resDict["AD"] = ["", ""]
                        else:
                            try:
                                resDict[formatHeader] = formatDict[formatHeader]
                            except Exception:
                                resDict[formatHeader] = self.gap_filer
                    self.variants.append(resDict)

if __name__ == "__main__":
    vcf = "/home/edi/ngs/storage/qnas/edi/test.vcf.gz"
    x = vep_vcf_input(vcf)

    with open("testVep.json", "w") as T:
        import json
        json.dump(x.variants, T, indent=1)