from pack.anne import Anne
import os 


mainout = "/home/edi/ngs/09-reports"
mainIn = "/home/edi/ngs/storage/qnas/04.annotations/819-G-MuErd_hg38"

clinic = ["Seizure", "Intellectual disability", "Global developmental delay"]
allfiles = [i for i in os.listdir(mainIn)]
vepFile = os.path.join(mainIn, [i for i in allfiles if i.endswith(".vep")][0])
Rohin = os.path.join(mainIn, [i for i in allfiles if i.endswith(".HomRegions.tsv")][0])
base = str(os.path.basename(vepFile)).split(".")[0]
output_clinic = os.path.join(mainout, f"{base}_clinic")
Anne(vepFile, Rohin, clinic, output_clinic)