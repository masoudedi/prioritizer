from pack.anne import Anne
import os 
file = "/home/edi/ngs/10-prioritizer/Result/Yavuz/Ilyome_patients_all.xlsx"

import pandas as pd 
df = pd.read_excel(file, sheet_name=0, dtype= str)
columns = df.columns.values

val = df.values.tolist()
for i in val:
    name = str(i[0]).strip("").split("\n")[0]
    clinic = [i.strip() for i in str(i[1]).strip("").split(",")]
    
    #start Analysis
    mainout = "/home/edi/ngs/10-prioritizer/Result/Yavuz/run02-03"
    output_clinic = os.path.join(mainout, f"{name}_clinic")
    output_genetic = os.path.join(mainout, f"{name}_genetic")
    
    mainIn = "/home/edi/ngs/10-prioritizer/Result/00-Vep_homo/run02-03"
    vepFile = os.path.join(mainIn, f"{name}.vep")
    Rohin = os.path.join(mainIn, f"{name}.HomRegions.tsv")
    
    if os.path.isfile(vepFile) and os.path.isfile(Rohin):
        Anne(vepFile, Rohin, clinic, output_clinic, output_genetic)
    else:
        print(f"{vepFile} or {Rohin} was not found")


