from readData import Read_input_vep
from filtration import Filter_write
import pandas as pd

listOfClinic = "Intellectual disability".split(",")

diseaseInheritance = "AR"

inn = "/home/edi/ngs/09-reports/733-G-OnGuz_S22_hg38/733-G-OnGuz_S22_hg38.vep"
inn = "/home/edi/ngs/10-prioritizer/feed/99012056_hg38.vep"

input_data = Read_input_vep(inn, listOfClinic, diseaseInheritance)
variants_result = input_data.Annotate_variant()
Filter_write(variants_result)