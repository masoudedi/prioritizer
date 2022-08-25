from .readData import Read_input_vep
from .filtration import Filter_write
from .readROH import ReadROH

class Anne:
    def __init__(self, VepFile, ROhFile, listoflinic, outputclinic):
        Roh_number = ReadROH(ROhFile)
        Di = "AR" if Roh_number >= 50 else "AD"
        #clinic
        input_data = Read_input_vep(VepFile, listoflinic, Di, Roh_number)
        variants_result = input_data.Annotate_variant()
        Filter_write(variants_result, outputclinic)
        
        #genetic
        # input_data = Read_input_vep(VepFile, [], Di, Roh_number)
        # variants_result = input_data.Annotate_variant()
        # Filter_write(variants_result, outputgenetic)