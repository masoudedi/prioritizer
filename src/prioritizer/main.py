import argparse
from cProfile import run
from .pack.config import Config
from .pack.BWA import Mapper
from .pack.Dedupe import DedupBam
from .pack.haplotypeCaller_scattering import HaplotypeCaller
from .pack.VCF import GetVCF
from .pack.Metrics import Metrics
import shutil
import timeit

parser = argparse.ArgumentParser(prog='ILYOME', description='Ilyome; a comprehensive pipeline for Exome sequencing analysis.')
parser.add_argument("--version", action="version", version="Ilyome v.0.3.1")
parser.add_argument('-i', metavar="--input", required=True, type=str, help="Path to the FASTQ files")
parser.add_argument("-o", metavar="--output", required=True, type=str, help="Output directory")
parser.add_argument('-l', metavar="--interval", default="/home/gnks/Desktop/edi/01.bedfiles/06.agilentXT-HS-V8_hg38/S33266340_Regions.interval_list",  type=str, help="Path to the interval file")
parser.add_argument('-r', metavar = '--reference', type = str, default="/home/gnks/Desktop/99.reference_genomes/hg38/Homo_sapiens_assembly38.fasta", help = 'Reference Genome in fasta format')
parser.add_argument('-f', metavar="--forks", default=45, type=int, help="Number of CPU cores to be implemented")
parser.add_argument('-pv',metavar='-padding', default=100, help="Padding value")
parser.add_argument('-m',metavar='-memory', default=200, help="maximum memory allocated")
parser.add_argument('-mx',metavar='-Max-records', default=100000000, help="Maximum recordes in RAM when Sorting and deduping BAM file")
args = parser.parse_args()

config = Config(args.i, args.o, args.l, args.r, args.pv, max_workers=args.f, memory=args.m, MAX_RECORDS_IN_RAM=args.mx)
scattered_list = config.list_scattered_intervals()

#start looping
times = []
for fastq_dict in config.fastqs:
    start = timeit.default_timer()
    if Mapper(**{"1":fastq_dict, "2":config.kwargs}).run() == 0:
        gvcf = HaplotypeCaller(ScatterIntevallist=scattered_list, **{"1":fastq_dict, "2":config.kwargs})
        vcf = GetVCF(**{"1":fastq_dict, "2":config.kwargs})
        metrics = Metrics(**{"1":fastq_dict, "2":config.kwargs}).MultiCore()
    
    stop = timeit.default_timer()
    timee = f"Time spent on {fastq_dict['rgsm']}: {((stop - start)/60)} Mins"
    times.append(timee)
    print(timee)

shutil.rmtree(config.scattered_path, ignore_errors=True)
#print time
for time in times:
    print(time)



