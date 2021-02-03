import sys,string
import os
from rpy2 import robjects as ro
r = ro.r

def STEP1_COUNT(ref):
	file_list=os.listdir('/home/hykim/FASTQ/')
	a=[file for file in file_list if file.endswith('PS')]
	for i in a:
       	 os.system('/home/hykim/SC/cellranger-4.0.0/cellranger count --id='+i+' --sample='+i+' --fastqs='+path+i+' --jobmode=local --transcriptome='+ref)

def MAKE_AGGRINPUT(path):
	fpout=open("aggr_N.csv","w")
	fpout.write("library_id"+","+"molecule_h5"+","+"Status""\n")
	fpout1=open("aggr_T.csv","w")
	fpout1.write("library_id"+","+"molecule_h5"+","+"Status""\n")

	file_list=os.listdir(path)
	sample=[file for file in file_list if file.startswith('PM-PS')]
	for i in sample:
       		if "N" in i:
               		fpout.write(i+","+path+i+"/outs/molecule_info.h5"+","+"Normal"+"\n")
	        if "T" in i:
        	        fpout1.write(i+","+path+i+"/outs/molecule_info.h5"+","+"Tumor"+"\n")
	fpout.close()
	fpout1.close()

def STEP2_AGGR():
	os.system('/home/hykim/cellranger-4.0.0/bin/cellranger aggr --id=colon_N --csv=aggr_N.csv --normalize=mapped')
	os.system('/home/hykim/cellranger-4.0.0/bin/cellranger aggr --id=colon_T --csv=aggr_T.csv --normalize=mapped')

def STEP3_SEURAT():
	r.source("2-1.seurat.r")
	r.source("2-2.SingleR.r")
	r.source("2-2.seurat.r")

def main():
	ref='/home/hykim/refdata-gex-GRCh38-2020-A'
	path='/home/hykim/'
	STEP1_COUNT(ref)
	MAKE_AGGRINPUT(path)
	STEP2_AGGR()
	STEP3_SEURAT()
main()
