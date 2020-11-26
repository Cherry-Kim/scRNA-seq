import string,sys,glob,os

path='/home/hykim/Project/Cancer/Colon/scRNA-seq/FASTQ/FASTQ/'
ref='/home/hykim/SC/refdata-gex-GRCh38-2020-A'
#file_list=os.listdir(path)
#a=[file for file in file_list if file.endswith('PM-PS')]
#for i in a:
#	os.system('/home/hykim/SC/cellranger-4.0.0/cellranger count --id='+i+' --sample='+i+' --fastqs='+path+i+' --jobmode=local --transcriptome='+ref)

os.system('/home/hykim/SC/cellranger-4.0.0/cellranger aggr --id=colon_T --csv=aggr_T.csv --normalize=mapped')
os.system('/home/hykim/SC/cellranger-4.0.0/cellranger aggr --id=colon_N --csv=aggr_N.csv --normalize=mapped')


