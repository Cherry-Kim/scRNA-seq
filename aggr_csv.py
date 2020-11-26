import string,sys,os

fpout=open("aggr_N.csv","w")
fpout.write("library_id"+","+"molecule_h5"+","+"Status""\n")
fpout1=open("aggr_T.csv","w")
fpout1.write("library_id"+","+"molecule_h5"+","+"Status""\n")
path='/Cancer/Colon/scRNA-seq/FASTQ/count/'
file_list=os.listdir(path)
sample=[file for file in file_list if file.startswith('PM-PS')]
for i in sample:
	if "N" in i:
		fpout.write(i+","+path+i+"/outs/molecule_info.h5"+","+"Normal"+"\n")
	if "T" in i:
		fpout1.write(i+","+path+i+"/outs/molecule_info.h5"+","+"Tumor"+"\n")

fpout.close() 
fpout1.close() 
