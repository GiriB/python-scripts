import argparse,csv
import numpy as np
import os
'''
Finds the norm_count in the formula :
	log2(norm_count + pseudo_count) = Y

Vectorized to work with numpy arrays
'''
@np.vectorize
def get_norm_count(Y,pseudo_count=1):
	# norm_count = 2^Y - 1
	norm_count = np.exp2(Y) - pseudo_count
	return norm_count


'''
Returns log2((X / UQ) * 1e6 + 1)
'''
@np.vectorize
def transform_log2(X,UQ,pseudo_count=1):
	return np.log2(X/UQ * 1e6 + pseudo_count)


'''
Reads a tabular separated file.
Returns row_names (the first column in the file),
		column_names (the first row in the file),
		data (as numpy.ndarray)
'''
def get_data_from_tsv(filename):
	data = []
	row_names = []
	column_names = []
	with open(filename,'r') as tsvf:
		# Read the tabular separated values 
		tsvf = csv.reader(tsvf, delimiter='\t')
		
		# Get the headers from the first row 
		# but skip the first column (sample)
		column_names = next(tsvf)[1:]
		
		for row in tsvf:
			# The first column has Gene names 
			row_names.append(row[0])
			# Read values as float, append to data
			data.append(list(map(float,row[1:])))

	# Make data numpy.ndarray
	data = np.asarray(data)
	return row_names,column_names,data


'''
Returns rows corresponding to the genes in genelist 
'''
def extract_for_genes(row_names, data, genelist):
	
	# row_map is a list of gene names that we found
	# The index of the gene represents the index of its 
	# row in mat
	row_map = []

	temp = []
	# Make a matrix only with genes in the genelist
	for i in range(len(row_names)):
		if row_names[i] in genelist:
			temp.append(i)
			row_map.append(row_names[i])
	
	mat = data[temp,:] 
	return mat,row_map

'''
Write UQ values for each column to a filename
'''
def write_uqfile(filename,uq,cols):
	assert len(uq) == len(cols)
	with open(filename,'w') as f:
		for i in range(len(cols)):
			f.write(str(cols[i])+"\t"+str(uq[i])+"\n")

if __name__ == '__main__':

	# Parse the inputs
	parser = argparse.ArgumentParser()
	# Required params
	parser.add_argument("--genmatfile", help="Xena Genomic Matrix File.",required=True)
	parser.add_argument("--genefile", help="Gene File. (Each row is a Gene name)",required=True)
	args = parser.parse_args()

	GENOMICMATRIX_FILE = args.genmatfile
	GENE_FILE = args.genefile

	# Get the data from xena genomicMatrix file
	rows,cols,data = get_data_from_tsv(GENOMICMATRIX_FILE)
	
	# Get the name of the genes from gene file (each row is a gene)
	genelist = []
	with open(GENE_FILE,'r') as genefile:
		for i in genefile:
			genelist.append(i.strip());
	# Converting to set results in 2x better performance
	genelist = set(genelist)

	# Extract rows for genes in genefile
	mat,row_map = extract_for_genes(rows, data, genelist)
	
	# UQ values for each column
	UQ_values = []

	# for each column in mat
	for i in range(mat.shape[1]):
		column = mat[:,i:i+1]

		# Get back norm_count
		column = get_norm_count(column)

		# Get 75 percentile value for the column
		UQ = np.percentile(column,75)
		UQ_values.append(UQ)
			
		# log2(X/UQ * 1e6 + 1)
		column = transform_log2(column,UQ)

		# Save the change back to the matrix
		mat[:,i:i+1] = column
		
	#Save the final matrix and UQ per column in a file 
	print("Writing to UQMatrix ...")
	np.savetxt('UQmatrix', mat)
	print("Writing to UQColumns ...")
	write_uqfile('UQcolumns',UQ_values,cols)