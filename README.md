## Upper Quantile Normalization Code 
[Take a look here](https://github.com/ucscXena/python-scripts/issues/1)

#### Installation 
    pip install -r requirements

#### Run
    python upper_quantile_normalize.py --genefile GENEFILE --genmatfile GENEMATRIXFLE

This generates two files `UQmatrix` and `UQcolumns` to the same directory.

#### Sample files 

For GENEFILE - use `genefile` from the code repo. It contains the genes from the first column of [this]( https://github.com/ucscXena/python-scripts/blob/master/geneLists/GTEX_TCGA_genes).
For GENEMATRIXFLE - download the dataset from  [here](http://ec2-52-23-185-93.compute-1.amazonaws.com/datapages/?dataset=TCGA.COAD.sampleMap/HiSeqV2&host=https://tcga.xenahubs.net)

**Note:** The script is Python3 compatible :)

### Using Docker 
#### Build 
    docker -t ucsc-python .
    
#### Run
    docker run ucsc-python
    

#### Getting the results 
The output files(`UQmatrix` and `UQcolumns`) are written to `/output` in the container. 
To get them in any directory on the host machine, use -

    docker run -v /absolute/path/to/directory:/output ucsc-python

For example, to get the results in the current directory -

    docker run -v $(pwd):/output ucsc-python
  
