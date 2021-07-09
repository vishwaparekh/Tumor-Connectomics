# Tumor Connectomics

The function getTCF computes full image tumor connectomics for the input single or multiparametric imaging data. 

The getTCF function requires the following input variables:
1. img: Input single (mxn) or multiparametric (mxnxo) radiological data
2. mask: The ROI to compute tumor connectomics
3. t: The threshold value (range: 0-1) to indicate the neighborhood size

The output to the getTCF function is a structure data type consisting of the following fields:
1. degreeCentrality
2. betweennessCentrality
3. eigenvectorCentrality
4. nodeStrength
5. clusteringCoefficient
6. NumConnComponents
7. avgPathLength

Each of the fields would be of size mxn pixel-wise colorcoded with the graph metrics.
