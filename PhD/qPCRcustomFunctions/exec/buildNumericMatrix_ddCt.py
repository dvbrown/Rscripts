#!/usr/bin/python2.7

import argparse, csv

def readAfile(filenameString):
    'Reads the input file into a dictionary object where the key is the first row'
    fileA = open(filenameString, 'U')
    inputA = csv.reader(fileA, delimiter='\t')
    data = []
    for row in inputA:
        data.append(row)
    fileA.close()
    return data

parser = argparse.ArgumentParser(description="Reads a file with ddCT values and changes it to a numeric matrix for statistical testing")
parser.add_argument('-i', '--inputData', required=True, help='The file containing ddCt data')
parser.add_argument('-g', '--geneColumn', required=True, help='The column number containing the gene names')
parser.add_argument('-o', '--originColumn', required=True, help='The column number containing the origin names')
parser.add_argument('-d', '--ddCtColumn', required=True, help='The column number containing the ddCt values')
args = parser.parse_args()

# Set the type of the incoming arguments
g = int(args.geneColumn)
o = int(args.originColumn)
d = int(args.ddCtColumn)

# Intialise a dictionary
dic = {}
inputData = readAfile(args.inputData)
# Remove the headers
firstRow = inputData[0]
inputData = inputData[1:]

# Get the gene names
header = ['origin']
for line in inputData:
    genes = line[g]
    if genes in header:
        pass
    else:
        header.append(genes)

for line in inputData:
    origin = line[o]
    ddCt = line[d]
    # Check membership in dictionary. If present then append, otherwise add new value
    if origin in dic.keys():
        dic[origin].append(ddCt)
    else:
        dic[origin] = []
        dic[origin].append(ddCt)
    
print '\t'.join(header)
    
for k, v in dic.items():
    print k + '\t' + '\t'.join(v)