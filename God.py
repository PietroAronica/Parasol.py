#!/usr/bin/python

import sys, getopt

# Define input files and output files 
def main(argv):
	inputfile= ''
	outputfile= ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print 'God.py -i <inputfile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Usage: God.py -i <inputfile> -o <outputfile>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
	return inputfile, outputfile
if __name__ == "__main__":
	main(sys.argv[1:])
# Define dictionary of input files
inf = {}
# Load up input file
f = open(main(sys.argv[1:])[0])
# Parse input, assign values to variables
data = f.readlines()
for line in data:
    key, value = line.split("=")
    inf[key.strip()] = value.strip()
f.close()
