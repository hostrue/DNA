# DNA
DNA search contest

# Step to convert a DNA to Time series,

1. Put the input file and the editChromosome.exe in the same folder
2. Run editChromosome.exe,
	command: editChromosome.exe DNAfile timeseriesfile
	example: editChromosome.exe chr.fa ts.txt

3. Sample input: chr.fa, contains first 1000bp of chimp chromosomose 2a

4. Time series is stored in binary format, data type is integer.

4. You can use readDNAts.m to read time series from the binary file 'ts.txt'.