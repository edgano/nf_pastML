#!/usr/bin/env python
import sys
from Bio import SeqIO

delimiter= ","
gap="-"

file_name = sys.argv[1]
tax = sys.argv[2]
outfile= tax+".csv"

resultFile= open(outfile,"w+")
commandLine= open("command.txt","w+")

fullFile = list(SeqIO.parse(file_name, "fasta"))
#   print header
    # id, ...loop to len(record)
#sys.stdout.write("id%s" % (delimiter))
resultFile.write("id%s" % (delimiter))

for x in range(1,len(fullFile[1].seq)):
    #sys.stdout.write("c%s%s" % (x,delimiter))
    resultFile.write("c%s%s" % (x,delimiter))
    commandLine.write("c%s " % (x))
#sys.stdout.write("\n")
resultFile.write("\n")
for record in fullFile:
    #sys.stdout.write("%s%s" % (record.id,delimiter))
    resultFile.write("%s%s" % (record.id,delimiter))
#   loop seq by seq to populate
#   loop throught seq to check if it's a residue or a gap
    for base in record.seq:
        # if base == gap
        #   print 0 + delimiter
        #else                   #go without indexes!! just printing column by column... idk if its to risky
        #   print 1 + delimiter
        if (base==gap):
            #sys.stdout.write("0%s" % (delimiter))
            resultFile.write("0%s" % (delimiter))
        else:
            #sys.stdout.write("1%s" % (delimiter))
            resultFile.write("1%s" % (delimiter))
    #print newLine
    #sys.stdout.write("\n")
    resultFile.write("\n")

resultFile.close() 
commandLine.close()
