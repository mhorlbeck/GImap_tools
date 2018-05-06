import os
import gzip
import sys
from Bio import Seq
import multiprocessing
import itertools
import argparse
import pandas as pd

def writeToCounts(fileTup):
	print '\t'.join(fileTup[:-2])
	sys.stdout.flush()
	
	mismatchDicts = fileTup[-2]
	testRun = fileTup[-1]
	
	statsCounts = {'A sgRNA not mapped': 0,
	'B sgRNA not mapped': 0, 
	'A barcode not mapped': 0,
	'B barcode not mapped': 0,
	'A sgRNA multiple mappings': 0,
	'B sgRNA multiple mappings': 0, 
	'A barcode multiple mappings': 0,
	'B barcode multiple mappings': 0,
	'All sgRNAs and barcodes uniquely map': 0,
	'A sgRNA and barcode do not match': 0,
	'B sgRNA and barcode do not match': 0,
	'A and B sgRNAs and barcodes match': 0}
	
	pairCounts_barcodes = dict()
	pairCounts_sgRNAs = dict()
	pairCounts_triple = dict()
	
	with gzip.open(fileTup[0]) as infile_r1:
		with gzip.open(fileTup[1]) as infile_r2:
			with gzip.open(fileTup[2]) as infile_r3:
					for i, (r1,r2,r3) in enumerate(itertools.izip(infile_r1, infile_r2, infile_r3)):
						if i%4 == 1:
							protospacer_a_r1 = matchBarcode(mismatchDicts['protospacer_a_r1'], r1[:19].strip())#r1[1:20].strip()) for earlier rapid runs
							
							protospacer_b_r2 = matchBarcode(mismatchDicts['protospacer_b_r2'], r2[:19].strip())#r2[1:20].strip())

							barcodes = str(Seq.Seq(r3[:38]).reverse_complement())#[1:38+1]).reverse_complement()) for earlier rapid runs
							barcode_a_r3 = matchBarcode(mismatchDicts['barcode_a_r3'], barcodes[:16])
							barcode_b_r3 = matchBarcode(mismatchDicts['barcode_b_r3'], barcodes[16+6:].strip())
							
							if protospacer_a_r1 == 'none' \
								or protospacer_b_r2 == 'none' \
								or barcode_a_r3 == 'none' \
								or barcode_b_r3 == 'none' \
								or protospacer_a_r1 == 'multiple' \
								or protospacer_b_r2 == 'multiple' \
								or barcode_a_r3 == 'multiple' \
								or barcode_b_r3 == 'multiple':
								
								if protospacer_a_r1 == 'none':
									statsCounts['A sgRNA not mapped'] += 1
								if protospacer_a_r1 == 'none':
									statsCounts['B sgRNA not mapped'] += 1
								if barcode_a_r3 == 'none':
									statsCounts['A barcode not mapped'] += 1
								if barcode_b_r3 == 'none':
									statsCounts['B barcode not mapped'] += 1
								if protospacer_a_r1 == 'multiple':
									statsCounts['A sgRNA multiple mappings'] += 1
								if protospacer_a_r1 == 'multiple':
									statsCounts['B sgRNA multiple mappings'] += 1
								if barcode_a_r3 == 'multiple':
									statsCounts['A barcode multiple mappings'] += 1
								if barcode_b_r3 == 'multiple':
									statsCounts['B barcode multiple mappings'] += 1

								#sgRNA sequencing only--moved from following else statement to include scenarios where the barcodes do not map but protospacers do
								if protospacer_a_r1 != 'none' and protospacer_a_r1 != 'multiple' and protospacer_b_r2 != 'none' and protospacer_b_r2 != 'multiple':
									combinedSgId = protospacer_a_r1 + '++' + protospacer_b_r2
									if combinedSgId not in pairCounts_sgRNAs:
										pairCounts_sgRNAs[combinedSgId] = 0

									pairCounts_sgRNAs[combinedSgId] += 1



								#barcode sequencing only--moved from following else statement to include scenarios where the barcodes do not map but protospacers do
								if barcode_a_r3 != 'none' and barcode_a_r3 != 'multiple' and barcode_b_r3 != 'none' and barcode_b_r3 != 'multiple':
									combinedSgId = barcode_a_r3 + '++' + barcode_b_r3
									if combinedSgId not in pairCounts_barcodes:
										pairCounts_barcodes[combinedSgId] = 0
										
									pairCounts_barcodes[combinedSgId] += 1
									
							else:
								statsCounts['All sgRNAs and barcodes uniquely map'] += 1

								combinedSgId = protospacer_a_r1 + '++' + protospacer_b_r2
								if combinedSgId not in pairCounts_sgRNAs:
									pairCounts_sgRNAs[combinedSgId] = 0

								pairCounts_sgRNAs[combinedSgId] += 1

								
								combinedSgId = barcode_a_r3 + '++' + barcode_b_r3
								if combinedSgId not in pairCounts_barcodes:
									pairCounts_barcodes[combinedSgId] = 0
									
								pairCounts_barcodes[combinedSgId] += 1

								if protospacer_a_r1 != barcode_a_r3 or protospacer_b_r2 != barcode_b_r3:
									if protospacer_a_r1 != barcode_a_r3:
										statsCounts['A sgRNA and barcode do not match'] += 1
									if protospacer_b_r2 != barcode_b_r3:
										statsCounts['B sgRNA and barcode do not match'] += 1
							
								else:
									statsCounts['A and B sgRNAs and barcodes match'] += 1
									
									combinedSgId = protospacer_a_r1 + '++' + protospacer_b_r2
									
									if combinedSgId not in pairCounts_triple:
										pairCounts_triple[combinedSgId] = 0
										
									pairCounts_triple[combinedSgId] += 1
									
						if testRun and i == testLines:
							break

	with open(fileTup[3] + '.sgRNA.counts.txt', 'w') as outfile:
		for pair in sorted(pairCounts_sgRNAs.keys()):
			outfile.write(pair + '\t' + str(pairCounts_sgRNAs[pair]) + '\n')

	with open(fileTup[3] + '.barcode.counts.txt', 'w') as outfile:
		for pair in sorted(pairCounts_barcodes.keys()):
			outfile.write(pair + '\t' + str(pairCounts_barcodes[pair]) + '\n')

	with open(fileTup[3] + '.tripleseq.counts.txt', 'w') as outfile:
		for pair in sorted(pairCounts_triple.keys()):
			outfile.write(pair + '\t' + str(pairCounts_triple[pair]) + '\n')

	numReads = (i+1)/4
	print fileTup[3], numReads, 'reads'
	
	print 'Percent A sgRNAs mapping', 100.0 - (statsCounts['A sgRNA not mapped'] * 100.0 / numReads)
	print 'Percent B sgRNAs mapping', 100.0 - (statsCounts['B sgRNA not mapped'] * 100.0 / numReads)
	print 'Percent A barcodes mapping', 100.0 - (statsCounts['A barcode not mapped'] * 100.0 / numReads)
	print 'Percent B barcodes mapping', 100.0 - (statsCounts['B barcode not mapped'] * 100.0 / numReads)
# 	'A sgRNA multiple mappings': 0,
# 	'B sgRNA multiple mappings': 0, 
# 	'A barcode multiple mappings': 0,
# 	'B barcode multiple mappings': 0,
	print 'Percent all sgRNAs mapping', statsCounts['All sgRNAs and barcodes uniquely map'] * 100.0 / numReads
	print '--percent A sgRNA and barcode mismatch', statsCounts['A sgRNA and barcode do not match'] * 100.0 / statsCounts['All sgRNAs and barcodes uniquely map']
	print '--percent B sgRNA and barcode mismatch', statsCounts['B sgRNA and barcode do not match'] * 100.0 / statsCounts['All sgRNAs and barcodes uniquely map']
	print '--percent both A and B match', statsCounts['A and B sgRNAs and barcodes match'] * 100.0 / statsCounts['All sgRNAs and barcodes uniquely map']

	print 'Total percent matching and mapping reads', statsCounts['A and B sgRNAs and barcodes match'] * 100.0 / numReads
	
	sys.stdout.flush()


def getMismatchDict(table, column, trim_range = None, allowOneMismatch = True):
    mismatch_dict = dict()

    if trim_range:
        col = table[column].apply(lambda seq: seq[trim_range[0]:trim_range[1]])
    else:
        col = table[column]
        
    for sgRNA, seq in col.iteritems():
#         print seq
        if seq in mismatch_dict:
            print 'clash with 0 mismatches', sgRNA, seq
            mismatch_dict[seq] = 'multiple'
        else:
            mismatch_dict[seq] = sgRNA

        if allowOneMismatch:
            for position in range(len(seq)):
                mismatchSeq = seq[:position] + 'N' + seq[position + 1:]

                if mismatchSeq in mismatch_dict:
                    print 'clash with 1 mismatch', sgRNA, mismatchSeq
                    mismatch_dict[seq] = 'multiple'
                else:
                    mismatch_dict[mismatchSeq] = sgRNA

    return mismatch_dict	
    
def matchBarcode(mismatch_dict, barcode, allowOneMismatch = True):
    if barcode in mismatch_dict:
        match = mismatch_dict[barcode]
        
    elif allowOneMismatch:
        match = 'none'
        for position in range(len(barcode)):
            mismatchSeq = barcode[:position] + 'N' + barcode[position + 1:]

            if mismatchSeq in mismatch_dict:
                if match == 'none':
                    match = mismatch_dict[mismatchSeq]
                else:
                    match = 'multiple'
                   
    else:
        match = 'none'
        
    return match					
	
testLines = 100000

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Process raw sequencing data from screens to counts files in parallel')
	parser.add_argument('Singles_Table', help='Table of individual sgRNAs in the library.')
	parser.add_argument('Out_File_Path', help='Directory where output files should be written.')
	parser.add_argument('Seq_File_Names', nargs='+', help='Name(s) of sequencing file(s). Unix wildcards can be used to select multiple files at once. The script will search for all *.fastq.gz, *.fastq, and *.fa(/fasta/fna) files with the given wildcard name.')
			
# 	parser.add_argument('-p','--processors', type=int, default = 1)
# 	parser.add_argument('--trim_start', type=int)
# 	parser.add_argument('--trim_end', type=int)
	parser.add_argument('--test', action='store_true', default=False, help='Run the entire script on only the first %d reads of each file. Be sure to delete or move all test files before re-running script as they will not be overwritten.' % testLines)

	args = parser.parse_args()

	outputDirectory = args.Out_File_Path

	inputFileList = sorted(args.Seq_File_Names) #assuming sorting will yield groups of R1-3 files
	#print inputFileList

	singlesTable = pd.read_csv(args.Singles_Table, sep='\t', header=None).iloc[:, [0,1,5,6,7]]
	singlesTable.columns = [
		'id',
		'sgId',
		'sgRNA sequence',
		'upstream barcode',
		'downstream barcode',
	]
	singlesTable = singlesTable.set_index('sgId').dropna()
	#print singlesTable.head()
	print 'sgRNAs in library', len(singlesTable)

	combinedMismatchDicts = {'protospacer_a_r1': getMismatchDict(singlesTable, 'sgRNA sequence', (1,20), allowOneMismatch=True),
	'protospacer_b_r2': getMismatchDict(singlesTable, 'sgRNA sequence', (1,20), allowOneMismatch=True),
	'barcode_a_r3': getMismatchDict(singlesTable, 'downstream barcode', allowOneMismatch=True),
	'barcode_b_r3': getMismatchDict(singlesTable, 'upstream barcode', allowOneMismatch=True)}

	fileTups = []
	for i, fastqfile in enumerate(inputFileList):
		if i%3 == 0:
			r1file = fastqfile
		elif i%3 == 1:
			r2file = fastqfile
		else:
			r3file = fastqfile
		
			outputfile = os.path.join(outputDirectory, os.path.split(fastqfile)[-1].split('_R')[0])
		
			fileTups.append((r1file, r2file, r3file, outputfile, combinedMismatchDicts, args.test))
		
	pool = multiprocessing.Pool(len(fileTups))

	# print fileTups
	pool.map(writeToCounts, fileTups)

	pool.close()
	pool.join()
