import os
import time
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import scrublet as scr
import numpy as np

def Running_Scrublet( matrixArg, labelArg, outputArg, seedArg, threArg=0.25, estRate=None ):

	# read the matrix
	f = open(matrixArg, 'r')
	fileLines = f.readlines()
	f.close()

      	# Transpose the matrix into format cells X genes
	genes = [x.split('\t')[0] for x in fileLines[1:]]
	barcodes = fileLines[0].rstrip().split('\t')
	counts_matrix = [x.rstrip().split('\t')[1:] for x in fileLines[1:]]
	counts_matrix = np.array(counts_matrix, dtype = 'float64')
	counts_matrix = counts_matrix.transpose()

      	# Initialize Scrublet object
	if estRate is None:
		scrub = scr.Scrublet(counts_matrix, random_state=seedArg)
	else:
		scrub = scr.Scrublet(counts_matrix, random_state=seedArg, expected_doublet_rate=estRate)
	# Run the default pipeline
	doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=3, min_cells=3, 
						min_gene_variability_pctl=85, get_doublet_neighbor_parents=True)

        # if manual threshold is given
	if threArg:
		print ('\nUsing manual threshold for doublet score\n')
		predicted_doublets = scrub.call_doublets(threshold=threArg)
	else:   pass
	print ("Doublet scores are : ", doublet_scores, '\n')
	print ("Predicted doubles are : ", predicted_doublets, "\n")

	# Plot histograms of observed transcriptomes and simulated dbls
	scrub.plot_histogram()
	plt.savefig(outputArg + 'plot_histogram.png')

	# Check Labels of data
	l = open(labelArg, 'r')
	labelLines = l.readlines()
	l.close()

	if len(labelLines)-1 != len(doublet_scores):
		print ('\n#######################\nlabel length is different\n#######################\n')
		sys.exit

      	# Make result table
	writeLines = [labelLines[0].rstrip() + '\tScrublet_score\tScrublet_label\n']
	for i in range(len(predicted_doublets)):
		writeLine = labelLines[i+1].rstrip() + '\t' + str(doublet_scores[i]) + '\t' + str(predicted_doublets[i]) + '\n'
		writeLines.append(writeLine)
	outwrite = open( outputArg + 'scrublet_result.txt', 'w')
	outwrite.write(''.join(writeLines))
	outwrite.close()



def Multiple_Running( expMtx, labelT, OutDir, nRun, estRate=None ):

	for i in range(int(nRun)): 
		seed = (i+1)*2022
		outdir = OutDir + "scr_" + str(i+1) + "."
		
		print("\n\n=====================> Iteration ", i+1," <=====================")
		print(" Seed : ", seed)
		print(" Output path : ", outdir, '\n')

		Running_Scrublet( expMtx, labelT, outdir, seedArg=seed, estRate=estRate )			


def Assembling_Data( setnumber, InputDir, OutputDir, Datatype, nRun, estRate=None ):
		
	if InputDir.endswith("/")==False: 
		InputDir = InputDir + '/'
	if OutputDir.endswith("/")==False:
		OutputDir = OutputDir + '/'

	# Algorithm to use
	print('\n===========================================================================================')
	print("                                Dataset number ", str(setnumber+1))
	print('   Detection algorithm : Scrublet')		

	# Data type
	if Datatype=='raw': 
		print('   Data type : Raw count matrix ')
		Matfilelist = [os.path.join(InputDir,file) for file in os.listdir(InputDir) if file.endswith('ExpMtx.nonZero.txt')]
		Labfilelist = [os.path.join(InputDir,file) for file in os.listdir(InputDir) if file.endswith('LabelTable.txt')]
		expMtx=Matfilelist[0]
		labelT=Labfilelist[0]
		
		OutDir = OutputDir + '1.Raw/'
		if not os.path.exists(OutDir):
			 os.mkdir(OutDir)
		
	else:
		print('   Data type : Simulated count matrix ')
		expMtx = InputDir + 'simul_' + str(setnumber+1) + '/simulDbls.ExpMtx.txt'
		labelT = InputDir + 'simul_' + str(setnumber+1) + '/simulDbls.LabelTable.txt' 	
	
		OutDir = OutputDir + 'simul_' + str(setnumber+1) + '/'
		if not os.path.exists(OutDir):
			 os.mkdir(OutDir)		
	
	print('   - Input file : ', expMtx)
	print('   - Output directory : ', OutDir)
	print('===========================================================================================')
	
	Multiple_Running( expMtx, labelT, OutDir, nRun, estRate) 


def timecheck(start, end):
	h, r = divmod(end-start, 3600)
	m, s = divmod(r, 60)
	print("The total execution time : {:0>2}:{:0>2}:{:05.2f}".format( int(h), int(m), int(s)))


def RunningScrublet(InputDir, OutputDir, Datatype, nSimul=100, nRun=10, estRate=None):
	if Datatype=='raw':
		start_time = time.time()
		if estRate is None:
			Assembling_Data( 0, InputDir, OutputDir, Datatype, nRun)
		else:
			Assembling_Data( 0, InputDir, OutputDir, Datatype, nRun, estRate)
		print("\n-------------------->  Finished <--------------------\n")
		end_time = time.time()
		timecheck(start_time, end_time)
	
	else :  #simulation data
		start_time = time.time()
		for i in range(int(nSimul)):
			if estRate is None:
				Assembling_Data( i, InputDir, OutputDir, Datatype, nRun)
			else:
				Assembling_Data( i, InputDir, OutputDir, Datatype, nRun, estRate=float(estRate))
			print("--------------------> Dataset ", str(i+1), " Finished ! <--------------------\n\n")
		end_time = time.time()
		timecheck(start_time, end_time)


