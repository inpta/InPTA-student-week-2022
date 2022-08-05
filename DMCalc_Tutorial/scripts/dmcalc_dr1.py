#!$HOME/anaconda3/envs/pulsar/bin/

################### WIDEBAND DM ESTIMATION SCRIPT ###################
# Script for estimating the DM of an observation using PSRCHIVE and 
# TEMPO2. This code is extracted from another one written for LOFAR 
# data analysis by Caterina (trimtim.py). 
#
# Dependencies: 
# PSRCHIVE python interface: 
# http://psrchive.sourceforge.net/manuals/python/
# SKLEARN: https://scikit-learn.org/stable/install.html
# TEMPO2: https://bitbucket.org/psrsoft/tempo2
# SCIPY
#
# Usage: 
# ./dmcalc.py -E test.par -M test.sm test.fits
#
# For more options and information, please check help section.
#
# If you have a directory with all the model files in one directory 
# (PWD/templates/) with '.sm' extension and all the parameter files 
# in 'ephemerides' directory with the name 'JXXXX-YYYY.par', you do 
# not need to give the -E and -M options and simply do
#
# ./dmcalc.py test.fits
#
#####################################################################

# import modules...
import os
import sys
import numpy as np
import psrchive
import argparse
import time
import warnings
warnings.filterwarnings("ignore")
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle


start = time.time()

def existing_file(filename):
	if not(os.path.exists(filename) and os.path.isfile(filename)):
		raise argparse.ArgumentTypeError("**"+filename+" does not exist**")
	else:
		return filename

parser = argparse.ArgumentParser(description='Code for measuring in-band '+ 
                                 'DM for pulsar data in psrfits format.')
parser.add_argument('files', nargs='+', type=existing_file, 
					help='The list of fits file(s) for processing')
parser.add_argument('-E', '--ephem', type=existing_file, 
					help='Ephemeris file to update the model. Exits if not '+
					      'given or is not available in "PWD/ephemerides" '+
					      'directory')
parser.add_argument('-M', '--model', nargs='+', type=existing_file,
					help='Model template for ToA generation. Exits if not '+ 
					     'given or is not available in "PWD/templates" '+
					     'directory')
parser.add_argument('-nch','--nchan', type=int, default=8,
					help='Number of frequency channels to use while '+
						 'estimating DM (Def: 8)')
parser.add_argument('-b3n','--b3nchan', type=int, default=16,
					help='Number of frequency channels to use in '+ 
					     'band3 uGMRT data (Def: 16)')
parser.add_argument('-b5n','--b5nchan', type=int, default=8,
					help='Number of frequency channels to use in '+ 
					     'band5 uGMRT data (Def: 8)')
parser.add_argument('-q', '--quiet', action='store_true', 
							help='Only print warnings')
parser.add_argument('-v5','--v5orold', action='store_true',
					help='Employs backend delay correction for old uGMRT processed '+
						 'data using the measured delays. Def: False')

"""
Estimates the Dispersion Measure (DM) from the data in psrfits file format.

Returns the value of DM with its uncertainty and reduced chi-square from
either tempo2 or an MCMC method based on fitting a second order polynomial
using scipy.optimize.

Parameters
----------
file(s) :  Input file(s) in psrfits format

tempo2  :  bool, optional. default: False. If True, performs DM calculation 
           using  tempo2. Otherwise, uses  the MCMC  method to  estimate DM.
          
ephem   :  Ephemeris (or parameter) file  of the  pulsar. This is  required 
           to update the model. It can be  given as a command line argument. 
           If it is available in "PWD/ephemerides" folder, one can use that.
           Giving the file with this option overrides the default one.

model   :  Template profile for cross-correlating  with the observation  to
           obtain DM. It can be given as a command line argument, otherwise
           it will look  for a matching one in  "PWD/ephemerides" directory
           and if found, will use that instead. One can use this  option to
           override the default selection.
           
nchan   : int, optional, default: 8. Number of frequency channels to use in
          the estimation of DM.

b3nchan : int, optional, default: 16. Number of frequency channels in band3
          of uGMRT data.
          
b5nchan : int, optional, default: 8. Number of frequency channels in band5
          of uGMRT data.

v5       : This is used for processing the uGMRT data folded with pinta v5 
           or earlier versions. The backend delays are measured and are kept
           based on the NCRA internal report.

quiet    : bool, optional,  default: False. Supresses all print  statements
           except warnings and errors.

Returns
-------
Dispersion Measure with uncertainty.

Notes
-----

Examples
--------
# (a) for a simple DM estimation with built-in fitting function and default
# directories:
#
dmcalc.py inputfile.fits
#
# (b) for using tempo2 to measure DM:
#
dmcalc.py inputfile.fits
#
# (c) to use different ephemeris and template files:
#
dmcalc.py -E aaaa.eph -M bbbb.fits inputfile.fits
#

"""
# Module that connects all the others and does the job.
def main():
	
	# parses the input arguments
	args = parser.parse_args()

	# checks status of quiet and tempo2
	quiet=False
	if args.quiet:
		quiet=True
	if args.ephem != None:
		ephemeris = args.ephem
	else:
		ephemeris = "ephemerides/"+ar_psr+".par"
		if not (os.path.exists(ephemeris)):
			print("\nError! No parameter file is given. Exiting!\n")
			sys.exit(1)
	if not quiet:
		print ("Pulsar Parameter file is:"+ephemeris+'\n')
		
	if not quiet:
		print("Loading the archive files for DM estimation... "),

	# loads the data and template file(s)
	pwd=os.getcwd()
	archives = []
	finalarchives = []
	model = []
	for filename in args.files:
		archives.append(psrchive.Archive_load(filename))
		finalarchives.append(psrchive.Archive_load(filename))
	narch = len(archives)
	for filename in args.model:
		model.append(psrchive.Archive_load(filename))
	nmod = len(model)	
	arfrq = np.zeros(narch)
	modfrq = np.zeros(nmod)
	# Reads the DM from ephemeris and applies it to both data and template files.
	# This step is no longer needed, but is kept as is for being safe than sorry.
	templ_dm = 0.0
	with open (args.ephem, 'r') as read_eph:
		for line in read_eph:
			if line.startswith('DM\t') or line.startswith('DM '):
				templ_dm = float (line.split()[1])
	for i in range (narch):
		arfrq[i] = archives[i].get_centre_frequency()
		archives[i].set_dispersion_measure(templ_dm)
		archives[i].set_ephemeris(ephemeris)
		archives[i].update_model()
		finalarchives[i].set_ephemeris(ephemeris)
		finalarchives[i].update_model()
	for j in range(nmod):
		modfrq[j] = model[j].get_centre_frequency()
		model[i].set_dispersion_measure(templ_dm)
		model[i].set_ephemeris(ephemeris)
		model[i].update_model()
	# If there are more than one input and template file, then the following part
	# makes the ordering of those files correct, if it was in a different order.
	# Additionally, for uGMRT data prior to pinta v6.2, it does the proper time 
	# corrections. It will also make the nchan of the data and template as per the
	# input. 
	if (narch == nmod and narch > 1):
		modpos = np.zeros(nmod)
		for i in range(narch):
			if (300. < archives[i].get_centre_frequency() < 500.):
				if args.v5orold:
					archives[i] = Correct_delay(archives[i])
				archives[i].fscrunch_to_nchan(args.b3nchan)
			if (1100. < archives[i].get_centre_frequency() < 1500.):				
				if args.v5orold:
					archives[i] = Correct_delay(archives[i])
				archives[i].fscrunch_to_nchan(args.b5nchan)
		for i in range(nmod):
			if (300. < model[i].get_centre_frequency() < 500.):
				model[i].fscrunch_to_nchan(args.b3nchan)
			if (1100. < model[i].get_centre_frequency() < 1500.):
				model[i].fscrunch_to_nchan(args.b5nchan)
		for i in range(narch):
			for j in range(nmod):
				if (np.around(arfrq[i]) == np.around(modfrq[j])):
					modpos[i] = j

	if not quiet:
		print(" done!")
	if (nmod == 2):
		model[0], model[1] = model[int(modpos[0])], model[int(modpos[1])]
	nchan = args.nchan
	if (narch == 1):
		if (args.b3nchan != 16):
			nchan = args.b3nchan
		if (args.b5nchan != 8):
			nchan = args.b5nchan
		archives[0].fscrunch_to_nchan(nchan)
		model[0].fscrunch_to_nchan(nchan)


	# Obtains the ToAs of each fits files and combines them. Up to 512 ToAs are 
	# possible per file. If required, this limit can be changed. 
	fr_o = np.zeros(shape=(narch,512))
	re_o = np.zeros(shape=(narch,512))
	tE_o = np.zeros(shape=(narch,512))
	freq = np.zeros(shape=(narch,512))
	toa  = np.zeros(shape=(narch,512))
	toaE = np.zeros(shape=(narch,512))
	resd = np.zeros(shape=(narch,512))
	for i in range(narch):
		#print (archives[i], model[i])
		tx1, tx2, tx3, td1, td2, td3, td4 = (get_TOAs(archives[i], model[i], ephemeris, quiet))
		for j in range(np.size(td1)):
			freq[i][j] = td1[j]
			toa[i][j]  = td2[j]
			resd[i][j] = td3[j]
			toaE[i][j] = td4[j]
		for k in range(np.size(tx1)):
			fr_o[i][k] = tx1[k]
			re_o[i][k] = tx2[k]
			tE_o[i][k] = tx3[k]
	tt1 = np.reshape(freq,np.size(freq))
	tt2 = toa.flatten()
	tt3 = resd.flatten()
	tt4 = toaE.flatten()
	condition = tt1 != 0
	init_freqs = np.extract(condition,tt1)
	init_toas  = np.extract(condition,tt2)
	init_resid = np.extract(condition,tt3)
	init_toasE = np.extract(condition,tt4)
	tt1 = np.reshape(fr_o,np.size(fr_o))
	tt2 = re_o.flatten()
	tt3 = tE_o.flatten() * 1e+6
	condition = tt1 != 0
	orig_freqs = np.extract(condition,tt1)
	orig_resid = np.extract(condition,tt2)
	orig_toasE = np.extract(condition,tt3)
	#print (orig_resid, orig_toasE)
	# Obtains the DM of all the files (together and separate) using tempo2 DM fit.
	if not quiet:
		print ("\n\nNow estimating DM...\n")
	dm, dmerror, chisq = get_DM(init_freqs, init_toas, init_toasE, archives, ephemeris)
	if not quiet:
		print ("\n\nUpdating the archives with new DM... "),

	mjd_start=archives[0].get_Integration(0).get_start_time().in_days()
	mjd_end=archives[0].get_Integration(0).get_end_time().in_days()
	ar_mjd = mjd_start + (mjd_end-mjd_start)/2.
	ar_psr = archives[0].get_source()
	ar_tel = archives[0].get_telescope()
	# Removing the DM and DMEPOCH from a copy of the ephemeris file given.
	oldpar = open(ephemeris,"r")
	partmp = ar_psr+'_tmp.par'
	newpar = open(partmp,"w+")
	for i, line in enumerate(oldpar):
		if not line.lstrip().startswith('DM'):
				if not line.lstrip().startswith('DMEPOCH'):
					newpar.write(line)
	oldpar.close()
	newpar.close()

	# updating the ephemeris file with measured DM
	dmline = "DM             "+str(dm[0])+"\t1\t"+str(dmerror[0])
	dmepochline  = "DMEPOCH	       "+str(round(ar_mjd,2))
	f = open(partmp,'a')
	f.write('%s\n%s\n' % (dmline, dmepochline))
	f.close()
	# Correcting the observed files with the obtained DM and getting the 
	# DM corrected ToAs for plotting and statistics.
	for i in range(narch):
		archives[i].set_ephemeris(partmp)
		archives[i].set_dispersion_measure(dm[0])
		archives[i].update_model()
		finalarchives[i].set_ephemeris(partmp)
		finalarchives[i].set_dispersion_measure(dm[0])
		finalarchives[i].update_model()
		finalarchives[i].tscrunch()
		finalarchives[i].dedisperse()
		finalarchives[i].remove_baseline()
		
	if not quiet:
		print(" done!\n")
		print("\n Getting the DM corrected ToAs... ")		
	freqf1 = np.zeros(shape=(narch,512))
	toaf1  = np.zeros(shape=(narch,512))
	toaEf1 = np.zeros(shape=(narch,512))
	resdf1 = np.zeros(shape=(narch,512))
	for i in range(narch):
		_, _, _, ty1, ty2, ty3, ty4 = (get_TOAs(archives[i], model[i], partmp, quiet))
		for j in range(np.size(ty1)):
			freqf1[i][j] = ty1[j]
			toaf1[i][j]  = ty2[j]
			resdf1[i][j] = ty3[j]
			toaEf1[i][j] = ty4[j]
	os.remove(partmp)
	ty1 = np.reshape(freqf1,np.size(freqf1))
	ty2 = toaf1.flatten()
	ty3 = resdf1.flatten()
	ty4 = toaEf1.flatten()
	condition = ty1 != 0
	final_freqs = np.extract(condition,ty1)
	final_toas  = np.extract(condition,ty2)
	final_resid = np.extract(condition,ty3)
	final_toasE = np.extract(condition,ty4)
	init_resid -= np.median(init_resid)
	final_resid -= np.median(final_resid)
	prefit_rms = np.zeros(np.size(dm)); postfit_rms = np.zeros(np.size(dm))
	med_toaE = np.zeros(np.size(dm)); centre_freqs = np.zeros(np.size(dm))
	bw = np.zeros(np.size(dm))
	prefit_rms[0] = np.sqrt(np.cov(init_resid, aweights=init_toasE))
	postfit_rms[0] = np.sqrt(np.cov(final_resid, aweights=final_toasE))
	med_toaE[0] = np.median(final_toasE)
	# Setting the correct band flags and obtaining the BW, centre freq of all the 
	# files for writing them out to a file.
	bandflag = [None] * 3
	if (len(archives) == 1):
		centre_freqs[0] = archives[0].get_centre_frequency()
		bw[0] = archives[0].get_bandwidth()
		if (300. < centre_freqs[0] < 500.0):
			bandflag[0] = 'band3'
		if (1160. < centre_freqs[0] < 1460.0):
			bandflag[0] = 'band5'
	if (len(archives) > 1):
		bandflag[0] = "band3+5"
		for i in range(len(archives)):
			cfrq = archives[i].get_centre_frequency()
			if (300. < cfrq < 500.):
				bw[1] = archives[i].get_bandwidth()
				centre_freqs[1] = archives[i].get_centre_frequency()
				bandflag[1] = "band3"
				condition = (init_freqs < 500.) & (init_freqs > 300.)
				tx = np.extract(condition,init_resid)
				tz = np.extract(condition,init_toasE)
				prefit_rms[1] = np.sqrt(np.cov(tx, aweights=tz))
				condition = (final_freqs < 500.) & (final_freqs > 300.)
				ty = np.extract(condition,final_resid)
				tz = np.extract(condition,final_toasE)
				postfit_rms[1] = np.sqrt(np.cov(ty, aweights=tz))
				med_toaE[1] = np.median(tz)
			if (1160. < cfrq < 1460.):
				bw[2] = archives[i].get_bandwidth()
				centre_freqs[2] = archives[i].get_centre_frequency()
				bandflag[2] = "band5"
				condition = (init_freqs < 1460.) & (init_freqs > 1160.)
				tx = np.extract(condition,init_resid)
				tz = np.extract(condition,init_toasE)
				prefit_rms[2] = np.sqrt(np.cov(tx, aweights=tz))
				condition = (final_freqs < 1460.) & (final_freqs > 1160.)
				ty2 = np.extract(condition,final_resid)
				tz2 = np.extract(condition,final_toasE)
				postfit_rms[2] = np.sqrt(np.cov(ty2, aweights=tz2))
				med_toaE[2] = np.median(tz2)
			bw[0] = bw[1] + bw[2]		
			centre_freqs[0] = (centre_freqs[1] + centre_freqs[2])/2.


	print (dmerror)
	if (str(dmerror[0]) == 'inf'):
		print ("\nERROR: The DM uncertainty is 'inf'! Please check the output plot and make sure data is correct. Not writing the DMs and ToAs out.")	
	if not (str(dmerror[0]) == 'inf'):

		# Getting the frequency resolved ToAs
		for i in range(narch):
			get_finalTOA(archives[i], model[i], ephemeris, templ_dm, quiet)
	
		# Creating a par file with DMMODEL parameters

		#dmmodelpar = ar_psr+"_"+str("%.f" % bw[0])+"MHz.DMMODEL.par"
		dmmodelpar = ar_psr+".DMMODEL.par"
		if (os.path.isfile(dmmodelpar)):
			oldpar = open(dmmodelpar,"r")
		if not (os.path.isfile(dmmodelpar)):
			f1 = open(ephemeris,"r")
			partmp = dmmodelpar
			dmmodelpar = open(partmp,"w+")
			for line in f1:
				dmmodelpar.write(line)
			dmmodelpar.write("DMMODEL DM 0\n")
			dmmodelpar.close()
			del dmmodelpar
			f1.close()
			#dmmodelpar = ar_psr+"_"+str("%.f" % bw[0])+"MHz.DMMODEL.par"
			dmmodelpar = ar_psr+".DMMODEL.par"
			oldpar = open(dmmodelpar,"r")
		#oldpar = dmmodelpar
		partmp = ar_psr+'_tmp.par'
		newpar = open(partmp,"w+")
		for i, line in enumerate(oldpar):
			dmomjdstr="DMOFF "+str("%.6f" % ar_mjd)
			#print(dmomjdstr)
			if line.lstrip().startswith(dmomjdstr):
				print("\nERROR: DM for MJD %.6f already exists. Maybe try moving the DMMODEL and DMX files and retry." % ar_mjd)
				sys.exit(1)
			if not line.lstrip().startswith('CONSTRAIN'):
						newpar.write(line)
		oldpar.close()
		newpar.close()
	
		ttmp = os.popen("mv %s %s" % (partmp,dmmodelpar)).read()
		f1 = open(dmmodelpar,"a")
		f1.write("DMOFF %.6f %.8f %.8f\n" %(ar_mjd,dm[0]-templ_dm,dmerror[0]))
		f1.write("CONSTRAIN DMMODEL")
		f1.close()
		# Creating a par file with DMX parameters
		#dmxpar = ar_psr+"_"+str("%.f" % bw[0])+".DMX.par"
		dmxpar = ar_psr+".DMX.par"
		if (os.path.isfile(dmxpar)):
			with open(dmxpar,'r') as f:
				last_line = f.readlines()[-1]
				last_dmx =int(last_line.strip().split()[0].split('_')[1])
				this_dmx = last_dmx+1
				dmx1 = f"DMX_{this_dmx:04d}"; dmx2 = f"DMXEP_{this_dmx:04d}"
				dmx3 = f"DMXR1_{this_dmx:04d}"; dmx4 = f"DMXR2_{this_dmx:04d}"
				dmx5 = f"DMXF1_{this_dmx:04d}"; dmx6 = f"DMXF2_{this_dmx:04d}"
			if (len(dm)>1):
				f1 = open(dmxpar,"a")
				f1.write("%s\t%.8e\t0\t%.8e\n" % (dmx1, dm[0]-templ_dm, dmerror[0]))
				f1.write("%s\t%.5f\n" % (dmx2, ar_mjd))
				f1.write("%s\t%.5f\n%s\t%.5f\n" % (dmx3, ar_mjd-0.5, dmx4, ar_mjd+0.5))
				f1.write("%s\t%.5f\n%s\t%.5f\n" % (dmx5, centre_freqs[1], dmx6, centre_freqs[2]))
				f1.close()
			if (len(dm)==1):
				f1 = open(dmxpar,"a")
				f1.write("%s\t%.8e\t0\t%.8e\n" % (dmx1, dm[0]-templ_dm, dmerror[0]))
				f1.write("%s\t%.5f\n" % (dmx2, ar_mjd))
				f1.write("%s\t%.5f\n%s\t%.5f\n" % (dmx3, ar_mjd-0.5, dmx4, ar_mjd+0.5))
				f1.write("%s\t%.5f\n%s\t%.5f\n" % (dmx5, centre_freqs[0]-bw[0]/2., dmx6, centre_freqs[0]+bw[0]/2))
				f1.close()
	
		#dmxpar = ar_psr+"_"+str("%.f" % bw[0])+".DMX.par"
		dmxpar = ar_psr+".DMX.par"
		if not (os.path.isfile(dmxpar)):
			if (len(dm) >1):
				f1 = open(ephemeris,"r")
				partmp = dmxpar
				dmxpar = open(partmp,"w+")
				for line in f1:
					dmxpar.write(line)
				#dmxpar.write("DMX\t\t6.500000\n")
				dmxpar.write("DMX_0001\t%.8e\t0\t%.8e\n" % (dm[0]-templ_dm,dmerror[0]))
				dmxpar.write("DMXEP_0001\t%.5f\n" % (ar_mjd))
				dmxpar.write("DMXR1_0001\t%.5f\n" % (ar_mjd-0.5))
				dmxpar.write("DMXR2_0001\t%.5f\n" % (ar_mjd+0.5))
				dmxpar.write("DMXF1_0001\t%.5f\n" % (centre_freqs[1]))
				dmxpar.write("DMXF2_0001\t%.5f\n" % (centre_freqs[2]))
				dmxpar.close()
				f1.close()
			if (len(dm) ==1):
				f1 = open(ephemeris,"r")
				partmp = dmxpar
				dmxpar = open(partmp,"w+")
				for line in f1:
					dmxpar.write(line)
				#dmxpar.write("DMX\t\t6.500000\n")
				dmxpar.write("DMX_0001\t%.8e\t0\t%.8e\n" % (dm[0]-templ_dm,dmerror[0]))
				dmxpar.write("DMXEP_0001\t%.5f\n" % (ar_mjd))
				dmxpar.write("DMXR1_0001\t%.5f\n" % (ar_mjd-0.5))
				dmxpar.write("DMXR2_0001\t%.5f\n" % (ar_mjd+0.5))
				dmxpar.write("DMXF1_0001\t%.5f\n" % (centre_freqs[0]-bw[0]/2))
				dmxpar.write("DMXF2_0001\t%.5f\n" % (centre_freqs[0]+bw[0]/2))
				dmxpar.close()
				f1.close()
	
		# Printing the results to the file and also in the terminal
		f= open(ar_psr+"_DM_timeseries.txt",'a')
		for i in range(np.size(dm)):
			f.write('%.6f %.6f %.6f %.2f %.4f %.4f %.4f %.2f %.2f %s %s\n' %(
				ar_mjd, dm[i], dmerror[i], chisq[i], prefit_rms[i], postfit_rms[i], med_toaE[i], centre_freqs[i], bw[i], ar_tel, bandflag[i]))
		f.close()
	
	
	# Plotting the pre- and post fit ToAs with all the details
	if (len(finalarchives) > 1):
		prof2Db3 = []
		profb3 = []
		prof2Db5 = []
		profb5 = []
		b3_bw = []; b5_bw = []; b3_freq = []; b5_freq = []; b3_nbin = []; b5_nbin = []
		for i in range(len(finalarchives)):
			if (finalarchives[i].get_centre_frequency() < 500.):
				ar_nchan = finalarchives[i].get_nchan()
				b3_nbin  = finalarchives[i].get_nbin()
				b3_bw = finalarchives[i].get_bandwidth()
				b3_freq = finalarchives[i].get_centre_frequency()
				prof2Db3 = finalarchives[i].get_data()[:,0,:,:].flatten().reshape(ar_nchan,b3_nbin)
				prof = finalarchives[i].clone()
				prof.fscrunch()
				profb3 = prof.get_data().flatten()
				profb3 /= np.max(profb3)
			if (finalarchives[i].get_centre_frequency() > 1000.):
				ar_nchan = finalarchives[i].get_nchan()
				b5_nbin  = finalarchives[i].get_nbin()
				b5_bw = finalarchives[i].get_bandwidth()
				b5_freq = finalarchives[i].get_centre_frequency()
				prof2Db5 = finalarchives[i].get_data()[:,0,:,:].flatten().reshape(ar_nchan,b5_nbin)
				prof = finalarchives[i].clone()
				prof.fscrunch()
				profb5 = prof.get_data().flatten()
				profb5 /= np.max(profb5)
		
		condition = orig_freqs < 500.
		orig_b3fr = np.extract(condition,orig_freqs)
		orig_b3re = np.extract(condition,orig_resid)
		orig_b3Er = np.extract(condition,orig_toasE)

		condition = orig_freqs > 1000.
		orig_b5fr = np.extract(condition,orig_freqs)
		orig_b5re = np.extract(condition,orig_resid)
		orig_b5Er = np.extract(condition,orig_toasE)

		condition = init_freqs < 500.
		init_b3fr = np.extract(condition,init_freqs)
		init_b3re = np.extract(condition,init_resid)
		init_b3Er = np.extract(condition,init_toasE)

		condition = init_freqs > 1000.
		init_b5fr = np.extract(condition,init_freqs)
		init_b5re = np.extract(condition,init_resid)
		init_b5Er = np.extract(condition,init_toasE)

		condition = final_freqs < 500.
		final_b3fr = np.extract(condition,final_freqs)
		final_b3re = np.extract(condition,final_resid)
		final_b3Er = np.extract(condition,final_toasE)

		condition = final_freqs > 1000.
		final_b5fr = np.extract(condition,final_freqs)
		final_b5re = np.extract(condition,final_resid)
		final_b5Er = np.extract(condition,final_toasE)

		fig = plt.figure(3, figsize=(8, 6))
		fig.subplots_adjust(hspace=0.05)
		ax0 = plt.subplot2grid((9, 8), (0,0), rowspan=3, colspan=3)
		ax1 = plt.subplot2grid((9, 8), (3,0), rowspan=1, colspan=3)
		ax2 = plt.subplot2grid((9, 8), (5,0), rowspan=3, colspan=3)
		ax3 = plt.subplot2grid((9, 8), (8,0), rowspan=1, colspan=3)
		
		ax4 = plt.subplot2grid((9, 8), (0,4), colspan=4, rowspan=3)
		ax5 = plt.subplot2grid((9, 8), (3,4), colspan=4, rowspan=3)
		ax6 = plt.subplot2grid((9, 8), (6,4), colspan=4, rowspan=3)
		leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)
		
		ax0.imshow((np.sqrt(prof2Db5**2))**0.5, origin='lower', extent=(0,b5_nbin-1,(np.around(b5_freq)-b5_bw/2),(np.around(b5_freq)+b5_bw/2)), aspect='auto', cmap='hot')
		ax0.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=8)
		ax0.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
		ax1.plot(np.arange(b5_nbin, dtype=float),profb5, color='black', linewidth=0.5)
		ax1.set_xlim(0,b5_nbin-1)
		ax1.set_ylabel('Intensity', fontweight='bold', fontsize=8)

		ax2.imshow((np.sqrt(prof2Db3**2))**0.5, origin='lower', extent=(0,b3_nbin-1,(np.around(b3_freq)-b3_bw/2),(np.around(b3_freq)+b3_bw/2)), aspect='auto', cmap='hot')
		ax2.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=8)
		ax2.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
		ax3.plot(np.arange(b3_nbin, dtype=float),profb3, color='black', linewidth=0.5)
		ax3.set_xlim(0,b3_nbin-1)
		ax3.set_xlabel('Pulse Phase (bins)', fontweight='bold', fontsize=8)
		ax3.set_ylabel('Intensity', fontweight='bold', fontsize=8)

		ax4.errorbar(orig_b3fr, orig_b3re, yerr=orig_b3Er, fmt='.', color='#D81B60', capsize=2)
		ax4.errorbar(orig_b5fr, orig_b5re, yerr=orig_b5Er, fmt='.', color='#1E88E5', capsize=2)
		ax4.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b5_freq)+b3_bw/0.8))
		ax4.grid()
		ax4.legend([leg1], ['Prefit: Unfiltered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax4.axes.xaxis.set_ticklabels([])

		ax5.errorbar(init_b3fr, init_b3re, yerr=init_b3Er, fmt='.', color='#D81B60', capsize=2)
		ax5.errorbar(init_b5fr, init_b5re, yerr=init_b5Er, fmt='.', color='#1E88E5', capsize=2)
		ax5.grid()
		ax5.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b5_freq)+b3_bw/0.8))
		ax5.legend([leg1], ['Prefit: Filtered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax5.axes.xaxis.set_ticklabels([])
		ax5.set_ylabel(r'ToA Residuals ($\mu$s)', fontweight='bold', fontsize=8)
		
		ax6.errorbar(final_b3fr, final_b3re, yerr=final_b3Er, fmt='.', color='#D81B60', capsize=2)
		ax6.errorbar(final_b5fr, final_b5re, yerr=final_b5Er, fmt='.', color='#1E88E5', capsize=2)
		ax6.grid()
		ax6.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b5_freq)+b3_bw/0.8))
		ax6.legend([leg1], ['Postfit'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax6.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=9)
		fig.suptitle('Source: PSR %s; MJD: %.4f; Prefit Wrms: %.2f $\mu$s; Postfit Wrms: %.2f $\mu$s\nMedian ToA Err: %.2f $\mu$s; DM: %.6f $\pm$ %.6f pc cm$^{-3}$;  Reduced $\chi^2$: %.2f'%(ar_psr, ar_mjd, prefit_rms[0], postfit_rms[0], med_toaE[0], dm[0], dmerror[0], chisq[0]), fontsize=10, fontweight='bold')


		dirplot=os.path.join(pwd,ar_psr+"_"+ar_tel+"_plots")
		if not os.path.exists(dirplot):
		   os.makedirs(dirplot)
		plotfile=dirplot+"/"+ar_psr+"_"+str(ar_mjd)+"_"+str(centre_freqs[0])+"_"+ar_tel+"_DMfitResid.pdf"
		plt.savefig(plotfile, format='pdf')
		plt.close()

	if (len(finalarchives) == 1):
		prof2D = []
		prof1D = []
		condition = orig_freqs < 500.
			
		ar_nchan = finalarchives[0].get_nchan()
		ar_nbin  = finalarchives[0].get_nbin()
		ar_bw = finalarchives[0].get_bandwidth()
		ar_freq = finalarchives[0].get_centre_frequency()
		prof2D = finalarchives[0].get_data()[:,0,:,:].flatten().reshape(ar_nchan,ar_nbin)
		prof = finalarchives[0].clone()
		prof.fscrunch()
		prof1D = prof.get_data().flatten()
		prof1D /= np.max(prof1D)

		fig,axs = plt.subplots(2, sharex=True, figsize=(8, 6))
		fig.subplots_adjust(hspace=0.05)
		ax0 = plt.subplot2grid((9, 8), (0,0), rowspan=7, colspan=3)
		ax1 = plt.subplot2grid((9, 8), (7,0), rowspan=2, colspan=3)

		ax2 = plt.subplot2grid((9, 8), (0,4), colspan=4, rowspan=3)
		ax3 = plt.subplot2grid((9, 8), (3,4), colspan=4, rowspan=3)
		ax4 = plt.subplot2grid((9, 8), (6,4), colspan=4, rowspan=3)
		leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)
		ax0.imshow((np.sqrt(prof2D**2))**0.5, origin='lower', extent=(0,ar_nbin-1,(np.around(ar_freq)-ar_bw/2),(np.around(ar_freq)+ar_bw/2)), aspect='auto', cmap='hot')
		ax0.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=8)
		ax0.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
		ax1.plot(np.arange(ar_nbin, dtype=float),prof1D, color='black', linewidth=0.5)
		ax1.set_xlim(0,ar_nbin-1)
		ax1.set_ylabel('Intensity', fontweight='bold', fontsize=8)

		if ar_freq < 500.:
			ax2.errorbar(orig_freqs, orig_resid, yerr=orig_toasE, fmt='.', color='#D81B60', capsize=2)
		if ar_freq > 1000.:
			ax2.errorbar(orig_freqs, orig_resid, yerr=orig_toasE, fmt='.', color='#1E88E5', capsize=2)
		ax2.set_xlim((np.around(ar_freq)-ar_bw/1.5), (np.around(ar_freq)+ar_bw/1.5))
		ax2.grid()
		ax2.legend([leg1], ['Prefit: Unfiltered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax2.axes.xaxis.set_ticklabels([])

		if ar_freq < 500.:
			ax3.errorbar(init_freqs, init_resid, yerr=init_toasE, fmt='.', color='#D81B60', capsize=2)
		if ar_freq > 1000.:
			ax3.errorbar(init_freqs, init_resid, yerr=init_toasE, fmt='.', color='#1E88E5', capsize=2)
		ax3.grid()
		ax3.set_xlim((np.around(ar_freq)-ar_bw/1.5), (np.around(ar_freq)+ar_bw/1.5))
		ax3.legend([leg1], ['Prefit: Filtered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax3.axes.xaxis.set_ticklabels([])
		ax3.set_ylabel(r'ToA Residuals ($\mu$s)', fontweight='bold', fontsize=8)
		
		if ar_freq < 500.:
			ax4.errorbar(final_freqs, final_resid, yerr=final_toasE, fmt='.', color='#D81B60', capsize=2)
		if ar_freq > 1000.:
			ax4.errorbar(final_freqs, final_resid, yerr=final_toasE, fmt='.', color='#1E88E5', capsize=2)
		ax4.grid()
		ax4.set_xlim((np.around(ar_freq)-ar_bw/1.5), (np.around(ar_freq)+ar_bw/1.5))
		ax4.legend([leg1], ['Postfit'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax4.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=9)
		fig.suptitle('Source: PSR %s; MJD: %.4f; Prefit Wrms: %.2f $\mu$s; Postfit Wrms: %.2f $\mu$s\nMedian ToA Err: %.2f $\mu$s; DM: %.6f $\pm$ %.6f pc cm$^{-3}$;  Reduced $\chi^2$: %.2f'%(ar_psr, ar_mjd, prefit_rms[0], postfit_rms[0], med_toaE[0], dm[0], dmerror[0], chisq[0]), fontsize=10, fontweight='bold')
		
		#axs[0].set_ylabel(r'Prefit residuals ($\mu$s)', fontweight='bold', fontsize=12)
		#axs[0].errorbar(init_freqs,init_resid,init_toasE,fmt='.k', capsize=2)	
		#axs[1].errorbar(final_freqs,final_resid,final_toasE,fmt='.k', capsize=2)
		#axs[1].set_ylabel(r'Postfit residuals ($\mu$s)', fontweight='bold', fontsize=12)
		#axs[1].set_xlabel(r'Frequency (MHz)', fontweight='bold', fontsize=12)
		fig.suptitle('Source: PSR %s; MJD: %.4f; Prefit Wrms: %.2f $\mu$s; Postfit Wrms: %.2f $\mu$s\nMedian ToA Err: %.2f $\mu$s; DM: %.6f $\pm$ %.6f pc cm$^{-3}$;  Reduced $\chi^2$: %.2f'%(ar_psr, ar_mjd, prefit_rms[0], postfit_rms[0], med_toaE[0], dm[0], dmerror[0], chisq[0]), fontsize=11, fontweight='bold')

		dirplot=os.path.join(pwd,ar_psr+"_"+ar_tel+"_plots")
		if not os.path.exists(dirplot):
		   os.makedirs(dirplot)
		plotfile=dirplot+"/"+ar_psr+"_"+str(ar_mjd)+"_"+str(centre_freqs[0])+"_"+ar_tel+"_DMfitResid.pdf"
		plt.savefig(plotfile, format='pdf')
		plt.close()

	import time
	end = time.time()
	total = end - start
	print ('\n-----------------------------------------------------------------------------')
	print ('MJD\t\tDM\t\tDMerr\t\tChisq\tC_Fr\tBW\tTel')
	print ('%.6f\t%.6f\t%.6f\t%.2f\t%.1f\t%.1f\t%s' % (ar_mjd, dm[0], dmerror[0], 
			chisq[0], centre_freqs[0], bw[0], ar_tel) )
	
	print ('-----------------------------------------------------------------------------')

	print("\nThe program took %.1f seconds to finish"%total)

# Function that obtains the ToAs for DM estimation
def get_TOAs(ar, std, ephemeris, quiet):
	
	init_dm = ar.get_dispersion_measure()
	if not quiet:
		print("Calculating DM of %s"%(ar.get_filename()))
		print("Using the ArrivalTime with FDM:mcmc=1 in Tempo2 format")
	arrtim = psrchive.ArrivalTime()
	arrtim.set_shift_estimator('FDM:mcmc=1')
	arrtim.set_format('Tempo2')
	arrtim.set_format_flags('IPTA')
	if not quiet:
		print("Loading the template file for processing... "),
	ar_psr = ar.get_source()
	ar_bw  = ar.get_bandwidth()
	ar_mjd = ar.get_Integration(0).get_start_time().in_days()
	ar_frq = ar.get_centre_frequency()
	pta_flag = '-pta InPTA'
	sys_flag = ''
	grp_flag = ''
	bnd_flag = ''
	if (ar_mjd < 58600.):
		if (300. < ar_frq < 500.):
			if (ar_bw == 200.):
				grp_flag = '-group GM_GWB_500_200_b1_pre36'
				sys_flag = '-sys GM_GWB_500_200_b1'
				bnd_flag = '-bandno 3'
			if (ar_bw == 100.):
				grp_flag = '-group GM_GWB_500_100_b1_pre36'
				sys_flag = '-sys GM_GWB_500_100_b1'
				bnd_flag = '-bandno 3'
		if (1260. < ar_frq < 1460.):
			if (ar_bw == 200.):
				grp_flag = '-group GM_GWB_1460_200_b0_pre36'
				sys_flag = '-sys GM_GWB_1460_200_b0'
				bnd_flag = '-bandno 5'
			if (ar_bw == 100.):
				grp_flag = '-group GM_GWB_1460_100_b1_pre36'
				sys_flag = '-sys GM_GWB_1460_100_b1'
				bnd_flag = '-bandno 5'
	if (ar_mjd > 58600.):
		if (300. < ar_frq < 500.):
			if (ar_bw == 200.):
				grp_flag = '-group GM_GWB_500_200_b1_post36'
				sys_flag = '-sys GM_GWB_500_200_b1'
				bnd_flag = '-bandno 3'
		if (1260. < ar_frq < 1460.):
			if (ar_bw == 200.):
				grp_flag = '-group GM_GWB_1460_200_b0_post36'
				sys_flag = '-sys GM_GWB_1460_200_b0'
				bnd_flag = '-bandno 5'
	ar_nbin  = ar.get_nbin()
	std_nbin = std.get_nbin()
	nbin = ar_nbin
	if (ar_nbin > std_nbin):
		nbin = std_nbin
	std.bscrunch_to_nbin(nbin)
	std.pscrunch()
	std.tscrunch()
	std.dedisperse()
	arrtim.set_standard(std)
	if not quiet:
		print(" done!")
	ar.bscrunch_to_nbin(nbin)
	ar.pscrunch()
	ar.tscrunch()
	arrtim.set_observation(ar)
	if not quiet:
		print("Finding the ToAs... "),
	# Finding the ToAs and reading it into numpy arrays
	toas = arrtim.get_toas()
	toas_filtered = [x.split()[:5] for x in toas] 
	str_filename,str_freq,str_mjd,str_toaErr,str_site = zip(*toas_filtered)
	freq = np.asarray(str_freq, dtype=np.float64)
	amjd = np.asarray(str_mjd, dtype=np.float64)
	terr = np.asarray(str_toaErr, dtype=np.float64)
	# removing the ToAs with zero errors
	condition1 = terr != 0.
	freqnew = np.extract(condition1,freq)
	amjdnew = np.extract(condition1,amjd)
	terrnew = np.extract(condition1,terr)
	if not quiet:
		print(" done!")
		print("Removing the bad ToAs using Huber Regression... "),
	# removing the ToAs with 3 sigma errorbars 
	condition1 = terrnew < 3*np.median(terrnew)
	freqnew = np.extract(condition1,freqnew)
	amjdnew = np.extract(condition1,amjdnew)
	terrnew = np.extract(condition1,terrnew)
	# writing the ToAs to a temporary file for getting the non-phase resolved ToAs using general2
	tempfile = ar.get_source()+"_tmp.txt"
	f = open(tempfile,"w+")
	head="FORMAT 1\n"
	f.write('%s' % head)
	for i in range(0,np.size(freqnew)):
		f.write('%s %.8f %.18f %.6f %s %s %s %s %s\n' % (str_filename[0], freqnew[i], 
				amjdnew[i], terrnew[i], str_site[0], pta_flag, sys_flag, grp_flag, bnd_flag))
	f.close()
	tmp = os.popen("tempo2 -output general2 -f %s %s -s \"1111111 {freq} {pre} {err}\n\" | grep '1111111'" 
					% (ephemeris,tempfile)).read()
	#os.remove(tempfile)

	# extracting the data from general2 output
	tmp1 = tmp.split('\n')
	freqtmp = np.zeros(np.size(amjdnew), dtype=np.float64)
	toastmp = np.zeros(np.size(amjdnew), dtype=np.float64)
	TErrtmp = np.zeros(np.size(amjdnew), dtype=np.float64)
	for i in range(np.size(amjdnew)):
		_,freqtmp[i],toastmp[i],TErrtmp[i] = (tmp1[i].split())
	TErrtmp /= 1e+6

	# importing libraries for outlier removal
	from sklearn import linear_model
	from sklearn.linear_model import HuberRegressor
	from sklearn.preprocessing import PolynomialFeatures
	from sklearn.pipeline import make_pipeline
	# changing the shape of frequency array
	freqarr = freqtmp.reshape(-1,1)
	# making a polynomial model	and fitting using Huber Regression
	toastmp *= 1e+6
	toashift = (np.min(toastmp)*-1.5)
	toastmp += toashift
	Terrtmp = TErrtmp*1e+6
	model = make_pipeline(PolynomialFeatures(2), HuberRegressor())
	model.fit(freqarr,toastmp,
			  huberregressor__sample_weight=np.ravel(1./Terrtmp))
	y_pred = model.predict(freqarr)
	residuals = toastmp - y_pred
	median = np.median(residuals)
	MAD = np.median(np.abs(residuals-np.median(residuals)))/0.6744897501960817
	# filtering the good ToAs
	condition2 = (residuals > median - 3*MAD) & (residuals < median + 3*MAD)
	freqf = np.extract(condition2,freqarr)
	if (len(freqf) == 0):
		freqf = freqarr
	if (len(freqf) != 0):
		amjdf = np.extract(condition2,amjdnew)
		toasf = np.extract(condition2,toastmp)
		terrf = np.extract(condition2,TErrtmp)
	terrf *= 1e+6
	if not quiet:
		print(" done!")
	return(freqtmp, toastmp, TErrtmp, freqf, amjdf, toasf, terrf)

# Function to get the final ToAs with proper uGMRT flags
def get_finalTOA(ar, std, ephemeris, templ_dm, quiet):
	
	if not quiet:
		print("Getting the final frequency resolved ToAs of %s..."%(ar.get_filename()))
	temppar = "tempnoefac.par"
	with open(ephemeris,"r") as ft:
		lines = ft.readlines()
	with open(temppar,"w") as ff:
		for line in lines:
			if not line.startswith("T2EFAC"):
				ff.write(line)
	ff.close()
	ar.set_ephemeris(temppar)
	ar.set_dispersion_measure(templ_dm)
	ar.update_model()
	ar_nchan = ar.get_nchan()
	snr = np.zeros(ar_nchan)
	artmp = ar.clone()
	for i in range(ar_nchan):
		snr[i] = artmp.get_Profile(0,0,i).snr()
	del (artmp)

	ar_psr = ar.get_source()
	ar_bw  = ar.get_bandwidth()
	ar_mjd = ar.get_Integration(0).get_start_time().in_days()
	ar_frq = ar.get_centre_frequency()
	pta_flag = '-pta InPTA'
	sys_flag = ''
	grp_flag = ''
	bnd_flag = ''
	if (ar_mjd < 58600.):
		if (300. < ar_frq < 500.):
			if (ar_bw == 200.):
				grp_flag = '-group GM_GWB_500_200_b1_pre36'
				sys_flag = '-sys GM_GWB_500_200_b1'
				bnd_flag = '-bandno 3'
			if (ar_bw == 100.):
				grp_flag = '-group GM_GWB_500_100_b1_pre36'
				sys_flag = '-sys GM_GWB_500_100_b1'
				bnd_flag = '-bandno 3'
		if (1260. < ar_frq < 1460.):
			if (ar_bw == 200.):
				grp_flag = '-group GM_GWB_1460_200_b0_pre36'
				sys_flag = '-sys GM_GWB_1460_200_b0'
				bnd_flag = '-bandno 5'
			if (ar_bw == 100.):
				grp_flag = '-group GM_GWB_1460_100_b1_pre36'
				sys_flag = '-sys GM_GWB_1460_100_b1'
				bnd_flag = '-bandno 5'
	if (ar_mjd > 58600.):
		if (300. < ar_frq < 500.):
			if (ar_bw == 200.):
				grp_flag = '-group GM_GWB_500_200_b1_post36'
				sys_flag = '-sys GM_GWB_500_200_b1'
				bnd_flag = '-bandno 3'
		if (1260. < ar_frq < 1460.):
			if (ar_bw == 200.):
				grp_flag = '-group GM_GWB_1460_200_b0_post36'
				sys_flag = '-sys GM_GWB_1460_200_b0'
				bnd_flag = '-bandno 5'

	tmp_filename = os.path.basename(std.get_filename())
	std.set_filename(tmp_filename)
	ar_nbin  = ar.get_nbin()
	std_nbin = std.get_nbin()
	nbin = ar_nbin
	if (ar_nbin > std_nbin):
		nbin = std_nbin

	arrtim = psrchive.ArrivalTime()
	arrtim.set_shift_estimator('FDM:mcmc=1')
	arrtim.set_format('Tempo2')
	arrtim.set_format_flags('IPTA')
	std.bscrunch_to_nbin(nbin)
	std.pscrunch()
	std.tscrunch()
	std.dedisperse()
	arrtim.set_standard(std)
	tmp_filename = os.path.basename(ar.get_filename())
	ar.set_filename(tmp_filename)
	ar.bscrunch_to_nbin(nbin)
	ar.pscrunch()
	ar.tscrunch()
	arrtim.set_observation(ar)
	# Finding the ToAs and reading it into numpy arrays
	toas = arrtim.get_toas()
	toas_filtered = [x.split()[:] for x in toas]
	finaltoasfile = ar_psr+"_allToAs.tim"
	if not (os.path.isfile(finaltoasfile)):
		ft = open(finaltoasfile,"a")
		head="FORMAT 1"
		ft.write('%s\n' % head)
		ft.close()
	if (os.path.isfile(finaltoasfile)):
		ft = open(finaltoasfile,"a")

	for i in range(len(toas)):
		if (snr[i]>=8):
			ft.write('%s -prof_snr %.2f %s %s %s %s\n' % (toas[i], snr[i], pta_flag, sys_flag, grp_flag, bnd_flag))
		if (snr[i]<8):
			ft.write('C %s -prof_snr %.2f %s %s %s %s\n' % (toas[i], snr[i], pta_flag, sys_flag, grp_flag, bnd_flag))

	if not quiet:
		print(" done!")

# Module that splits the ToAs of different bands for obtaining DM
def get_DM(freqs, toas, toaE, archives, ephemeris):
	bands = np.zeros(len(archives))
	bandw = np.zeros(len(archives))
	files = []
	site = archives[0].get_telescope()
	dm  = []
	dmE = []
	chi = []
	ta1,ta2,ta3 = dmfit("temp.ar", freqs, toas, toaE, site, ephemeris)
	dm.append(ta1); dmE.append(ta2); chi.append(ta3)	
	if ((site == 'GMRT' or site == 'gmrt') and (len(archives) > 1)):
		band3_file=[]; band5_file=[]
		for i in range(len(archives)):
			if (300. < archives[i].get_centre_frequency() < 500.):
				band3_file = archives[i].get_filename()
			if (1160. < archives[i].get_centre_frequency() < 1460.):
				band5_file = archives[i].get_filename()
		condition = (freqs < 500.) & (freqs > 300.)
		band3_freqs = np.extract(condition,freqs)
		band3_toas  = np.extract(condition,toas)
		band3_toasE = np.extract(condition,toaE)
		if (np.size(band3_freqs) != 0):
			tb1,tb2,tb3 = dmfit(band3_file, band3_freqs, band3_toas, band3_toasE, site, ephemeris)
			dm.append(tb1); dmE.append(tb2); chi.append(tb3)
		condition = (freqs > 1160.) & (freqs < 1460.)
		band5_freqs = np.extract(condition,freqs)
		band5_toas  = np.extract(condition,toas)
		band5_toasE = np.extract(condition,toaE)
		if (np.size(band5_freqs) != 0):
			tc1,tc2,tc3 = dmfit(band5_file, band5_freqs, band5_toas, band5_toasE, site, ephemeris)
			dm.append(tc1); dmE.append(tc2); chi.append(tc3)	
	return(dm,dmE,chi)

# Module that estimates DM with tempo2
def dmfit(filename, freqs, toas, toaE, tel, ephemeris):
	tempfile = filename+"_"+str(freqs[0])+"_toas.txt"
	f = open(tempfile,"w+")
	head="FORMAT 1\n"
	f.write('%s' % head)
	for i in range(np.size(freqs)):
		f.write('%s %.8f %.18f %.6f %s\n' % (filename, freqs[i], toas[i], toaE[i], tel))
	f.close()
	awkstr = "-nofit -fit dm | grep 'DM (cm^-3 pc)'| awk \'{print $5,$6}\'"
	dmstr=os.popen("tempo2 -f %s %s %s" % (ephemeris, tempfile, awkstr)).read()
	(dm, dmerr) = dmstr.split()
	dmval = float(dm)
	dmverr = float(dmerr)
	# doing the fit again to read the chisquare
	chisqstr=os.popen("tempo2 -f %s %s -nofit -fit dm | grep 'Fit Chisq'| awk \'{print $9}\'" % (ephemeris, tempfile)).read()
	fitchisq = float(chisqstr)
	os.remove(tempfile)
	return(dmval,dmverr,fitchisq)

# Function to correct the backend delay for pre pinta V6.2 processed data
def Correct_delay(ar):
	if (ar.get_telescope() == 'GMRT' or ar.get_telescope() == 'gmrt'):
		ar_mjd = ar.get_Integration(0).get_start_time().in_days()
		ar_frq  = ar.get_centre_frequency()
		period = (ar.get_Integration(0).get_folding_period())
		# cycle 34-35
		if (ar_mjd >= 58230. and ar_mjd < 58550.):
			if (ar_frq> 300. and ar_frq < 500.):
				ar.unload('temp1.ar')
				tt=os.popen('psredit -c be:delay=-4.02653184 -m temp1.ar').read()
				ar = psrchive.Archive_load('temp1.ar')
				ar.update_model()
				os.remove('temp1.ar')
			if (ar_frq> 1200. and ar_frq < 1500.):
				ar.unload('temp2.ar')
				tt=os.popen('psredit -c be:delay=-4.02653184 -m temp2.ar').read()
				ar = psrchive.Archive_load('temp2.ar')
				ar.update_model()
				os.remove('temp2.ar')
		# cycle 37
		if (ar_mjd >= 58810. and ar_mjd < 58991.):
			if (ar_frq> 300. and ar_frq < 500.):
				ar.unload('temp1.ar')
				tt=os.popen('psredit -c be:delay=-2.01326592 -m temp1.ar').read()
				ar = psrchive.Archive_load('temp1.ar')
				ar.update_model()
				os.remove('temp1.ar')
			if (ar_frq> 1200. and ar_frq < 1500.):
				ar.unload('temp2.ar')
				tt=os.popen('psredit -c be:delay=-1.34217728 -m temp2.ar').read()
				ar = psrchive.Archive_load('temp2.ar')
				ar.update_model()
				os.remove('temp2.ar')
	return(ar)


main()
