#!/usr/bin/env python


'''
Function: Smooth the TDOA file from KFST in one of the following modes:

	1. every given X seconds
	2. every given utterance / segment (segmentation provided externally)

Example folder:
	/share/spandh.ami1/usr/yulan/dev/btk/tools/KFST/dbg/cir_test/bigwav/

Example command:
	smoothTDOA.py  --scp=smooth_vad.scp  --mode=0

	The scp file is expected to have 3 columns when mode is 0: 

	<input .tdoa file>  <output .tdoa file>  <corresponding segmatation>

	The scp file is expected to have at least 2 columns when mode is a
	positive number X:

	<input .tdoa file>  <output .tdoa file>

Yulan Liu, 2 May 2016
'''

import sys, subprocess, optparse, math
from optparse import OptionParser
import numpy



def load_tdoa(fname):
   '''
   Load the information of one TDOA file.
   '''
# [yulan@snarl out3]$ head -1 AMI-N1005A_f4041.tdoa 
#  0.00 4 7 -5.000000e-04 4 6 -1.250000e-04 4 5 -6.250000e-05 5 6 0.000000e+00 5 7 -3.750000e-04 2 3 1.875000e-04 2 5 1.250000e-04 2 4 2.500000e-04 2 7 -1.875000e-04 2 6 1.250000e-04 3 4 6.250000e-05 3 5 6.250000e-05 3 6 1.875000e-04 3 7 -3.125000e-04 0 3 3.750000e-04 0 2 1.250000e-04 0 1 0.000000e+00 0 7 -6.250000e-05 0 6 1.875000e-04 0 5 3.750000e-04 0 4 4.375000e-04 1 2 0.000000e+00 1 3 3.750000e-04 1 6 3.750000e-04 1 7 -6.250000e-05 1 4 5.000000e-04 1 5 3.750000e-04 6 7 -1.875000e-04 

   fid = open(fname, 'r')
   infodict = dict()
   for line in fid:
	line = line.strip()
	if line=='' or line[0]=='#':
	    continue
	t = float(line.split()[0])
	tmp = line.split()[1:]

	try:
	    infodict[t]
	    print '\n[WARNING] Duplicated record for time: ', t, ' in file: ', fname
	    print '          Existing record: ', infodict[t]
	    print '          Overwite with line: ', line
	except KeyError:
	    infodict[t] = dict()

	for i in xrange(0, len(tmp), 3):
	    pair = tmp[i]+' '+tmp[i+1]
	    tdoa = tmp[i+2]
	    try:
		infodict[t][pair]
	    except KeyError:
		infodict[t][pair] = []
	    infodict[t][pair].append(tdoa)
   fid.close()
   return infodict


def smooth_vad(infname, outfname, vadfname):
    '''
    Smooth one TDOA file according to given segmenation information.
    Note that it only smooth the postion record within one utterance.
    '''
    print infname, ' -> ', outfname, '\n\twith segmentation from ', vadfname
   
    # Buffer the segmentation information
    buff = []
    buff_str = []
    fid = open(vadfname, 'r')
    for line in fid:
	line = line.strip()
	if line=='' or line[0]=='#':
	    continue
	data = line.split()
	assert(len(data)==2)
	buff.append([float(x) for x in data])
   	buff_str.append(data)
    fid.close()
     
    # Load TDOA file
    infodict = load_tdoa(infname)

    # Smoothing and saving results
    buff_idx = 0
    tdoa_buff = dict()
    outdir = '/'.join(outfname.split('/')[:-1])
    subprocess.call('mkdir -p '+outdir, shell=True)
    outfid = open(outfname, 'w')
    WriteFlag = 0
    aver_tdoa = dict()
    for t in sorted(infodict.keys()):
	tdoa_data_tmp = infodict[t]
	t = float(t)
	try:
	    if t<buff[buff_idx][0]:		# One utterance hasn't started yet, copy the codes.
		outfid.write('%.2f '%t)
		for pair in sorted(tdoa_data_tmp.keys()):
		    tdoa_tmp = tdoa_data_tmp[pair][0]
		    outfid.write(' %s %s' % (pair, tdoa_tmp))
		outfid.write('\n')
		aver_tdoa = {}
		tdoa_buff = {}
	    elif t<buff[buff_idx][1]:		# In the middle of one utterance
		for pair in sorted(tdoa_data_tmp.keys()):
		    try:
			tdoa_buff[pair]
		    except KeyError:
			tdoa_buff[pair] = []
		    tdoa_buff[pair] += tdoa_data_tmp[pair]
	    	continue
	    elif t>=buff[buff_idx][1]:
	    	if t==buff[buff_idx][1] or len(tdoa_buff)==0:
                    for pair in sorted(tdoa_data_tmp.keys()):
                    	try:
                            tdoa_buff[pair]
                    	except KeyError:
                            tdoa_buff[pair] = []
                    	tdoa_buff[pair] += tdoa_data_tmp[pair]
		for pair in tdoa_buff.keys():
		    aver_tdoa[pair] = numpy.mean([float(x) for x in tdoa_buff[pair]], axis=0)
		outfid.write('%.2f '%t)
                for pair in sorted(aver_tdoa.keys()):
                    outfid.write(' %s %f' % (pair, aver_tdoa[pair]))
                outfid.write('\n')
                aver_tdoa = {}
		tdoa_buff = {}

	    	if t==buff[buff_idx][1]:
		    buff_idx += 1
	    	    continue
	    	else:
		    buff_idx += 1
		    try:
		    	buff[buff_idx]
                    	for pair in sorted(tdoa_data_tmp.keys()):
                            try:
                            	tdoa_buff[pair]
                            except KeyError:
                            	tdoa_buff[pair] = []
                            tdoa_buff[pair] += tdoa_data_tmp[pair]
		    	continue
		    except IndexError:
		    	WriteFlag = 1
        except IndexError:
            # Already finished all segments, going to copy the raw position as it is.
            for pair in sorted(tdoa_data_tmp.keys()):
                try:
                    tdoa_buff[pair]
                except KeyError:
                    tdoa_buff[pair] = []
                tdoa_buff[pair] += tdoa_data_tmp[pair]

            if len(tdoa_buff)>0:
                for pair in tdoa_buff.keys():
                    aver_tdoa[pair] = numpy.mean([float(x) for x in tdoa_buff[pair]], axis=0)
            WriteFlag = 1

        if WriteFlag==1:
            outfid.write('%.2f '%t)
            for pair in sorted(aver_tdoa.keys()):
                outfid.write(' %s %f' % (pair, aver_tdoa[pair]))
            outfid.write('\n')
            aver_tdoa = {}
            tdoa_buff = {}

	    WriteFlag = 0
	    continue

    if len(tdoa_buff)>0:
        for pair in tdoa_buff.keys():
            aver_tdoa[pair] = numpy.mean([float(x) for x in tdoa_buff[pair]], axis=0)
        outfid.write('%.2f '%t)
        for pair in sorted(aver_tdoa.keys()):
            outfid.write(' %s %f' % (pair, aver_tdoa[pair]))
        outfid.write('\n')
        aver_tdoa = {}
        tdoa_buff = {}
    outfid.close()




def smooth_fixed(infname, outfname, win):
   '''
   Smooth one TDOA file with a given local window size.
   '''
   print infname, ' -> ', outfname
   # Load TDOA file
   infodict = load_tdoa(infname)
   tdoa_buff = {}; aver_tdoa = {}
   outdir = '/'.join(outfname.split('/')[:-1])
   subprocess.call('mkdir -p '+outdir, shell=True)
   outfid = open(outfname, 'w')
   i = 1
   for t in sorted(infodict.keys()):
        tdoa_data_tmp = infodict[t]
	t = float(t)
	if t>=i*win and len(tdoa_buff)>0:
            for pair in tdoa_buff.keys():
		try:
                    aver_tdoa[pair] = numpy.mean([float(x) for x in tdoa_buff[pair]], axis=0)
		except:
		    print '[DEBUG] tdoa_buff[', pair, '] = ', tdoa_buff[pair]
		    exit()
            outfid.write('%.2f '%((i-1)*win))
            for pair in sorted(aver_tdoa.keys()):
            	outfid.write(' %s %f' % (pair, aver_tdoa[pair]))
            outfid.write('\n')
	    tdoa_buff = {}; tdoa_buff = {}

	    i += 1
	    continue
	else:
            for pair in sorted(tdoa_data_tmp.keys()):
                try:
                    tdoa_buff[pair]
                except KeyError:
                    tdoa_buff[pair] = []
                tdoa_buff[pair] += tdoa_data_tmp[pair]

   if len(tdoa_buff)>0:
        for pair in tdoa_buff.keys():
            aver_tdoa[pair] = numpy.mean([float(x) for x in tdoa_buff[pair]], axis=0)

        outfid.write('%.2f '%t)
        for pair in sorted(aver_tdoa.keys()):
            outfid.write(' %s %f' % (pair, aver_tdoa[pair]))
        outfid.write('\n')
        tdoa_buff = {}; tdoa_buff = {}
   outfid.close()


def smooth_all(options):
    if abs(float(options.mode))<1e-5:		# Almost 0
	print 'Going to smooth file with external segmenation information.'
	fid = open(options.scp, 'r')
	for line in fid:
	    line = line.strip()
	    if line=='' or line[0]=='#':
		continue
	    data = line.split()
	    assert(len(data)>=3)
	    smooth_vad(data[0], data[1], data[2])
	fid.close()
    elif float(options.mode)>0:
	win = float(options.mode)
	print 'Going to smooth file every ', win, ' seconds.'
	fid = open(options.scp, 'r')
	for line in fid:
	    line = line.strip()
	    if line=='' or line[0]=='#':
		continue
	    data = line.split()
	    assert(len(data)>=2)
	    smooth_fixed(data[0], data[1], win)
	fid.close()
    else:
	print '[ERROR] Unrecognized mode: ', options.mode
	print '        It has to be either positive value or 0. Try --help for more information.'
	exit()


def parse_option(cmd):
    parser = OptionParser(usage='smoothTDOA.py  --scp=<scp file>  --mode=X', version='0.1')
    parser.add_option('--scp', dest='scp', help='List of input TDOA files, output TDOA files and segmentation files (when in mode 0).', default=None)
    parser.add_option('--mode', dest='mode', help='Smoothing mode: 0 - smoothing by external segmentation; ANY POSITIVE VALUE X - smoothing every X seconds..', default=0)

    (options, args) = parser.parse_args(cmd)
    assert(options.scp is not None)
    return options


def main():
    usage = 'smoothTDOA.py  --scp=<scp>  --mode=X\n'
    if len(sys.argv)<3 and not (sys.argv[1] in ['-h', '--help']):
	print usage
        print 'Please try option "--help" to get help information.'
	exit()

    options = parse_option(sys.argv)
    smooth_all(options)


if __name__ == '__main__':
    main()


