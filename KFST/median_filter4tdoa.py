#!/usr/bin/env python


'''
Function: Run moving median filter over BTK TDOA values.

Example folder:
	/Users/yulan/tmp/btk/KFST/11May2016/

Example command:
	median_filter4tdoa.py  $LISTOFTDOA  $NWIN  $OUTDIR

Yulan Liu, 24 May 2016
'''

import sys, os, subprocess, numpy
import scipy
import scipy.signal
from scipy.signal import medfilt 


usage = 'median_filter4tdoa.py  <file listing path to TDOA files>  <size of window>  <output folder>'

if len(sys.argv)<4:
    print usage
    exit()


def load_tdoa(tdoa_fname):
    '''
    Read TDOA from files.
    '''
    fid = open(tdoa_fname, 'r')
    T = []
    TDOAs = dict()
    NR = None
    for line in fid:
	line = line.strip()
	if line=='' or line[0]=='#':
	    continue
#   0.00 4 7 -5.000000e-04 4 6 -1.250000e-04 4 5 -6.250000e-05 5 6 0.000000e+00 5 7 -3.750000e-04 2 3 1.875000e-04 2 5 1.250000e-04 2 4 2.500000e-04 2 7 -1.875000e-04 2 6 1.250000e-04 3 4 6.250000e-05 3 5 6.250000e-05 3 6 1.875000e-04 3 7 -3.125000e-04 0 3 3.750000e-04 0 2 1.250000e-04 0 1 0.000000e+00 0 7 -6.250000e-05 0 6 1.875000e-04 0 5 3.750000e-04 0 4 4.375000e-04 1 2 0.000000e+00 1 3 3.750000e-04 1 6 3.750000e-04 1 7 -6.250000e-05 1 4 5.000000e-04 1 5 3.750000e-04 6 7 -1.875000e-04 

	data = line.split()
	if NR is None:
	    NR = len(data)
	else:
	    if not NR==len(data):
		print '[WARNING] Current line differs in number of fields, supposed to be %d, but having %d ' % (NR, len(data))
		print '		 To avoid the problem caused to median filter, skipping line: \n', line
		continue

	t = float(data[0])
	T.append(t)
	tmp = data[1:]
	for i in xrange(0, len(tmp), 3):
	    tag = '%s %s' % (tmp[i], tmp[i+1])
	    tdoa = float(tmp[i+2])
	    try:
		TDOAs[tag]
	    except KeyError:
		TDOAs[tag] = []
	    TDOAs[tag].append(tdoa)
	
    # A quick double check whether the number of samples of different tags equals
    tag_list = sorted(TDOAs.keys())
    try:
    	N = len(TDOAs[tag_list[0]])
    except IndexError:
	print '[WARNING] IndexError in "N = len(TDOAs[tag_list[0]])" when processing file: ', tdoa_fname
	print '          This might be caused by an empty input TDOA file. Skipping file ', tdoa_fname
	print '[DEBUG] tag_list = ', tag_list
	print '[DEBUG] TDOAs = ', TDOAs
	return {}, []
    for i in xrange(len(tag_list)):
	assert(N==len(TDOAs[tag_list[i]]))

    fid.close()
    return TDOAs, T


def write_tdoa(TDOAs, T, tdoa_fname):
    '''
    Save TDOA in files.
    '''
    if len(T)==0 or len(TDOAs.keys())==0:
	print '[WARNING] Skipping writing file: ', tdoa_fname
    else:
    	fid = open(tdoa_fname, 'w')
    	tag_list = sorted(TDOAs.keys())
    	for i in xrange(len(T)):
	    t = T[i]
	    fid.write('%.2f' % t)
	    for tag in tag_list:
	    	tdoa = TDOAs[tag][i]
	    	fid.write(' %s %g' %(tag, tdoa))
	    fid.write('\n')
    	fid.close()


def median_filter_seq(in_seq, Nwin):
    '''
    Apply median filter to a sequence / list.
    '''
    if numpy.mod(Nwin, 2)==0:
	print '[WARNING] The size of median filter window should be odd number. Changing window size from %d to %d.' % (Nwin, Nwin+1)
    	Nwin += 1
#    print '[DEBUG] in_seq = ', in_seq
#    print '[DEBUG] Nwin = ', Nwin
    out_seq = medfilt(numpy.array(in_seq), Nwin)
    return out_seq


def median_filter_tdoa(tdoa_fname, Nwin, out_fname):
    '''
    Apply median filter to one TDOA file, and save the output.
    '''
    print '%s -> %s' % (tdoa_fname, out_fname)
    TDOAs, T = load_tdoa(tdoa_fname)
    tag_list = sorted(TDOAs.keys())
    TDOAs_out = dict()
    for tag in tag_list:
	TDOAs_out[tag] = median_filter_seq(TDOAs[tag], Nwin)

    write_tdoa(TDOAs_out, T, out_fname)




def main():
    scp_fname = sys.argv[1]
    Nwin = int(sys.argv[2])
    outdir = sys.argv[3]
    subprocess.call('mkdir -p '+outdir, shell=True)
    fid = open(scp_fname, 'r')
    for line in fid:
	line = line.strip()
	if line=='' or line[0]=='#':
	    continue
    	tdoa_fname = line
	out_fname = outdir + '/'+tdoa_fname.split('/')[-1]

	median_filter_tdoa(tdoa_fname, Nwin, out_fname)

    fid.close()


if __name__=='__main__':
    main()

