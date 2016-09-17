#!/usr/bin/env python

'''
Function: convert the raw location output from KFST into TDOA, to see how 
	much difference there is. This script only supports circular array.

Yulan Liu, 18 May 2016
'''

import os, sys, subprocess
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../')
from utils import load_geometry, pos2tdoa_circular

usage = 'txt2tdoa.py  <list of input tdoa file full path>  <microphone geometry file>  <output folder>'
if len(sys.argv)<4:
    print usage
    exit()


def txt2tdoa_circular(infile, outfile, geometry, refX=0):
    '''
    Convert one .txt file (with arrival angles) into .tdoa file, for circular array only.
    By default using the first channel as reference channel
    '''
    Nch = len(geometry)
    fid = open(infile, 'r')
    out_fid = open(outfile, 'w')
    for line in fid:
        line = line.strip()
        if line=='' or line[0]=='#':
            continue
	else:
            # Time   0.00 : (theta, phi) = (8.480370e-01, -5.294148e-02)
            data = line.split()
            t = float(data[1])
            theta, phi = [float(x) for x in line.split('(')[-1].split(')')[0].split(',')]
            #   0.00 4 7 -5.000000e-04 4 6 -1.250000e-04 4 5 -6.250000e-05 5 6 0.000000e+00 5 7 -3.750000e-04 2 3 1.875000e-04 2 5 1.250000e-04 2 4 2.500000e-04 2 7 -1.875000e-04 2 6 1.250000e-04 3 4 6.250000e-05 3 5 6.250000e-05 3 6 1.875000e-04 3 7 -3.125000e-04 0 3 3.750000e-04 0 2 1.250000e-04 0 1 0.000000e+00 0 7 -6.250000e-05 0 6 1.875000e-04 0 5 3.750000e-04 0 4 4.375000e-04 1 2 0.000000e+00 1 3 3.750000e-04 1 6 3.750000e-04 1 7 -6.250000e-05 1 4 5.000000e-04 1 5 3.750000e-04 6 7 -1.875000e-04 
            out_fid.write('%.2f  '% t)

            for ich in xrange(Nch):
		if ich==refX:	# skip
		    continue
                t = pos2tdoa_circular(theta, phi, geometry, refX, ich)
                out_fid.write(' %d %d %e'%(refX, ich, t))
            out_fid.write('\n')

    fid.close()
    out_fid.close()




def main():
    tdoa_fname = sys.argv[1]
    micgeo_fname = sys.argv[2]
    outdir = sys.argv[3]

    if outdir=='/'.join(tdoa_fname.split('/')[:-1]):
	print '[ERROR] Input and output TDOA files seems to be in the same folder! Abort now to avoid over-writing.'
	exit()

    subprocess.call('mkdir -p '+outdir, shell=True)
    arraytype, geometry = load_geometry(micgeo_fname)
    
    out_fname = tdoa_fname+'.tdoa'
    fid = open(tdoa_fname, 'r')
    out_fid = open(out_fname, 'w')
    for line in fid:
	line = line.strip()
	if line=='' or line[0]=='#':
	    continue
    	if arraytype==1:    # circular array
	    infile = line
	    outfile = outdir+'/'+'.'.join(infile.split('/')[-1].split('.')[:-1])+'.tdoa'
	    print infile, ' -> ', outfile
	    txt2tdoa_circular(infile, outfile, geometry)
	else:
	    print '[ERROR] Unsupported microphone array type. This script only supports circular array (1) currently.'
	    exit()

    fid.close()
    out_fid.close()


if __name__=='__main__':
    main()


