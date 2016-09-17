#!/usr/bin/env python


'''
Function: Smooth the position file from KFST in one of the following modes:

	1. every given X seconds
	2. every given utterance / segment (segmentation provided externally)

Example folder:
	/share/spandh.ami1/usr/yulan/dev/btk/tools/KFST/dbg/cir_test/bigwav/

Example command:
	smoothPos.py  --scp=smooth_vad.scp  --mode=0

	The scp file is expected to have 3 columns when mode is 0: 

	<input .txt file>  <output .txt file>  <corresponding segmatation>

	The scp file is expected to have at least 2 columns when mode is a
	positive number X:

	<input .txt file>  <output .txt file>

Yulan Liu, 17 Mar 2016
'''

import sys, subprocess, optparse, math
from optparse import OptionParser
import numpy


def unwrap(inlist, L):
    '''
    Unwrap a sequence of angles.
    '''
    buff = []
    pre = ''
    for i, a in enumerate(inlist):
	if pre=='':
	    buff.append(a)
	    pre = a
	    continue
	else:
	    if abs(a-pre)<=max(abs(a+L-pre), abs(a-L-pre)):
		pass
	    elif abs(a+L-pre)>abs(a-pre):
		a += L
	    elif abs(a-L-pre)>abs(a-pre): 
   		a -= L
	    buff.append(a)
	    pre = a
	    continue
    return buff


def load_pos(fname):
   '''
   Load the information of one position file.
   '''
## Linear array
# Time   1.28 : 1.214301e+00
# Time   1.54 : 1.243904e+00
# Time   1.79 : 1.253513e+00
# Time   2.05 : 1.258232e+00


## Circular array
# Time   0.00 : (theta, phi) = (1.277241e+00, -3.088651e+00)
# Time   0.26 : (theta, phi) = (1.194526e+00, -2.827238e+00)
# Time   0.51 : (theta, phi) = (2.417945e-01, -2.376973e+00)
# Time   0.77 : (theta, phi) = (5.149969e-01, -2.339190e+00)
# Time   1.02 : (theta, phi) = (5.203482e-01, -2.280806e+00)
# Time   1.28 : (theta, phi) = (3.667310e-01, -2.195042e+00)
# Time   1.54 : (theta, phi) = (3.238956e-01, -2.191035e+00)

   fid = open(fname, 'r')
   info = []
   for line in fid:
	line = line.strip()
	if line=='' or line[0]=='#':
	    continue
	T = line.split(':')[0].strip().split()[-1]
	if len(line.split(':')[1].split('='))==2:	# Probably circular array
	    K = line.split(':')[1].split('=')[0].strip()	# (theta, phi)
	else:	# Probably linear array
	    K = ''
	if len(line.split('=')[1].split(','))>1:	# Probably circular array
	    pos = [x.strip() for x in line.split('=')[1].split('(')[-1].split(')')[0].split(',')]
	else:	# Probably linear array
	    pos = [line.split(':')[-1].strip()]
	info.append([T, pos])
   fid.close()
   return info, K


def smooth_vad(infname, outfname, vadfname):
    '''
    Smooth one position file according to given segmenation information.
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
     
    # Load position file
    info, K = load_pos(infname)

    # Smoothing and saving results
    buff_idx = 0
    pos_buff = []
    outdir = '/'.join(outfname.split('/')[:-1])
    subprocess.call('mkdir -p '+outdir, shell=True)
    outfid = open(outfname, 'w')
    CopyFlag = 0
    for i, tmp in enumerate(info):
	t = float(tmp[0])
	pos = tmp[1]
	try:
	    if t<buff[buff_idx][0]:		# One utterance hasn't started yet, copy the codes.
	    	CopyFlag = 1
	    elif t<buff[buff_idx][1]:	# In the middle of one utterance
	    	pos_buff.append([float(x) for x in pos])
	    	continue
	    elif t>=buff[buff_idx][1]:
	    	if t==buff[buff_idx][1] or len(pos_buff)==0:
	    	    pos_buff.append([float(x) for x in pos])
            	tmp = numpy.transpose(pos_buff)
            	tmp[0] = unwrap(tmp[0], math.pi)
            	tmp[1] = unwrap(tmp[1], math.pi*2)
            	pos_buff = numpy.transpose(tmp)
	    	aver_pos = numpy.array(pos_buff).mean(axis=0)
	    	if K=='':	# Probably linear array
		    outfid.write('Time   '+buff_str[buff_idx][0]+' : '+str(aver_pos[0])+'\n')
	    	else:	# Probably circular array
		    outfid.write('Time   '+buff_str[buff_idx][0]+' : '+K+' = ('+', '.join([str(x) for x in aver_pos])+')\n')

	    	pos_buff = []
	    	if t==buff[buff_idx][1]:
		    buff_idx += 1
	    	    continue
	    	else:
		    buff_idx += 1
		    try:
		    	buff[buff_idx]
		    	pos_buff.append([float(x) for x in pos])
		    	continue
		    except IndexError:
		    	CopyFlag = 1
        except IndexError:
            # Already finished all segments, going to copy the raw position as it is.
            if len(pos_buff)>0:
            	tmp = numpy.transpose(pos_buff)
            	tmp[0] = unwrap(tmp[0], math.pi)
            	tmp[1] = unwrap(tmp[1], math.pi*2)
            	pos_buff = numpy.transpose(tmp)
                aver_pos = numpy.array(pos_buff).mean(axis=0)
                if K=='':   # Probably linear array
                    outfid.write('Time   '+buff_str[buff_idx][0]+' : '+str(aver_pos[0])+'\n')
                else:       # Probably circular array
                    outfid.write('Time   '+buff_str[buff_idx][0]+' : '+K+' = ('+', '.join([str(x) for x in aver_pos])+')\n')
                pos_buff = []
            CopyFlag = 1

        if CopyFlag==1:
            if K=='':   # Probably linear array
                # Time   1.28 : 1.214301e+00
                outfid.write('Time   '+tmp[0]+' : '+pos[0]+'\n')
            else:       # Probably circular array
                # Time   0.00 : (theta, phi) = (1.277241e+00, -3.088651e+00)
                outfid.write('Time   '+tmp[0]+' : '+K+' = ('+', '.join(pos)+')\n')
            continue

    if len(pos_buff)>0:
        tmp = numpy.transpose(pos_buff)
        tmp[0] = unwrap(tmp[0], math.pi)
        tmp[1] = unwrap(tmp[1], math.pi*2)
        pos_buff = numpy.transpose(tmp)
        aver_pos = numpy.array(pos_buff).mean(axis=0)
        if K=='':   # Probably linear array
            outfid.write('Time   '+buff_str[buff_idx][0]+' : '+str(aver_pos[0])+'\n')
        else:       # Probably circular array
            outfid.write('Time   '+buff_str[buff_idx][0]+' : '+K+' = ('+', '.join([str(x) for x in aver_pos])+')\n')
    pos_buff = []
    outfid.close()




def smooth_fixed(infname, outfname, win):
   '''
   Smooth one position file with a given local window size.
   '''
   print infname, ' -> ', outfname
   # Load position file
   info, K = load_pos(infname)
   pos_buff = []
   outdir = '/'.join(outfname.split('/')[:-1])
   subprocess.call('mkdir -p '+outdir, shell=True)
   outfid = open(outfname, 'w')
   i = 1
   for tmp in info:
	t = float(tmp[0])
	pos = tmp[1]
	if t>=i*win and len(pos_buff)>0:
#	    # DEBUG
#	    if i==1:
#		print '[DEBUG] pos_buff:', pos_buff
	    tmp = numpy.transpose(pos_buff)
	    tmp[0] = unwrap(tmp[0], math.pi)
	    tmp[1] = unwrap(tmp[1], math.pi*2)
	    pos_buff = numpy.transpose(tmp)
	    aver_pos = numpy.array(pos_buff).mean(axis=0)
            if K=='':   # Probably linear array
                outfid.write('Time   '+'{:.2f}'.format((i-1)*win)+' : '+str(aver_pos[0])+'\n')
            else:       # Probably circular array
                outfid.write('Time   '+'{:.2f}'.format((i-1)*win)+' : '+K+' = ('+', '.join([str(x) for x in aver_pos])+')\n')
	    pos_buff = []
	    pos_buff.append([float(x) for x in pos])
	    i += 1
	    continue
	else:
	    pos_buff.append([float(x) for x in pos])

   if len(pos_buff)>0:
        tmp = numpy.transpose(pos_buff)
        tmp[0] = unwrap(tmp[0], math.pi)
        tmp[1] = unwrap(tmp[1], math.pi*2)
        pos_buff = numpy.transpose(tmp)
        aver_pos = numpy.array(pos_buff).mean(axis=0)
        if K=='':   # Probably linear array
            outfid.write('Time   '+'{:.2f}'.format(i*win)+' : '+str(aver_pos[0])+'\n')
        else:       # Probably circular array
            outfid.write('Time   '+'{:.2f}'.format(i*win)+' : '+K+' = ('+', '.join([str(x) for x in aver_pos])+')\n')
   pos_buff = []
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
    parser = OptionParser(usage='smoothPos.py  --scp=<scp file>  --mode=X', version='0.1')
    parser.add_option('--scp', dest='scp', help='List of input position files, output position files and segmentation files (when in mode 0).', default=None)
    parser.add_option('--mode', dest='mode', help='Smoothing mode: 0 - smoothing by external segmentation; ANY POSITIVE VALUE X - smoothing every X seconds..', default=0)

    (options, args) = parser.parse_args(cmd)
    assert(options.scp is not None)
    return options


def main():
    usage = 'smoothPos.py  --scp=<scp>  --mode=X\n'
    if len(sys.argv)<3 and not (sys.argv[1] in ['-h', '--help']):
	print usage
        print 'Please try option "--help" to get help information.'
	exit()

    options = parse_option(sys.argv)
    smooth_all(options)


if __name__ == '__main__':
    main()


