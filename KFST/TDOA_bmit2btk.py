#!/usr/bin/env python

'''
Function: Convert the TDOA from BeamformIt to the format supported by BTK.
	It only takes one argument: the file listing full paths of input 
	TDOA file and output TDOA file in two columns.


Yulan Liu, 18 May 2016
'''

import sys, subprocess

usage = 'TDOA_bmit2btk.py   <scp>   <sampling frequency>'
if len(sys.argv)<3:
    print usage
    exit()

def convert_tdoa(infile, outfile, fs=16000):
    '''
    Convert TDOA file from BeamformIt format to BTK format.
    '''
    print infile, ' -> ', outfile
    infid = open(infile, 'r')
    outfid = open(outfile, 'w')
    for line in infid:
	line = line.strip()
	if line=='' or line[0]=='#':
	    continue
# 0 -> 2 0.208784  0 1.000000  -2 0.177302  0 1.000000  0 1.000000  0 1.000000  0 1.000000  0 1.000000 
# 250 -> -2 0.156027  0 1.000000  2 0.178206  -4 0.109189  0 0.076296  0 0.075005  0 0.084063  2 0.115262 
# 500 -> 2 0.146953  0 1.000000  -2 0.177980  -4 0.107577  0 0.086221  0 0.090553  0 0.099456  2 0.111682 
# 750 -> 3 0.245355  0 1.000000  -2 0.390337  -2 0.385181  1 0.391067  4 0.234785  6 0.211617  6 0.272227 
# 1000 -> 3 0.162097  0 1.000000  -2 0.277204  -2 0.225740  1 0.242711  5 0.151510  6 0.121896  6 0.132545 
# 1250 -> 2 0.140391  0 1.000000  -2 0.275900  -2 0.174146  1 0.204948  5 0.139635  6 0.114704  6 0.117395 
# 1500 -> -2 0.186325  0 1.000000  2 0.177865  -2 0.092477  1 0.074652  5 0.071161  6 0.072501  6 0.091982 

	t = float(line.split()[0])/1000.0		# Convert to ms to s
	outfid.write('%.2f ' % t)
	tmp = line.split()[2:]
	tdoa_cnt_list = []
	for i in xrange(0, len(tmp), 2):
	    tdoa_cnt = int(tmp[i])
	    tdoa_cnt_list.append(tdoa_cnt)
	ref_bias = tdoa_cnt_list[0]			# Regularize to use channel 0 as reference channel
	for i, tdoa_cnt in enumerate(tdoa_cnt_list):
	    tdoa_cnt = tdoa_cnt - ref_bias
	    tdoa = - tdoa_cnt*1.0/fs			# NOTE!!! This minus sign is very important!
	    if i==0:
		continue
	    outfid.write('%d %d %e ' %(0, i, tdoa))
	outfid.write('\n')

    outfid.close()
    infid.close()


def main():
    scp = sys.argv[1]
    fs = int(sys.argv[2])
    fid = open(scp, 'r')
    for line in fid:
	line = line.strip()
	if line=='' or line[0]=='#':
	    continue
	data = line.split()
	assert(len(data)==2)
	infile = data[0]
	outfile = data[1]
	outdir = '/'.join(outfile.split('/')[:-1])
	subprocess.call('mkdir -p '+outdir, shell=True)
	convert_tdoa(infile, outfile, fs)
    fid.close()

if __name__=='__main__':
    main()


