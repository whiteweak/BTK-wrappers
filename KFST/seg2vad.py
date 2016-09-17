#!/usr/bin/env python

'''
Function: Convert the segment lines in mlf file into vad format that is accepted by
	other BTK tool scripts written by Yulan.

	The script will perform segment boundary validation by comparing the starting
	time and ending time of each segment with processing block size in KFST, thus
	change the segment boundary slightly so that there won't be too little size
	of data to process, as that will cause some unexpected error.

Example folder:
	/share/spandh.ami1/asr/dev/mtg/swc/exp/mdm/KFST/TBL1/manseg/

Example command:
	seg2vad.py  swc1.lst  vad_manseg

Yulan Liu, 7 Apr 2016
'''

import sys, subprocess

usage = 'seg2vad.py  <input segment list>  <output folder>'

if len(sys.argv)<3:
    print usage
    exit()

print '\nExecute:\n\t'+' '.join(sys.argv)+'\n'


def read_segs(infile):
    '''
    Read the segments and save the information.
    '''
    fid = open(infile, 'r')
    infodict = dict()
    for line in fid:
	line = line.strip()
	if line=='' or line[0] in ['#', '!']:
	    continue
	# Correct segment lines should look like
	# "*/SWC1-00001_mn0001_000110_000180.lab"
	# <ses> _ <spkr> _ <st in 10ms> _ <et in 10ms> 
	if line[0]=='"':
	    seg = line.split('/')[-1].split('.')[0].split('"')[-1]
	    ses, spk, st, et = seg.split('_')
	    st = int(st)/100.0
	    et = int(et)/100.0
	if st>et:
	    print '[WARNING] Illegal segment: starting time (',st,') larger than ending time (',et,')'
	    print '          Skipping line: ', line
	    continue
	try:
	    infodict[ses]
	except KeyError:
	    infodict[ses] = dict()
	try:
	    infodict[ses][spk]
	except KeyError:
	    infodict[ses][spk] = dict()
	try:
	    infodict[ses][spk][st]
	    print '[WARNING] Duplicated record for speaker "'+spk+'" in session "'+ses+'", time starting at:', st, ' (second)'
	    print '          Skipping line: ', line
	except KeyError:
	    infodict[ses][spk][st] = et

    fid.close()
    return infodict



def write_vad(outdir, infodict):
    '''
    Write information into VAD format that are accepted by other BTK tool scripts
    written by Yulan.
    '''
    subprocess.call('mkdir -p '+outdir, shell=True)   
    for ses in sorted(infodict.keys()):
	for spk in sorted(infodict[ses].keys()):
	    outfile = outdir+'/'+ses+'_'+spk
	    fid = open(outfile, 'w')
	    fid.write('# '+ses+'_'+spk+'\n')
	    for st in sorted(infodict[ses][spk].keys()):
		et = infodict[ses][spk][st]
		fid.write('{:.2f}'.format(st)+' '+'{:.2f}'.format(et)+'\n')
	    fid.close()


def main():
    seg_infile = sys.argv[1]
    outdir = sys.argv[2]
    infodict = read_segs(seg_infile)
    write_vad(outdir, infodict)


if __name__ == '__main__':
    main()

