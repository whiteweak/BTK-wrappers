#!/usr/bin/env python


'''
Function: Delay and sum beamforming. This version is adapted from the
	version by John McDonough and Kenichi Kumatani for the convenience
	of large scale experiments.

 	This version assumes that you have smoothed the postion file from 
	KFST with external segmenation information already. Thus it only 
	have two running modes:

	  - Smooth the postion for the whole audio before doing beamforming
 	  - Using the raw position information and re-calculate the beamforming
	    parameters each time the postion changes

	It has removed the support for the following mode which was found
	buggy in previous version:

	  - Read external segmentation file and smooth position within any
	    utterance before doing beamforming.

Example folder:

	/share/spandh.ami1/usr/yulan/dev/btk/tools/DSBF/dbg/

Example command:

	/share/spandh.ami1/sw/std/python/v2.7.6/x86_64/bin/python2.7 -u ../DSBeamforming.py  --cfg=dbg.cfg

Yulan Liu, 6 Feb 2016
'''


import sys, os, os.path, pickle, numpy, subprocess
import string, time, wave

import optparse
from optparse import OptionParser

from numpy import *
from copy import *

from btk.common import *
from btk.stream import *
from btk.feature import *
from btk.utils import *

from btk import dbase
from btk.modulated import *
from btk.subbandBeamforming import *
from btk.beamformer import *
from btk.postfilter import *

sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../')
from utils import pos2tdoa_linear
from utils import pos2tdoa_circular
from utils import readPositionFile
from utils import readTDOAFile
from utils import load_geometry


class DSBeamformer():
    '''
    Delay and sum beamforming.
    '''

    def __init__(self, options):
	'''
	Prepare beamforming environments so that we could easily loop over 
	each wav files later.
	'''

    	self.postfilterType = int(options.postfilterType)
    	self.protoPath      = options.protoPath
    	self.M		    = int(options.Nsubband)
    	self.m 		    = int(options.FL)
    	self.r		    = int(options.ED)
    	self.alpha	    = float(options.alpha)
    	self.chanN 	    = int(options.chanN)
	self.chanNwav	    = int(options.chanNwav)
	self.Nbit	    = int(options.Nbit)
	self.mode	    = int(options.mode)

    	self.arraytype, self.geometries = load_geometry(options.geometryFile)
#	print '[DEBUG] arraytype:', self.arraytype, '\ngeometries:', self.geometries

    	assert(self.chanN<=len(self.geometries))
    	if(self.chanN<len(self.geometries)):
	    print '[WARNING] %d microphones are involved in beamforming while the geometry information indicates there are %d microphones. Using the geometries of the first %d microphones only.' % (self.chanN, len(self.geometries), self.chanN)


    	self.D              = int(self.M / 2 ** self.r) 	# frame shift
    	self.fftLen  	    = int(self.M)
    	self.NC 	    = 1
    	self.sampleRate     = int(options.sampleRate)
    	self.outSampleRate  = self.sampleRate		# Output wav file use the same sampling rate as the input wav file
    	self.timeIncrement  = (1.0 * self.D) / float(self.sampleRate) # *must* be float division    


    	# To distinguish different configurations, save the output wav files into subfolder named with configuration details.
    	if self.postfilterType > 0 :
            self.outputID = 'doa.M-%d_m-%d_r-%d_a-%0.2f' %(self.M, self.m, self.r, self.alpha)
    	else:
            self.outputID = 'doa.M-%d_m-%d_r-%d.woPF'    %(self.M, self.m, self.r)
    	self.outdir = options.OutDir+'/'+self.outputID
    	subprocess.call('mkdir -p '+self.outdir, shell=True)


    	# Load filter bank prototype
    	# Read analysis prototype 'h'
    	protoFile = '%s/h-M=%d-m=%d-r=%d.txt' %(self.protoPath, self.M, self.m, self.r)
    	print 'Loading analysis prototype from \'%s\'' %protoFile
    	fp = open(protoFile, 'r')
    	self.h_fb = pickle.load(fp)
    	fp.close()


    	# Read synthesis prototype 'g'
    	protoFile = '%s/g-M=%d-m=%d-r=%d.txt' %(self.protoPath, self.M, self.m, self.r)
    	print 'Loading synthesis prototype from \'%s\'' %protoFile
    	fp = open(protoFile, 'r')
    	self.g_fb = pickle.load(fp)
    	fp.close()


    	# Build the analysis chain
    	self.sampleFeats = []
    	self.analysisFBs = []
    	for chX in range(self.chanN):
             #sampleFeature = IterativeSampleFeaturePtr(blockLen = D, chanX = chX)
             #sampleFeature = IterativeSingleChannelSampleFeaturePtr( blockLen = D )
             sampleFeature = SampleFeaturePtr( blockLen=self.D, shiftLen=self.D, padZeros=False )
             self.sampleFeats.append(sampleFeature)
             analysisFB = OverSampledDFTAnalysisBankPtr(sampleFeature, prototype=self.h_fb, M=self.M, m=self.m, r=self.r, delayCompensationType=2 )
             self.analysisFBs.append(analysisFB)

    	# Initialize beamformer
    	self.pBeamformer = SubbandGSCPtr( fftLen=self.fftLen, halfBandShift=False)

    	# Add the post-filter
    	if self.postfilterType > 0 :
            self.output = ZelinskiPostFilterPtr( self.pBeamformer, self.fftLen, self.alpha, self.postfilterType )
    	else:
            self.output = self.pBeamformer

    	# Set the subband analysis filter bank to the beamformer as the input
    	for chX in range(self.chanN):
            self.pBeamformer.setChannel(self.analysisFBs[chX])

   	# Set the beamformer to the synthesizer to transform subbands into the time domain
    	self.synthesisFB = OverSampledDFTSynthesisBankPtr(self.output, prototype=self.g_fb, M=self.M, m=self.m, r=self.r, delayCompensationType=2, gainFactor=self.D )
			# delayCompensationType:	1 - compensate delays in the synthesis filter bank
			#				2 - compensate delays in the analythesis and synthesis filter banks



    def run_beamforming(self, inwav, posFile=None, TDOAFile=None, mode=None, outaudio_prefix=None, overlapSample=1024):
	'''
	Do actual beamforming over each wav file.

	At the margin of each process block, you can set an amount of overlapping samples 
	so that they will be smoothed. By default this is 1024 samples.
	'''
	if mode is None:
	    mode = self.mode

        # Output wav file full path
        self.prefix = inwav.split('/')[-1]
	if outaudio_prefix is None:
            self.outwav = self.outdir+'/'+ self.prefix + '.wav'
	else:
	    self.outwav = self.outdir+'/'+ outaudio_prefix.split('/')[-1] + '.wav'
	
	if posFile is None:
	    assert(TDOAFile is not None)
	    self.delaysL, Time = readTDOAFile( TDOAFile, mode=mode )
	else:
	    if TDOAFile is not None:
		print '[WARNING] Both postion file (', posFile, ') and TDOA file (', TDOAFile, ') are given, read position file by default.'
	    self.delaysL, Time, Pos = readPositionFile( posFile, self.geometries, self.arraytype, mode=mode)


	# Prepare to write audio samples processed with the beamformer
        wavefp = wave.open( self.outwav, 'w' )
        wavefp.setnchannels(1)
        wavefp.setsampwidth(2)
        wavefp.setframerate(int(self.outSampleRate))


	# Loop over the chunk
	startSample = 0
	endSample = 0
	overlapBuff = []
	writeBuff = []
	for i, tmp in enumerate(Time):
	    if i>0:
		st = et		# Using the ending time of previous iteration as the starting time
	    else:
		st = tmp[0]
	    et = tmp[1]

	    delaysL = - numpy.array(self.delaysL[i])	# NOTE!!! The minus sign here is very important!
#	    tmp_pos = numpy.array(Pos[i])	# [DEBUG]
	    self.pBeamformer.reset()
	    self.pBeamformer.calcGSCWeights( self.sampleRate, delaysL )	
	    
	    for fbinX in range(self.fftLen/2+1):
                self.pBeamformer.setActiveWeights_f( fbinX, zeros( 2*(self.chanN-self.NC), numpy.float ) )

	    # Read audio files
	    startT       = st
	    startSample  = endSample		# int(startT * self.sampleRate)
	    endT 	 = et
	    if endT==-1:
		endSample = -1
	    else:
            	endSample    = int(endT * self.sampleRate)
		# IMPORTANT: I have to do this otherwise the output sample size might not match the input sample size, resulting in a big alignment problem for large audio files
		durSample = int((endSample-startSample)/64)*64	# This is because the systhesis function works on such chunks
		endSample = startSample + durSample

	    if i>0:
		startSample = startSample-overlapSample

            print '\n[DEBUG] i / len(Time): ', i, '/', len(Time)
            print '[DEBUG] st -- et: ', st, ' -- ', et
#	    print '[DEBUG] average position (tmp_pos): \n', tmp_pos
	    print '[DEBUG] average delay (delaysL): ', delaysL


	    if i==(len(Time)-1):
	  	endSample = -1			# Process all rest data

	    for chX in range(self.chanN):
                if self.chanNwav==1:
                    sampleFile = inwav+'_'+str(chX+1).zfill(2)+'.wav'
                    if not os.path.exists(sampleFile):
                        print '[WARNING] Can not find file: %s skipping it now.' % sampleFile
			sampleFile = inwav+'-'+str(chX+1).zfill(2)+'.wav'
			if not os.path.exists(sampleFile):
			     print '[WARNING] Can not find file: %s either skipping it now.' % sampleFile
	                     continue
                    print 'Reading file %s as Ch %d, starting time: %f, ending time: %f; starting sample: %d, ending sample: %d ' %(sampleFile, chX, startT, endT, startSample, endSample)
                    self.sampleFeats[chX].read(sampleFile, cfrom=startSample, to=endSample, chX=1, chN=1, samplerate=self.sampleRate)
                elif self.chanNwav>1:
                    sampleFile = inwav+'.wav'
                    if not os.path.exists(sampleFile):
                        print 'Can not find file: %s skipping it now.' % sampleFile
                        return -1
                    print 'Reading file %s: Ch %d, starting time: %f, ending time: %f; starting sample: %d, ending sample: %d ' %(sampleFile, chX, startT, endT, startSample, endSample)
                    self.sampleFeats[chX].read(sampleFile, cfrom = startSample, to = endSample, chX = (chX + 1), chN = self.chanNwav, samplerate = self.sampleRate)
                self.analysisFBs[chX].reset()

            if self.postfilterType > 0 :
                self.output.setBeamformer( self.pBeamformer )


	    # Buffer the raw output
	    writeBuff = []
	    currentTime = 0.0
	    tmp_L = 0
            for b in self.synthesisFB:
		writeBuff+=list(b)
                currentTime += self.timeIncrement
		tmp_L += len(b)

	    if not endSample==-1:	# Just a lazy check whether the input output length match or not
		if not (len(writeBuff)==endSample-startSample):
		    print '[ERROR] len(writeBuff) != endSample-startSample '
		    print '[DEBUG] len(writeBuff) = ', len(writeBuff)
		    print '[DEBUG] endSample-startSample = ', endSample-startSample
	    	assert(len(writeBuff)==endSample-startSample)


	    # Smooth overlapping part
	    if i==0:
		samples2write = numpy.array(writeBuff[:-overlapSample])
		wavefp.writeframes( (samples2write/(2**(self.Nbit-16))).astype(numpy.short).tostring() )
		overlapBuff = writeBuff[-overlapSample:]
		continue

	    if i>0:
		for l in range(overlapSample):
		    writeBuff[l] = (overlapBuff[l]+writeBuff[l])/2
	    if i<len(Time)-1:
		samples2write = numpy.array(writeBuff[:-overlapSample])

                wavefp.writeframes( (samples2write/(2**(self.Nbit-16))).astype(numpy.short).tostring() )
		overlapBuff = writeBuff[-overlapSample:]
	    else:
		wavefp.writeframes( (numpy.array(writeBuff)/(2**(self.Nbit-16))).astype(numpy.short).tostring() )

		
#                # 24bit -> 16bit
#                wavefp.writeframes( (b/(2**(self.Nbit-16))).astype(numpy.short).tostring() )
#                # 24bit per sample -> 32 bit
#                #soundFileWriter.writeInt( numpy.array( b, numpy.int ) * 256.0 )
#                currentTime += self.timeIncrement
#	    print '[DEBUG] type(b)', type(b)
#	    print '[DEBUG] b:', b

        if self.mode==1:        # Smoothed over the whole file
            samples2write = numpy.array(writeBuff[-overlapSample:])
            wavefp.writeframes( (samples2write/(2**(self.Nbit-16))).astype(numpy.short).tostring() )

        wavefp.close()




def parse_option(cmd):
    '''
    Parse the option from command directly, or from configuration file.
    It also does some basical paramter check and throws errors if some
    very fundamental options are missing or in wrong format.
    '''

    parser = OptionParser(usage='DSBeamforming.py  --cfg=<configuration file>', version='0.1')

    # Parameters via configuration file
    parser.add_option('--cfg', dest='cfg', help='Import configuration file.', default=None)

    parser.add_option('--utteranceList', dest='utteranceList', help='List of full path prefix (without _01.wav or .wav) of utterances to be processed in the first column, and correponding position file prefix in the second column (e.g. output .txt files of Kalman speaker tracking), and output audio file prefix in the third column.', default='./uttList.txt')
    parser.add_option('--protoPath', dest='protoPath', help='The full path of the folder storing propertype configurations.', default='../DesignPrototype/prototypes/Nyquist')
    parser.add_option('--postfilterType', dest='postfilterType', help='If integer value >0: Zelinski with abs() real operator; else: use beamformer output.', default=0)
    parser.add_option('--geometryFile', dest='geometryFile', help='Microphone array geometry information: array type and coordinates of each microphone. Usually this file is saved automatically when running speaker tracking.', default=None)
    parser.add_option('--alpha', dest='alpha', help='Smoothing factor when calculating short-term cross spectral density (CSD).', default=0.6)
    parser.add_option('--Nbit', dest='Nbit', help='How many bit per sample per channel in wav file.', default=16)
    parser.add_option('--mode', dest='mode', help='0: using raw position without smoothing; 1: smooth position over whole auido.', default=0)

    parser.add_option('--Nsubband', dest='Nsubband', help='Filter bank parameters: number of subbands.', default=512)
    parser.add_option('--FL', dest='FL', help='Filter bank parameters: filter prototype length factor.', default=2)
    parser.add_option('--ED', dest='ED', help='Filter bank parameters: exponential decimation factor.', default=3)

    parser.add_option('--OutDir', dest='OutDir', help='Full path of the folder where beamforming output wav files will be saved.', default='./out')

    parser.add_option('--sampleRate', dest='sampleRate', help='Signal sample rate (Hz).', default=16000)

    parser.add_option('--chanN', dest='chanN', help='Number of channels for speaker tracking / beamforming.', default=4)
    parser.add_option('--chanNwav', dest='chanNwav', help='Number of channels in each wav file.', default=1)

    (options, args) = parser.parse_args(cmd)

    # !!! NOTE: The values in configuration file always overwrites those in the command line in the end.

    # Load the parameters from configuration file
    if options.cfg is not None:
        fid = open(options.cfg, 'r')
        for line in fid:
            line = line.strip()
            if line=='' or line[0]=='#':
                continue
            if '#' in line:
                line = line.split('#')[0].strip()

            data = line.split('=')
            cmd_tmp = 'options.'+data[0].strip()+' = '+data[1].strip()

            try:
                exec cmd_tmp
            except:
                print '[WARNING] Unrecognized configuration option: ', data[0]
                print '          Ignoring line in configuration file (', options.cfg, '): ', line
                continue
        fid.close()

    return options



def main():
    if len(sys.argv)<2:
        print 'Please try option "--help" to get help information.'
        exit()

    options = parse_option(sys.argv)

    print '\nConfiguration used:\n', options

    sys.stdout.write( '\nBegan beamforming on workstation  process %d on %s\n\n' %(os.getpid(), time.asctime()) )

    beamformer = DSBeamformer(options)

    fid = open(options.utteranceList, 'r')
    for line in fid:
	line = line.strip()
	if line=='' or line[0]=='#':
	    continue
	data = line.split()
	wav = data[0]
	outaudio_prefix = data[2]
	tmpFile = data[1]
	if '.tdoa' in tmpFile:	# Assume it is TDOA file
	    TDOAFile = tmpFile
	    posFile = None
	else:			# Assume it is position file
	    posFile = tmpFile
	    TDOAFile = None

	beamformer.run_beamforming(wav, posFile=posFile, TDOAFile=TDOAFile, mode=int(options.mode), outaudio_prefix=outaudio_prefix)
    fid.close()
    sys.stdout.write( '\nFinished beamforming on workstation process %d on %s\n\n' %(os.getpid(), time.asctime()) )



if __name__ == '__main__':
    main()



