#!/usr/bin/env python



'''
Function:       Speaker tracking.

Example folder:
        /share/spandh.ami1/asr/dev/mtg/wsj5k/exp.mdm/dev/KFST/

Example command:
        ./speakerTrackingEKF.py  --cfg=dbg/dbg.cfg

About the algorithms used, please cite the following paper:

        U. Klee, T. Gehrig and J. W. McDonough, "Kalman Filters for Time Delay of Arrival-Based Source Localization," EURASIP J. Adv. Sig., 2006. 

The following website also provides some information:

        http://mlsp.cs.cmu.edu/people/johnmcd/papers_7.html


Modified from John McDonough and Kenichi Kumatani's BTK experiment scripts by Yulan Liu to enable a better user interface.

22 Jan 2016
'''


import sys, os, string, copy, time
import glob, socket, shutil, itertools

from subprocess import check_output
import optparse
from optparse import OptionParser


import numpy, math
from numpy import *

from btk.utils import *
from btk.feature import *
from btk.tdoa import *
from btk.kalman import *


def unpack_pair(s, t=0):
    '''
    Unpack the string that looks like

	"['1 2', '1 3', '1 4', '2 3', '2 4', '3 4']"
	"['0 0 0', '2 0 0', '2 2 0', '0 2 0']"
	
    into python lists int (t=1), or float (t=2), or strings.
    '''

    tmp = eval('[' + s.split('[')[-1].split(']')[0] + ']')

    l = []
    for x in tmp:
	if t==1:
	    l.append([int(y) for y in x.split()])
	elif t==2:
	    l.append([float(y) for y in x.split()])
	else:
	    l.append(x.split())
    return l

def bound_angle(x, min_x, max_x):
   '''
   Wrap the angle "x" into a range of [min_x, max_x].
   '''
   if min_x<=x<=max_x:
	return x
   else:
	L = max_x-min_x
	if x>max_x:
	    while x>max_x:
		x -= L
	elif x<min_x:
	    while x<min_x:
		x += L
	return x


def getIndices(Pairs):
    '''
    Pair    		==>    Indices

    [[1, 2], [1, 3]]           {1:1, 2:1, 3:1}

    '''

    Indices = {}
    for pair in Pairs:
        firstX                      = int(pair[0])
        secondX                     = int(pair[1])
        Indices[firstX]   = 1
        Indices[secondX]  = 1

    indexList = []
    for index in Indices:
        indexList.append(index)

    return sorted(indexList)



def get_pairs(N, m):
    '''
    Generate all possible combination of microphone pairs, 
    from N channels and with each pair composed of m 
    microphones.
    '''
    pairs = "[ "
    a = range(1, N+1, 1)
    for tmp in itertools.combinations(a, m):
	pairs += "'"
	pairs += ' '.join([str(x) for x in tmp])
	pairs += "', "

    pairs += "]"
    return pairs



def parse_option(cmd):
    '''
    Parse the option from command directly, or from configuration file.
    It also does some basical paramter check and throws errors if some
    very fundamental options are missing or in wrong format.
    '''

    parser = OptionParser(usage='speakerTrackingEKF.py  --chanN=<#channel>  --sampleRate=<sampling rate>  --samplePath=<audio folder full path>\n\n\tOr\n\nspeakerTrackingEKF.py  --cfg=<configuration file>', version='0.1')


    # Parameters via configuration file
    parser.add_option('--cfg', dest='cfg', help='Import configuration file.', default=None)

    # Parameters via command
    parser.add_option('--utteranceList', dest='utteranceList', help='List of full path prefix (without _01.wav or -01.wav) of utterances to be processed as left column, and the correponding output file prefix as right column.', default='./uttList.txt')
#    parser.add_option('--delay', dest='delay', help='', default=0.0)
#    parser.add_option('-F', '--sampleRate', dest='sampleRate', help='Signal sample rate (Hz).', default=16000)
    parser.add_option('--sampleRate', dest='sampleRate', help='Signal sample rate (Hz).', default=16000)
    parser.add_option('--D', dest='D', help='Half of FFT length.', default=4096)
    parser.add_option('--threshold', dest='threshold', help='Threshold for GCC-PHAT.', default=0.1)
    parser.add_option('--chanN', dest='chanN', help='Number of channels for speaker tracking / beamforming.', default=4)
    parser.add_option('--minimumPairs', dest='minimumPairs', help='Minimal pair size (number of channels in a pair) when calculating TDOA.', default=2)
#    parser.add_option('--interpolate', dest='interpolate', help='Interpolate or not.', default=False)
    parser.add_option('--microphonePairs', dest='microphonePairs', help='Pairs of microphone that we are going to calculate TDOA.', default=[])
    parser.add_option('--diameter', dest='diameter', help='The diameter of circular array or sphere array (in mm).', default='')
    parser.add_option('--dis2mic', dest='dis2mic', help='The distance of two neighbouring microphones (equally distributed) in linear array (in mm).', default='')

    parser.add_option('--arraytype', dest='arraytype', help='Microphone array type, 0 for linear array, 1 for circular array, 2 for sphere array.', default=1)
    parser.add_option('--arrayCoordinates', dest='arrayCoordinates', help='The 3D coordinants for microphones in the array.', default=[])

    parser.add_option('--chanNwav', dest='chanNwav', help='Number of channels in each wav file.', default=1)

    parser.add_option('--batchMode', dest='batchMode', help='-1: Process each audio file in one batch; 0: Use given segmentation information (either from manual segmentation or from external VAD) via option "--vadList"; Any positive value: Re-initialize the Kalman Filter every X seconds, X is the value you give here.', default=-1)
    parser.add_option('--vadList', dest='vadList', help='The full path to the file listing the full path for VAD files. The order should match those you gave to option "--utteranceList", and each VAD file is composed of two columns, the first column indicates the starting time of one speech segment, and the second column the ending time, unit is in second.', default='./vadList.txt')


    (options, args) = parser.parse_args(cmd)


    print '!!! NOTE: The values in configuration file always overwrites those in the command line in the end.\n'

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

#    v = optparse.Values(cmd)

    # Convert string to list: microphone pairs and geometry
    if type(options.microphonePairs)==str:      # e.g. in bash command --microphonePairs="['1 2', '1 3', '1 4', '2 3', '2 4', '3 4']"
	options.microphonePairs = unpack_pair(options.microphonePairs, 1)
    elif (type(options.microphonePairs)==list) and (len(options.microphonePairs)==0):
        options.microphonePairs = unpack_pair(get_pairs(options.chanN, options.minimumPairs), 1)
	print 'Using automatically generated microphone pairs: ', options.microphonePairs
    else:
        print 'WARNING! Unrecognized type for paramter "--microphonePairs"! If you are not sure how to use it, please just use the default setup.'


    if type(options.arrayCoordinates)==str:	# e.g. in bash command --arrayCoordinates="['0 0 0', '2 0 0', '2 2 0', '0 2 0']"
	options.arrayCoordinates = unpack_pair(options.arrayCoordinates, 2)

    print '[LOG] Configuration used: ', options

    return options


class SingleSpeakerTracker:
    '''
    Single speaker tracking with Kalman filter for linear array or circular array.
    '''

    def __init__(self, options):
	self.sampleRate = options.sampleRate
	self.f_s = float(self.sampleRate)
	self.chanNwav = int(options.chanNwav)

	self.arraytype_dict = {0:'linear', 1:'circular', 2:'sphere'}


	# FFT
	self.D = int(options.D)
	self.deltaT = float(self.D) / self.f_s
	self.fftLen= self.D * 2
	
	# Default sound speed (mm/s)
	self.c = 343000.0

	# cross-correlation threshold for source detection
	self.threshold = options.threshold

	# detect a source if # mic. pairs with the CC threshold exceeds this number
        self.minimumPairs= int(options.minimumPairs)
	self.microphonePairs = options.microphonePairs


	# set the default array geometry
	self.channelsN = int(options.chanN)
	self.arraytype = int(options.arraytype)

	if self.arraytype == 0:		# Linear array use coordinates
	    if len(options.arrayCoordinates)==self.channelsN:
	    	self.geometries = numpy.array([[x[0]] for x in options.arrayCoordinates])
	    else:
		self.geometries = []
		self.dis2mic = float(options.dis2mic)
		self.configLinearArray()
	elif self.arraytype == 1:	# Circular array could use both coordinates (higher priority) or diameter
	    if len(options.arrayCoordinates)==self.channelsN:
		self.geometries = numpy.array(options.arrayCoordinates)
	    else:
            	self.geometries = []
	    	self.diameter = float(options.diameter)	# in mm
            	self.configCircularArray()
	elif self.arraytype == 2:		# Sphere array could use both coordinates (higher priority) or diameter
	    if len(options.arrayCoordinates)==self.channelsN:
		self.geometries = numpy.array(options.arrayCoordinates)
	    else:
		self.geometries = []
                self.diameter = float(options.diameter) # in mm

	else:
	    print '[ERROR] Unrecognized array type, try 0 (linear array), 1 (circular array), or 2 (sphere array) for "--arraytype"'
	    exit()


    def configSphereArray(self):
	'''
	Automatically generate the coordinates of each microphone given the 
	number of channels and the array diameter.
	'''

	if not (self.arraytype==2):
	    print '[ERROR] Trying to configure as a sphere array while the array type is: ', self.arraytype_dict[self.arraytype]
	    exit()

	print '[ERROR] Currently there is no support for sphere array yet, however it will come soon. Please contact the support for more details.'
        exit()


    def configCircularArray(self):
	'''
	Automatically generate the coordinates of each microphone given the 
        number of channels and the array diameter.
	'''

        if not (self.arraytype==1):
            print '[ERROR] Trying to configure as a circular array while the array type is: ', self.arraytype_dict[self.arraytype]
            exit()

	radiusArray = self.diameter / 2
	azimuthDelta            = - 2.0 * numpy.pi / self.channelsN
        arrgeom                 = []
        for n in range(self.channelsN):
            m_x                 = radiusArray * math.cos(n * azimuthDelta)
            m_y                 = radiusArray * math.sin(n * azimuthDelta)
            m_z                 = 0.0
            arrgeom.append([m_x, m_y, m_z])
	self.geometries = numpy.array(arrgeom)


    def configLinearArray(self):
	'''
	Automatically generate the coordinates of each microphone given the
	number of channels and the distance between each two neighbouring
	microphones.
	'''

        if not (self.arraytype==0):
            print '[ERROR] Trying to configure as a linear array while the array type is: ', self.arraytype_dict[self.arraytype]
            exit()

	n = int(math.floor(self.channelsN/2))	# Index of microphone (starting from 0) that is assumed to be [0, 0, 0]
	arrgeom = []
	for i in range(self.channelsN):
	   m_x = (i-n) * self.dis2mic
	   arrgeom.append([m_x])
	self.geometries = numpy.array(arrgeom)


    def makeFrontEnd(self, sources, microphonePositions):
	'''
	Pair the features (sources) from different channel according to the
	microphone pairs.
	'''

        microphonePairsSource   = []
        pairX = 0

        for pair in self.microphonePairs:
            firstX              = int(pair[0])
            secondX             = int(pair[1])
            firstSource         = sources[firstX - 1]
            secondSource        = sources[secondX - 1]
            phat                = PHATFeature(firstSource, secondSource, self.fftLen)
            tdoa                = TDOAFeature(phat, self.fftLen, self.f_s)
            newPair             = MicrophonePairSource(pairX, firstX-1, secondX-1, tdoa)
            microphonePairsSource.append(newPair)

            pairX += 1
	
	if self.arraytype==0:
            return FarfieldLinearArrayTDOAFeatureVector(microphonePairsSource, microphonePositions, self.minimumPairs, self.threshold, self.c)

	elif self.arraytype==1:
            return FarfieldCircularArrayTDOAFeatureVector(microphonePairsSource, microphonePositions, self.minimumPairs, self.threshold, self.c)

	elif self.arraytype==2:		# TODO
	    print '[ERROR] Not supporting sphere array for now.'
	    exit()
	    # return FarfieldSphereArrayTDOAFeatureVector(microphonePairsSource, microphonePositions, self.minimumPairs, self.threshold, self.c)

	else:
	    print '[ERROR] Unrecognized array type: ', self.arraytype
	    exit()

	return 1


    def buildExtKFSourceTracker(self, initialXk, boundaries):
	'''
	Extend Kalman Filter source tracker and build the processing chain for
	all channels. 
	'''

        self.sampleFeats        = []
        self.spectra            = []

        for chanX in range(self.channelsN):
            sampleFeature       = SampleFeaturePtr(blockLen = self.D, shiftLen = self.D, padZeros = True)
            self.sampleFeats.append(sampleFeature)
            hammingFeature      = HammingFeaturePtr(sampleFeature)
            fft         	= FFTFeaturePtr(hammingFeature, self.fftLen)
            self.spectra.append(fft)

        self.frontEnd   = self.makeFrontEnd(self.spectra, self.geometries )
        self.microphoneList     = getIndices(self.microphonePairs)


        # Miscellaneous parameters
	if self.arraytype==0:	# Linear array
	    sizeXk		= 1
            sigmaK2             = 10.0
            sigmaV2             = 4.0E-04
	elif self.arraytype==1:	# Circular array
	    sizeXk              = 2
            sigmaK2             = 0.1
            sigmaV2             = 3.35E-04

        sizeYk                  = len(self.microphonePairs)
        sigmaK2initial          = 1.0E+10
        gateProbability         = 0.95

	# Save for running time re-initialization
	self.sizeXk		= sizeXk
	self.sizeYk		= sizeYk
	self.sigmaK2initial	= sigmaK2initial
	self.sigmaK2		= sigmaK2
	self.sigmaV2		= sigmaV2
	self.gateProbability	= gateProbability

        # Transition matrix is the identity matrix
        F                       = numpy.identity(sizeXk, numpy.float)
	self.F			= F

        # Observation and process noise are diagonal with diagonal terms sigma2
        U                       = sigmaK2 * numpy.identity(sizeXk)
	self.U			= U

        self.kalmanFilter       = ExtendedKalmanFilter(self.frontEnd, F, U, sigmaV2, sigmaK2initial, self.deltaT,
                                                       initialXk = initialXk,
                                                       gateProbability = gateProbability, boundaries = boundaries)
        return self.kalmanFilter


    def updateExtKFSourceTracker(self, initialXk, boundaries):
	'''
	Update KF with new initialXk and boundaries, without initializing the front-end
	'''

	self.kalmanFilter       = ExtendedKalmanFilter(self.frontEnd, self.F, self.U, self.sigmaV2, self.sigmaK2initial, self.deltaT,
                                                       initialXk = initialXk,
                                                       gateProbability = self.gateProbability, boundaries = boundaries)
	return self.kalmanFilter


    def init_kalmanFilter(self, initialXk, boundaries):
	'''
	Initiate Kalman Filter at running time.
	'''
	self.kalmanFilter.F	= numpy.identity(self.sizeXk)
	self.kalmanFilter.H	= None
	self.kalmanFilter.U	= self.sigmaK2 * numpy.identity(self.sizeXk)
	self.kalmanFilter.sigmaV2 		= self.sigmaV2
        self.kalmanFilter.stateLength           = self.kalmanFilter.F.shape[0]
        self.kalmanFilter.I                     = numpy.identity(self.kalmanFilter.stateLength, numpy.float)
	self.kalmanFilter.deltaT		= self.deltaT
	self.kalmanFilter.gateProbability	= self.gateProbability
	self.kalmanFilter.boundaries		= boundaries
        self.kalmanFilter.observed              = False
	

        print 'Gate probability = %10.4f' % self.gateProbability

        if self.kalmanFilter.gateProbability == 0.0:
            self.kalmanFilter.innovationFilter       = False
        else:
            self.kalmanFilter.innovationFilter       = True
            self.kalmanFilter.chi                    = scipy.stats.distributions.chi_gen(a=0.0,name='chi',shapes='df')

        self.kalmanFilter.K_filter                   = self.sigmaK2initial * numpy.identity(self.kalmanFilter.stateLength, numpy.float)
        self.kalmanFilter.K_predict                  = self.sigmaK2initial * numpy.identity(self.kalmanFilter.stateLength, numpy.float)

        self.kalmanFilter.lastUpdateT                = -1		# ??
        self.kalmanFilter.time                       = -1		# ??


        if initialXk is None:
            self.kalmanFilter.xk_filter = zeros(self.kalmanFilter.stateLength, numpy.float)
        else:
            self.kalmanFilter.xk_filter = initialXk


    def loadSamples(self, samplePath, startSample=0, endSample=-1):
	'''
	If the input audio file is single channel, then the assumed file name in full path is:

		ch1: samplePath + '_01.wav' / '-01.wav'
		ch2: samplePath + '_02.wav' / '-02.wav'
		...
		ch8: samplePath + '_08.wav' / '-08.wav'
		...
		ch10: samplePath + '_10.wav' / '-10.wav'
		...

	If the input audio file is multi-channel, then the assumed file name is:

		ch1 - chN: samplePath + '.wav'
	'''

 	startT = startSample*1.0 / self.sampleRate
	if endSample==-1:
	    endT = -1
	else:
	    endT   = endSample*1.0 / self.sampleRate
	if self.chanNwav==1:		# .wav files of single channel, assuming the filename in "*-??.wav" or "*_??.wav" pattern
            for chanX in range(self.channelsN):
		sampleFile = samplePath+'_'+str(chanX+1).zfill(2)+'.wav'
            	if os.path.exists(sampleFile):
		    print 'Reading file %s as Channel %d' %(sampleFile, chanX + 1)
		else:
		    print '[WARNING] Unable to find sample file %s' %sampleFile
		    sampleFile = samplePath+'-'+str(chanX+1).zfill(2)+'.wav'
		    if os.path.exists(sampleFile):
			print 'Reading file %s as Channel %d' %(sampleFile, chanX + 1)
		    else:
		    	print '[ERROR] Unable to find sample file %s either.' %sampleFile
		    	exit()

                self.spectra[chanX].reset()
		self.sampleFeats[chanX].read(sampleFile, cfrom=startSample, to=endSample, chX=1, chN=1, samplerate=self.sampleRate)

	elif self.chanNwav>1:		# .wav files of multiple channels
	    sampleFile = samplePath + '.wav'
            if os.path.exists(sampleFile):
	    	for chanX in range(self.channelsN):
                    print 'Reading file %s from %3.2f to %3.2f : Ch %d' %(sampleFile, startT, endT, chanX)
		    self.spectra[chanX].reset()
                    self.sampleFeats[chanX].read(sampleFile, cfrom = startSample, to = endSample, chX = (chanX + 1), chN = self.chanNwav, samplerate = self.sampleRate)

	    else:
		print '[ERROR] Unable to find sample file %s' %sampleFile
                exit()
	else:				# Negative number of channels in .wav file
	    print '[ERROR] The number of channels in .wav file can only be positive interger, however it is now ', self.chanNwav
	    print '        Check the parameters for "--chanNwav" option.'
	    exit()

    def save_geometry(self, file):
	'''
	Save the geometry details.
	'''
	fid = open(file, 'w')
	fid.write('#Arraytype : %d\n' %self.arraytype)	# Array type
	print '\n[LOG] microphone geometries: ', self.geometries, '\n'
	for tmp in self.geometries:
	    try:
	    	fid.write(' '.join([str(x) for x in tmp])+'\n')
	    except TypeError:
	    	fid.write(str(tmp)+'\n')
	fid.close()


    def load_geometry(self, file):
   	'''
   	Load the geometry information from file.
   	'''
   	fid = open(file, 'r')
   	self.geometries = []
   	for line in fid:
            line = line.strip()
            if line=='':
            	continue
            if line[0]=='#':        # Arraytype information
            	self.arraytype = int(line.split(':')[-1].strip())
            else:
            	self.geometries.append([float(x) for x in line.split()])
   	fid.close()
   	return self.arraytype, self.geometries




def track(in_prefix, out_prefix, options, vadFile=''):
    '''
    Kalman Filter based sound source tracking for one utterance of microphone 
    array recordings.
 

	options.arraytype:	0	Linear array
				1	Circular array
				2	Sphere array

    '''
    arraytype = int(options.arraytype)
    sampleRate = int(options.sampleRate)

    # define the search space for sound sources
    if arraytype==0:	# Linear array
	originalXk = numpy.array([0], numpy.float)
	originalBoundaries = []
	originalBoundaries.append([-numpy.pi, numpy.pi])
    elif arraytype==1:	# Circular array
	originalXk = numpy.array([numpy.pi / 4.0, 0.0])
    	originalBoundaries = []
#    	originalBoundaries.append([0.0, numpy.pi / 2.0])	# This might not be true
	originalBoundaries.append([0.0, numpy.pi])
    	originalBoundaries.append([-numpy.pi, numpy.pi])
    elif arraytype==2:	# Sphere array
	print '[ERROR] Sphere array not supported currently.'
	exit()
    else:
	print '[ERROR] Unrecognized array type, try 0 for linear array, 1 for circular array and 2 for sphere array.'
	exit()

    initialXk       = originalXk
    boundaries      = originalBoundaries


    # Prepare filename and folder structure
    outputDir = '/'.join(out_prefix.split('/')[:-1])
    positionfile = out_prefix + '.txt'
    tdoafile = out_prefix + '.tdoa'

    if not os.path.exists(outputDir):
        try:
            os.makedirs(outputDir)
        except:
            print 'could not make %s' % outputDir
    print 'Output files: ', positionfile, tdoafile
    posfp           = open(positionfile, 'w')
    tdoafp          = open(tdoafile,     'w')

    sst = SingleSpeakerTracker(options)

    # Save microphone geometry for future use
    sst.save_geometry(outputDir+'/mic_geometry')


    # Check whether we want to re-initialize the KF regularly or not.
    batchMode = float(options.batchMode)
    TIME = []	
    flag = -1		# Default: process each file in one batch
    if (abs(batchMode) > 1e-5) and (batchMode>0):               # Re-initialize the KF at a constant rate
        print 'Re-initialize the KF at every ', batchMode, ' seconds.'
	flag = 1
    elif abs(batchMode) < 1e-5:                                 # Read external segmentation information
        print 'Reading external segmentation information from file: ', vadFile
	flag = 0
	fid = open(vadFile, 'r')
	for line in fid:
	    line = line.strip()
	    if line=='' or line[0]=='#':
		continue
	    tmp = line.split()
	    assert(len(tmp)==2)
	    assert(float(tmp[1])>=float(tmp[0]))
	    TIME.append([float(x) for x in tmp])
	fid.close()
    elif abs(batchMode-(-1)) < 1e-5:                            # Process each file in one batch (the default mode)
        print 'Process each audio file in one batch.'
	flag = -1

    # Control the update rate of KF tracking
    startT = 0.0
    deltaT = sst.deltaT		

    # Control the initialilzation rate of KF
    SegIdx = 0


    # Align the boundary between segment and processing block
    ChunkSize = int(options.D)		# Feature processing block size
    startSample = 0
    average = zeros(len(initialXk), numpy.float)

    kalmanFilter = sst.buildExtKFSourceTracker(initialXk, boundaries)
    if flag==1 or flag==-1:
	endSample = -1
    elif flag==0:
	tmp_et = TIME[SegIdx][0]
        endSample = int(tmp_et*sampleRate)
    sst.loadSamples(in_prefix, startSample, endSample)
    LoadFlag = 0
    InitFlag = 0
    frameX = 0
    wasSource = False


#   NOTE: you can also use following codes to compare raw location tracking and the one with KF.
#          This is one example for Linear array.
    if options.arraytype==0:    # Linear array, need to update the initialization
    	while True:
            initialXk = sst.frontEnd.getInstantaneousPosition( frameX )
            if initialXk > -1e10:
            	print 'Time %6.2f : theta = %10.2f' %((startT + frameX * deltaT), initialXk*180.0/pi)
            	break
            frameX += 1
    	print 'Update initialXk for linear array to: ', initialXk
	# We implicitly skipped some samples at the beginning which does not pass the initialXk>-1e10 test
	kalmanFilter = sst.updateExtKFSourceTracker(initialXk, boundaries)

    while True:
	# Check whether we need to re-initialize the KF again, do so if needed
	if flag==1:
	    if (startT+(frameX-1)*deltaT > (SegIdx+1)*batchMode):
		InitFlag = 1
		SegIdx += 1
	    else:
		InitFlag = 0
	elif flag==0:
	    try:
                if (not startSample==-1) and ((startT+(frameX-1)*deltaT - TIME[SegIdx][0])*sampleRate > -ChunkSize):
		    LoadFlag = 1
		    SegIdx += 1
		else:
		    InitFlag = 0
		    LoadFlag = 0
	    except IndexError:	# Reaching the last bit in this file
		pass
	else:
	    InitFlag = 0

	if InitFlag==1 and LoadFlag==0:		# Initialize the KF without re-loading features
            print '[LOG] Initializing the KF tracker without re-loading features.'
	    sst.init_kalmanFilter(initialXk, boundaries)
	    average = zeros(len(initialXk), numpy.float)
	    wasSource = False
	    InitFlag = 0


	# KF
	try:
            if options.arraytype==0:    # Linear array
		try:
                    azimuth = sst.kalmanFilter.next(frameX)
		    azimuth = bound_angle(azimuth, boundaries[0][0], boundaries[0][1])

		    # Added by Yulan Liu, 28 Jan 2016
		    if not sst.kalmanFilter.isSourceObserved():
		    	print 'Time %6.2f : Nothing detected' %(startT + frameX * deltaT)
		    	frameX += 1
		    	average += azimuth

		    else:
                    	posfp.write('Time %6.2f : %e\n' %((startT + frameX * deltaT), azimuth[0] ))
                    	print 'Time %6.2f : theta = %10.2f'   %((startT + frameX * deltaT), azimuth[0]*180.0/pi)
                    	tdoafp.write('%6.2f ' %(startT + frameX * deltaT))
                    	for k,v in sst.frontEnd.getTDOAs().items():
                    	    tdoafp.write('%s %e ' %(k,v) )
                    	tdoafp.write('\n')
		    	wasSource = True
                    	frameX += 1
                    	average += azimuth
		except TypeError:	# Probably raised by "sst.kalmanFilter.next(frameX)", you should rarely reach this.
		    LoadFlag = 1
		    pass

    	    elif options.arraytype==1:	# Circular array
		try:
                    polarX          = sst.kalmanFilter.next(frameX)

		    # Bound the angle as sometimes it really go wild for some reason
		    polarX[0] = bound_angle(polarX[0], boundaries[0][0], boundaries[0][1])
                    polarX[1] = bound_angle(polarX[1], boundaries[1][0], boundaries[1][1])

                    if not sst.kalmanFilter.isSourceObserved():
		    	print       'Time %6.2f : Nothing detected' %(startT + frameX * deltaT + TimeBias)
		    	frameX          += 1
		    	average         += polarX
		    	continue
		    else:
                    	posfp.write('Time %6.2f : (theta, phi) = (%e, %e)\n' %((startT + frameX * deltaT), polarX[0], polarX[1]))
                    	print       'Time %6.2f : (theta, phi) = (%6.2f, %6.2f)'   %((startT + frameX * deltaT),
                                                                             polarX[0] * 180.0 / numpy.pi, polarX[1] * 180.0 / numpy.pi)

                   	tdoafp.write('%6.2f ' %(startT + frameX * deltaT))
                    	for k,v in sst.frontEnd.getTDOAs().items():
                            tdoafp.write('%s %e ' %(k, v))
                    	tdoafp.write('\n')
                    	wasSource = True
		    	frameX   += 1
		    	average  += polarX

                except TypeError:       # Probably raised by "sst.kalmanFilter.next(frameX)", you should rarely reach this.
                    LoadFlag = 1
                    pass


            elif options.arraytype==2:  # Sphere array
            	print '[ERROR] No support for sphere array yet.'
            	exit()
            else:
            	print '[ERROR] Unsupported array type: ', options.arraytype
            	exit()

        except StopIteration:
	    if (flag==1) or (flag==-1) or (flag==0 and endSample==-1):
                print 'end'
                break;
        except:
	    print '[WARNING] There might be an un-addressed issue, probably caused by trying to process zero size block. Lazily skip it for now. You could avoid this by slightly adjust your segment boundaries in the external VAD output, which is not really necessary once the bug is fixed properly.'
	    print sys.exc_info()
            continue



	# Re-load features when needed
	if flag==0 and LoadFlag==1:
            tmp_st = tmp_et                 
            startSample = endSample         
	    if startSample==-1:
		break
	    startT = startSample*1.0/sampleRate
            try:
                tmp_et = TIME[SegIdx][0]
                endSample = int(tmp_et*sampleRate)
            except IndexError:
                tmp_et = -1
                endSample = -1
	    if startSample==endSample:
		print '[WARNING] Trying to load audio of zero duration, skipping it.'
		SegIdx += 1
		continue
            print '[LOG] Reload features from ', tmp_st, ' to ', tmp_et, '(in sec), or ', startSample, ' to ', endSample, ' (in sample index)'
	    print '[LOG] KF re-construction...'

	    # Rebuild the object
            del sst
	    del kalmanFilter
            sst = SingleSpeakerTracker(options)
	    kalmanFilter = sst.buildExtKFSourceTracker(initialXk, boundaries)
            sst.loadSamples(in_prefix, startSample, endSample)
            InitFlag = 0
	    LoadFlag = 0
	    frameX = 0
#           average   = numpy.zeros(len(initialXk), numpy.float)
	    wasSource = False

    	    if options.arraytype==0:    # Linear array: update the initialization
            	while True:
		    # -----------------
		    # NOTE:
		    # 	"getInstantaneousPosition" is implemented for linear array only, in
		    # 	class "FarfieldLinearArrayTDOAFeatureVector". That is why we are not
		    # 	doing the same thing for other arrays.
		    # -----------------
                    initialXk = sst.frontEnd.getInstantaneousPosition( frameX )
            	    if initialXk > -1e10:
                    	print 'Time %6.2f : theta = %10.2f' %((startT + frameX * deltaT), initialXk*180.0/pi)
                    	break
            	    frameX += 1
            	print 'Update initialXk for linear array to: ', initialXk
        	# We implicitly skipped some samples at the beginning which does not pass the initialXk>-1e10 test
		kalmanFilter = sst.updateExtKFSourceTracker(initialXk, boundaries)

    if wasSource == False:
        print 'Missed speech'

    posfp.close()
    tdoafp.close()



def main():
    if len(sys.argv)<2:
	print 'Please try option "--help" to get help information.'
	exit()

    print '\nExecute: \n\t' + ' '.join(sys.argv) + '\n'

    options = parse_option(sys.argv)
    vadfiles = []

    batchMode = float(options.batchMode)
    if (abs(batchMode) > 1e-5) and (batchMode>0):               # Re-initialize the KF at a constant rate
        print 'Re-initialize the KF at every ', batchMode, ' seconds.'
    elif abs(batchMode) < 1e-5:                                 # Read external segmentation information
        print 'Reading external segmentation information from files whose path is listed in file: ', options.vadList
	fid = open(options.vadList, 'r')
	for line in fid:
	    line = line.strip()
	    if line=='' or line[0] =='#':
		continue
	    vadfiles.append(line)
	fid.close()
    elif abs(batchMode-(-1)) < 1e-5:                                                    # Process each file in one batch
        print 'Process each audio file in one batch.'


    fid = open(options.utteranceList, 'r')
    i = 0
    for line in fid:
	line = line.strip()
	if line=='' or line[0]=='#':
	    continue
	in_prefix, out_prefix, = line.split()
	if len(vadfiles)>0:
	    track(in_prefix, out_prefix, options, vadfiles[i])
	else:
	    track(in_prefix, out_prefix, options)
	i += 1
    fid.close()

if __name__ == '__main__':
    main()


