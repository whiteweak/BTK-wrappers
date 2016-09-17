#!/usr/bin/env python

#
#                           Beamforming Toolkit
#                                  (btk)
#
#   Purpose: Collection of frequently used functions by beamforming wrappers
#   Authors: Yulan Liu, John McDonough and Kenichi Kumatani
#   Contact: yulan.liu.wings@foxmail.com
#   Date:    June 16, 2016


import numpy, pickle, gzip

def load_geometry(file):
   '''
   Load the geometry information from file.
   '''
   fid = open(file, 'r')
   geometry = []
   for line in fid:
        line = line.strip()
        if line=='':
            continue
        if line[0]=='#':        # Arraytype information
            arraytype = int(line.split(':')[-1].strip())
        else:
            geometry.append([float(x) for x in line.split()])
   fid.close()
   return arraytype, geometry


def load_geometry_force3D(file):
   '''
   Load the geometry information from file.
   '''
   fid = open(file, 'r')
   geometry = []
   for line in fid:
        line = line.strip()
        if line=='':
            continue
        if line[0]=='#':        # Arraytype information
            arraytype = int(line.split(':')[-1].strip())
        else:
            tmp = [float(x) for x in line.split()]
            if len(tmp)<3:      # NOTE: This is different from "load_geometry"
                tmp += [0]*(3-len(tmp))
            geometry.append(tmp)
   fid.close()
   return arraytype, geometry



def pos2tdoa_circular(avgTheta, avgPhi, arraygeometry, refX, chX):
    '''
    Calculate the TDOA according to the geometry information of microphone array and the
    arriving angle of the sound signal, assuming that the sound wave has reached far-field.
    This function is particularly for circular array.
    '''
    sspeed   = 343000.0         # Sound speed (mm/s)

    if chX==refX:       # Reference channel and target channel are the same
        return 0.0
    else:
        v1 = - numpy.array([numpy.sin(avgTheta)*numpy.cos(avgPhi), numpy.sin(avgTheta)*numpy.sin(avgPhi), numpy.cos(avgTheta)])     # Unit vector of the sound arrival direction
        v2 = numpy.array(arraygeometry[refX]) - numpy.array(arraygeometry[chX])     # Vector in the direction from reference channel to current channel
        v1_l = 1
        v2_l = numpy.sqrt(sum(numpy.multiply(v2,v2)))                                               # Distance between two channels
        dist = v2_l
        v1_dot_v2 = sum(numpy.multiply(v1,v2))
        cos_avgDOAs = v1_dot_v2/(v1_l*v2_l) # Cosine of the angle between the two vectors
        timedelay = - numpy.sign(v1_dot_v2) * dist * cos_avgDOAs / sspeed		# The minus sign is needed because that is how other modules in BTK works! Confusing? Yes!
        return timedelay



def pos2tdoa_linear(avgDOA, arraygeometry, refX, chX):
    '''
    Calculate the TDOA according to the geometry information of microphone array and the
    arriving angle of the sound signal, assuming that the sound wave has reached far-field.
    This function is particularly for circular array.
    '''
    sspeed   = 343000.0         # Sound speed (mm/s)

    if chX==refX:       # Reference channel and target channel are the same
        return 0.0
    else:
        dist = abs( arraygeometry[chX][0] - arraygeometry[refX][0] )		# This is from original codes but it is probably wrong if you are not careful with avgDOA
	timedelay = dist * numpy.cos( avgDOA ) / sspeed
        return timedelay


def readTDOAFile( tdoaFile, mode=0, refX=0 ):
    '''
    Read TDOA file (.tdoa files), with 0-th channel as reference channel
    by default. It only supports two modes:

	mode = 0	Read as it is.
	mode = 1	Return TDOA smoothed over the whole file.

    '''
    print 'Reading %s' %(tdoaFile)

    TIME = []
    TDOA = []
    fid = open(tdoaFile)
    pre_t = 0.0
    for line in fid:
 	line = line.strip()
	if line=='' or line[0]=='#':
	    continue
	data = line.split()
	t = float(data[0])
	tdoa_dict = dict()
	for i in xrange(1, len(data), 3):
	    pair = [int(data[i]), int(data[i+1]), float(data[i+2])]
	    if refX==pair[0]:
		tdoa = pair[-1]
		try:
		    tdoa_dict[pair[1]]
		    print '[WARNING] Duplicated record for channel pair: (', pair[0], pair[1], ')'
		    print '          Skipping line: ', line
		    continue
		except KeyError:
		    tdoa_dict[pair[1]] = tdoa
	    elif refX==pair[1]:
		tdoa = -pair[-1]
		try:
		    tdoa_dict[pair[0]]
                    print '[WARNING] Duplicated record for channel pair: (', pair[0], pair[1], ')'
                    print '          Skipping line: ', line
                    continue
		except KeyError:
		    tdoa_dict[pair[0]] = tdoa
	tdoa_list = []
	for ch in sorted(tdoa_dict.keys()):
	    if refX==(ch-1):
		tdoa_list.append(0.0)
	    tdoa_list.append(tdoa_dict[ch])
	if refX==ch+1:
	    tdoa_list.append(0.0)
	TDOA.append(tdoa_list)

	if pre_t==t:	# t=pre_t=0.0
	    continue
	else:
	    TIME.append([pre_t, t])
	    pre_t = t
    TIME.append([t, -1])
    fid.close()
    if mode==0:
	return TDOA, TIME
    elif mode==1:
	aver_TDOA = []
	for i in xrange(len(TDOA[0])):
	    aver_TDOA.append(mean([TDOA[j][i] for j in xrange(len(TDOA))]))
	return [aver_TDOA], [[0, -1]]


def readPositionFile( posFile, arraygeometry, arraytype, mode=0, refX=0):
    '''
    Read speaker position file (.txt files from Kalman filter speaker tracking),
    with 0-th channel as reference channel by default.

        mode = 0        Read as it is.
        mode = 1        Return position smoothed over the whole file.

    '''

    print 'Reading %s' %(posFile)
    chanN    = len(arraygeometry)

    # Linear array
    AZIMUTH  = []

    # Circular array
    THETA    = []
    PHI      = []

    sspeed   = 343000.0         # Sound speed (mm/s)
    TIME     = []               # Collection of raw time index
    TDOA     = []               # Collection of raw TDOAs

    fp = open( posFile, 'r' )
    for line in fp:
        line = line.strip()
        if line=='':
            continue
        elems1  = line.split(':')
        timestamp = float(elems1[0].split()[1])
        TIME.append(timestamp)

        if arraytype==0:        # Linear array
            # Time   1.28 : 1.214301e+00
            azimuth = float(elems1[1])
            AZIMUTH.append(azimuth)

            tmp = []
            for chX in range(chanN):
                tmp.append( pos2tdoa_linear(azimuth, arraygeometry, refX, chX) )
            TDOA.append(tmp)

        elif arraytype==1:      # Circular array
            # Time   0.77 : (theta, phi) = (1.466826e+00, -4.734720e-01)
            data = [float(x) for x in elems1[1].split('(')[-1].split(')')[0].split(',')]
            theta = data[0]
            phi = data[1]
            THETA.append(theta)
            PHI.append(phi)

            tmp=[]
            for chX in range(chanN):
                tmp.append( pos2tdoa_circular(theta, phi, arraygeometry, refX, chX) )
            TDOA.append(tmp)

        elif arraytype==2:      # Sphere array
            print '[ERROR] Currently not supporting sphere array yet.'
            exit()
        else:
            print '[ERROR] Unrecognized arraytype: ', arraytype
            exit()
    fp.close()

    # Calculate TDOA from averaged angle when required
    if mode==0:
        print 'Using the position file as it is.'
        if arraytype==0:        # Linear array
            return numpy.array(TDOA), [ [TIME[i], TIME[i+1]] for i in range(len(TIME)-1)], [AZIMUTH]
        elif arraytype==1:      # Circular array
            return numpy.array(TDOA), [ [TIME[i], TIME[i+1]] for i in range(len(TIME)-1)], zip(THETA, PHI)
    elif mode==1:
        print 'Average the position output over the whole file.'
        if arraytype==0:        # Linear array
            aver_azimuth = numpy.array(AZIMUTH).mean(axis=0)
            tmp = []
            for chX in range(chanN):
                tmp.append( pos2tdoa_linear(azimuth, arraygeometry, refX, chX) )
            aver_TDOA = tmp
            return numpy.array([aver_TDOA]), [[0, -1]], [[aver_azimuth]]
        elif arraytype==1:      # Circular array
            aver_theta = numpy.array(THETA).mean(axis=0)
            aver_phi = numpy.array(PHI).mean(axis=0)
            tmp=[]
            for chX in range(chanN):
                tmp.append( pos2tdoa_circular(theta, phi, arraygeometry, refX, chX) )
            aver_TDOA = tmp
            return numpy.array([aver_TDOA]), [[0, -1]], [[aver_theta, aver_phi]]
        elif arraytype==2:      # Sphere array
            print '[ERROR] Currently not supporting sphere array yet.'
            exit()
        else:
            print '[ERROR] Unrecognized arraytype: ', arraytype
            exit()
    else:
        print '[ERROR] Unrecognized mode: ', mode
        print '        Try "--help" to get more information.'
        exit()


def myLoadWeights( fileName, fftLen, chanN, NC=1, idx=0):
    '''
    Load the weights file generated by script "estimateWeights.py".

    Difference to the original version:

     - To support large audio files with a lot of utterances, the input
    weight file might have the statistics for many utterances. Thus idx
    is used to find the correct place to load.
    '''
    fp = gzip.open(fileName, 'rb',1)
    waH = []
    packWa = []
    try:
        counter = 0
        while True:
            if counter==idx:
                nFrame4Adaptation = pickle.load(fp)
                waHK = pickle.load(fp)
                for fbinX in range(fftLen/2+1):
                    packWa.append( numpy.zeros(2 * (chanN - NC), numpy.float) )
                    for chX in range(chanN-NC):
                        packWa[fbinX][2*chX]     = waHK[fbinX][0][chX].real
                        packWa[fbinX][2*chX + 1] = waHK[fbinX][0][chX].imag
                waH.append( packWa )
                print nFrame4Adaptation, chanN, NC
                break
            else:
                pickle.load(fp)
                pickle.load(fp)
                counter += 1
    except EOFError:
        print 'Loaded the weights of %d frames' %counter
    except:
        print 'cannot read %s' %fileName

    fp.close()

    return [nFrame4Adaptation,waH]

