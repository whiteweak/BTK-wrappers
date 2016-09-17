# BTK-wrappers
The wrappers created by Yulan Liu are based on an original version from CMU.
BTK and its orignal wrappers from CMU could be accessed via the "DOWNLOAD"
link in the home page of BTK:

	http://distantspeechrecognition.sourceforge.net

# Major changes made compared to the CMU version of BTK wrappers

	1. Enables a portable python interface. Original wrappers embed all
	   configuration and parameters in the python wrappers. This makes
	   it quite difficult to setup a new experiments for a new user, 
	   as it requires a decent knowledge of BTK itself and a decent level
	   of expertise in beamforming and dereverberation. The wrappers here
	   separates parameters from scripts, and all python wrapper scripts
	   accept parameters in both command format and configuration file 
	   format. Thus it allows various experiments configurations that 
	   could be transferred and replicated easily.

	2. Avoids configuration headache in large scale experiments. Original
	   wrappers embed all parameters in the script, as well as a pre-
	   assumed folder structure for audio files which are segmented per
	   utterance. This version of wrapper reads list files (similar to
	   .scp file in HTK) so that parallel jobs running in cluster could
	   share the same configuration while with different input/output
	   file paths.


# A few extra features of the wrappers

	1. Compatibility with HTK. Some file format is borrowed from HTK
	   equivalents, like ".scp" file and ".mlf" file. Please refer to
	   the comments in each python wrapper script for more details.

	2. Compatibility with BeamformIt. Since BeamformIt has gained a wide
	   range of popularity due to its simplicity in use and robustness
	   in performance, there are some scripts that convert the TDOA
	   from BeamformIt to BTK format for cross comparison, potential
	   combination and improvement.

# More information
For more details about BTK, please visit the following page:

	http://distantspeechrecognition.sourceforge.net

Please note that BTK and its wrappers are all based on GPL license. 

# Contact
Yulan Liu (yulan.liu.wings@foxmail.com)

