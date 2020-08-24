"""
Calculate the mobility transfer function 
run as python3 mobTrans.py [path] [file name] [number of frames] [0: all particles or A or B]
"""
import numpy as np
import pylab as pl
import sys                                  # to be able to load the file from the command line
path='../facilitation/'
sys.path.append(path)
import singPartDist as nf


cutoff = 1.3 			# cutoff for particles being 'close' to each other (~ first min of g(r))
rho = 1.4                       # number density for periodic boundaries
numFrames = int(sys.argv[3])	# number of frames
path = sys.argv[1]		# path to input file
filexyz = sys.argv[2] 		# coordinate file 
partType = sys.argv[4]	# type of particles to compute: 0: all particles, 'A': type A particle: 'B': type B particles
if partType == '0':
	numPart = 10002
	partLabel = ['A','B'] # for loading coordinates of A and B particles
elif partType =='A': 
	numPart = 6668
	partLabel = ['A','A']
else: 
	numPart = 3334
	partLabel = ['B','B']
side  = (numPart / rho)**(1/3)	# lenght of box
numFast = int(numPart*0.1)	# number of particles that are considered 'fast'


def readCoords(filexyz, numFrames, numPart):
        """
        Reads data from an xyz file
        Args:
                filexyz(string): name of the xyz file to read
                numFrames (int): number of frames in file
                numPart (int): number of particles
        Return:
                allCoords (list of list) for each frame a list of all 
                particles consisting of list of all three coordinates
                for each particle (x,y,z)
        
        """

        frame = -1
        allCoords = np.zeros((numFrames,numPart,3))
        with open(filexyz, 'r') as readFile:
                for line in readFile:
                        splitL = line.split()
                        if len(splitL) ==1:
                                frame +=1
                                if frame == numFrames:
                                        break
                                particleCounter = 0
                        elif (not splitL[0] == 'Atoms.') and (splitL[0] == partLabel[0] or splitL[0] == partLabel[1]):
                                allCoords[frame][particleCounter,0] =splitL[1]
                                allCoords[frame][particleCounter,1] =splitL[2]
                                allCoords[frame][particleCounter,2] =splitL[3]
                                particleCounter+=1
        return allCoords
allCoords = readCoords(path+filexyz,numFrames,numPart)

def selectRand(allCoords,delta,numFast,numFrames,fastPart):
	"""
	Randomly reassignes particle types in the same ratio as before (for chosen particle type)
	"""
	randFast = []
	fastPart = np.array(fastPart)
	for i,frame in enumerate(range(delta,numFrames,delta)):
		Ids = np.array(range(numPart))
		Ids = Ids[np.in1d(Ids,fastPart[i,:])==False]	
		randType = np.random.random_sample(size=len(Ids))#numPart-fastPart.shape[1])
		randType= (randType<(numFast/(numPart-numFast))).astype(int)
		randFast.append(Ids[np.where(randType == 1)[0]])
	return randFast

def minDistPart(allCoords,ID,IDlist,frame,L):
	"""
	Calculate the distances of one particle to a set of particle in the previous frame
	and compute the minimum. 
	"""
	distances = allCoords[frame][ID,:]-allCoords[frame][IDlist,:]
	for dist in range(len(distances)):
		distances[dist,:] = nf.periodic_boundary(distances[dist,:],L)
	dist1D = np.sqrt(distances[:,0]**2 + distances[:,1]**2 + distances[:,2]**2)
	return min(dist1D)


def selectFast(allCoords,delta,numFast,numFrames,side):
	"""
		split trajectory into bits of length delta and identify the numFast fastest particles 
	"""
	fastPart = []
	for frame in range(delta,numFrames,delta):
		dist =np.array([nf.squareDist(allCoords[:,particle,:],frame-delta,frame,side) for particle in range(numPart)])
		fastPart.append(dist.argsort()[-numFast:])
	return fastPart


def distFast(fastPart,refPart,allCoords,delta):
	"""
	Calculate the distance between minimum distance between sets of particles
	"""
	minDists = []
	for frame in range(1,len(fastPart)):
		for partID in fastPart[frame]: #do I compare with the correct sample for random???
			if (not  partID in fastPart[frame-1]) and (not partID in refPart[frame-1]):
				distMin = minDistPart(allCoords,partID,refPart[frame-1],frame*delta,side)
				minDists.append(distMin)
	return minDists
	
for delta in [200,400]:
	pl.figure()
	print('length of trajcectory section delta is %.3d frames'%(delta))
	fastPart = selectFast(allCoords,delta,numFast,numFrames,side)	# select all fast particles for all frames
	randPart = selectRand(allCoords,delta,numFast,numFrames,fastPart)	# randomly select particles that are not fast for each frame
	minDistsFast = np.array(distFast(fastPart,fastPart, allCoords,delta))
	minDistsRand = np.array(distFast(fastPart,randPart, allCoords,delta))	# Calculate distance between particles in different frames
	pl.hist(minDistsFast,bins = 20,density = True,histtype = 'step', color = 'red',label = 'fast')	#plot distribution of minimum distance
	pl.hist(minDistsRand,bins = 20,density= True,histtype = 'step', color = 'blue', label  = 'rand')
	pl.legend(frameon=False)
	pl.title('delta= '+str(delta))
	print('Fast particles with previous fast neighbours: %.2f percent'%(len(minDistsFast[minDistsFast<cutoff])/len(minDistsFast)))		#compute ratio of close particles
	print('Fast particles with previous selection of random neighbours: %.2f percent'%(len(minDistsRand[minDistsRand<cutoff])/len(minDistsRand)))		#compute ratio of close particles
pl.show()
