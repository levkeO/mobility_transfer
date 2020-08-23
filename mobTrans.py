"""
Calculate the mobility transfer function 
run as python3 mobTrans_v2.py [path] [file name] [number of frames] [number of particles]
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
numPart = 10002			# number of particles
numFast = int(numPart*0.1)	# number of particles that are considered 'fast'
path = sys.argv[1]		# path to input file
filexyz = sys.argv[2] 		# coordinate file 
side  = (numPart / rho)**(1/3)	# lenght of box
numA =3334			# number of A particles to compute particle types separately


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
                        elif (not splitL[0] == 'Atoms.') and (splitL[0] == 'B' or splitL[0] == 'B'):
                                allCoords[frame][particleCounter,0] =splitL[1]
                                allCoords[frame][particleCounter,1] =splitL[2]
                                allCoords[frame][particleCounter,2] =splitL[3]
                                particleCounter+=1
        print(particleCounter)
        return allCoords
allCoords = readCoords(path+filexyz,numFrames,numA)
print(allCoords.shape)

def selectRand(allCoords,delta,numFast,numFrames,fastPart):
	"""
	Randomly reassignes particle types in the same ratio as before (for chosen particle type)
	SOMETHIng is wrong there should be an evaluattion by frame!
	"""
	randFast = []
	fastPart = np.array(fastPart)
	for i,frame in enumerate(range(delta,numFrames,delta)):
		Ids = np.array(range(numPart))
		Ids = Ids[np.in1d(Ids,fastPart[i,:])==False]	
		#print(Ids.shape, Ids.max())
		#print(Ids)
		randType = np.random.random_sample(size=len(Ids))#numPart-fastPart.shape[1])
		#print(len(randType),'len rand')
		#print(randType.shape)
		randType= (randType<(numFast/(numPart-numFast))).astype(int)
		randFast.append(Ids[np.where(randType == 1)[0]])
		#print(np.array(randFast[i]))
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
	print('numFrames',numFrames)
	for frame in range(delta,numFrames,delta):
		if frame%100 ==0:
			print(frame)
		dist =np.array([nf.squareDist(allCoords[:,particle,:],frame-delta,frame,side) for particle in range(numPart)])
		fastPart.append(dist.argsort()[-numFast:])
	return fastPart


def distFast(fastPart,refPart,allCoords,delta):
	"""
	Calculate the distance between minimum distance between sets of particles
	"""
	minDists = []
	print(np.array(refPart).shape)
	for frame in range(1,len(fastPart)):
		counter =  0
		for partID in fastPart[frame]: #do I compare with the correct sample for random???
			if (not  partID in fastPart[frame-1]) and (not partID in refPart[frame-1]):
				distMin = minDistPart(allCoords,partID,refPart[frame-1],frame*delta,side)
				minDists.append(distMin)
				counter+=1
		#print(counter)		
	return minDists
	
numPart = numA
numFast = int(numA*0.1)
for delta in [200]:
	pl.figure()
	fastPart = selectFast(allCoords,delta,numFast,numFrames,side)	# select all fast particles for all frames
	randPart = selectRand(allCoords,delta,numFast,numFrames,fastPart)	# randomly select particles that are not fast for each frame
	minDistsFast = np.array(distFast(fastPart,fastPart, allCoords,delta))
	minDistsRand = np.array(distFast(fastPart,randPart, allCoords,delta))	# Calculate distance between particles in different frames
	pl.hist(minDistsFast,bins = 20,normed = True,histtype = 'step', color = 'red',label = 'fast')	#plot distribution of minimum distance
	pl.hist(minDistsRand,bins = 20,normed= True,histtype = 'step', color = 'blue', label  = 'rand')
	pl.title('delta= '+str(delta))
	pl.legend(frameon=False)
	print(np.array(minDistsFast).mean())
	print(np.array(minDistsRand).mean())
	print('delta:',delta,'fast: ',len(minDistsFast[minDistsFast<cutoff])/len(minDistsFast))		#compute ratio of close particles
	print('delta',delta,'random: ',len(minDistsRand[minDistsRand<cutoff])/len(minDistsRand))
pl.show()
