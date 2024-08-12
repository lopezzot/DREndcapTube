from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import cm, mm, GeV, MeV
SIM = DD4hepSimulation()

## Add sensitive detector action to calorimeter
SIM.action.calo = "DREndcapTubesSDAction"
SIM.action.calorimeterSDTypes = [u'calorimeter']
SIM.action.mapActions['DREndcapTubes'] = "DREndcapTubesSDAction"
SIM.filter.calo = ""
## end of calorimeter sensitive detector action

## The compact XML file, or multiple compact files, if the last one is the closer.
SIM.compactFile = ["DRdetector/DRcalo/compact/DREndcapTubes.xml"]
SIM.enableGun = True
## Macro file to execute for runType 'run' or 'vis'
SIM.macroFile = ""
## number of events to simulate, used in batch mode
SIM.numberOfEvents = 100
## Outputfile from the simulation: .slcio, edm4hep.root and .root output files are supported
SIM.outputFile = "edm4hep_test.root"
## Physics list to use in simulation
SIM.physicsList = "QGSP_BERT"
## Verbosity use integers from 1(most) to 7(least) verbose
## or strings: VERBOSE, DEBUG, INFO, WARNING, ERROR, FATAL, ALWAYS
SIM.printLevel = 3
## The type of action to do in this invocation
## batch: just simulate some events, needs numberOfEvents, and input file or gun
## vis: enable visualisation, run the macroFile if it is set
## qt: enable visualisation in Qt shell, run the macroFile if it is set
## run: run the macroFile and exit
## shell: enable interactive session
SIM.runType = "batch"

################################################################################
## Configuration for the Detector Construction. 
################################################################################
#SIM.geometry.dumpGDML = "DREndcapTubes.gdml"
#SIM.geometry.dumpHierarchy = 0

################################################################################
## Configuration for the DDG4 ParticleGun 
################################################################################

##  direction of the particle gun, 3 vector 
SIM.gun.direction = (0, 0, 1)

## choose the distribution of the random direction for theta
## 
##     Options for random distributions:
## 
##     'uniform' is the default distribution, flat in theta
##     'cos(theta)' is flat in cos(theta)
##     'eta', or 'pseudorapidity' is flat in pseudorapity
##     'ffbar' is distributed according to 1+cos^2(theta)
## 
##     Setting a distribution will set isotrop = True
##     
SIM.gun.distribution = None

## Total energy (including mass) for the particle gun.
## 
## If not None, it will overwrite the setting of momentumMin and momentumMax
SIM.gun.energy = None

## Maximal pseudorapidity for random distibution (overrides thetaMin)
SIM.gun.etaMax = None

## Minimal pseudorapidity for random distibution (overrides thetaMax)
SIM.gun.etaMin = None

##  isotropic distribution for the particle gun
## 
##     use the options phiMin, phiMax, thetaMin, and thetaMax to limit the range of randomly distributed directions
##     if one of these options is not None the random distribution will be set to True and cannot be turned off!
##     
SIM.gun.isotrop = False

## Maximal momentum when using distribution (default = 0.0)
SIM.gun.momentumMax = 5.0*GeV

## Minimal momentum when using distribution (default = 0.0)
SIM.gun.momentumMin = 0.0
SIM.gun.multiplicity = 1
SIM.gun.particle = "mu-"

## Maximal azimuthal angle for random distribution
SIM.gun.phiMax = None

## Minimal azimuthal angle for random distribution
SIM.gun.phiMin = None

##  position of the particle gun, 3 vector 
SIM.gun.position = (0.0, 200.0*cm, 0.0)

## Maximal polar angle for random distribution
SIM.gun.thetaMax = None

## Minimal polar angle for random distribution
SIM.gun.thetaMin = None
