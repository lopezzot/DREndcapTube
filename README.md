# DREndcapTube
DD4hep description of a dual-readout endcap calorimeter using the capillary tube technology

![DREndcapTube](https://github.com/user-attachments/assets/f1acdfd0-afd3-4279-b557-055692a43477)

## Presentations
- FCC Full Simulation Meeting, 24/7/2024, [DD4hep implementation of the IDEA Endcap Calo with the capillary tubes technology](https://indico.cern.ch/event/1439207/contributions/6056623/attachments/2903299/5092292/lopezzot_fccsim_2472024.pdf) (errata corrige: slide 5 "~1.1 Gb" -> "~2.1 Gb")

## How to
### Build and compile using the key4hep stack
```sh
git clone https://github.com/lopezzot/DREndcapTube.git
source /cvmfs/sw.hsf.org/key4hep/setup.sh
cd DREndcapTube
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../install/ ..
make install
cd ../install
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib
```

### Visualize the geometry
```sh
geoDisplay DRdetector/DRcalo/compact/DREndcapTubes.xml
```
or
```sh
geoWebDisplay DRdetector/DRcalo/compact/DREndcapTubes.xml
```

### Convert geometry to root and browse it
Convert the geometry to root
```sh
./scripts/dd4hep2root -c DRdetector/DRcalo/compact/DREndcapTubes.xml -o DREndcapTubes.root
```
browse the root file
```sh
rootbrowse DREndcapTubes.root
```

### Check geometry overlaps
The following command converts the geometry to geant4 and checks for overlaps
```sh
ddsim --compactFile DRdetector/DRcalo/compact/DREndcapTubes.xml --runType run --macroFile scripts/overlap.mac --part.userParticleHandler=''
```

### Scan material along a line
The following command performs a material scan along a line from (0,0,0) to (10,200,500), unit must be cm.
```sh
materialScan ./DRdetector/DRcalo/compact/DREndcapTubes.xml 0 0 0 10 200 500
```

### Run simulations
The following command does a simple simulation with the particle gun.
```sh
ddsim --compactFile DRdetector/DRcalo/compact/DREndcapTubes.xml --enableGun --gun.particle geantino --gun.energy 1000*MeV --gun.direction "0 0 -1" --gun.position "0 200*cm 0" --outputFile out_edm4hep.root -N 100 --part.userParticleHandler=""
```
The following command does a simulation according to a steering file.
```sh
ddsim --steeringFile scripts/steer_mu_0_5GeV.py
```
a steering file template can be obtained with `ddsim --dumpSteeringFile`.
A geometry description `gdml` file can be created using the dedicated entry in the steering file.\
The following command does a simulation with optical photons and a DD4hep regex-sensitive-detector.
```sh
ddsim --steeringFile scripts/steer_e-_10GeV_regexSDAction.py
```

### Alternative approach to build, compile and visualize geometry using a local DD4hep installation
I experienced crashes while visualizing the geometry over ssh connection to lxplus-alma9 machines.
The following are instructions on how to visualize the geometry using a local DD4hep installation.
```sh
git clone https://github.com/lopezzot/DREndcapTube.git
cd DREndcapTube/DRdetector/DRcalo/
mkdir build && cd build
cmake ../onlygeocmake/
make
source bin/thisDREndcapTubes.sh
geoDisplay ../compact/DREndcapTubes.xml
```
