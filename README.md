# DREndcapTube
DD4hep description of a dual-readout endcap calorimeter using the capillary tube technology

![DREndcapTube](https://github.com/user-attachments/assets/f1acdfd0-afd3-4279-b557-055692a43477)

## Presentations
- FCC Full Simulation Meeting, 24/7/2024, [DD4hep implementation of the IDEA Endcap Calo with the capillary tubes technology](https://indico.cern.ch/event/1439207/contributions/6056623/attachments/2903299/5092292/lopezzot_fccsim_2472024.pdf)

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
