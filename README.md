# README: The Extended Finite Element Method (XFEM) plugin in C++#
The Extended Finite Element Method (XFEM), allows for discontinuities - appearing for example after a fracture or the cutting of an object - in simulations based on the Finite Element Method (FEM). This repository contains an implementation of XFEM in C++, that can be added to an existing FEM framework. For an introduction into the XFEM and detailed information on our implementation, we refer to Paulus et al - ["Simulation of Complex Cuts in Soft Tissue with the Extended Finite Element Method (X-FEM)"](https://journals.ub.uni-heidelberg.de/index.php/emcl-pp/article/view/17635).

In this README, we explain our XFEM implementation by means of a simple example: a low level implementation of the Finite Element Method (FEM), that can be replaced by your FEM framework.

The python script simplifiedUsage.py simplifies the (first) usage of the code in the git repository.
For more information see the steps below and/or command: python simplifiedUsage.py -h

### 1 - Download ###
* install git: sudo apt-get install git
* clone the repository: git clone https://chrijopa@bitbucket.org/chrijopa/xfem-in-cpp.git

### 2 - Installation of third party libraries ###
Using the python script (prerequisit: python installed and 1): python simplifiedUsage.py -l

Without using the python script:

* install necessary third party libraries: sudo apt-get install libeigen3-dev libvtk5-dev libcgal-dev python-matplotlib libtinyxml-dev paraview cmake make

### 3 - Compilation ###
The compilation has been tested on Ubuntu 14.04 using g++ 4.8.2. If you encounter any problems compiling on your system, please contact ChristophJPaulus@gmail.com

Using the python script (prerequisit: python installed and 1): python simplifiedUsage.py -lc

Without using the python script:

* prerequisits: 1 and 2
* go into the folder where you cloned the xfem folder and make a new directory: mkdir build
* go into the directory: cd build
* run cmake (if not installed, use: sudo apt-get install cmake) : cmake ..
* compile the plugin: make

### 4 - Launch the example Liver.xml ###
Using the python script (prerequisit: python installed and 1): python simplifiedUsage.py -lcpoe

Without using the python script:

* prerequisits: 1, 2 and 3
* update the global paths in scenes/Liver.xml (replacing PATH)
* introduce the output folders: mkdir -p output/Liver/displacedObject output/Liver/stiffness
* Launch the example: build/XFEMExec scenes/Liver.xml

The result of the calculation can be opened using paraview.

### 5 - Start the evaluation/analysis ###
The evaluation can take a while to calculate. If you need to accelerate the computation, please use -l in evaluateConvergence.py to allow for parallization.

Using the python script (prerequisit: python installed and 1): python simplifiedUsage.py -lca

Without using the python script:

* prerequisits: 1, 2 and 3
* lance the evaluation: python evaluateConvergence.py FEM.json XFEM.json