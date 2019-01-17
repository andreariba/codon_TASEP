# Simulation of the inhomogeneous TASEP model

Implementation of inhomogeneous TASEP model, the user can input ORF sequences and the elongation rate of each codon to simulate the ribosome dynamics. 
Supplementary material related to the paper https://www.biorxiv.org/content/early/2018/11/08/465914

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

g++ compiler is needed. The current version has been tested on GCC v4.8.4 and operating system CentOS v7.5.1804


### Installing

To compile the program you have to run the following commands:

```
cd simulator

make

cd ..

```

## Running

Run the simulation for a specific ORF typing the following command

```

simulator/TASEP_simulator.exe -c <codonspeedfile> -f <elongationspeed> -n <ORFname> -s <ORFsequence>

```
the output includes a file \"\<ORFname\>.dat\" containing six columns:

1. factor multiplying the codon speeds (default=6)
2. initiation rate
3. protein synthesis rate
4. ribosome density
5. fraction of queuing events
6. vector of ribosome density per codon


## Authors

* **Andrea Riba** - https://scholar.google.ch/citations?user=aZ9tkZgAAAAJ&hl=en


