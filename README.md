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

Explain how to run the automated tests for this system

```

simulator/TASEP_simulator.exe -c <codonspeedfile> -f <elongationspeed> -n <ORFname> -s <ORFsequence>

```
the output includes a file \"\<ORFname\>.dat\" containing six columns

| factor | initiation rate | protein synthesis rate | ribosome density | fraction of queuing events | vector of ribosome density per codon |
| 6 | 0.01 | 0.00975261 | 0.00187338 | 0.00150791 | [0.000401433,0.00284735,0.00206447, ... ] |
| --- | --- | --- | --- | --- | --- |

## Authors

* **Andrea Riba** - https://scholar.google.ch/citations?user=aZ9tkZgAAAAJ&hl=en


