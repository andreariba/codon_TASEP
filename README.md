# Simulation of the inhomogeneous TASEP model

Implementation of inhomogeneous TASEP model, the user can input ORF sequences and the elongation rate of each codon to simulate the ribosome dynamics.

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
the output 

## Authors

* **Andrea Riba** - *Initial work*

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


