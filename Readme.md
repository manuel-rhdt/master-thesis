# Gillespie

This program allows you to simulate trajectories given a chemical reaction network. It also allows to have time-varying
propensities.

## Usage

### Input file format

The input file format is human readable and has to be processed by the script `parse.pl` in this directory. A common
invocation would be `$ ./parse.pl gene-expr.txt > gene-expr.inp` which will in total generate three files from a given
input file `gene-expr.txt`:
- The first output file is `gene-expr.inp` which provides general information about the simulation run
- the second output file is `gene-expr.reactions` which contains information about all reactions
- the third output file is `gene-expr.components` which contains the names and initial concentrations of all components

A very simple input file for a gene expression reaction looks as follows:

```
# gene-expr.txt

gene-expr : name

S X: components

S == 1000 : initial
X == 1000 : initial

S --> X : k = 1.0
X --> NULL : k = 0.5
```

This describes a simulation named "gene-expr" with two reaction components, "S" and "X". Both reaction components contain
1000 copies at the beginning of the simulation. 2 reactions are simulated:
- The first reaction converts one "S" to a "X" with a reaction constant of k = 1.0.
- The second reaction describes the decay of one "X" to nothing with rate k = 0.5.

A more detailed and complete description of the input file format is provided as a comment in `parse.pl`.

### Command line usage

You can display the following help message by running `./Gillespie -h`:

```
Usage: Gillespie [options] INPUT

Options:
	-o FILE		    Set the output trajectory file
	--overwrite X=Y	Overwrite the initial value of component X with Y
	-s SEED		    Specify which seed to use
	-t TRAJECTORY	Load a trajectory to use for the simulation
```

The `Gillespie` program expects one input file as argument (This file typically has file extension `.inp`. If `-` is
provided as the input file name, `Gillespie` will read the input from standard input). The program then expects there
to exist a `[name].reactions` and a `[name].components` files to exist in the working directory where `[name]` refers
to the name specified in the input file.

If `-o` is not used to override it, the program will write a file `[name].traj` containing the simulated trajectory.
