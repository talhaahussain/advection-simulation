# advection-simulation

Continuous Assessment for ECM3446 - High Performance Computing, set by Prof. Luo Man (Year 3, Semester 2). Simulates the advection of material from a chimney by wind in the planetary boundary layer. The main program is optimised for high performance computers by parallelisation, using OpenMP. This repository includes this file, a single C program, three gnuplot scripts, and four PNG images (located in the results subdirectory).

This work received a final mark of 100/100. 

Please see `specification.pdf` for specification.

### Prerequisites

The following instructions assume that you are trying to run the program and view the results with a Linux machine, with gcc, OpenMP, gnuplot.

### Usage

The compiled program has NOT been included - please compile the serial program with

```
gcc -o advection2D -std=c99 advection2D.c -lm 
```

or the parallelised program on an OpenMP-enabled machine with

```
gcc -fopenmp -o advection2D -std=c99 advection2D.c -lm
```

Please run the compiled program with 

```
./advection2D
```

As aforementioned, gnuplot scripts have been provided. Please note that you must first compile and run the program, and generate initial.dat, final.dat and vertical\_averages.dat. To use these to view the output of the program, please use
```
gnuplot plot_initial
```
```
gnuplot plot_final
```
```
gnuplot plot_average
```

or, alternatively (assuming the scripts are in the current directory)

```
gnuplot plot_*
```

These will generate 3 images - initial.png, final.png and vertical\_averages.png. You can view these with an image viewer of your choice, such as feh. If you are unable to generate these images by executing the program, please refer to the results directory.

The solution plots to each section are provided in the results directory. You can view these with an image viewer of your choice.
