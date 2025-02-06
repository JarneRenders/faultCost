# FaultCost
This program can be used to determine the fault cost, minimum leaf number and or leaf-guaranteedness of a given graph. See "J. Goedgebeur, J. Renders, G. Wiener, and C.T. Zamfirescu, Network fault costs based on minimum leaf spanning trees, manuscript"

The latest version of this program can be obtained from <https://github.com/JarneRenders/FaultCost>.

### Installation

This requires a working shell and `make`. On Windows an easy way to simulate this is by using Windows Subsystem for Linux (WSL).

- Download and extract [`nauty`](https://pallini.di.uniroma1.it/).
- Copy the file `splay.c` to the `utilities` folder.
- Compile the program using: 
  * `make` to create a binary for the 64 bit version
  * `make 128bit` to create a binary for the 128 bit version
  * `make all` to create all of the above

The 64 bit version can handle graphs up to 64 vertices, the 128 bit version up to 128 vertices.
The lower bit version is faster than the higher bit one, hence it is recommended to use the version which is strictly higher, but closest to the order of the graphs you want to inspect.


### Usage of faultCost

This helptext can be found by executing `./faultCost -h`.

Usage: `./faultCost [-m|-l#|-L] [-aco#O#v] [-h]`

Compute fault cost or properties related to number of leaves in a
spanning tree.

Graphs are read from stdin in graph6 format. Graphs are sent to stdout
in graph6 format.

If -m, -l# or -L is absent the fault cost will be computed.

```
    -a, --all
            send distinct degree sequences of every computed
            ml-subgraph to stderr
    -c, --count-branches
            count the number of branches in every computed ml-subgraph;
            introduces overhead
    -h, --help
            print this help message
    -l#, --k-leaf-guaranteed=#
            send all #-leaf-guaranteed graphs to stdout
    -L, --leaf-guaranteed
            send all leaf-guaranteed graphs to stdout
    -m, --ml-number
            compute the minimum leaf number of the input graphs;
            combine with -o# or -O# to send graphs with certain min
            leaf numbers to stdout
    -o#, --output=#
            combine with no arguments or -m to send graphs with fault
            cost # or ml number # to stdout respectively; can be used
            multiple times to output for more values; can be combined
            with -O#
    -O#, --output-from=#
            combine with no arguments or -m to send graphs with fault
            cost at least # or ml number at least # to stdout
            respectively; can be combined with -o#
            
    -v, --verbose
            send more verbose output to stderr
```