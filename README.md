# Stochastic Progressive Photon Mapping

Implementation of the Stochastic Progressive Photon Mapping algorithm.

## Requirements

* CMake
* OpenMP support (optional)

## Usage

To compile:

``` 
cmake .
make
```

To run:

```
cd bin
./sppm <num_rounds> <num_photons_per_round>
```

Recommened setting for the sample scene is 2500 rounds and 200000 photons per round.

## Reference

Hachisuka, T., & Jensen, H. W. (2009). Stochastic progressive photon mapping. *ACM Transactions on Graphics (TOG)*, *28*(5), 141.

