# DiceForge

## Overview
This is a simple C++ library for *Pseudo Random Number Generation*. It provides a functionality similar to that of the python standard library *random*. 

## Features

- **Ease of Use**: Straightforward API for quick integration into your C++ projects.
- **Versatility**: Supports different PRNG algorithms, allowing you to choose based on your requirements.
- **Convenience**: Provides various standard distributions used in scientific computation and also has rudimentary support for custom probability distributions. 

### Upcoming Features
- To be decided

## Available PRNGs and Distributions

**Pseudorandom Number Generators**
1. Mersenne Twister (MT)
2. Linear Feedback Shift Register (LFSR)
3. XOR-Shift (XOR)
4. Blum Blum Shub (BBS)
5. Naor-Reingold (NR)

**Distrbutions**

1. Cauchy
2. Exponential
3. Gaussian
4. Maxwell
5. Weibull
6. Bernoulli
7. Binomial
8. Gibbs
9. Hypergeometric
10. Negative-Hypergeometric
11. Poisson

## Benchmarks

### Diehard Tests

Add results here

### Uniformity Analysis

Using DiceForge's pseudo random number generators to create a uniformly distributed random variable in [0,1), the statistical analysis for *100000000* random samples is:

| Generator | Mean | Variance |
| --------- | ---- | -------- |
| BBS32  | 0.501565 | 0.0830039 |
| BBS64  | 0.501565 | 0.0830038 |
| XORShift32  | 0.499978 | 0.08334   |
| XORShift64  | 0.499931 | 0.0833327 |
| MT32   | 0.499956 | 0.0833425 |
| MT64   | 0.500034 | 0.0833241 |
| LFSR32 | 0.500006 | 0.0833379 |
| LFSR64 | 0.499971 | 0.0833363 |
| NaorReingold | 0.49824  | 0.0827696 |

The theoretical mean and variance for a uniformly distributed random variable in [0,1) are \
mean = 1/2 (0.5) \
variance = 1/12 (0.08333...)

### Time performance

Time taken for generating *100000000* random numbers:

#### DiceForge:
| Generator | Time Taken |
| --------- | ---------- |
| BBS32 | 15210.2 ms |
| BBS64 | 24677.3 ms |
| XORShift32 | 484.018 ms |
| XORShift64 | 441.947 ms |
| MT32 | 1742.37 ms |
| MT64 | 1751.81 ms |
| LFSR32 | 10000.9 ms |
| LFSR64 | 18665.6 ms |
| NaorReingold | 195678 ms  |

For comparison, benchmarking other existing standard libraries for the same test.

| Generator | Time Taken |
| --------- | ---------- |
| C++'s Mersenne Twister | 4032.7 ms |
| C rand() function| 1388.03 ms |
| python's random | 175175.89 ms (~ 3min) |
| numpy's randint | 165600.48 ms (~ 3min) |


## Documentation

Check out the [Documentation](https://www.overleaf.com/project/65d9eea60dbb4690fe6ff8be) for detailed information on library usage, supported algorithms, and more!

## Getting Started

### Prerequisites

- C++ compiler (minimum C++17 or later, note that C++20 is required for using 2D probability functions)
- CMake (>=3.28) [If building library]

### Installation

#### Using a prebuilt library
1. Download the compiled binaries from [out](out) folder along with the include files [include](include)
2. If you want compact compilation commands later on, add the path to *include/* to the environment variable PATH
3. Start generating pseudo-random numbers!

#### Building it yourself
1. Clone the repository: `git clone https://github.com/yourusername/diceforge-library-clone.git`
2. Create the build folder by using the CMakeLists.txt
2. Build library using the given CMake configurations (```cmake --build <build-folder>```)

#### After successfully building you can
* Install the library to *usr/local/* (```cmake --install <build-folder>```) 
* Use it without installation by following the steps mentioned in [Using a prebuilt library](#using-a-prebuilt-library)

## Usage

Here's a quick example to get you started:

*example.cpp*
```cpp
#include "diceforge.h"

int main() {
    // Create a PRNG object
    DiceForge::XORShift32 prng = DiceForge::XORShift32(123);

    // Generate and print a random number
    std::cout << "Random Number: " << prng.next() << std::endl;

    return 0;
}
```

Now to compile it with a compiler of your choice (we use g++ in the following example, also it is assumed that g++ has access toand searches usr/local for installed libraries), enter the following command in your terminal

1. If you installed diceforge
```bash
g++ example.cpp -ldiceforge -o example.out
```

2. If you did not install diceforge and you included diceforge to PATH
```bash
g++ example.cpp libdiceforge.a -o example.out
```
where *libdiceforge.a* is the built static library

3. If you did not install diceforge and you did not include diceforge to PATH 
```bash
g++ example.cpp libdiceforge.a -I "<include>" -o example.out
```
where *\<include\>* is the path to the folder corresponding to [include](include) and *libdiceforge.a* is the built static library

Feel free to explore the library and experiment with different algorithms!

## Contributing

Currently, this is an IIT-M Math Club exclusive project and we aren't accepting outside contributions.

### For Project Members
- View the [Tasks folder](Contributors_Only/Tasks) to view the tasks to be done.

- Reference materials are provided in the [Contributors_Only](Contributors_Only) folder.

Current TODO

Final pending tasks
- [ ] Proofreading and finalizing documentation
- [ ] Update README with latest results

Good to haves
- [ ] Get bigger primes for BBS

Completed
- [x] Bernoulli docuumentation
- [x] Custom pdf distribution
- [x] Geometric distribution
- [x] Fix bug in Naor Reingold RNG, also decide what to do when seed is zero and mention it in code documentation
- [x] Exception handling
- [x] 2d rv
- [x] Testing the RNGs
- [x] More examples on using DiceForge in Documentation
- [x] Non-linear transformation in LFSR
- [x] A defualt random singleton for convenience