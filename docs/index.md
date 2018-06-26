# SisterBayes -- Simultaneous Divergence Time Testing Using the Joint Site Frequency Spectrum

*SisterBayes* is a package for simultaneous divergence time testing in an ABC (Approximate Bayesian Computation) framework, using the entire joint site frequency spectrum as summary statistic suite.
By using the joint site frequency spectrum as a summary statistic, we are maximizing the information available to the ABC framework, and thus maximizing power.

## Introduction

Given a collection of molecular sequence alignments sampled from multiple pairs of diverging sister taxa, *SisterBayes* answers estimates the posterior probabilities of series of *divergence models* that describe the timing of the divergences of each of these pairs relative to each other, from complete synchrony or simultaneous divergence (all pairs diverged together) to complete independence (all pairs diverged at different times) and everything in between.

*SisterBayes* thus builds upon previous programs such as:
-   [msBayes](http://msbayes.sourceforge.net/)
-   [dpp-Bayes](https://github.com/joaks1/dpp-msbayes)
-   [PyMsBayes](https://joaks1.github.io/PyMsBayes/)

Like the above programs, *SisterBayes* uses an [Approximate Bayesian Computation]() or ABC framework for inference, which involves reducing the observed data to summary statistics, simulating summary statistics from the prior under various different models, and using rejection sampling based on variously weighted distance measures between the observed and simulated summar statistics to obtain an approximation of samples from the posterior of the model/parameter space.

*SisterBayes* distinguishes itself from the above by incorporating the following features:

-   *SisterBayes* uses the *entire* Joint Summary Statistic Spectrum, thus maximizing the amount of information that could possibly be used, instead of just a small handful of relatively highly-lossy summary statistics (e.g., Tajima's D, Pi, etc.),
-   *SisterBayes* allows the full use of multi-locus data as well as SNP's; in contrast, *msBayes* only supports multi-locus data by averaging summary statistics across loci, resulting in extreme and extremely undesirable jettisoning of power.
-   *SisterBayes* supports multiple different types of prior distributions on $\tau$, the time of the root split between the two populations of each species pair: a gamma distribution, a uniform plus beta "buffer", as well as a simple uniform.
-   *SisterBayes* is written from scratch in modern object-oriented Python, natively supports multi-threading, and, in addition to incorporating all the excellent model enhancements in [PyMsBayes](https://joaks1.github.io/PyMsBayes/), also incorporates features from other approaches (e.g., the "Beta" buffer parameter) as well as some of its own (e.g., the post-analytic clustering heuristic).
-   *SisterBayes* provides special options useful in theoretical studies, including fixing divergence times, numbers of divergence classes, as well as output of raw alignments and trees.

## Installation

### Prerequisites

-   [Python](https://www.python.org/downloads/) (2.7 or greater)
-   [pip](https://pip.pypa.io/en/latest/installing/)
-   [DendroPy](https://www.dendropy.org/downloading.html) (4.3.0 or greater)
-   [FastSimCoal2](http://cmpg.unibe.ch/software/fastsimcoal2/)

### Installing

You can install directly from the main GitHub repository using:

```
$ pip install git+https://github.com/jeetsukumaran/SisterBayes.git
```

or

```
$ pip install git+git://github.com/jeetsukumaran/SisterBayes.git
```

If you do not have administrative privileges on the machine, you could use the ``--user`` flag to carry out local (user-directory) installation:

```
$ pip install --user git+https://github.com/jeetsukumaran/SisterBayes.git
```

or

```
$ pip install --user git+git://github.com/jeetsukumaran/SisterBayes.git
```

## Quick Start: Workflow Overview

A full *SisterBayes* analysis consists of four steps.
Each of these steps corresponds to a distinct program in the *SisterBayes* package.
The output of each step typically becomes the input for the subsequent step (with some steps combining the results of multiple steps).
Using *SisterBayes* can thus be seen as analogous to a pipeline, with data/information flowing from one program to another.
The four major steps and associated programs are:

1.  Generate samples from the prior using ``sisterbayes-simulate.py``.
2.  Calculate Joint Frequency Spectrum summary statistics from observed or empirical data using ``sisterbayes-sumstats.py``.
3.  Carry out rejection sampling on samples from the prior to approximate samples from the posterior based on observed data using ``sisterbayes-reject.py``.
4.  Summarize posterior using ``sisterbayes-summarize.py``.

## Usage


A basic understanding of the following concepts are required to effectively and sensibly use *SisterBayes*:

-   Approximate Bayesian Computation
-   The Simultaneous Divergence Time testing problem

<!--- A great overview of this tradition of analysis and programs are given in the [PyMsBayes](https://joaks1.github.io/PyMsBayes/) documentation, and we strongly encourage the user to refer to this site before proceeding if the user is not already familiar with it. -->

As noted in the "Quick Start" section above, a *SisterBayes* analysis consists of four major steps, each executed with a different program in the package:

1.  ``sisterbayes-simulate.py``
2.  ``sisterbayes-sumstats.py``
3.  ``sisterbayes-reject.py``
4.  ``sisterbayes-summarize.py``

All these programs (as well as others in the *SisterBayes* package) are command-line programs that need to be executed from the shell.
They all take a "--help" argument, that provides detailed listing of all the options, arguments, flags, and subcommands they support, as well as what they do.


### Simulation Data from the Prior: ``sisterbayes-simulate.py``

The first step in an ABC analysis is to draw samples from the prior distribution of the join space of the models using simulations.
The ``sisterbayes-simulate.py`` program achieves this.
You can see a full list of options and arguments from this program by running:

~~~
$ sisterbayes-simulate.py --help
~~~

#### Configuration File

At a minimum, you need to provide ``sisterbayes-simulate.py`` a *configuration* file, which describes in detail the simulation model and the priors for all its parameters.
*SisterBayes* supports two configuration file formats.
The first is a legacy format, which conforms to that of [PyMsBayes](https://joaks1.github.io/PyMsBayes/tutorials/configuration.html#config), with some extensions to support *SisterBayes*-specific extensions and features.
The second is a more modern (and robust) JSON format.

The legacy format is specifically designed to be 100% compatible with  [PyMsBayes](https://joaks1.github.io/PyMsBayes/tutorials/configuration.html#config) so as to allow PyMsBayes analytical setups to be re-used with little or no modification with *SisterBayes*.
A full description of the legacy format is given here: https://joaks1.github.io/PyMsBayes/tutorials/configuration.html#config.
Basically, it consists of a *plain text* file with two parts: a preamble and a sample table.
The preamble describes the model and the sample table describes the data.
An example of the legacy file format is:

~~~
concentrationShape = 1000.0
concentrationScale = 0.00437
thetaShape = 4.0
thetaScale = 0.001
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 000
tauShape = 1.0
tauScale = 0.02
timeInSubsPerSite = 1
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 0
migrationShape = 0
migrationScale = 0
numTauClasses = 0

BEGIN SAMPLE_TBL
species-1   locus-1 1.0     1.0     10      8       32.42   389     0.27    0.24    0.26    species-1-locus-1.fasta
species-1   locus-2 1.0     1.0     8       6       5.51    500     0.25    0.22    0.24    species-1-locus-2.fasta
species-1   locus-3 1.0     1.0     6       8       8.38    524     0.26    0.23    0.26    species-1-locus-3.fasta
species-1   locus-4 1.0     1.0     8       10      5.20    345     0.25    0.23    0.24    species-1-locus-4.fasta
species-1   locus-5 1.0     1.0     8       8       29.59   417     0.27    0.23    0.21    species-1-locus-5.fasta
species-1   mito-1  0.25    4.0     5       5       8.15    600     0.22    0.24    0.27    species-1-mito-1.fasta
species-2   locus-1 1.0     1.0     6       10      7.53    400     0.25    0.24    0.26    species-2-locus-1.fasta
species-2   locus-3 1.0     1.0     10      8       11.14   550     0.27    0.22    0.24    species-2-locus-3.fasta
species-2   locus-4 1.0     1.0     8       8       9.39    350     0.24    0.24    0.23    species-2-locus-4.fasta
species-2   locus-5 1.0     1.0     10      10      13.32   450     0.26    0.24    0.22    species-2-locus-5.fasta
species-2   mito-1  0.25    4.0     4       5       7.59    549     0.23    0.26    0.23    species-2-mito-1.fasta
species-3   locus-1 1.0     1.0     10      6       17.03   367     0.25    0.23    0.27    species-3-locus-1.fasta
species-3   locus-3 1.0     1.0     8       10      59.17   541     0.26    0.22    0.25    species-3-locus-3.fasta
species-3   locus-4 1.0     1.0     6       8       6.90    333     0.28    0.23    0.21    species-3-locus-4.fasta
species-3   mito-1  0.25    4.0     5       4       11.42   587     0.22    0.22    0.25    species-3-mito-1.fasta
END SAMPLE_TBL
~~~

The preamble section consists of the following key-value pairs:

-   concentrationShape
-   concentrationScale
-   thetaShape
-   thetaScale
-   ancestralThetaShape
-   ancestralThetaScale
-   thetaParameters
-   tauShape
-   tauScale
-   timeInSubsPerSite
-   bottleProportionShapeA
-   bottleProportionShapeB
-   bottleProportionShared
-   migrationShape
-   migrationScale
-   numTauClasses

