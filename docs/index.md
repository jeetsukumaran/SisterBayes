# SisterBayes -- Simultaneous Divergence Time Testing Using the Joint Site Frequency Spectrum

*SisterBayes* is a package for simultaneous divergence time testing in an ABC (Approximate Bayesian Computation) framework, using the entire joint site frequency spectrum as summary statistic suite.
By using the joint site frequency spectrum as a summary statistic, we are maximizing the information available to the ABC framework, and thus maximizing power.

### Introduction

Given a collection of molecular sequence alignments sampled from multiple pairs of diverging sister taxa, *SisterBayes* answers estimates the posterior probabilities of series of *divergence models* that describe the timing of the divergences of each of these pairs relative to each other, from complete synchrony or simultaneous divergence (all pairs diverged together) to complete independence (all pairs diverged at different times) and everything in between.

*SisterBayes* thus builds upon previous programs such as:
-   [msBayes](http://msbayes.sourceforge.net/)
-   [dpp-Bayes](https://github.com/joaks1/dpp-msbayes)
-   [PyMsBayes](https://joaks1.github.io/PyMsBayes/)
Like the above programs, *SisterBayes* uses an [Approximate Bayesian Computation]() or ABC framework for inference, which involves reducing the observed data to summary statistics, simulating summary statistics from the prior under various different models, and using rejection sampling based on variously weighted distance measures between the observed and simulated summar statistics to obtain an approximation of samples from the posterior of the model/parameter space.

*SisterBayes* distinguishes itself from the above by incorporating the following features:

-   Instead of just a small handful of relatively highly-lossy summary statistics (e.g., Tajima's D, Pi, etc.), *SisterBayes* uses the *entire* Joint Summary Statistic Spectrum, thus maximizing the amount of information that could possibly be used.
-   *SisterBayes* allows the full use of multi-locus data (while *msBayes* supports multi-locus data by averaging summary statistics across loci, resulting in extreme and extremely undesirable jettisoning of power) as well as SNP's.
-   *SisterBayes* is written from scratch in modern object-oriented Python, natively supports multi-threading, and, in addition to incorporating all the excellent model enhancements in [PyMsBayes](https://joaks1.github.io/PyMsBayes/), also incorporates features from other approaches (e.g., the "Beta" buffer parameter) as well as some of its own (e.g., the post-analytic clustering heuristic).

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

