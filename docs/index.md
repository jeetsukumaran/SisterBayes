# SisterBayes -- Simultaneous Divergence Time Testing Using the Joint Site Frequency Spectrum

*SisterBayes* is a package for simultaneous divergence time testing in an ABC (Approximate Bayesian Computation) framework, using the entire joint site frequency spectrum as summary statistic suite.
By using the joint site frequency spectrum as a summary statistic, we are maximizing the information available to the ABC framework, and thus maximizing power.

## Getting Started

### Prerequisites

-   [Python](https://www.python.org/downloads/) (2.7 or greater)
-   [pip](https://pip.pypa.io/en/latest/installing/)
-   [DendroPy](https://www.dendropy.org/downloading.html) (4.3.0 or greater)

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

