#! /usr/bin/env python

##############################################################################
## Copyright (c) 2017 Jeet Sukumaran.
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * The names of its contributors may not be used to endorse or promote
##       products derived from this software without specific prior written
##       permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
## IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
## THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
## PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
## AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

import os
import collections
import random
from sisterbayes import utility
from sisterbayes import fsc2

def compose_concatenated_locus_label(lineage_pair):
    # concatenated_locus_label = "concat{}".format(len(lineage_pair.locus_definitions))
    concatenated_locus_label = "concat"
    return concatenated_locus_label

def weighted_choice(seq, weights, rng):
    """
    Selects an element out of seq, with probabilities of each element
    given by the list `weights` (which must be at least as long as the
    length of `seq` - 1).
    """
    if weights is None:
        weights = [1.0/len(seq) for count in range(len(seq))]
    else:
        weights = list(weights)
    if len(weights) < len(seq) - 1:
        raise Exception("Insufficient number of weights specified")
    sow = sum(weights)
    if len(weights) == len(seq) - 1:
        weights.append(1 - sow)
    return seq[weighted_index_choice(weights, sow, rng)]

def weighted_index_choice(weights, sum_of_weights, rng):
    """
    (From: http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/)
    The following is a simple function to implement weighted random choice in
    Python. Given a list of weights, it returns an index randomly, according
    to these weights [1].
    For example, given [2, 3, 5] it returns 0 (the index of the first element)
    with probability 0.2, 1 with probability 0.3 and 2 with probability 0.5.
    The weights need not sum up to anything in particular, and can actually be
    arbitrary Python floating point numbers.
    If we manage to sort the weights in descending order before passing them
    to weighted_choice_sub, it will run even faster, since the random call
    returns a uniformly distributed value and larger chunks of the total
    weight will be skipped in the beginning.
    """
    rnd = rng.uniform(0, 1) * sum_of_weights
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

def sample_partition(
        number_of_elements,
        scaling_parameter,
        rng,):
    groups = []
    element_ids = [i for i in range(number_of_elements)] # decouple actual element index values, from indexes used to run Dirichlet process
    rng.shuffle(element_ids)                             # ... thus allowing for us to ensure the "first" index is randomized
    a = scaling_parameter
    for i, element_id in enumerate(element_ids):
        probs = []
        n = i + 1
        if i == 0:
            groups.append([element_id])
            continue
        p_new = a/(a + n - 1.0)
        probs.append(p_new)
        for group in groups:
            p = len(group)/(a + n - 1.0)
            probs.append(p)
        assert abs(sum(probs) - 1.0) <= 1e-5
        selected_idx = weighted_index_choice(
                weights=probs,
                sum_of_weights=1.0,
                rng=rng)
        if selected_idx == 0:
            groups.append([element_id])
        else:
            groups[selected_idx-1].append(element_id)
    return groups

def compose_lineage_pair_label(lineage_pair_idx):
    return "spp{}".format(lineage_pair_idx)

def compose_deme_label(deme_idx):
    return "deme{}".format(deme_idx)
_DEME0_LABEL = compose_deme_label(0)
_DEME1_LABEL = compose_deme_label(1)
_ANCESTOR_DEME_LABEL = compose_deme_label("A")

class LocusDefinition(object):

    def __init__(self, locus_d):
        self.configure(locus_d)

    def configure(self, locus_d):
        # Doc/comments for parameters from, and following, PyMsBayes (Jamie Oaks; https://github.com/joaks1/PyMsBayes)
        # label for this locus
        self.locus_label = locus_d.pop("locus_label")
        # The number in this column is used to scale for differences in ploidy among loci
        # or for differences in generation-times among taxa. In our example configuration
        # file 1.0 is used for loci from a diploid nuclear genome, whereas 0.25 is used
        # for a region of the mitochondrial genome (because its haploid and maternally
        # inherited). However, if a taxon "species-3" had 1/4 the generation times of
        # the other two taxa, we would specify "1.0" for the third column for its
        # mitochondrial locus, and "4.0" for the third column for its nuclear loci.
        self.ploidy_factor = float(locus_d.pop("ploidy_factor"))
        # The number in this column is used to scale for differences in mutation rates
        # among taxa and/or loci.
        self.mutation_rate_factor = float(locus_d.pop("mutation_rate_factor"))
        # Number of genes/sequences/taxa from first population
        self.num_genes_deme0 = int(locus_d.pop("num_genes_deme0"))
        # Number of genes/sequences/taxa from second population
        self.num_genes_deme1 = int(locus_d.pop("num_genes_deme1"))
        # This is the transition/transversion rate ratio ("Kappa") of the HKY85
        # model of nucleotide substitution [4] for this alignment. NOTE: This is
        # the transition/transversion rate ratio, not the "count" ratio. I.e.,
        # Kappa = 1 is equal to the Jukes-Cantor model.
        self.ti_tv_rate_ratio = float(locus_d.pop("ti_tv_rate_ratio"))
        # Number of sites
        self.num_sites = int(locus_d.pop("num_sites"))
        # Equilibrium frequency of nucleotide
        self.freq_a = float(locus_d.pop("freq_a"))
        # Equilibrium frequency of nucleotide
        self.freq_c = float(locus_d.pop("freq_c"))
        # Equilibrium frequency of nucleotide
        self.freq_g = float(locus_d.pop("freq_g"))
        # Path to alignment file (optional)
        self.alignment_filepath = locus_d.pop("alignment_filepath", None)
        # Done!
        if locus_d:
            raise Exception("Unrecognized locus definition entries: {}".format(locus_d))

class LineagePair(object):

    def __init__(self, sister_pair_label):
        self.label = sister_pair_label
        self.locus_definitions = []

    def add_locus_definition(self, locus_d):
        locus = LocusDefinition(locus_d)
        self.locus_definitions.append(locus)
        return locus

class SisterBayesModel(object):

    def __init__(self, params_d, locus_info,):
        if locus_info is not None:
            self.configure_loci(locus_info) # do this first, so we know the number of sister pairs and loci before working on params
        if params_d is not None:
            self.configure_params(params_d)
        self.fsc2_params_adjustment_hack = 1E8

    def configure_loci(self, locus_info):
        self.label_to_lineage_pair_map = {}
        self.lineage_pairs = [] # list to maintain order for indexing during Dirichlet process partitioning
        self.lineage_pairs_loci_labels = {}
        for locus_d in locus_info:
            sister_pair_label = locus_d.pop("sister_pair_label")
            try:
                lineage_pair = self.label_to_lineage_pair_map[sister_pair_label]
            except KeyError:
                lineage_pair = LineagePair(
                        sister_pair_label=sister_pair_label)
                self.label_to_lineage_pair_map[sister_pair_label] = lineage_pair
                self.lineage_pairs.append(lineage_pair)
                self.lineage_pairs_loci_labels[lineage_pair] = set()
            locus_definition = lineage_pair.add_locus_definition(locus_d)
            if locus_definition.locus_label in self.lineage_pairs_loci_labels[lineage_pair]:
                raise ValueError("Lineage pair '{}': locus with label '{}' has already been defined".format(lineage_pair.label, locus_definition.locus_label))
            self.lineage_pairs_loci_labels[lineage_pair].add(locus_definition.locus_label)

    def configure_params(self, params_d):
        # Doc/comments for parameters from, and following, PyMsBayes (Jamie Oaks; https://github.com/joaks1/PyMsBayes)
        params_d = utility.CaseInsensitiveDict(params_d)
        # Shape and scale of Gamma hyperprior on
        # concentration parameter of Dirichlet process to partition pairs
        if "concentrationShape" in params_d or "concentrationScale" in params_d:
            self.prior_concentration = (
                    float(params_d.pop("concentrationShape")),
                    float(params_d.pop("concentrationScale"))
                    )
        elif "numTauClasses" not in params_d:
            raise ValueError("Concentration shape ('concentrationShape') and scale ('concentrationScale') parameters missing and 'numTauClasses' or 'fixedDivTimeModel' not specified")

        # Special param for simulation/testing: fixed theta to a particular
        # value. Either '0' (ignore) or a three comma-separate values.
        # Values are for deme0, deme1, and ancestral deme, respectively.
        # Any value of '0' will be ignored (allowed to vary according to
        # conditions specified previously, e.g., 'thetaParameters', 'theta', or
        # 'ancestralTheta'.
        fixed_thetas_str = params_d.pop("fixedThetas", "")
        if fixed_thetas_str and fixed_thetas_str != "0":
            fdt = [float(i.strip()) for i in fixed_thetas_str.split(",")]
            if len(fdt) != 3:
                raise ValueError("Expecting 3 values for 'fixedThetas' entry but only found {}".format(len(fdt)))
            self.fixed_thetas = fdt
            if "thetaShape" in params_d or "thetaScale" in params_d:
                raise ValueError("Cannot specify theta prior shape or scale if 'fixedThetas' are specified")
        else:
            self.fixed_thetas = None
            # # Shape and scale of Gamma hyperprior on theta PyMsBayes does not
            # # independently model N and \mu, but FastsimCoal2 does.
            self.prior_theta = (
                    float(params_d.pop("thetaShape")),
                    float(params_d.pop("thetaScale"))
                    )

        # Shape and scale of Gamma hyperprior on population size.
        # PyMsBayes does not independently model N and \mu. Here we do.
        # BEAST* uses a gamma distribution with a mean 2\psi
        # and a shape of 2, with user specifying a (hyper-)prior on \psi.
        # Here, for now, we just use Gamma directly, with default shape
        # parameter of 2.
        # self.prior_popsize = (
        #         float(params_d.pop("popsizeShape", 2))
        #         float(params_d.pop("popsizeScale"))
        #         )
        # # Shape and scale of Gamma hyperprior on mutation rate.
        # # PyMsBayes does not independently model N and \mu. Here we do.
        # self.prior_mutRate = (
        #         float(params_d.pop("mutRateShape"))
        #         float(params_d.pop("mutRateScale"))
        #         )

        # Shape and scale of Gamma hyperprior on
        # theta (population) parameters for ancestral deme
        self.prior_ancestral_theta = (
                float(params_d.pop("ancestralThetaShape")),
                float(params_d.pop("ancestralThetaScale"))
                )
        # specification of fixed/free parameters
        self.theta_constraints = str(params_d.pop("thetaParameters", "000"))
        if len(self.theta_constraints) != 3:
            raise ValueError("Incorrectly specified 'thetaParameters' constraints: '{}'".format(self.theta_constraints))
        for idx, i in enumerate(self.theta_constraints):
            if i not in ["0", "1", "2"]:
                raise ValueError("Incorrectly specified 'thetaParameters' constraints: '{}'".format(self.theta_constraints))
        # Shape and scale of Gamma hyperprior on
        # divergence times
        self.prior_tau = (
                float(params_d.pop("tauShape", 0)),
                float(params_d.pop("tauScale", 0))
                )
        # Alternate time prior parameterization: uniform distribution + beta
        # parameter to control pulses
        pbb = params_d.pop("pulseBufferBeta", None)
        try:
            self.pulse_buffer_beta = float(pbb)
        except TypeError:
            self.pulse_buffer_beta = None
        self.tau_min = float(params_d.pop("tauMin", 0))
        self.tau_max = float(params_d.pop("tauMax", 0))
        if self.tau_min != 0 or self.tau_max != 0:
            if self.prior_tau[0] != 0 or self.prior_tau[1] != 0:
                raise ValueError("Cannot specify gamma-distributed prior and uniform prior on tau at the same time")
            if self.pulse_buffer_beta is None:
                self.tau_parameterization = "uniform"
            else:
                self.tau_parameterization = "uniform+beta"
            # self.tau_parameterization = "uniform"
        else:
            if self.prior_tau[0] == 0 or self.prior_tau[1] == 0:
                if "fixedTaus" not in params_d or params_d["fixedTaus"] == "0":
                    raise ValueError("Must specify either gamma-distributed prior or uniform prior on tau")
            self.tau_parameterization = "gamma"
        # Shape and scale of Gamma hyperprior on
        # divergence times
        self.prior_migration = (
                float(params_d.pop("migrationShape", 0)),
                float(params_d.pop("migrationScale", 0))
                )
        if self.prior_migration[0] != 0 or self.prior_migration[1] != 0:
            raise NotImplementedError("Migration is not yet supported")
        # 1: Time units are in expected substitutions per site. For example, a
        #    divergence of 0.05 means that, on average, 5% of sites have changed
        #    since the populations diverged (so you expect 10% divergence between
        #    the populations since the population divergence). Thus, you can
        #    convert these units to the number of generations or years by dividing
        #    by the mutation rate.
        tss = int(params_d.pop("timeInSubsPerSite", 1))
        self.time_in_subs_per_site = bool(tss)
        if not self.time_in_subs_per_site:
            raise NotImplementedError("Time not in expected substitutions per site not supported")
        # If both are positive, these settings define a beta prior on the magnitude of a
        # post-divergence bottleneck in each of the descendant populations.
        # bottleProportionShapeA and bottleProportionShapeB correspond to the shape
        # parameters alpha and beta, respectively, of the beta prior.
        # The bottleneck magnitude is the proportion of the effective population size
        # that remains during the bottleneck. For example, a value of 0.95 would mean
        # that bottleneck reduces the effective population size by 5%.
        # If either or both are zero or less, there is no post-divergence population
        # bottleneck in the descendant populations (i.e., the bottleneck-magnitude
        # parameters, along with the timing of each bottleneck, are removed from
        # the model).
        self.prior_bottleneck_proportion = (
                float(params_d.pop("bottleProportionShapeA", 0)),
                float(params_d.pop("bottleProportionShapeB", 0)),
                )
        # If bottleProportionShared = 0, then there are two free
        # bottleneck-magnitude parameters for each population pair (one for
        # each descendant population). If bottleProportionShared = 1, then
        # there is one bottleneck-magnitude parameter for each population pair
        # (i.e., the descendant populations of each pair share the same
        # bottleneck magnitude; the bottleneck magnitude still varies among the
        # pairs).
        self.bottle_proportion_shared = bool(int(params_d.pop("bottleProportionShared", 0)))
        # Fixed the divergence time model.
        self.fixed_divergence_time_model = params_d.pop("fixedDivTimeModel", None)
        if self.fixed_divergence_time_model is not None:
            original_fdtm = self.fixed_divergence_time_model
            if self.fixed_divergence_time_model.upper().startswith("M"):
                self.fixed_divergence_time_model = self.fixed_divergence_time_model[1:]
            if len(self.fixed_divergence_time_model) != self.num_lineage_pairs:
                raise ValueError("Expecting {} elements in divergence time model description, but only found {}: '{}'".format(
                    self.num_lineage_pairs,
                    len(self.fixed_divergence_time_model),
                    original_fdtm))
        # If this setting is zero (the default), the number of divergence
        # events is free to vary according to the Dirichlet process prior on
        # divergence models. If it is greater than zero, then the model is
        # constrained to numTauClasses divergence events. This is useful for
        # simulation-based power analyses, but should not be used for empirical
        # analyses.
        self.num_tau_classes = int(params_d.pop("numTauClasses", 0))
        if self.num_tau_classes > 0 and self.fixed_divergence_time_model is not None:
            num_grps = len(set(self.fixed_divergence_time_model))
            if num_grps != self.num_tau_classes:
                raise ValueError("Fixed divergence time model 'M{}' requires {} 'numTauClasses', but {} specified".format(
                    self.fixed_divergence_time_model, num_grps, self.num_tau_classes))
        # Special param for simulation/testing: fixed tau to a set of
        # particular values. Either '0' (ignored) or a comma-separated list of
        # values. For each potential cluster (= number of species pairs if
        # 'numTauClasses' is not specified or 0, or 'numTauClasses' otherwise),
        # there needs to be a divergence time given. These times will be used
        # in order (e.g., the first cluster of species pairs will get the first
        # divergence time, the second the next, etc.).
        fixed_divergence_times_str = params_d.pop("fixedTaus", "")
        if fixed_divergence_times_str and fixed_divergence_times_str != "0":
            fdt = [float(i.strip()) for i in fixed_divergence_times_str.split(",")]
            if self.num_tau_classes == 0:
                self.num_tau_classes = len(fdt)
            if len(fdt) != self.num_lineage_pairs and (self.num_tau_classes > 0 and len(fdt) != self.num_tau_classes):
                raise ValueError("Expecting {} values for 'fixedTaus' entry but only found {}".format(self.num_lineage_pairs, len(fdt)))
            self.fixed_divergence_times = fdt
        else:
            self.fixed_divergence_times = None
        # Done!
        if params_d:
            raise ValueError("Unrecognized parameter configuration entries: {}".format(params_d))

    def _get_num_lineage_pairs(self):
        return len(self.lineage_pairs)
    num_lineage_pairs = property(_get_num_lineage_pairs)

    def sample_parameter_values_from_prior(self, rng):
        params = collections.OrderedDict()

        ## div time
        params["param.divTimeModel"] = "NA" # initialize here, so first column
        if self.fixed_divergence_time_model:
            num_grps = len(set(self.fixed_divergence_time_model))
            groups = [[] for idx in range(num_grps)]
            group_id_group_idx_map = {}
            assert len(self.fixed_divergence_time_model) == self.num_lineage_pairs
            for group_id, element_id in zip( self.fixed_divergence_time_model, range(self.num_lineage_pairs) ):
                try:
                    group_idx = group_id_group_idx_map[group_id]
                except KeyError:
                    group_idx = len(group_id_group_idx_map)
                    group_id_group_idx_map[group_id] = group_idx
                groups[group_idx].append(element_id)
            assert len(groups) == num_grps
            assert len(groups) == len(group_id_group_idx_map)
        elif self.num_tau_classes:
            if self.num_tau_classes >= self.num_lineage_pairs:
                groups = [[idx] for idx in range(self.num_lineage_pairs)]
            else:
                element_ids = [i for i in range(self.num_lineage_pairs)] # decouple actual element index values, from indexes used to run Dirichlet process
                rng.shuffle(element_ids)
                groups = [[] for idx in range(self.num_tau_classes)]
                for group in groups:
                    group.append(element_ids.pop(0))
                for element_id in element_ids:
                    rng.choice(groups).append(element_id)
        else:
            concentration_v = rng.gammavariate(*self.prior_concentration)
            # params["param.concentration"] = concentration_v
            groups = sample_partition(
                    number_of_elements=self.num_lineage_pairs,
                    scaling_parameter=concentration_v, # sample from prior
                    rng=rng,
                    )
        params["param.numDivTimes"] = len(groups)
        if self.fixed_divergence_times:
            div_time_values = [self.fixed_divergence_times[i] for i in range(len(groups))]
        elif self.tau_parameterization == "gamma":
            div_time_values = [rng.gammavariate(*self.prior_tau) for i in groups]
        elif self.tau_parameterization == "uniform":
            div_time_values = [rng.uniform(self.tau_min, self.tau_max) for i in groups]
        elif self.tau_parameterization == "uniform+beta":
            div_time_values = []
            available_time_ranges = [ [self.tau_min, self.tau_max] ]
            available_time_range_sizes = [ self.tau_max-self.tau_min ]
            for group_idx, group in enumerate(groups):
                assert len(available_time_ranges) == len(available_time_range_sizes)
                # print("---")
                # print(div_time_values)
                # print("ranges: {}".format(available_time_ranges))
                # print("sizes: {}".format(available_time_range_sizes))
                time_range = weighted_choice(
                        seq=available_time_ranges,
                        weights=available_time_range_sizes,
                        rng=rng)
                selected_div_time = rng.uniform(time_range[0], time_range[1])
                # print(selected_div_time)
                div_time_values.append(selected_div_time)
                new_available_time_ranges = []
                new_available_time_range_sizes = []
                prev_break = self.tau_min
                for dt in div_time_values:
                    t0 = dt - self.pulse_buffer_beta
                    if t0 <= prev_break:
                        continue
                    new_available_time_ranges.append( [prev_break, t0] )
                    new_available_time_range_sizes.append(t0 - prev_break)
                    prev_break += self.pulse_buffer_beta
                    if prev_break >= self.tau_max:
                        break
                if new_available_time_ranges:
                    available_time_ranges = new_available_time_ranges
                    available_time_range_sizes = new_available_time_range_sizes
                else:
                    raise ValueError("Unable to sample divergence times: pulse buffer parameter value too large relative to time range")
            # print("\n\n\n---\n{}\n---\n\n\n".format(div_time_values))
        else:
            raise ValueError("Unsupported tau parameterization: '{}'".format(self.tau_parameterization))
        fsc2_run_configurations = collections.OrderedDict()
        div_time_model_desc = [None for i in range(self.num_lineage_pairs)]

        expected_lineage_pair_idxs = set([i for i in range(self.num_lineage_pairs)])
        lineage_pair_div_times = {}
        # groups sorted by earliest occurring lineage pair index, to retain consistent div time coding
        for group_id, group in enumerate(sorted(groups, key=lambda group: min(lineage_pair_idx for lineage_pair_idx in group))):
            for lineage_pair_idx in group:
                assert lineage_pair_idx in expected_lineage_pair_idxs
                assert lineage_pair_idx not in fsc2_run_configurations
                lineage_pair = self.lineage_pairs[lineage_pair_idx]
                ## divergence time
                div_time_model_desc[lineage_pair_idx] = str(group_id+1) # divergence time model description
                div_time = div_time_values[group_id]
                lineage_pair_div_times[lineage_pair] = div_time
        del group_id
        del group
        del div_time
        for lineage_pair in self.lineage_pairs:
            div_time = lineage_pair_div_times[lineage_pair]
            params["param.divTime.{}".format(lineage_pair.label)] = div_time
            ## population parameters --- theta parameterization
            if self.fixed_thetas and self.fixed_thetas[0] > 0:
                deme0_theta = self.fixed_thetas[0]
            else:
                deme0_theta = rng.gammavariate(*self.prior_theta)
            if self.fixed_thetas and self.fixed_thetas[1] > 0:
                deme1_theta = self.fixed_thetas[1]
            elif self.theta_constraints[1] == self.theta_constraints[0]:
                deme1_theta = deme0_theta
            else:
                deme1_theta = rng.gammavariate(*self.prior_theta)
            if self.fixed_thetas and self.fixed_thetas[2] > 0:
                deme2_theta = self.fixed_thetas[2]
            elif self.theta_constraints[2] == self.theta_constraints[0]:
                deme2_theta = deme0_theta
            elif self.theta_constraints[2] == self.theta_constraints[1]:
                deme2_theta = deme1_theta
            elif self.prior_ancestral_theta[0] != 0 and self.prior_ancestral_theta[1] != 0:
                deme2_theta = rng.gammavariate(*self.prior_ancestral_theta)
            else:
                deme2_theta = rng.gammavariate(*self.prior_theta)
            params["param.theta.{}.{}".format(lineage_pair.label, _DEME0_LABEL)] = deme0_theta
            params["param.theta.{}.{}".format(lineage_pair.label, _DEME1_LABEL)] = deme1_theta
            params["param.theta.{}.{}".format(lineage_pair.label, _ANCESTOR_DEME_LABEL)] = deme2_theta

            for locus_id, locus_definition in enumerate(lineage_pair.locus_definitions):
                #
                #   Fastsimecoal2 separates pop size and mutation rate, but
                #   the msBayes/PyMsBayes model does not separate the two,
                #   using theta.
                #
                #   We could just reparameterize the PyMsBayes model here,
                #   sampling over N and mu independently. But say we want to
                #   stick to the theta parameterization.
                #
                #   ## Scheme 1
                #
                #   We could simply scale everything by mutation rate --
                #   i.e., population size and div time specified in units of
                #   N * mu, and in the sequence generation assume a base
                #   mutation rate of 1.0. Problem with this is that
                #   Fastsimcoal coerces the population size to an integer,
                #   and so anything less than 1 becomes zero. So we multiply
                #   the population size time by a large number and adjust
                #   this in the actual mutation rate so N mu remains the
                #   same:
                #
                #     theta = 4 N mu = 4 * (N * C) * (mu/C)
                #
                #   Of course, VERY important to also apply the adjustment
                #   factor to the divergence time, or, indeed, any other time
                #   variable!
                #
                #   In detail:
                #   - We "know" theta (because that is how we have parameterized our model and we are sampling theta from the prior)
                #   - We do not know N and mu -- these are part of theta:
                #         theta = 4 N mu
                #   - For fastsimcoal, we need to know N and mu separately
                #   - We scale everything to a real mutation rate, so *our* mu is an arbitrary 1.0
                #   - For the fastsimcoal population size, we divide theta by (4 * mu) = 4 to convert
                #         fsc_N = (theta / (4 * 1.0) * C
                #     where C is the adjustment factor so that fsc_N > 1.
                #   - Our tau is expressed in the expected number of substitutions per site =
                #         = mutation rate (expected number of substitutions per site) * abs_time
                #   - fastsimcoal tau is expressed in generations, so we divide our tau by mu to convert
                #     - fsc_tau = (tau / 1.0) * C
                #
                mu_factor = (locus_definition.mutation_rate_factor * locus_definition.ploidy_factor)

                fsc2_config_d = {
                    "d0_population_size": (deme0_theta/4.0 * mu_factor) * self.fsc2_params_adjustment_hack,
                    "d1_population_size": (deme1_theta/4.0 * mu_factor) * self.fsc2_params_adjustment_hack,
                    "d0_sample_size": locus_definition.num_genes_deme0,
                    "d1_sample_size": locus_definition.num_genes_deme1,
                    "div_time": div_time * self.fsc2_params_adjustment_hack, # ditto
                    "num_sites": locus_definition.num_sites,
                    "recombination_rate": 0,
                    "mutation_rate":  mu_factor / self.fsc2_params_adjustment_hack,
                    "ti_proportional_bias": (1.0 * locus_definition.ti_tv_rate_ratio)/3.0,
                    }
                fsc2_run_configurations[locus_definition] = fsc2_config_d
        params["param.divTimeModel"] = "M{}".format("".join(div_time_model_desc))
        return params, fsc2_run_configurations

