#! /usr/bin/env python

from __future__ import division
import math
import bisect
import csv
import os
import sys
import argparse
import collections
import re
from sisterbayes import utility
from sisterbayes import calclib

class SisterBayesSummarizer(object):

    def __init__(self,
            cluster_criteria=None,
            cluster_criteria_value=None,
            field_delimiter="\t",
            exclude_field_patterns=None,
            include_only_field_patterns=None,
            ):
        self.cluster_criteria = cluster_criteria
        self.cluster_criteria_value = cluster_criteria_value
        self.field_delimiter = field_delimiter
        self.all_fieldnames = None
        self.other_fieldnames = None
        self.stat_fieldnames = None
        self.stat_fieldnames_check = None
        self.other_fieldname_check = None
        self.stat_values = []
        self.other_values = []

    def cluster_by_absolute_difference_threshold(self,
            sp_labels,
            realized_div_time_samples,
            absolute_difference_threshold):
        model_counter = collections.Counter()
        for row_idx, realized_div_time_sample in enumerate(realized_div_time_samples):
            div_time_model_desc = [None for i in sp_labels]
            group_idx = 1
            for sp_idx, sp_label in enumerate(sp_labels):
                if sp_idx == 0:
                    div_time_model_desc[sp_idx] = str(group_idx)
                    group_idx += 1
                else:
                    current_dt = realized_div_time_sample[sp_label]
                    for prev_sp_idx, prev_sp_label in enumerate(sp_labels[:sp_idx-1]):
                        ref_dt = realized_div_time_sample[prev_sp_label]
                        if abs(current_dt - ref_dt) < (absolute_difference_threshold):
                            div_time_model_desc[sp_idx] = div_time_model_desc[prev_sp_idx]
                            break
                    else:
                        div_time_model_desc[sp_idx] = str(group_idx)
                        group_idx += 1
            model_name = "M"+"".join(div_time_model_desc)
            try:
                model_counter[model_name] += 1
            except KeyError:
                model_counter[model_name] = 1
        return model_counter

    # post-hoc derivation of div time model, clustering species that have
    # population splits +/- t time units from each other.
    def cluster_by_relative_difference_threshold(self,
            sp_labels,
            realized_div_time_samples,
            all_div_times,
            relative_difference_threshold):
        v0 = min(all_div_times)
        v1 = max(all_div_times)
        absolute_difference_threshold = (v1 - v0) * relative_difference_threshold
        return self.cluster_by_absolute_difference_threshold(
                sp_labels=sp_labels,
                realized_div_time_samples=realized_div_time_samples,
                absolute_difference_threshold=absolute_difference_threshold)

    def cluster_by_bin_size(self,
            sp_labels,
            realized_div_time_samples,
            all_div_times,
            bin_size):
        bin_size = float(bin_size)
        model_counter = collections.Counter()
        for row_idx, realized_div_time_sample in enumerate(realized_div_time_samples):
            div_time_model_desc = [None for i in sp_labels]
            group_idx = 1
            bin_idx_group_idx_map = {}
            for sp_idx, sp_label in enumerate(sp_labels):
                current_dt = realized_div_time_sample[sp_label]
                assigned_bin_idx = calclib.bin_index(current_dt, bin_size)
                try:
                    group_idx = bin_idx_group_idx_map[assigned_bin_idx]
                except KeyError:
                    group_idx = len(bin_idx_group_idx_map) + 1
                    bin_idx_group_idx_map[assigned_bin_idx] =  group_idx
                div_time_model_desc[sp_idx] = str(group_idx)
            model_name = "M"+"".join(div_time_model_desc)
            try:
                model_counter[model_name] += 1
            except KeyError:
                model_counter[model_name] = 1
        return model_counter

    def cluster_by_num_bins(self,
            sp_labels,
            realized_div_time_samples,
            all_div_times,
            num_bins):
        v0 = min(all_div_times)
        v1 = max(all_div_times)
        bin_size = (v1-v0)/float(num_bins)
        return self.cluster_by_bin_size(
            sp_labels=sp_labels,
            realized_div_time_samples=realized_div_time_samples,
            all_div_times=all_div_times,
            bin_size=bin_size)

    def summarize(self, target_data_filepath,):
        with utility.universal_open(target_data_filepath) as src:
            reader = csv.DictReader(
                    src,
                    delimiter=self.field_delimiter,
                    quoting=csv.QUOTE_NONE)
            categorical_params = collections.OrderedDict()
            continuous_params = collections.OrderedDict()
            realized_div_time_samples = []
            all_div_times = []
            sp_labels = []
            for row_idx, row in enumerate(reader):
                realized_div_time_samples.append({})
                for key_idx, key in enumerate(reader.fieldnames):
                    if key in categorical_params:
                        categorical_params[key][row[key]] += 1
                    elif key in continuous_params:
                        continuous_params[key].append(float(row[key]))
                    else:
                        if key in ("param.divTimeModel", "param.numDivTimes"):
                            val = row[key]
                            is_categorical = True
                        else:
                            try:
                                val = float(row[key])
                                is_categorical = False
                            except ValueError:
                                val = row[key]
                                is_categorical = True
                        if is_categorical:
                            categorical_params[key] = collections.Counter()
                            categorical_params[key][val] += 1
                        else:
                            continuous_params[key] = [val]
                    if key.startswith("param.divTime."):
                        sp_label = key.replace("param.divTime.", "")
                        realized_div_time_samples[-1][sp_label] = continuous_params[key][-1]
                        all_div_times.append(val)
                        if row_idx == 0:
                            sp_labels.append(sp_label)
            ### EXPERIMENTAL ###
            # categorical_params["param.effectiveDivTimeModel"] = self.cluster_by_relative_difference_threshold(
            #         sp_labels=sp_labels,
            #         realized_div_time_samples=realized_div_time_samples,
            #         all_div_times=all_div_times,
            #         relative_difference_threshold=0.01)
            if self.cluster_criteria is not None:
                if self.cluster_criteria == "bin_size":
                    cluster_results = self.cluster_by_bin_size(
                        sp_labels=sp_labels,
                        realized_div_time_samples=realized_div_time_samples,
                        all_div_times=all_div_times,
                        bin_size=self.cluster_criteria_value)
                elif self.cluster_criteria == "num_bins":
                    cluster_results = self.cluster_by_num_bins(
                        sp_labels=sp_labels,
                        realized_div_time_samples=realized_div_time_samples,
                        all_div_times=all_div_times,
                        num_bins=self.cluster_criteria_value)
                elif self.cluster_criteria == "absolute_difference_threshold":
                    cluster_results = self.cluster_by_absolute_difference_threshold(
                        sp_labels=sp_labels,
                        realized_div_time_samples=realized_div_time_samples,
                        absolute_difference_threshold=self.cluster_criteria_value)
                elif self.cluster_criteria == "relative_difference_threshold":
                    cluster_results = self.cluster_by_relative_difference_threshold(
                        sp_labels=sp_labels,
                        realized_div_time_samples=realized_div_time_samples,
                        all_div_times=all_div_times,
                        relative_difference_threshold=self.cluster_criteria_value)
                else:
                    raise ValueError("Unrecognized cluster criteria: '{}'".format(self.cluster_criteria))
                categorical_params["param.effectiveDivTimeModel"] = cluster_results
            ### EXPERIMENTAL ###

            output_prefix = os.path.splitext(os.path.basename(target_data_filepath))[0]
            with utility.universal_open(output_prefix + ".summary.continuous.tsv", "w") as dest:
                row_results = collections.OrderedDict()
                for param_idx, param_name in enumerate(continuous_params):
                    values = continuous_params[param_name]
                    row_results["param"] = param_name
                    summary = calclib.summarize(values)
                    row_results["mean"] = summary["mean"]
                    row_results["var"] = summary["var"]
                    row_results["sd"] = summary["sd"]
                    row_results["min"] = summary["range"][0]
                    row_results["max"] = summary["range"][1]
                    row_results["hpd5"] = summary["hpd95"][0]
                    row_results["hpd95"] = summary["hpd95"][1]
                    try:
                        row_results["quant5"] = summary["quant_5_95"][0]
                        row_results["quant95"] = summary["quant_5_95"][1]
                    except TypeError:
                        row_results["quant5"] = "NA"
                        row_results["quant95"] = "NA"
                    if param_idx == 0:
                        dest.write(self.field_delimiter.join(row_results.keys()) + "\n")
                    dest.write(self.field_delimiter.join("{}".format(v) for v in row_results.values()) + "\n")
            for param_idx, param_name in enumerate(categorical_params):
                with utility.universal_open(output_prefix + ".summary.{:02d}.{}.tsv".format(param_idx+1, param_name), "w") as dest:
                    param_counter = categorical_params[param_name]
                    total = float(sum(param_counter.values()))
                    for category_idx, (category_name, category_count) in enumerate(param_counter.most_common()):
                        row_results = collections.OrderedDict()
                        row_results["label"] = category_name
                        row_results["freq"] = category_count/total
                        row_results["count"] = category_count
                        if category_idx == 0:
                            dest.write(self.field_delimiter.join(row_results.keys()) + "\n")
                        dest.write(self.field_delimiter.join("{}".format(v) for v in row_results.values()) + "\n")

def main():
    parser = argparse.ArgumentParser(
            description="SISTERBAYES Summarizer",
            )
    parser.add_argument(
            "posteriors_filepath",
            help="Path to posteriors parameter file.")
    clustering_options = parser.add_argument_group("Clustering Options",
            "Calculate an 'effective' divergence time model based on clustering lineage pairs into simultaneous diverging groups using actual divergence times rather than labeled generating model.")
    # clustering_options =  parser.add_mutually_exclusive_group(required=True)
    clustering_options.add_argument("--bin-size",
        type=float,
        default=None,
        help="Cluster using divergence times bins of the specified size.")
    clustering_options.add_argument("--num-bins",
        type=float,
        default=None,
        help="Cluster by splitting the range of divergence times uniformly into the specified number of bins.")
    clustering_options.add_argument("--absolute-difference-threshold",
        type=float,
        default=None,
        help="Cluster by grouping together divergence times less than the specified threshold.")
    clustering_options.add_argument("--relative-difference-threshold",
        type=float,
        default=None,
        help="Cluster by grouping together divergence times less than the specified threshold, expressed as a proportion of the range of divergence times.")
    processing_options = parser.add_argument_group("Processing Options")
    processing_options.add_argument("--field-delimiter",
        type=str,
        default="\t",
        help="Field delimiter (default: <TAB>).")
    run_options = parser.add_argument_group("Run Options")
    run_options.add_argument(
            "-q", "--quiet",
            action="store_true",
            help="Work silently.")
    args = parser.parse_args()

    num_selected_clustering_options = 0
    cluster_criteria = None
    cluster_criteria_value = None
    for a in ("bin_size", "num_bins", "absolute_difference_threshold", "relative_difference_threshold"):
        if getattr(args, a) is not None:
            num_selected_clustering_options += 1
            cluster_criteria = a
            cluster_criteria_value = getattr(args, a)
    if num_selected_clustering_options > 1:
        sys.exit("Multiple clustering strategies selected: plese select just one")

    summarizer = SisterBayesSummarizer(
            cluster_criteria=cluster_criteria,
            cluster_criteria_value=cluster_criteria_value,
            field_delimiter=args.field_delimiter,
            )
    summarizer.summarize(args.posteriors_filepath)

if __name__ == "__main__":
    main()




