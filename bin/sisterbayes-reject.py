#! /usr/bin/env python

import math
import csv
import os
import sys
import argparse
import heapq
import collections
import sisterbayes
from sisterbayes import utility

class SisterBayesRejector(object):

    def __init__(self,
            rejection_criteria_type,
            rejection_criteria_value,
            run_logger,
            stats_field_prefix="stat",
            logging_frequency=1000,
            field_delimiter="\t",
            is_write_summary_stats=False,
            is_write_rejection_score=False,
            is_output_target_params=False,
            is_suppress_checks=False,
            ):
        self.rejection_criteria_type = rejection_criteria_type
        self.rejection_criteria_value = rejection_criteria_value
        self.run_logger = run_logger
        self.stats_field_prefix = stats_field_prefix
        self.logging_frequency = logging_frequency
        self.field_delimiter = field_delimiter
        self.is_write_summary_stats = is_write_summary_stats
        self.is_write_rejection_score = is_write_rejection_score
        self.is_output_target_params = is_output_target_params
        self.is_suppress_checks = is_suppress_checks
        self.stat_fieldnames = None
        self.stat_fieldnames_check = None
        self.non_stat_fieldnames = None
        self.distance_score_fieldname = "rejection.score"

    def euclidean_distance(self, vector1, vector2):
        assert len(vector1) == len(vector2)
        dist = [(a - b)**2 for a, b in zip(vector1, vector2)]
        dist = math.sqrt(sum(dist))
        return dist

    def extract_stat_fieldnames(self, fieldnames):
        stat_fieldnames = []
        non_stat_fieldnames = []
        for fieldname in fieldnames:
            if fieldname.startswith(self.stats_field_prefix):
                stat_fieldnames.append(fieldname)
            else:
                non_stat_fieldnames.append(fieldname)
        assert len(stat_fieldnames) == len(set(stat_fieldnames))
        assert len(non_stat_fieldnames) == len(set(non_stat_fieldnames))
        return stat_fieldnames, non_stat_fieldnames

    def extract_stats_data_vector_from_csv_row(self, row):
        data_vector = [float(v) for v in (row[k] for k in self.stat_fieldnames)]
        return data_vector

    def process(self,
            target_data_filepath,
            priors_data_filepaths,
            output_prefix,
            output_suffix):
        if output_prefix is None:
            output_prefix = os.path.splitext(os.path.basename(target_data_filepath))[0]
        if output_suffix is None:
            output_suffix = ""
        else:
            output_suffix = "." + output_suffix
        with utility.universal_open(target_data_filepath) as src:
            target_data_reader = csv.DictReader(
                    src,
                    delimiter=self.field_delimiter,
                    quoting=csv.QUOTE_NONE)
            for target_row_idx, target_row in enumerate(target_data_reader):
                if target_row_idx == 0:
                    self.stat_fieldnames, self.non_stat_fieldnames = self.extract_stat_fieldnames(target_data_reader.fieldnames)
                    self.stat_fieldnames_set = set(self.stat_fieldnames)
                self.run_logger.info("Scoring target data {}".format(target_row_idx+1))
                target_data_vector = self.extract_stats_data_vector_from_csv_row(target_row)
                posteriors_filepath = "{}.posterior.{:03d}{}.tsv".format(output_prefix, target_row_idx+1, output_suffix)
                self.accept_reject(
                        target_data_vector=target_data_vector,
                        priors_data_filepaths=priors_data_filepaths,
                        output_filepath=posteriors_filepath)
                if self.is_output_target_params:
                    target_params_filepath = "{}.posterior.{:03d}.target{}.tsv".format(output_prefix, target_row_idx+1, output_suffix)
                    with open(target_params_filepath, "w") as target_params_f:
                        target_params_f.write(self.field_delimiter.join(self.non_stat_fieldnames))
                        target_params_f.write("\n")
                        target_params_f.write(self.field_delimiter.join(str(target_row[k]) for k in self.non_stat_fieldnames))
                        target_params_f.write("\n")

    def accept_reject(self,
            target_data_vector,
            priors_data_filepaths,
            output_filepath):
        if self.rejection_criteria_type == "num":
            num_to_retain = self.rejection_criteria_value
        else:
            num_to_retain = None
        dest = utility.universal_open(output_filepath, "w")
        all_prior_fieldnames = []
        all_prior_fieldnames_set = None
        accepted_heap = []
        for fidx, priors_data_filepath in enumerate(priors_data_filepaths):
            self.run_logger.info("Reading simulation file {} of {}: '{}'".format(fidx+1, len(priors_data_filepaths), priors_data_filepath))
            with utility.universal_open(priors_data_filepath) as src:
                priors_data_reader = csv.DictReader(
                        src,
                        delimiter=self.field_delimiter,
                        quoting=csv.QUOTE_NONE)
                for row_idx, row in enumerate(priors_data_reader):
                    if self.logging_frequency and row_idx > 0 and row_idx % self.logging_frequency == 0:
                        self.run_logger.info("Reading simulation file {} of {}, row {}".format(fidx+1, len(priors_data_filepaths), row_idx+1))
                    if row_idx == 0:
                        if fidx == 0:
                            all_prior_fieldnames = list(priors_data_reader.fieldnames)
                            all_prior_fieldnames_set = set(all_prior_fieldnames)
                            current_file_stat_fieldnames = set(self.extract_stat_fieldnames(priors_data_reader.fieldnames)[0])
                            s1 = current_file_stat_fieldnames - self.stat_fieldnames_set
                            if s1:
                                raise ValueError("File '{}': Following summary statistics fields not found in target: {}".format(
                                    priors_data_filepath, ", ".join(s1)))
                            s2 =  self.stat_fieldnames_set - current_file_stat_fieldnames
                            if s2:
                                raise ValueError("File '{}': Following summary statistics fields given in target but not found here: {}".format(
                                    priors_data_filepath, ", ".join(s2)))
                            header_row = []
                            for fnidx, fn in enumerate(all_prior_fieldnames):
                                if self.is_write_summary_stats or fn not in self.stat_fieldnames_set:
                                    header_row.append(fn)
                            if self.is_write_rejection_score:
                                header_row.append(self.distance_score_fieldname)
                            dest.write("{}\n".format(self.field_delimiter.join(header_row)))
                        else:
                            current_file_fieldnames = set(priors_data_reader.fieldnames)
                            s1 = current_file_fieldnames - all_prior_fieldnames_set
                            if s1:
                                raise ValueError("File '{}': Following fields found, but not found in previous files: {}".format(
                                    priors_data_filepath, ", ".join(s1)))
                            s2 =  all_prior_fieldnames_set - current_file_fieldnames
                            if s2:
                                raise ValueError("File '{}': Following fields found in previous files, but not found here: {}".format(
                                    priors_data_filepath, ", ".join(s2)))
                    prior_data_vector = self.extract_stats_data_vector_from_csv_row(row)
                    distance_score = self.euclidean_distance(target_data_vector, prior_data_vector)
                    row_values = self.field_delimiter.join(row[fn] for fn in priors_data_reader.fieldnames if self.is_write_summary_stats or fn not in self.stat_fieldnames_set)
                    if self.is_write_rejection_score:
                        row_values = "{}{}{}".format(row_values, self.field_delimiter, distance_score)
                    heap_score = -1 * (distance_score)
                    heap_entry = (heap_score, row_values)
                    if self.rejection_criteria_type == "distance":
                        if distance_score <= self.rejection_criteria_value:
                            accepted_heap.append(heap_entry)
                    elif self.rejection_criteria_type == "num":
                        if len(accepted_heap) < num_to_retain:
                            accepted_heap.append(heap_entry)
                            if len(accepted_heap) == num_to_retain:
                                heapq.heapify(accepted_heap)
                        else:
                            heapq.heappushpop(accepted_heap, heap_entry)
                    else:
                        raise NotImplementedError(self.rejection_criteria_type)
                    # for fnidx, fn in enumerate(all_prior_fieldnames):
                    #     value = row[fn]
                    #     if self.is_write_summary_stats or fn not in self.stat_fieldnames_set:
                    #         dest.write("{}{}".format(value, self.field_delimiter))
                    # dest.write("{}\n".format(distance))
        accepted_heap.sort(reverse=True)
        for hidx, heap_entry in enumerate(accepted_heap):
            heap_entry = accepted_heap[hidx]
            dest.write(heap_entry[1])
            dest.write("\n")
        dest.flush()
        dest.close()

    # def accept_posteriors(self, distanced_scored_params_filepath, num_samples, output_filepath):
    #     dest = utility.universal_open(output_filepath, "w")
    #     if self.rejection_criteria_type == "num":
    #         num_to_retain = self.rejection_criteria_value
    #     elif self.rejection_criteria_type == "proportion":
    #         num_to_retain = int(self.rejection_criteria_value * num_samples)
    #     accepted_heap = []
    #     self.run_logger.info("Accepting/rejecting simulations from the prior ...")
    #     with utility.universal_open(distanced_scored_params_filepath) as src:
    #         priors_data_reader = csv.DictReader(
    #                 src,
    #                 delimiter=self.field_delimiter,
    #                 quoting=csv.QUOTE_NONE)
    #         for row_idx, row in enumerate(priors_data_reader):
    #             if self.logging_frequency and row_idx > 0 and row_idx % self.logging_frequency == 0:
    #                 self.run_logger.info("Accepting/rejecting: row {}".format(row_idx+1))
    #             if row_idx == 0:
    #                 dest.write(self.field_delimiter.join(priors_data_reader.fieldnames))
    #                 dest.write("\n")
    #             distance_score = float(row[self.distance_score_fieldname])
    #             row_values = self.field_delimiter.join(row[k] for k in priors_data_reader.fieldnames)
    #             if self.rejection_criteria_type == "distance":
    #                 if float(distance_score) <= self.rejection_criteria_value:
    #                     dest.write(row_values)
    #                     dest.write("\n")
    #             else:
    #                 # heap_score = (1.0/(distance_score + 1))
    #                 # heap_score = -1 * (distance_score + 1)
    #                 heap_score = -1 * (distance_score)
    #                 heap_entry = (heap_score, row_values)
    #                 if len(accepted_heap) < num_to_retain:
    #                     accepted_heap.append( heap_entry  )
    #                     if len(accepted_heap) == num_to_retain:
    #                         heapq.heapify(accepted_heap)
    #                 else:
    #                     heapq.heappushpop(accepted_heap, heap_entry)
    #     if self.rejection_criteria_type != "distance":
    #         accepted_heap.sort(reverse=True)
    #         for hidx, heap_entry in enumerate(accepted_heap):
    #             heap_entry = accepted_heap[hidx]
    #             dest.write(heap_entry[1])
    #             dest.write("\n")
    #     dest.flush()
    #     dest.close()

def main():
    package_id = sisterbayes.package_id()
    parser = argparse.ArgumentParser(
            description="SISTERBAYES Rejection Sampler",
            )
    parser.add_argument("--version", action="version", version=package_id)
    parser.add_argument(
            "target_data_filepath",
            help="Path to target or observed data file.")
    parser.add_argument(
            "simulations_data_filepaths",
            nargs="+",
            help="Path to samples from the prior data files.")
    rejection_criteria = parser.add_argument_group("Rejection Criteria")
    rejection_criteria.add_argument(
            "-n", "--max-num",
            type=int,
            metavar="#",
            default=None,
            help="Retain this number of samples from the prior into the posterior.")
    # rejection_criteria.add_argument(
    #         "-p", "--max-proportion",
    #         type=float,
    #         metavar="0.##",
    #         default=None,
    #         help="Retain this proportion (0 > 'p' > 1.0) of samples from the prior into the posterior.")
    rejection_criteria.add_argument(
            "-d", "--max-distance",
            type=float,
            metavar="#.##",
            default=None,
            help="Retain samples this distance or lower from the prior into the posterior.")
    processing_options = parser.add_argument_group("Processing Options")
    processing_options.add_argument("--field-delimiter",
        type=str,
        default="\t",
        help="Field delimiter (default: <TAB>).")
    processing_options.add_argument("--stats-field-prefix",
        type=str,
        default="stat",
        help="Prefix identifying summary statistic fields (default: '%(default)s').")
    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument('-o', '--output-prefix',
            action='store',
            type=str,
            default=None,
            metavar='NAME-PREFIX',
            help="Prefix for output filename.")
    output_options.add_argument('-O', '--output-suffix',
            action='store',
            type=str,
            default=None,
            metavar='NAME-SUFFIX',
            help="Suffix for output filename.")
    output_options.add_argument(
            "--write-summary-stats",
            action="store_true",
            help="Include summary stats in the samples from the posterior.")
    output_options.add_argument(
            "--write-rejection-score",
            action="store_true",
            help="Include rejection score in the output.")
    output_options.add_argument(
            "--output-target-params",
            action="store_true",
            help="For each row in the target data file processed, output a file of non-summary stats fields found.")
    run_options = parser.add_argument_group("Run Options")
    # run_options.add_argument(
    #         "-L", "--large-file",
    #         dest="limit_memory",
    #         action="store_true",
    #         default=False,
    #         help="Use two-pass processing that reduces memory footprint by not requiring entire simulations/priors file(s) to be read into memory at once.")
    run_options.add_argument(
            "-q", "--quiet",
            action="store_true",
            help="Work silently.")
    run_options.add_argument('--log-to-file',
            action='store_true',
            dest='log_to_file',
            default=None,
            help="Save log to file.")
    args = parser.parse_args()
    num_non_Nones = sum([1 for i in (args.max_num, args.max_distance) if i is not None])
    if num_non_Nones == 0:
        sys.exit("Require exactly one of '-n'/'--max-num', or '-d'/'--max-distance' to be specified.")
    elif num_non_Nones > 1:
        sys.exit("Require only one of '-n'/'--max-num', or '-d'/'--max-distance' to be specified.")
    if args.max_num:
        rejection_criteria_type = "num"
        rejection_criteria_value = args.max_num
    elif args.max_proportion:
        rejection_criteria_type = "proportion"
        rejection_criteria_value = args.max_proportion
    elif args.max_distance:
        rejection_criteria_type = "distance"
        rejection_criteria_value = args.max_distance
    run_logger = utility.RunLogger(
            name="sisterbayes-estimate",
            stderr_logging_level="info",
            log_to_stderr=not args.quiet,
            log_to_file=args.log_to_file,
            )
    run_logger.info("Running: {}".format(package_id))
    rejector = SisterBayesRejector(
            rejection_criteria_type=rejection_criteria_type,
            rejection_criteria_value=rejection_criteria_value,
            run_logger=run_logger,
            stats_field_prefix=args.stats_field_prefix,
            field_delimiter=args.field_delimiter,
            is_write_summary_stats=args.write_summary_stats,
            is_write_rejection_score=args.write_rejection_score,
            is_output_target_params=args.output_target_params,
            )
    rejector.process(
            target_data_filepath=args.target_data_filepath,
            priors_data_filepaths=args.simulations_data_filepaths,
            output_prefix=args.output_prefix,
            output_suffix=args.output_suffix)

if __name__ == "__main__":
    main()


