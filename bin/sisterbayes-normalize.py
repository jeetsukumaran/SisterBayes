#! /usr/bin/env python

import math
import csv
import os
import sys
import argparse
import collections
from sisterbayes import utility

class FileInfo(object):

    def __init__(self, filepath):
        self.filepath = filepath
        self.fieldnames = None
        self.data_row_idx_range = {}

class SisterBayesNormalizer(object):

    def __init__(self,
            compose_output_path_f,
            run_logger,
            stats_field_prefix="stat",
            field_delimiter="\t",
            logging_frequency=1000,
            missing_data_value="NA",
            ):
        self.compose_output_path_f = compose_output_path_f
        self.run_logger = run_logger
        self.stats_field_prefix = stats_field_prefix
        self.logging_frequency = logging_frequency
        self.field_delimiter = field_delimiter
        self.missing_data_value = missing_data_value
        self.fields = collections.OrderedDict()
        self.stat_field_names = collections.OrderedDict()
        self.file_infos = []
        self.current_data_row_idx = 0

    def _read_file(self, src):
        reader = csv.DictReader(
                src,
                delimiter=self.field_delimiter,
                quoting=csv.QUOTE_NONE)
        file_info = FileInfo(filepath=src.name)
        self.file_infos.append(file_info)
        file_info.fieldnames = reader.fieldnames
        file_data_row_start = self.current_data_row_idx
        for row_idx, row in enumerate(reader):
            if self.logging_frequency and row_idx > 0 and (row_idx % self.logging_frequency) == 0:
                self.run_logger.info("- Reading row {}".format(row_idx+1))
            assert len(row) == len(reader.fieldnames)
            for field_idx, field_name in enumerate(row):
                if field_name not in self.fields:
                    self.fields[field_name] = [None for idx in range(self.current_data_row_idx)]
                if field_name.startswith(self.stats_field_prefix):
                    self.stat_field_names[field_name] = 1
                    self.fields[field_name].append(float(row[field_name]))
                else:
                    self.fields[field_name].append(row[field_name])
            self.current_data_row_idx += 1
        file_data_row_end = self.current_data_row_idx
        file_info.data_row_idx_range = (file_data_row_start, file_data_row_end)

    def read_files(self, filepaths):
        for file_idx, filepath in enumerate(filepaths):
            self.run_logger.info("Reading file {} of {}: '{}'".format(file_idx+1, len(filepaths), filepath))
            with utility.universal_open(filepath) as src:
                self._read_file(src)

    def normalize(self):
        for field_name in self.stat_field_names:
            self.run_logger.info("- Normalizing field '{}'".format(field_name))
            values = [v for v in self.fields[field_name] if v is not None]
            min_value, max_value = min(values), max(values)
            normalization_factor = max_value - min_value
            if normalization_factor == 0:
                continue
            for vidx, v in enumerate(self.fields[field_name]):
                if v is not None:
                    self.fields[field_name][vidx] = (v - min_value)/normalization_factor

    def write_results(self):
        for file_idx, file_info in enumerate(self.file_infos):
            output_filepath = self.compose_output_path_f(file_info.filepath, file_idx)
            self.run_logger.info("Writing file {} of {}: '{}'".format(file_idx+1, len(self.file_infos), output_filepath))
            with utility.universal_open(output_filepath, "w") as dest:
                writer = utility.get_csv_writer(dest=dest,
                        fieldnames=file_info.fieldnames,
                        delimiter=self.field_delimiter,
                        restval=self.missing_data_value,
                        )
                writer.writeheader()
                for data_row_idx in range(*file_info.data_row_idx_range):
                    if self.logging_frequency and data_row_idx > 0 and (data_row_idx % self.logging_frequency) == 0:
                        self.run_logger.info("- Writing row {}".format(data_row_idx+1))
                    row = {}
                    for field_name in file_info.fieldnames:
                        row[field_name] = self.fields[field_name][data_row_idx]
                    writer.writerow(row)

    def process_files(self, filepaths):
        self.run_logger.info("Reading data ...")
        self.read_files(filepaths)
        self.run_logger.info("Normalizing ...")
        self.normalize()
        self.run_logger.info("Writing data ...")
        self.write_results()

def main():
    parser = argparse.ArgumentParser(
            description="SISTERBAYES Normalizer",
            )
    parser.add_argument(
            "data_filepaths",
            nargs="+",
            help="Path to data to be normalized.")
    parser.add_argument("--field-delimiter",
        type=str,
        default="\t",
        help="Field delimiter (default: <TAB>).")
    parser.add_argument("--stats-field-prefix",
        type=str,
        default="stat",
        help="Prefix identifying summary statistic fields to be normalized (default: '%(default)s').")
    parser.add_argument(
            "-q", "--quiet",
            action="store_true",
            help="Work silently.")
    args = parser.parse_args()
    run_logger = utility.RunLogger(
            name="sisterbayes-estimate",
            stderr_logging_level="info",
            log_to_stderr=not args.quiet,
            log_to_file=False,
            )
    compose_output_path_f = lambda x, row_idx: os.path.join(os.path.splitext(x)[0] + ".normalized.tsv")
    normalizer = SisterBayesNormalizer(
            compose_output_path_f=compose_output_path_f,
            run_logger=run_logger,
            stats_field_prefix=args.stats_field_prefix,
            field_delimiter=args.field_delimiter,
            )
    normalizer.process_files(filepaths=args.data_filepaths,)

if __name__ == "__main__":
    main()




