#! /usr/bin/env python

import csv
import os
import sys
import argparse
import collections
from sisterbayes import utility

def main():
    parser = argparse.ArgumentParser(
            description="SISTERBAYES Divergence Time Model Extractor",
            )
    parser.add_argument(
            "filepaths",
            nargs='+',
            help="Path to files.")
    processing_options = parser.add_argument_group("Processing Options")
    processing_options.add_argument(
            "-t", "--truncate-field-names",
            action="store_true",
            default=False,
            help="Truncate field names.")
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
    truncate_table = {}
    results = []
    out = sys.stdout
    for file_idx, filepath in enumerate(args.filepaths):
        src = open(filepath)
        result = collections.OrderedDict()
        if not args.quiet:
            sys.stderr.write("-- {}\n".format(filepath))
        reader = csv.DictReader(
                src,
                delimiter=args.field_delimiter,
                quoting=csv.QUOTE_NONE)
        for row_idx, row in enumerate(reader):
            for key_idx, key in enumerate(reader.fieldnames):
                normalized_case_key = key.lower()
                if not(
                        normalized_case_key.startswith("param.divtime")
                        or normalized_case_key.startswith("param.numdivtimes") ):
                    continue
                if args.truncate_field_names:
                    try:
                        truncated_key = truncate_table[key]
                    except KeyError:
                        truncated_key = key.replace("param.", "").replace("divTime.", "").replace("numDivTimes", "psi")
                        truncate_table[key] = truncated_key
                    result[truncated_key] = row[key]
                else:
                    result[key] = row[key]
            result["source"] = filepath
            results.append(result)
    utility.write_dict_csv(
            dest=out,
            list_of_dicts=results,
            delimiter=args.field_delimiter,
            is_no_header_row=False)

if __name__ == "__main__":
    main()





