#! /usr/bin/env python

import os
import sys
import argparse
import traceback
import time
import sisterbayes
from sisterbayes import simulate
from sisterbayes import utility

def main():
    parser = argparse.ArgumentParser()
    package_id = sisterbayes.package_id()
    parser.add_argument("--version", action="version", version=package_id)

    simulator_options = parser.add_argument_group("Simulation Configuration")
    simulator_options.add_argument("configuration_filepath",
            metavar="CONFIGURATION-FILE",
            help="Path to file defining the simulation model and parameters.")
    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument('-o', '--output-name-prefix',
            action='store',
            dest='output_name_prefix',
            type=str,
            default=None,
            metavar='NAME-PREFIX',
            help="Prefix for output filenames (default: same as configuration filename stem).")
    output_options.add_argument('-O', '--output-directory',
            action='store',
            dest='output_directory',
            type=str,
            default=None,
            metavar='DIRECTORY',
            help="Directory for output files (default: current working directory).")
    output_options.add_argument(
            "-U",
            "--unfolded-site-frequency-spectrum",
            "--derived-site-frequency-spectrum",
            action="store_true",
            default=False,
            help="Calculate the unfolded or derived site frequency spectrum."
            " Otherwise, defaults to the folded or minor site frequency"
            " spectrum."
            )
    output_options.add_argument(
            "--infinite-sites-model",
            action="store_true",
            default=False,
            help="Use infinite sites model instead of finite sites."
            )
    output_options.add_argument(
            "--calculate-single-population-site-frequency-spectrum",
            action="store_true",
            default=False,
            help="Calculate the single (within) population site frequency"
            " spectrum in addition to the joint."
            )
    output_options.add_argument("-l", "--labels",
            action="append",
            help="Addition field/value pairs to add to the output (in format <FIELD-NAME>:value;)")
    output_options.add_argument('--field-delimiter',
            type=str,
            default='\t',
            help="Delimiter string separating fields in output (default: <TAB>').")
    output_options.add_argument('--summary-stats-label-prefix',
            type=str,
            default='stat',
            metavar='PREFIX',
            help="Prefix for summar statistic field labels (default: '%(default)s').")
    output_options.add_argument( "--include-model-id-field",
            action="store_true",
            default=False,
            help="Include a 'model.id' field (with same value as 'param.divTimeModel' field) in output.")
    output_options.add_argument( "--append",
            action="store_true",
            default=False,
            help="Append instead of overwriting output file(s).")
    output_options.add_argument( "--no-write-header",
            action="store_true",
            default=False,
            help="Do not writer header row.")

    run_options = parser.add_argument_group("Run Options")
    run_options.add_argument("-n", "--num-reps",
            type=int,
            default=1,
            help="Number of replicates (default: %(default)s).")
    run_options.add_argument("-m", "--num-processes",
            default=1,
            type=int,
            help="Number of processes/CPU to run (default: %(default)s).")
    run_options.add_argument("-z", "--random-seed",
            default=None,
            help="Seed for random number generator engine.")
    run_options.add_argument("--log-frequency",
            default=None,
            type=int,
            help="Frequency that background progress messages get written to the log (0: do not log informational messages).")
    run_options.add_argument("--file-logging-level",
            default="none",
            choices=["debug", "info", "warning", "error", "critical", "none", ],
            help="Message level threshold for screen logs (default: %(default)s).")
    run_options.add_argument("--stderr-logging-level",
            default="info",
            choices=["debug", "info", "warning", "error", "critical", "none", ],
            help="Message level threshold for screen logs (default: %(default)s).")
    run_options.add_argument('-w', '--working-directory-parent',
            action='store',
            type=str,
            default=None,
            help="Directory within which to create temporary directories and files.")
    run_options.add_argument("--no-cleanup",
            action="store_true",
            default=False,
            help="Do not clean-up temporary files.")
    run_options.add_argument("--debug-mode",
            action="store_true",
            default=False,
            help="Run in debugging mode.")

    fsc2_options = parser.add_argument_group("FastSimCoal2 Options")
    fsc2_options.add_argument("--fsc2-path",
            metavar="FSC2-PATH",
            default="fsc25",
            help="Path to FastsimCoal2 application (default: %(default)s).")

    args = parser.parse_args()

    config_d = {}
    if not os.path.exists(args.configuration_filepath):
        sys.exit("ERROR: Configuration file '{}' not found.".format(args.configuration_filepath))
    utility.parse_legacy_configuration(
            filepath=args.configuration_filepath,
            config_d=config_d)
    config_d["output_prefix"] = utility.output_prefix(
            primary_source_filepath=args.configuration_filepath,
            output_name_prefix=args.output_name_prefix,
            output_directory=args.output_directory)
    if args.log_frequency is None:
        config_d["logging_frequency"] = int(args.num_reps/10.0)
    elif args.log_frequency == 0:
        config_d["logging_frequency"] = None
    else:
        config_d["logging_frequency"] = args.log_frequency
    config_d["fsc2_path"] = args.fsc2_path
    if utility.which(config_d["fsc2_path"]) is None:
        sys.exit("ERROR: FastSimCoal2 executable '{}' not found.\n"
                    "Install FastSimCoal2 and specify path to the executable\n"
                    "using the '--fsc2-path' argument.".format(config_d["fsc2_path"]))
    config_d["file_logging_level"] = args.file_logging_level
    config_d["standard_error_logging_level"] = args.stderr_logging_level
    # config_d["log_to_file"] = args.log_to_file
    # config_d["log_to_stderr"] = args.log_to_stderr
    config_d["is_unfolded_site_frequency_spectrum"] = args.unfolded_site_frequency_spectrum
    config_d["is_calculate_single_population_sfs"] = args.calculate_single_population_site_frequency_spectrum
    config_d["is_calculate_joint_population_sfs"] = True
    config_d["is_infinite_sites_model"] = args.infinite_sites_model
    config_d["stat_label_prefix"] = args.summary_stats_label_prefix
    config_d["supplemental_labels"] = utility.parse_fieldname_and_value(args.labels)
    config_d["field_delimiter"] = args.field_delimiter
    config_d["is_include_model_id_field"] = args.include_model_id_field
    with utility.TemporaryDirectory(
            prefix="sisterbayes-",
            parent_dir=args.working_directory_parent,
            is_suppress_cleanup=args.no_cleanup) as working_directory:
        config_d["working_directory"] = working_directory
        simulator = simulate.SisterBayesSimulator(
                config_d=config_d,
                num_processes=args.num_processes,
                is_verbose_setup=True,
                package_id=package_id,
                )
        filepath = config_d["output_prefix"] + ".tsv"
        dest = utility.universal_open(filepath, "a" if args.append else "w")
        # dest = utility.open_destput_file_for_csv_writer(
        #         filepath=filepath,
        #         is_append=args.append)
        if args.append or args.no_write_header:
            is_write_header = False
        else:
            is_write_header = True
        with dest:
            # writer = utility.get_csv_writer(
            #         dest=dest,
            #         delimiter=args.field_delimiter)
            try:
                results = simulator.execute(
                        nreps=args.num_reps,
                        dest=dest,
                        results_store=None,
                        is_write_header=is_write_header)
            except Exception as e:
                sys.stderr.write("Traceback (most recent call last):\n  {}{}\n".format(
                    "  ".join(traceback.format_tb(sys.exc_info()[2])),
                    e))
                sys.exit(1)

if __name__ == "__main__":
    main()

