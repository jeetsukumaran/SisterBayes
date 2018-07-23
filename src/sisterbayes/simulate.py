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

import collections
import random
import sys
import os
import time
try:
    # Python 3
    import queue
except ImportError:
    # Python 2.7
    import Queue as queue
import multiprocessing
import traceback

import sisterbayes
from sisterbayes import utility
from sisterbayes import model
from sisterbayes import fsc2

class SimulationWorker(multiprocessing.Process):

    def __init__(self,
            name,
            model,
            task_queue,
            results_queue,
            fsc2_path,
            working_directory,
            run_logger,
            logging_frequency,
            messenger_lock,
            random_seed,
            is_calculate_single_population_sfs,
            is_calculate_joint_population_sfs,
            is_unfolded_site_frequency_spectrum,
            is_infinite_sites_model,
            is_concatenate_loci,
            is_normalize_by_site_counts,
            stat_label_prefix,
            is_include_model_id_field,
            supplemental_labels,
            is_debug_mode,
            is_store_raw_alignment,
            is_store_raw_mutation_tree,
            is_store_raw_true_tree,
            raw_data_output_prefix,
            raw_data_alignment_format,
            raw_data_tree_format,
            ):
        multiprocessing.Process.__init__(self, name=name)
        self.model = model
        self.rng = random.Random(random_seed)
        self.task_queue = task_queue
        self.results_queue = results_queue
        self.run_logger = run_logger
        self.logging_frequency = logging_frequency
        self.messenger_lock = messenger_lock
        self.is_unfolded_site_frequency_spectrum = is_unfolded_site_frequency_spectrum
        self.stat_label_prefix = stat_label_prefix
        self.is_include_model_id_field = is_include_model_id_field
        self.supplemental_labels = supplemental_labels
        self.is_debug_mode = is_debug_mode
        self.kill_received = False
        self.num_tasks_received = 0
        self.num_tasks_completed = 0
        self.is_store_raw_alignment = is_store_raw_alignment
        # self.is_store_concatenated_alignment = True # TODO
        self.is_concatenate_loci = is_concatenate_loci
        self.is_store_raw_mutation_tree = is_store_raw_mutation_tree
        self.is_store_raw_true_tree = is_store_raw_true_tree
        self.raw_data_output_prefix = raw_data_output_prefix
        self.raw_data_alignment_format = raw_data_alignment_format
        self.raw_data_tree_format = raw_data_tree_format
        self.is_normalize_by_site_counts = is_normalize_by_site_counts,
        self.fsc2_handler = fsc2.Fsc2Handler(
                name=name,
                fsc2_path=fsc2_path,
                working_directory=working_directory,
                is_calculate_single_population_sfs=is_calculate_single_population_sfs,
                is_calculate_joint_population_sfs=is_calculate_joint_population_sfs,
                is_unfolded_site_frequency_spectrum=is_unfolded_site_frequency_spectrum,
                is_infinite_sites_model=is_infinite_sites_model,
                is_store_raw_alignment=self.is_store_raw_alignment,
                is_store_raw_mutation_tree=self.is_store_raw_mutation_tree,
                is_store_raw_true_tree=self.is_store_raw_true_tree,
                raw_data_alignment_format=self.raw_data_alignment_format,
                raw_data_tree_format=self.raw_data_tree_format,
                is_debug_mode=self.is_debug_mode,
                fsc2_params_adjustment_hack=model.fsc2_params_adjustment_hack
                )

    def send_worker_message(self, msg, level):
        if self.run_logger is None:
            return
        # if self.run_logger.messaging_level > level or self.messenger.silent:
        #     return
        msg = "{}: {}".format(self.name, msg)
        self.messenger_lock.acquire()
        try:
            self.run_logger.log(msg, level=level)
        finally:
            self.messenger_lock.release()

    def send_worker_critical(self, msg):
        self.send_worker_message(msg, utility.RunLogger.CRITICAL_MESSAGING_LEVEL)

    def send_worker_debug(self, msg):
        self.send_worker_message(msg, utility.RunLogger.DEBUG_MESSAGING_LEVEL)

    def send_worker_info(self, msg):
        self.send_worker_message(msg, utility.RunLogger.INFO_MESSAGING_LEVEL)

    def send_worker_warning(self, msg):
        self.send_worker_message(msg, utility.RunLogger.WARNING_MESSAGING_LEVEL)

    def send_worker_error(self, msg):
        self.send_worker_message(msg, utility.RunLogger.ERROR_MESSAGING_LEVEL)

    def run(self):
        result = None
        while not self.kill_received:
            try:
                rep_idx = self.task_queue.get_nowait()
            except queue.Empty:
                break
            if rep_idx is None:
                # poison pill
                break
            self.num_tasks_received += 1
            # self.send_worker_critical("Received task: '{task_name}'".format(
            #     _num_assigned_tasks=self.num_tasks_received,
            #     task_name=rep_idx))
            # rng = random.Random(random_seed)
            try:
                result = self.simulate(rep_idx)
            except (KeyboardInterrupt, Exception) as e:
                # traceback.print_exc()
                e.worker_name = self.name
                e.traceback_exc = traceback.format_exc()
                self.results_queue.put(e)
                break
            if self.kill_received:
                break
            self.results_queue.put(result)
            self.num_tasks_completed += 1
            # self.send_info("Completed task {_num_assigned_tasks}: '{task_name}'".format(
        if self.kill_received:
            self.send_worker_warning("Terminating in response to kill request")

    def simulate(self, rep_idx):
        results_d = collections.OrderedDict()
        if self.is_include_model_id_field:
            results_d["model.id"] = None
        if self.supplemental_labels:
            for key in self.supplemental_labels:
                results_d[key] = self.supplemental_labels[key]
        params, fsc2_run_configurations = self.model.sample_parameter_values_from_prior(rng=self.rng)
        results_d.update(params)
        for lineage_pair_idx, lineage_pair in enumerate(self.model.lineage_pairs):
            concatenated_locus_label = model.compose_concatenated_locus_label(lineage_pair)
            for locus_definition_idx, locus_definition in enumerate(lineage_pair.locus_definitions):
                if self.is_concatenate_loci:
                    field_name_prefix="{}.{}.{}".format(
                            self.stat_label_prefix,
                            lineage_pair.label,
                            concatenated_locus_label,
                            )
                    locus_results_store = collections.OrderedDict()
                else:
                    field_name_prefix="{}.{}.{}".format(
                            self.stat_label_prefix,
                            lineage_pair.label,
                            locus_definition.locus_label)
                    locus_results_store = results_d
                self.fsc2_handler.run(
                        field_name_prefix=field_name_prefix,
                        fsc2_config_d=fsc2_run_configurations[locus_definition],
                        random_seed=self.rng.randint(1, 1E6),
                        results_d=locus_results_store,
                        is_normalize_by_site_counts=self.is_normalize_by_site_counts,
                        raw_data_output_prefix="{}.{:04d}".format(self.raw_data_output_prefix, rep_idx+1),
                        lineage_pair=lineage_pair, # only needed for normalization or raw data output path composition
                        locus_definition=locus_definition, # only needed for normalization or raw data output path composition
                        )
                if self.is_concatenate_loci:
                    for key in locus_results_store:
                        try:
                            results_d[key] += locus_results_store[key]
                        except KeyError:
                            results_d[key] = locus_results_store[key]
        if self.is_include_model_id_field:
            results_d["model.id"] = results_d["param.divTimeModel"]
        return results_d

class TaskQueue(object):
    # If we have a large number of replicates, just pre-loading them on the
    # queue seems to stall. We could just load them one-by-one as each
    # worker completes, but as the master cycle also does other processing
    # of results as each worker completes this leaves open the
    # possibilities that the task queue will be empty waiting for the
    # master cycle to complete its results processing before loading the
    # task queue, which means also the possibilty of processes waiting idly
    # until they get assigned a tasks. Not very probable for any moderately
    # complex tasks, but still ...
    # Here we implement a subclass of Queue that allows for loading
    # of tasks in bite-sized chunks, that also supports some
    # housekeeping/tracking.

    def __init__(self, num_processes, job_size):
        self.num_processes = num_processes
        self.job_size = job_size
        self._num_assigned_tasks = 0
        self._num_completed_tasks = 0
        self._task_load_chunk = self.num_processes * 10
        self._tasks = multiprocessing.Queue()
        self.max_active_tasks = self.num_processes * 10

    def load_tasks(self):
        current_load = len(self)
        if current_load > self.max_active_tasks:
            return
        if self.job_size - (self._num_assigned_tasks + self._task_load_chunk) < 0:
            num_to_load = self.job_size - self._num_assigned_tasks
        else:
            num_to_load = min(self.max_active_tasks - current_load, self._task_load_chunk)
        for idx in range(num_to_load):
            self._tasks.put(self._num_assigned_tasks)
            self._num_assigned_tasks += 1
        if self._num_assigned_tasks >= self.job_size:
            for idx in range(self.num_processes):
                self._tasks.put(None) # poison pill

    def task_done(self):
        self._num_completed_tasks += 1

    def get(self, *args, **kwargs):
        return self._tasks.get(*args, **kwargs)

    def get_nowait(self, *args, **kwargs):
        return self._tasks.get_nowait(*args, **kwargs)

    def __len__(self):
        return self._num_assigned_tasks - self._num_completed_tasks

class SisterBayesSimulator(object):

    def __init__(self,
            config_d,
            num_processes=None,
            logging_frequency=1000,
            package_id=None,
            is_verbose_setup=True,
            is_store_raw_alignment=False,
            is_store_raw_mutation_tree=False,
            is_store_raw_true_tree=False,
            raw_data_output_prefix=None,
            raw_data_alignment_format="fasta",
            raw_data_tree_format="nexus",
            is_debug_mode=False,
            ):
        # configure
        if package_id is None:
            self.package_id = sisterbayes.package_id()
        else:
            self.package_id = package_id
        self.elapsed_time = 0.0 # need to be here for logging
        self.is_store_raw_alignment = is_store_raw_alignment
        self.is_store_raw_mutation_tree = is_store_raw_mutation_tree
        self.is_store_raw_true_tree = is_store_raw_true_tree
        self.raw_data_output_prefix = raw_data_output_prefix
        self.raw_data_alignment_format = raw_data_alignment_format
        self.raw_data_tree_format = raw_data_tree_format
        self.is_debug_mode = is_debug_mode
        config_d = dict(config_d) # make copy so we can pop items
        self.is_verbose_setup = is_verbose_setup
        self.configure_simulator(config_d)
        self.num_cpus = utility.cpu_count()
        if num_processes is None or num_processes <= 0:
            self.num_processes = self.num_cpus # problem is hyperthreading results in 2x "real" CPU's :(
        elif num_processes == 1 and self.num_cpus > 1 and self.is_verbose_setup:
            self.run_logger.info(
                    ("Multiple processors ({num_cpus}) available:"
                    " consider using the '-M' or '-m' options to"
                    " parallelize processing"
                    ).format(num_cpus=self.num_cpus))
            self.num_processes = 1
        else:
            self.num_processes = num_processes
        if self.is_verbose_setup:
            self.run_logger.info("Will run up to {} processes in parallel".format(self.num_processes))
            self.run_logger.info("{} lineage pairs in analysis:".format(self.model.num_lineage_pairs))
            for lineage_pair_idx, lineage_pair in enumerate(self.model.lineage_pairs):
                self.run_logger.info("  - '{}': {:>2d} loci (Samples: {})".format(
                        lineage_pair.label,
                        len(lineage_pair.locus_definitions),
                        ", ".join("{}/{}".format(locus.num_genes_deme0, locus.num_genes_deme1) for locus in lineage_pair.locus_definitions),
                        ))
        self.worker_class = SimulationWorker
        self._task_queue = None
        self._num_assigned_tasks = None
        self._task_load_chunk = None

    def configure_simulator(self, config_d, verbose=True):
        self.title = config_d.pop("title", "sisterbayes-{}-{}".format(time.strftime("%Y%m%d%H%M%S"), id(self)))
        self.output_prefix = config_d.pop("output_prefix", self.title)
        self.working_directory = config_d.pop("working_directory", self.title)
        self.run_logger = config_d.pop("run_logger", None)
        if self.run_logger is None:
            self.run_logger = utility.RunLogger(
                    name="sisterbayes-simulate",
                    stderr_logging_level=config_d.pop("standard_error_logging_level", "info"),
                    log_to_file=config_d.pop("log_to_file", True),
                    log_to_stderr=config_d.pop("log_to_stderr", True),
                    log_path=self.output_prefix + ".log",
                    file_logging_level=config_d.pop("file_logging_level", "info"),
                    )
        self.run_logger.system = self
        self.logging_frequency = config_d.pop("logging_frequency", 1000)
        if self.is_verbose_setup:
            self.run_logger.info("Running: {}".format(self.package_id))
            self.run_logger.info("Configuring simulation: '{}'".format(self.title))
            self.run_logger.info("Current directory: '{}'".format(os.path.abspath(os.curdir)))
            self.run_logger.info("Working directory: '{}'".format(self.working_directory))
            self.run_logger.info("Output directory: '{}'".format(os.path.dirname(os.path.abspath(self.output_prefix))))
            self.run_logger.info("Output filename prefix: '{}'".format(os.path.basename(self.output_prefix)))
        self.fsc2_path = config_d.pop("fsc2_path", "fsc")
        if self.is_verbose_setup:
            self.run_logger.info("FastSimCoal2 path: '{}'".format(self.fsc2_path))
        self.rng = config_d.pop("rng", None)
        if self.rng is None:
            self.random_seed = config_d.pop("random_seed", None)
            if self.random_seed is None:
                self.random_seed = random.randint(1, sys.maxsize)
            if self.is_verbose_setup:
                self.run_logger.info("Initializing with random seed: {}".format(self.random_seed))
            self.rng = random.Random(self.random_seed)
        else:
            if "random_seed" in config_d:
                raise TypeError("Cannot specify both 'rng' and 'random_seed'")
            if self.is_verbose_setup:
                self.run_logger.info("Using existing random number generator")
        if self.is_verbose_setup and self.is_debug_mode:
            self.run_logger.info("Running in DEBUG mode")
        self.site_frequency_spectrum_type = config_d.pop("site_frequency_spectrum_type", "unfolded").lower()
        self.is_unfolded_site_frequency_spectrum = config_d.pop("is_unfolded_site_frequency_spectrum", False)
        self.is_calculate_single_population_sfs = config_d.pop("is_calculate_single_population_sfs", False)
        self.is_calculate_joint_population_sfs = config_d.pop("is_calculate_joint_population_sfs", True)
        self.is_infinite_sites_model = config_d.pop("is_infinite_sites_model", False)
        self.is_concatenate_loci = config_d.pop("is_concatenate_loci", False)
        self.is_normalize_by_site_counts = config_d.pop("is_normalize_by_site_counts", False)
        if not self.is_calculate_single_population_sfs and not self.is_calculate_joint_population_sfs:
            raise ValueError("Neither single-population nor joint site frequency spectrum will be calculated!")
        self.stat_label_prefix = config_d.pop("stat_label_prefix", "stat")
        self.supplemental_labels = config_d.pop("supplemental_labels", None)
        self.field_delimiter = config_d.pop("field_delimiter", "\t")
        self.is_include_model_id_field = config_d.pop("is_include_model_id_field", False)
        if "params" not in config_d:
            raise ValueError("Missing 'params' entry in configuration")
        params_d = config_d.pop("params")
        if "locus_info" not in config_d:
            raise ValueError("Missing 'locus_info' entry in configuration")
        locus_info = config_d.pop("locus_info")
        self.model = model.SisterBayesModel(params_d=params_d, locus_info=locus_info,)
        if config_d:
            raise Exception("Unrecognized configuration entries: {}".format(config_d))

    def get_eta(self, start_time, result_count, nreps):
        rate = (time.time() - start_time) / float(result_count)
        return (nreps - result_count) * rate
        # if result_count < self.num_processes:
        #     return (nreps/float(self.num_processes)) * rate
        # else:
        #     return (nreps - result_count) * rate

    def formatted_eta(self, start_time, result_count, nreps):
        return utility.format_elapsed_time(self.get_eta(
            start_time=start_time,
            result_count=result_count,
            nreps=nreps))

    def execute(self,
            nreps,
            dest=None,
            results_store=None,
            params_only_dest=None,
            is_write_header=True,
            ):
        # load up queue
        self.run_logger.info("Priming task queue")
        # Does not work very well (stalls) with large
        # number of tasks (e.g., > 100K)
        # for rep_idx in range(nreps):
        #     _task_queue.put( rep_idx )
        # So we use our own task queue class, which abstracts out
        # the job of chunking tasks
        task_queue = TaskQueue(
                num_processes=self.num_processes,
                job_size=nreps)
        task_queue.load_tasks()
        # time.sleep(0.1) # to avoid: 'IOError: [Errno 32] Broken pipe'; https://stackoverflow.com/questions/36359528/broken-pipe-error-with-multiprocessing-queue
        self.run_logger.info("Launching {} worker processes".format(self.num_processes))
        results_queue = multiprocessing.Queue()
        messenger_lock = multiprocessing.Lock()
        workers = []
        main_time_start = time.time()
        for pidx in range(self.num_processes):
            worker = self.worker_class(
                    # name=str(pidx+1),
                    name="{}-{}".format(self.title, pidx+1),
                    model=self.model,
                    task_queue=task_queue,
                    results_queue=results_queue,
                    fsc2_path=self.fsc2_path,
                    working_directory=self.working_directory,
                    run_logger=self.run_logger,
                    logging_frequency=self.logging_frequency,
                    messenger_lock=messenger_lock,
                    random_seed=self.rng.randint(1, sys.maxsize),
                    is_calculate_single_population_sfs=self.is_calculate_single_population_sfs,
                    is_calculate_joint_population_sfs=self.is_calculate_joint_population_sfs,
                    is_unfolded_site_frequency_spectrum=self.is_unfolded_site_frequency_spectrum,
                    is_infinite_sites_model=self.is_infinite_sites_model,
                    is_concatenate_loci=self.is_concatenate_loci,
                    is_normalize_by_site_counts=self.is_normalize_by_site_counts,
                    stat_label_prefix=self.stat_label_prefix,
                    is_include_model_id_field=self.is_include_model_id_field,
                    supplemental_labels=self.supplemental_labels,
                    is_debug_mode=self.is_debug_mode,
                    is_store_raw_alignment=self.is_store_raw_alignment,
                    is_store_raw_mutation_tree=self.is_store_raw_mutation_tree,
                    is_store_raw_true_tree=self.is_store_raw_true_tree,
                    raw_data_output_prefix=self.raw_data_output_prefix,
                    raw_data_alignment_format=self.raw_data_alignment_format,
                    raw_data_tree_format=self.raw_data_tree_format,
                    )
            worker.start()
            workers.append(worker)

        # collate results
        result_count = 0
        try:
            while result_count < nreps:
                try:
                    result = results_queue.get_nowait()
                except queue.Empty:
                    continue
                if isinstance(result, KeyboardInterrupt):
                    raise result
                elif isinstance(result, Exception):
                    self.run_logger.error("Exception raised in worker process '{}'"
                                          "\n>>>\n{}<<<\n".format(
                                              result.worker_name,
                                              result.traceback_exc))
                    raise result
                task_queue.task_done()
                task_queue.load_tasks() # add another chunk of tasks
                if results_store is not None:
                    results_store.append(result)
                if dest is not None:
                    if result_count == 0 and is_write_header:
                        dest.write(self.field_delimiter.join(result.keys()))
                        dest.write("\n")
                    dest.write(self.field_delimiter.join("{}".format(v) for v in result.values()))
                    dest.write("\n")
                if params_only_dest is not None:
                    # result_params = {k:result[k] for k in result if k.startswith("param")}
                    result_params = collections.OrderedDict()
                    for k in result:
                        if not k.startswith("param"):
                            continue
                        result_params[k] = result[k]
                    if result_count == 0 and is_write_header:
                        params_only_dest.write(self.field_delimiter.join(result_params.keys()))
                        params_only_dest.write("\n")
                    params_only_dest.write(self.field_delimiter.join("{}".format(v) for v in result_params.values()))
                    params_only_dest.write("\n")
                result_count += 1
                if self.logging_frequency and result_count % self.logging_frequency == 0:
                    self.run_logger.info("Completed replicate {} (remaining job time: {})".format(
                        result_count,
                        self.formatted_eta(start_time=main_time_start, result_count=result_count, nreps=nreps),
                        ))
                elif result_count == self.num_processes and (self.logging_frequency is not None and self.num_processes < self.logging_frequency):
                    self.run_logger.info("Preliminary estimate of job time: {}".format(
                        self.formatted_eta(start_time=main_time_start, result_count=result_count, nreps=nreps),))
        except (Exception, KeyboardInterrupt) as e:
            for worker in workers:
                worker.terminate()
            raise
        self.run_logger.info("All {} worker processes terminated".format(self.num_processes))
        self.run_logger.info("Job done: {} replicates completed in {}".format(
            nreps,
            utility.format_elapsed_time(time.time() - main_time_start)))
        return results_store

