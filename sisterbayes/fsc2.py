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

import dendropy
import subprocess
import collections
import os
from sisterbayes import utility

FSC2_CONFIG_TEMPLATE = """\
// Number of population samples (demes)
2
// Population effective sizes (number of genes)
{d0_population_size}
{d1_population_size}
// Sample sizes
{d0_sample_size}
{d1_sample_size}
// Growth rates: negative growth implies population expansion
0
0
// Number of migration matrices : 0 implies no migration between demes
0
// Historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 4 historical event
1  historical event
{div_time} 0 1 1 2 0 0
// Number of independent loci [chromosome]; '0' => same structure for all loci
1 0
// Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
// Per Block:data type, number of loci, per generation recombination rate, per generation mutation rate and optional parameters
DNA {num_sites} {recombination_rate} {mutation_rate} {ti_proportional_bias}
// Command: {fsc2_command}
"""

class Fsc2RuntimeError(RuntimeError):
    def __init__(self, msg):
        RuntimeError.__init__(self, msg)

class Fsc2Handler(object):

    def __init__(self,
            name,
            fsc2_path,
            working_directory,
            is_calculate_single_population_sfs,
            is_calculate_joint_population_sfs,
            is_unfolded_site_frequency_spectrum,
            is_infinite_sites_model,
            ):
        self.name = name
        self.fsc2_path = fsc2_path
        self.working_directory = working_directory
        self.is_unfolded_site_frequency_spectrum = is_unfolded_site_frequency_spectrum
        if self.is_unfolded_site_frequency_spectrum:
            self.sfs_file_prefix = "DAF"
            self.fsc2_sfs_generation_command = "-d"
        else:
            self.sfs_file_prefix = "MAF"
            self.fsc2_sfs_generation_command = "-m"
        self.is_calculate_single_population_sfs = is_calculate_single_population_sfs
        self.is_calculate_joint_population_sfs = is_calculate_joint_population_sfs
        self.is_infinite_sites_model = is_infinite_sites_model
        self.is_output_dna_as_snp = False
        self._is_file_system_staged = False
        self._num_executions = 0
        self._current_execution_id = None
        self._parameter_filepath = None
        self._results_dirpath = None
        self._deme0_site_frequency_filepath = None
        self._deme1_site_frequency_filepath = None
        self._arlequin_filepath = None
        self._joint_site_frequency_filepath = None

    def _get_parameter_filepath(self):
        if self._parameter_filepath is None:
            # self._parameter_filepath = os.path.join(self.working_directory, ".".join([self.name, "par"]))
            self._parameter_filepath = ".".join([self.name, "par"])
        return self._parameter_filepath
    parameter_filepath = property(_get_parameter_filepath)

    def _get_results_dirpath(self):
        if self._results_dirpath is None:
            self._results_dirpath = os.path.join(self.working_directory, os.path.splitext(self.parameter_filepath)[0])
        return self._results_dirpath
    results_dirpath = property(_get_results_dirpath)

    def _get_result_deme0_site_frequency_filepath(self):
        if self._deme0_site_frequency_filepath is None:
            self._deme0_site_frequency_filepath = os.path.join(self.results_dirpath, "{}_{}pop0.obs".format(self.name, self.sfs_file_prefix))
        return self._deme0_site_frequency_filepath
    deme0_site_frequency_filepath = property(_get_result_deme0_site_frequency_filepath)

    def _get_result_deme1_site_frequency_filepath(self):
        if self._deme1_site_frequency_filepath is None:
            self._deme1_site_frequency_filepath = os.path.join(self.results_dirpath, "{}_{}pop1.obs".format(self.name, self.sfs_file_prefix))
        return self._deme1_site_frequency_filepath
    deme1_site_frequency_filepath = property(_get_result_deme1_site_frequency_filepath)

    def _get_result_joint_site_frequency_filepath(self):
        if self._joint_site_frequency_filepath is None:
            self._joint_site_frequency_filepath = os.path.join(self.results_dirpath, "{}_joint{}pop1_0.obs".format(self.name, self.sfs_file_prefix))
        return self._joint_site_frequency_filepath
    joint_site_frequency_filepath = property(_get_result_joint_site_frequency_filepath)

    def _get_result_arlequin_filepath(self):
        if self._arlequin_filepath is None:
            self._arlequin_filepath = os.path.join(self.results_dirpath, "{}_1_1.arp".format(self.name))
        return self._arlequin_filepath
    arlequin_filepath = property(_get_result_arlequin_filepath)

    def _new_execution_reset(self):
        self._current_execution_id = None
        self._parameter_filepath = None

    def _setup_for_execution(self):
        self._new_execution_reset()
        if not self._is_file_system_staged:
            self._stage_filesystem()

    def _stage_filesystem(self):
        if not os.path.exists(self.working_directory):
            os.makedirs(self.working_directory)
        self._is_file_system_staged = True

    def _generate_parameter_file(self,
            fsc2_config_d,):
        assert self.parameter_filepath
        with utility.universal_open(os.path.join(self.working_directory, self.parameter_filepath), "w") as dest:
            self._write_parameter_configuration(
                    dest=dest,
                    fsc2_config_d=fsc2_config_d,
                    )

    def _write_parameter_configuration(self, dest, fsc2_config_d):
            config = FSC2_CONFIG_TEMPLATE.format(**fsc2_config_d)
            dest.write(config)

    def _parse_deme_site_frequencies(self,
            filepath,
            field_name_prefix,
            results_d):
        with utility.universal_open(filepath) as src:
            lines = src.read().split("\n")
            assert len(lines) == 4 and lines[3] == ""
            header_row = lines[1].split("\t")
            results_d_row = lines[2].split("\t")
            assert len(header_row) == len(results_d_row)
            for key, val in zip(header_row, results_d_row):
                if not val:
                    continue
                results_d["{}.{}".format(field_name_prefix, key)] = float(val)
        return results_d

    def _parse_joint_site_frequencies(self,
            filepath,
            field_name_prefix,
            results_d):
        with utility.universal_open(filepath) as src:
            lines = src.read().split("\n")
            col_keys = [c for c in lines[1].split("\t")[1:] if "sites with multiple" not in c]
            row_idx = 0
            for line in lines[2:]:
                if not line:
                    continue
                cols = line.split("\t")
                if len(cols) - 1 != len(col_keys):
                    raise ValueError("Row {}: Expecting {} columns, but found {}: {}".format(
                        row_idx + 1,
                        len(col_keys),
                        len(cols) - 1,
                        cols))
                assert len(cols) - 1 == len(col_keys)
                row_key = cols[0]
                col_idx = 0
                for col_key, val in zip(col_keys, cols[1:]):
                    # results_d["{}.{}.{}".format(field_name_prefix, row_key, col_key)] = float(val)
                    results_d["{}.{}.{}".format(field_name_prefix, row_idx, col_idx)] = float(val)
                    col_idx += 1
                row_idx += 1
        return results_d

    def _harvest_run_results(self, field_name_prefix, results_d):
        if self.is_calculate_single_population_sfs:
            self._parse_deme_site_frequencies(
                    filepath=self.deme0_site_frequency_filepath,
                    field_name_prefix="{}.{}.sfs".format(field_name_prefix, compose_deme_label(0)),
                    results_d=results_d)
            self._parse_deme_site_frequencies(
                    filepath=self.deme1_site_frequency_filepath,
                    field_name_prefix="{}.{}.sfs".format(field_name_prefix, compose_deme_label(1)),
                    results_d=results_d)
        if self.is_calculate_joint_population_sfs:
            self._parse_joint_site_frequencies(
                    filepath=self.joint_site_frequency_filepath,
                    field_name_prefix="{}.joint.sfs".format(field_name_prefix),
                    results_d=results_d)
        return results_d

    def _harvest_raw_results(self,
            output_prefix,
            lineage_pair_label,
            locus_label):
        data_dicts = []
        with utility.universal_open(self.arlequin_filepath) as src:
            idx = 0
            in_alignment = False
            for row in src:
                row = row[:-1] # chomp terminating \n
                if in_alignment:
                    if row == "":
                        if len(data_dicts) == 2:
                            break
                        else:
                            in_alignment = False
                    else:
                        try:
                            x1, x2, data = row.split("\t")
                            data = data.replace("0", "A").replace("1", "C").replace("2", "G").replace("3", "T")
                            data_dicts[-1]["T"+x1] = data
                        except IndexError:
                            raise
                elif "SampleData=" in row:
                    in_alignment = True
                    data_dicts.append(collections.OrderedDict())
        assert len(data_dicts) == 2
        for idx, data_dict in enumerate(data_dicts):
            dna = dendropy.DnaCharacterMatrix.from_dict(data_dict)
            dna.write(
                    path="{}.{}.{}.{}.fasta".format(output_prefix, lineage_pair_label, idx+1, locus_label),
                    schema="fasta",
                    wrap=False)

    def _post_execution_cleanup(self):
        pass

    def run(self,
            field_name_prefix,
            fsc2_config_d,
            random_seed,
            results_d,
            is_store_raw_data=False,
            raw_data_output_prefix=None,
            lineage_pair_label=None,
            locus_label=None,
            ):
        self._setup_for_execution()
        cmds = []
        cmds.append(self.fsc2_path)
        cmds.extend(["-i", self.parameter_filepath])
        cmds.append(self.fsc2_sfs_generation_command)
        cmds.extend(["-n", "1"])                # number of simulations to perform
        cmds.extend(["-r", str(random_seed)])   # seed for random number generator (positive integer <= 1E6)
        if self.is_infinite_sites_model:
            cmds.append("-I")                    # -I  --inf               : generates DNA mutations according to an infinite site (IS) mutation model
        cmds.append("-S")                       # -S  --allsites          : output the whole DNA sequence, incl. monomorphic sites
        cmds.append("-s0")                      # -s  --dnatosnp 2000     : output DNA as SNP data, and specify maximum no. of SNPs to output (use 0 to output all SNPs). (required to calculate SFS)
        if not is_store_raw_data:
            cmds.append("-x")                       # -x  --noarloutput       : does not generate Arlequin output
        fsc2_config_d["fsc2_command"] = " ".join(cmds)
        self._generate_parameter_file(fsc2_config_d)
        p = subprocess.Popen(cmds,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=self.working_directory,
                )
        stdout, stderr = utility.communicate_process(p)
        if p.returncode != 0:
            raise Fsc2RuntimeError("FastSimCoal2 execution failure: {}".format(stderr))
        self._num_executions += 1
        if results_d is None:
            results_d = collections.OrderedDict()
        self._harvest_run_results(
                field_name_prefix=field_name_prefix,
                results_d=results_d)
        if is_store_raw_data:
            self._harvest_raw_results(
                    output_prefix=raw_data_output_prefix,
                    lineage_pair_label=lineage_pair_label,
                    locus_label=locus_label,
                    )
        self._post_execution_cleanup()

