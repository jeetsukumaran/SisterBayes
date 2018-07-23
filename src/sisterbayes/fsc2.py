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
import itertools
import random
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
            is_store_raw_alignment=False,
            is_store_raw_mutation_tree=False,
            is_store_raw_true_tree=False,
            raw_data_alignment_format="fasta",
            raw_data_tree_format="nexus",
            fsc2_params_adjustment_hack=None,
            is_debug_mode=False,
            rng=None
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
        self.is_store_raw_alignment=is_store_raw_alignment,
        self.is_store_raw_mutation_tree=is_store_raw_mutation_tree,
        self.is_store_raw_true_tree=is_store_raw_true_tree,
        if self.is_store_raw_alignment or self.is_store_raw_mutation_tree or self.is_store_raw_true_tree:
            self.is_store_raw_data = True
        self.raw_data_alignment_format = raw_data_alignment_format
        self.raw_data_tree_format = raw_data_tree_format
        self.fsc2_params_adjustment_hack = fsc2_params_adjustment_hack
        self.is_debug_mode = is_debug_mode
        self.is_infinite_sites_model = is_infinite_sites_model
        self.is_output_dna_as_snp = False
        self._is_file_system_staged = False
        self._num_executions = 0
        self._current_execution_id = None
        self._parameter_filepath = None
        self._results_dirpath = None
        self._deme0_site_frequency_filepath = None
        self._deme1_site_frequency_filepath = None
        self._joint_site_frequency_filepath = None
        self._arlequin_filepath = None
        self._mut_trees_filepath = None
        self._true_trees_filepath = None
        self.is_compose_raw_data_output_paths_de_novo = False
        if rng is None:
            self.rng = random.Random()
        else:
            self.rng = rng

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

    def _get_result_mut_trees_filepath(self):
        if self._mut_trees_filepath is None:
            self._mut_trees_filepath = os.path.join(self.results_dirpath, "{}_1_mut_trees.trees".format(self.name))
        return self._mut_trees_filepath
    mut_trees_filepath = property(_get_result_mut_trees_filepath)

    def _get_result_true_trees_filepath(self):
        if self._true_trees_filepath is None:
            self._true_trees_filepath = os.path.join(self.results_dirpath, "{}_1_true_trees.trees".format(self.name))
        return self._true_trees_filepath
    true_trees_filepath = property(_get_result_true_trees_filepath)

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
            is_normalize_by_site_counts,
            lineage_pair,
            locus_definition,
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
                key = "{}.{}".format(field_name_prefix, key)
                val = float(val)
                if is_normalize_by_site_counts:
                    val = val / locus_definition.num_sites
                results_d[key] = val
        return results_d

    def _parse_joint_site_frequencies(self,
            filepath,
            field_name_prefix,
            is_normalize_by_site_counts,
            lineage_pair,
            locus_definition,
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
                    val = float(val)
                    if is_normalize_by_site_counts:
                        val = val / locus_definition.num_sites
                    results_d["{}.{}.{}".format(field_name_prefix, row_idx, col_idx)] = val
                    col_idx += 1
                row_idx += 1
        return results_d

    def _harvest_run_results(self,
            field_name_prefix,
            is_normalize_by_site_counts,
            lineage_pair,
            locus_definition,
            results_d):
        if self.is_calculate_single_population_sfs:
            self._parse_deme_site_frequencies(
                    filepath=self.deme0_site_frequency_filepath,
                    field_name_prefix="{}.{}.sfs".format(field_name_prefix, compose_deme_label(0)),
                    is_normalize_by_site_counts=is_normalize_by_site_counts,
                    lineage_pair=lineage_pair,
                    locus_definition=locus_definition,
                    results_d=results_d)
            self._parse_deme_site_frequencies(
                    filepath=self.deme1_site_frequency_filepath,
                    field_name_prefix="{}.{}.sfs".format(field_name_prefix, compose_deme_label(1)),
                    is_normalize_by_site_counts=is_normalize_by_site_counts,
                    lineage_pair=lineage_pair,
                    locus_definition=locus_definition,
                    results_d=results_d)
        if self.is_calculate_joint_population_sfs:
            self._parse_joint_site_frequencies(
                    filepath=self.joint_site_frequency_filepath,
                    field_name_prefix="{}.joint.sfs".format(field_name_prefix),
                    is_normalize_by_site_counts=is_normalize_by_site_counts,
                    lineage_pair=lineage_pair,
                    locus_definition=locus_definition,
                    results_d=results_d)
        return results_d

    def _postprocess_raw_tree_taxa(self, tree):
        for taxon in tree.taxon_namespace:
            individual_id, population_id = taxon.label.split(".")
            taxon.label = self._compose_raw_data_taxon_label(
                    population_id=int(population_id),
                    individual_id=int(individual_id))
        tree.taxon_namespace.sort(key=lambda x:x.label)

    def _compose_raw_data_taxon_label(self, population_id, individual_id):
        return "T.{:02d}.{:>03}".format(population_id, individual_id)

    def _parse_raw_results_dna(self):
        data_dict = collections.OrderedDict()
        decodings = list(itertools.permutations("ACGT", 4))
        column_decodings = {}
        with utility.universal_open(self.arlequin_filepath) as src:
            idx = 0
            first_pop_max_ind_id = 0
            in_alignment = False
            for row in src:
                row = row[:-1] # chomp terminating \n
                if in_alignment:
                    if row == "":
                        if idx == 2:
                            break
                        else:
                            in_alignment = False
                    else:
                        try:
                            x1, x2, encoded_data = row.split("\t")
                            decoded_data = []
                            for col_idx, ch in enumerate(encoded_data):
                                try:
                                    decoding_lookup = column_decodings[col_idx]
                                except KeyError:
                                    decoding_lookup = self.rng.choice(decodings)
                                    column_decodings[col_idx] = decoding_lookup
                                decoded_data.append( decoding_lookup[ int(ch) ] )
                            pop_id, ind_id = [int(label_part) for label_part in x1.split("_")]
                            # this ugly hack is because fastsimcoal inconsistently numbers
                            # individuals on trees vs sequences ... even
                            # assuming there is a correspondence
                            if pop_id == 1:
                                first_pop_max_ind_id = max(first_pop_max_ind_id, ind_id)
                            else:
                                ind_id += first_pop_max_ind_id
                            taxon_label = self._compose_raw_data_taxon_label(
                                    population_id=pop_id,
                                    individual_id=ind_id)
                            # taxon_label = self._compose_raw_data_taxon_label(
                            #         population_id=idx,
                            #         individual_id=x1)
                            data_dict[taxon_label] = "".join(decoded_data)
                        except IndexError:
                            raise
                elif "SampleData=" in row:
                    idx += 1
                    in_alignment = True
        return dendropy.DnaCharacterMatrix.from_dict(data_dict)

    def _store_dna_results(self, path_stem):
        dna = self._parse_raw_results_dna()
        writer_kwargs = {}
        if self.raw_data_alignment_format == "fasta":
            writer_kwargs["wrap"] = False
        dna.write(
                path=path_stem + ".seqs.{}".format(self.raw_data_alignment_format),
                schema=self.raw_data_alignment_format,
                **writer_kwargs)

    def _store_true_tree_results(self, path_stem):
        # branch lengths in units of generations
        tree = dendropy.Tree.get(path=self.true_trees_filepath, schema="nexus")
        self._postprocess_raw_tree_taxa(tree)
        if self.is_debug_mode:
            tree.write(path=path_stem + ".tree.gens-artificial.{}".format(self.raw_data_tree_format), schema=self.raw_data_tree_format)
        for nd in tree:
            if nd.edge.length is None:
                nd.edge.length = 0.0
            else:
                nd.edge.length = float(nd.edge.length)
            nd.edge.length = nd.edge.length / self.fsc2_params_adjustment_hack
        tree.write(path=path_stem + ".tree.time.{}".format(self.raw_data_tree_format), schema=self.raw_data_tree_format)
        return tree

    def _store_mut_tree_results(self, path_stem, locus_definition):
        # branch lengths in (absolute) number of mutations
        tree = dendropy.Tree.get(path=self.mut_trees_filepath, schema="nexus")
        self._postprocess_raw_tree_taxa(tree)
        if self.is_debug_mode:
            tree.write(path=path_stem + ".tree.mut-counts.{}".format(self.raw_data_tree_format), schema=self.raw_data_tree_format)
        # convert branch lengths to per site number of mutations
        for nd in tree:
            if nd.edge.length is None:
                nd.edge.length = 0.0
            else:
                nd.edge.length = float(nd.edge.length)
            nd.edge.length = nd.edge.length / locus_definition.num_sites
        tree.write(path=path_stem + ".tree.mut-props.{}".format(self.raw_data_tree_format), schema=self.raw_data_tree_format)
        return tree

    def _store_raw_results(self,
            output_prefix,
            lineage_pair,
            locus_definition,
            ):
        if self.is_compose_raw_data_output_paths_de_novo:
            path_stem = "{}.{}.{}".format(output_prefix, lineage_pair.label, locus_definition.locus_label)
        else:
            path_stem = "{}.{}".format(output_prefix, os.path.splitext(os.path.basename(locus_definition.alignment_filepath))[0])
        if self.is_store_raw_alignment:
            self._store_dna_results(path_stem)
        if self.is_store_raw_true_tree:
            self._store_true_tree_results(path_stem)
        if self.is_store_raw_mutation_tree:
            self._store_mut_tree_results(path_stem, locus_definition)

    def _post_execution_cleanup(self):
        pass

    def run(self,
            field_name_prefix,
            fsc2_config_d,
            random_seed,
            results_d,
            is_normalize_by_site_counts=False,
            raw_data_output_prefix=None,
            lineage_pair=None,
            locus_definition=None,
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
        if not self.is_store_raw_data:
            cmds.append("-x")                       # -x  --noarloutput       : does not generate Arlequin output
        else:
            cmds.append("-T")
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
                is_normalize_by_site_counts=is_normalize_by_site_counts,
                lineage_pair=lineage_pair,
                locus_definition=locus_definition,
                results_d=results_d)
        if self.is_store_raw_data:
            self._store_raw_results(
                    output_prefix=raw_data_output_prefix,
                    lineage_pair=lineage_pair,
                    locus_definition=locus_definition,
                    )
        self._post_execution_cleanup()

