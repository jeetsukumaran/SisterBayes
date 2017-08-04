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
import dendropy
from sisterbayes import model

class SisterBayesSummaryStatsCalculator(object):

    def __init__(self, **kwargs):
        self.output_prefix = kwargs.pop("output_prefix", "sisterbayes")
        self.is_unfolded_site_frequency_spectrum = kwargs.pop("is_unfolded_site_frequency_spectrum", False)
        self.is_calculate_single_population_sfs = kwargs.pop("is_calculate_single_population_sfs", False)
        self.is_calculate_joint_population_sfs = kwargs.pop("is_calculate_joint_population_sfs", True)
        self.stat_label_prefix = kwargs.pop("stat_label_prefix", "stat")
        self.supplemental_labels = kwargs.pop("supplemental_labels", None)
        self.alignment_directory_head = kwargs.pop("alignment_directory_head", None)
        self.field_delimiter = kwargs.pop("field_delimiter", "\t")
        locus_info = kwargs.pop("locus_info", None)
        params = kwargs.pop("params", None) # ignore
        if locus_info:
            self.model = model.SisterBayesModel(params_d=None, locus_info=locus_info,)
        else:
            self.model = None
        if kwargs:
            raise Exception("Unrecognized configuration entries: {}".format(kwargs))
        self.default_state_alphabet = dendropy.new_standard_state_alphabet("0123456789ACGTU", case_sensitive=False)

    def read_data(self, filepath, datatype, schema):
        if not os.path.isabs(filepath) and self.alignment_directory_head is not None:
            filepath = os.path.join(self.alignment_directory_head, filepath)
        if datatype == "dna":
            data = dendropy.DnaCharacterMatrix.get(path=filepath, schema=schema)
        elif datatype == "standard" or datatype == "snp":
            data = dendropy.StandardCharacterMatrix.get(
                    path=filepath,
                    schema=schema,
                    default_state_alphabet=self.default_state_alphabet)
        return data

    def write_summary_stats(self,
            dest=None,
            results_store=None,
            is_write_header=True,
            ):
        results_d = collections.OrderedDict()
        if self.supplemental_labels:
            for key in self.supplemental_labels:
                results_d[key] = self.supplemental_labels[key]
        for lineage_pair_idx, lineage_pair in enumerate(self.model.lineage_pairs):
            for locus_definition in lineage_pair.locus_definitions:
                field_name_prefix="{}.{}.{}.joint.sfs".format(
                        self.stat_label_prefix,
                        lineage_pair.label,
                        locus_definition.locus_label)
                data = self.read_data(
                        filepath=locus_definition.alignment_filepath,
                        datatype="standard",
                        schema="fasta")
                sequences = data.sequences()
                d0_sequences = sequences[:locus_definition.num_genes_deme0]
                d1_sequences = sequences[locus_definition.num_genes_deme0:]
                assert len(d0_sequences) == locus_definition.num_genes_deme0
                assert len(d1_sequences) == locus_definition.num_genes_deme1
                assert len(sequences) == locus_definition.num_genes_deme0 + locus_definition.num_genes_deme1
                jsfs = self.folded_joint_site_frequency_spectrum(
                        d0_sequences=d0_sequences,
                        d1_sequences=d1_sequences,)
                for row_idx in range(len(jsfs)):
                    for col_idx in range(len(jsfs[row_idx])):
                        results_d["{}.{}.{}".format(field_name_prefix, row_idx, col_idx)] = float(jsfs[row_idx][col_idx])
        # dest.fieldnames = results_d.keys()
        if is_write_header:
            dest.write(self.field_delimiter.join(results_d.keys()))
            dest.write("\n")
        dest.write(self.field_delimiter.join("{}".format(v) for v in results_d.values()))
        dest.write("\n")
                # self.fsc2_handler.run(
                #         field_name_prefix="{}.{}.{}".format(
                #                 self.stat_label_prefix,
                #                 lineage_pair.label,
                #                 locus_definition.locus_label),
                #         fsc2_config_d=fsc2_run_configurations[locus_definition],
                #         random_seed=self.rng.randint(1, 1E6),
                #         results_d=results_d,
                #         )
        return results_d


    def folded_joint_site_frequency_spectrum(self,
            d0_sequences,
            d1_sequences,
            is_discard_multiple_mutation_site=True):
        deme_sequences = (d0_sequences, d1_sequences)
        # weirdly, FastsimCoal2 puts first deme second axis, i.e. columns,
        # while second deme gets put on rows
        jsfs = [[0 for i in range(len(d0_sequences)+1)] for j in range(len(d1_sequences)+1)]
        num_demes = 2
        nsites = None
        deme_site_columns = []
        for deme_idx in range(num_demes):
            deme_sites = zip(*(s.symbols_as_list() for s in deme_sequences[deme_idx]))
            if nsites is None:
                nsites = len(deme_sites)
            else:
                assert len(deme_sites) == nsites
            deme_site_columns.append(deme_sites)
        for site_idx in range(len(deme_site_columns[0])):
            deme_counters = []
            pooled_counter = collections.Counter()
            for deme_idx in range(num_demes):
                deme_counter = collections.Counter(deme_site_columns[deme_idx][site_idx])
                deme_counters.append(deme_counter)
                pooled_counter.update(deme_counter)
            if len(pooled_counter) == 1:
                jsfs[0][0] += 1
                continue
            majority_allele = pooled_counter.most_common(1)[0][0]
            del pooled_counter[majority_allele]
            if is_discard_multiple_mutation_site and len(pooled_counter) > 1:
                continue
            for deme_idx in range(num_demes):
                del deme_counters[deme_idx][majority_allele]
            jsfs[sum(deme_counters[1].values())][sum(deme_counters[0].values())] += 1
        return jsfs
