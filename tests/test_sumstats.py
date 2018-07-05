#! /usr/bin/env python

import os
import unittest
import dendropy
import sys
from sisterbayes import utility
from sisterbayes import sumstats
TESTS_DIR = os.path.realpath(os.path.dirname(__file__))
TESTS_DATA_DIR = os.path.join(TESTS_DIR, "data")

def _pathmap(subpath):
    return os.path.join(TESTS_DATA_DIR, "genomic", subpath)

# class SfsTestDataSet(object):

#     def __init__(self, filename_prefix):
#         self.filename_prefix = filename_prefix
#         self.dirpath = os.path.join(TESTS_DATA_DIR, "genomic")
#         self.read_source_data()

#     def read_source_data(self):
#         self.dna = []
#         self.snps = []
#         self.dna_fasta_filepaths = []
#         self.snps_fasta_filepaths = []
#         for idx in range(2):
#             self.dna.append([])
#             src_path = os.path.join(self.dirpath, "{}.data.dna.deme{}.txt".format(self.filename_prefix, idx))
#             with utility.universal_open(src_path) as src:
#                 for row in src:
#                     self.dna[idx].append(row.strip())
#             self.snps.append([])
#             src_path = os.path.join(self.dirpath, "{}.data.snps.deme{}.txt".format(self.filename_prefix, idx))
#             with utility.universal_open(src_path) as src:
#                 for row in src:
#                     self.snps[idx].append(row.strip())
#             self.dna_fasta_filepaths.append(os.path.join(self.dirpath, "testdataset01.data.dna.deme{}.fasta".format(idx)))
#             self.snps_fasta_filepaths.append(os.path.join(self.dirpath, "testdataset01.data.snps.deme{}.fasta".format(idx)))

# _TEST_DATASETS = [
#         SfsTestDataSet("testdataset01"),
#         ]

# class DataReadingTestCases(unittest.TestCase):

#     def test_read_dna(self):
#         sc = sumstats.SisterBayesSummaryStatsCalculator()
#         for test_ds in _TEST_DATASETS:
#             for idx in range(2):
#                 data = sc.read_data(
#                         filepath=test_ds.dna_fasta_filepaths[idx],
#                         datatype="dna",
#                         schema="fasta")
#                 obs_sequences = [s.symbols_as_string() for s in data.sequences()]
#                 self.assertEqual(obs_sequences, test_ds.dna[idx])

#     def test_read_snps(self):
#         sc = sumstats.SisterBayesSummaryStatsCalculator()
#         for test_ds in _TEST_DATASETS:
#             for idx in range(2):
#                 data = sc.read_data(
#                         filepath=test_ds.snps_fasta_filepaths[idx],
#                         datatype="standard",
#                         schema="fasta")
#                 obs_sequences = [s.symbols_as_string() for s in data.sequences()]
#                 self.assertEqual(obs_sequences, test_ds.snps[idx])

# class SiteFrequencySpectrumCalcTestCase(unittest.TestCase):

#     def test_unfolded_sfs(self):
#         print("ok")
#         sc = sumstats.SisterBayesSummaryStatsCalculator()
#         for test_ds in _TEST_DATASETS:
#             for datatype, source_filepaths in (
#                     ("dna", test_ds.dna_fasta_filepaths,),
#                     ("standard", test_ds.snps_fasta_filepaths),
#                     ):
#                 for idx in range(2):
#                     data = sc.read_data(
#                             filepath=source_filepaths[idx],
#                             datatype=datatype,
#                             schema="fasta")
#                     print("> {}: {}: {}".format(datatype, idx, data.folded_site_frequency_spectrum()))
#                     print("--")


class TwoPopSfsTests(unittest.TestCase):

    def read_expected_sfs(self, filename):
        filepath = _pathmap(filename)
        with open(filepath) as src:
            return [float(v) for v in src.read().strip().split(",")]

    def read_expected_jsfs(self, filename):
        filepath = _pathmap(filename)
        vectors = []
        with open(filepath) as src:
            for row in src:
                vectors.append( [float(v) for v in row.strip().split(",")] )
        return vectors

    def read_obs_data(self, test_data_name, datatype):
        obs_data = []
        for deme_idx in (0, 1):
            obs_data_path = _pathmap(test_data_name + ".data.{}.deme{}.fasta".format(datatype, deme_idx))
            if datatype == "dna":
                obs_data.append( dendropy.DnaCharacterMatrix.get(path=obs_data_path, schema="fasta") )
            else:
                obs_data.append( dendropy.StandardCharacterMatrix.get(path=obs_data_path, schema="fasta") )
        return obs_data

    def test_two_pop_sfs(self):
        ss = sumstats.SisterBayesSummaryStatsCalculator()
        for test_data_name in (
                "testdataset01",
                "testdataset02",
                "testdataset03",
                ):
            for datatype in ("dna", "std"):
                d0, d1 = self.read_obs_data(test_data_name=test_data_name, datatype=datatype)
                obs_folded_jsfs = ss.folded_joint_site_frequency_spectrum(d0.sequences(), d1.sequences())
                exp_folded_jsfs = self.read_expected_jsfs("{}.jsfs.folded.txt".format(test_data_name))
                self.assertEqual(obs_folded_jsfs, exp_folded_jsfs)

if __name__ == "__main__":
    unittest.main()
