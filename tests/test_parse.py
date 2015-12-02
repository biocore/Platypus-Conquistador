# ----------------------------------------------------------------------------
# Copyright (c) 2015--, platypus development team.
#
# Distributed under the terms of the BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import division

from copy import copy
from unittest import TestCase, main
from os.path import join, dirname

from platypus.parse import (
    parse_first_database, parse_second_database, process_results,
    M9, parse_m9)


class TopLevelTests(TestCase):
    """ Tests of compare blast databases code """

    def setUp(self):
        self.base = join(dirname(__file__), 'support_files')
        self.db1 = open(join(self.base, 'first_db.txt'), 'r')
        self.db2 = open(join(self.base, 'second_db.txt'), 'r')
        self.blasttest = open(join(self.base, 'blast_test_output.txt'))
        self.smrtest = open(join(self.base, 'sortmerna_test_output.txt'))
        self.smrtest_bad = open(join(self.base,
                                     'sortmerna_test_output_bad.txt'))
        self.m9_empty = M9(*[None] * 12)

    def test_parse_m9_bad(self):
        """Raise on an incomplete record"""
        with self.assertRaises(ValueError):
            list(parse_m9(self.smrtest_bad))

    def test_parse_m9_sortmerna(self):
        """Parse sortmerna output as expected"""
        fna_161_hits = [M9(query='QiimeExactMatch.4502806.3.fna_161',
                           subject='NZ_CAAQ01005666|642979345',
                           percent_id=87.8, aln_length=188, mismatches=23,
                           gapopenings=0, q_start=3, q_end=190, s_start=360,
                           s_end=547, evalue=1.05e-59, bitscore=234)]
        fna_282_hits = [M9(query='QiimeExactMatch.4502806.3.fna_282',
                           subject='NZ_CAAX01000347|642979348',
                           percent_id=75, aln_length=164, mismatches=41,
                           gapopenings=0, q_start=3, q_end=166, s_start=1745,
                           s_end=1908, evalue=1.19e-22, bitscore=111)]
        fna_504_hits = [M9(query='QiimeExactMatch.4502806.3.fna_504',
                           subject='NZ_CAAX01000347|642979348',
                           percent_id=96.4, aln_length=28, mismatches=1,
                           gapopenings=0, q_start=170, q_end=197, s_start=2575,
                           s_end=2602, evalue=0.00257, bitscore=47)]
        fna_562_hits = [M9(query='QiimeExactMatch.4502806.3.fna_562',
                           subject='NZ_CAAQ01004522|642979345_hit_a',
                           percent_id=91.6, aln_length=143, mismatches=12,
                           gapopenings=0, q_start=3, q_end=145, s_start=1369,
                           s_end=1511, evalue=2.63e-50, bitscore=203.0),
                        M9(query='QiimeExactMatch.4502806.3.fna_562',
                           subject='NZ_CAAQ01004522|642979345_hit_b',
                           percent_id=91.6, aln_length=143, mismatches=12,
                           gapopenings=0, q_start=4, q_end=145, s_start=1369,
                           s_end=1511, evalue=2.63e-50, bitscore=203.0)]
        fna_828_hits = [M9(query='QiimeExactMatch.4502806.3.fna_828',
                           subject='NZ_CAAU01001667|642979341',
                           percent_id=76, aln_length=122, mismatches=23,
                           gapopenings=7, q_start=10, q_end=131, s_start=88,
                           s_end=208, evalue=1.03e-12, bitscore=78.0)]

        exp = [('QiimeExactMatch.4502806.3.fna_161', fna_161_hits),
               ('QiimeExactMatch.4502806.3.fna_282', fna_282_hits),
               ('QiimeExactMatch.4502806.3.fna_504', fna_504_hits),
               ('QiimeExactMatch.4502806.3.fna_562', fna_562_hits),
               ('QiimeExactMatch.4502806.3.fna_828', fna_828_hits)]

        obs = list(parse_m9(self.smrtest))
        self.assertEqual(obs, exp)

    def test_parse_m9_blast(self):
        """Parse blast output as expected"""
        fna_6_hits = [M9(query='4502804.3.fna_6', subject='28_interest.fna',
                         percent_id=100.00, aln_length=14, mismatches=0,
                         gapopenings=0, q_start=42, q_end=55, s_start=98,
                         s_end=111, evalue=0.006, bitscore=28.2)]
        fna_8_hits = [M9(query='4502804.3.fna_8',
                         subject='20_interest_hita.fna', percent_id=100.00,
                         aln_length=13, mismatches=0, gapopenings=0,
                         q_start=19, q_end=31, s_start=529, s_end=541,
                         evalue=0.059, bitscore=26.3),
                      M9(query='4502804.3.fna_8',
                         subject='20_interest_hitb.fna', percent_id=100.00,
                         aln_length=13, mismatches=0, gapopenings=0,
                         q_start=19, q_end=31, s_start=529, s_end=541,
                         evalue=0.059, bitscore=26.3),
                      M9(query='4502804.3.fna_8',
                         subject='20_interest_hitc.fna', percent_id=100.00,
                         aln_length=13, mismatches=0, gapopenings=0,
                         q_start=19, q_end=31, s_start=529, s_end=541,
                         evalue=0.059, bitscore=26.3)]
        fna_9_hits = [M9(query='4502804.3.fna_9', subject='26_interest.fna',
                         percent_id=100.00, aln_length=14, mismatches=0,
                         gapopenings=0, q_start=114, q_end=127, s_start=3039,
                         s_end=3026, evalue=0.020, bitscore=28.2)]
        exp = [(None, [copy(self.m9_empty)]),
               ('4502804.3.fna_6', fna_6_hits),
               (None, [copy(self.m9_empty)]),
               ('4502804.3.fna_8', fna_8_hits),
               ('4502804.3.fna_9', fna_9_hits),
               (None, [copy(self.m9_empty)])]

        obs = list(parse_m9(self.blasttest))
        self.assertEqual(obs, exp)

    def test_parse_first_database(self):
        """Parse first db should build best_hits correctly"""

        out_total_queries, out_best_hits = parse_first_database(self.db1,
                                                                [.80], [50])
        self.assertEquals(out_total_queries, 4)
        self.assertEquals(
            out_best_hits,
            {'HABJ36W02EXF44': [
                {'a': {'evalue': 0.0,
                       'subject_id': 'NZ_ABEH01000018_641736102',
                       'bit_score': 1005.0, 'percentage_id': 99.42,
                       'alg_length': 519},
                 'b': {'subject_id': None, 'bit_score': -1}}],
             'BLANK-TEST-NOT-IN-SECOND': [
                {'a': {'evalue': 4e-133,
                       'subject_id': 'NZ_ACZD01000120_647000262',
                       'bit_score': 482.0, 'percentage_id': 88.79,
                       'alg_length': 455},
                 'b': {'subject_id': None, 'bit_score': -1}}],
             'HABJ36W02DLDSY': [
                {'a': {'evalue': 0.0,
                       'subject_id': 'NZ_ABEH01000005_641736102',
                       'bit_score': 959.0, 'percentage_id': 99.22,
                       'alg_length': 512},
                 'b': {'subject_id': None, 'bit_score': -1}}]})

    def test_parse_second_database(self):
        """Check the second database is parsed correctly & updates best_hits"""
        best_hits = {
            'HABJ36W02EXF44': [{
                'a': {'evalue': 0.0,
                      'subject_id': 'NZ_ABEH01000018_641736102',
                      'bit_score': 1005.0, 'percentage_id': 99.42,
                      'alg_length': 519},
                'b': {'subject_id': None,
                      'bit_score': -1}}],
            'HABJ36W02DLDSY': [{
                'a': {'evalue': 0.0,
                      'subject_id': 'NZ_ABEH01000005_641736102',
                      'bit_score': 959.0,
                      'percentage_id': 99.22, 'alg_length': 512},
                'b': {'subject_id': None, 'bit_score': -1}}]}

        parse_second_database(self.db2, best_hits, [0.30], [30])
        self.assertEquals(best_hits, {
            'HABJ36W02EXF44': [{
                'a': {'evalue': 0.0, 'bit_score': 1005.0,
                      'subject_id': 'NZ_ABEH01000018_641736102',
                      'alg_length': 519, 'percentage_id': 99.42},
                'b': {'evalue': 4.6, 'subject_id': 'NC_013421_646311947',
                      'bit_score': 40.1, 'percentage_id': 90.62,
                      'alg_length': 32}}],
            'HABJ36W02DLDSY': [{
                'a': {'evalue': 0.0, 'bit_score': 959.0,
                      'subject_id': 'NZ_ABEH01000005_641736102',
                      'alg_length': 512, 'percentage_id': 99.22},
                'b': {'evalue': 4e-133,
                      'subject_id': 'NZ_ACZD01000120_647000262',
                      'bit_score': 482.0, 'percentage_id': 88.79,
                      'alg_length': 455}}]})

    def test_process_results(self):
        """Check results are processed and summarized correctly"""
        best_hits = {
            'HABJ36W02EXF44': [{
                'a': {'evalue': 0.0, 'subject_id': 'NZ_ABEH01000018_641736102',
                      'bit_score': 1005.0, 'percentage_id': 99.42,
                      'alg_length': 519},
                'b': {'subject_id': None, 'bit_score': -1}}],
            'HABJ36W02DLDSY': [{
                'a': {'evalue': 0.0, 'subject_id': 'NZ_ABEH01000005_641736102',
                      'bit_score': 959.0, 'percentage_id': 99.22,
                      'alg_length': 512},
                'b': {'subject_id': None, 'bit_score': -1}}],
            'SAME-VALUES': [{
                'a': {'evalue': 0.0, 'subject_id': 'RESULT-A',
                      'bit_score': 959.0, 'percentage_id': 99.22,
                      'alg_length': 512},
                'b': {'evalue': 0.0, 'subject_id': 'RESULT-B',
                      'bit_score': 959.0, 'percentage_id': 99.22,
                      'alg_length': 512}}],
            'OTHER-BETTER': [{
                'a': {'evalue': 0.0, 'subject_id': 'RESULT-A',
                      'bit_score': 10.0, 'percentage_id': 10.0,
                      'alg_length': 10},
                'b': {'evalue': 0.0, 'subject_id': 'RESULT-B',
                      'bit_score': 959.0, 'percentage_id': 100,
                      'alg_length': 900}}]
        }

        out_results = process_results([0.80], [50], [0.30], [30], best_hits,
                                      self.base)
        # removing the sumout_resultsmary pointer so we don't need to test
        out_results[0].pop('summary_fh')
        self.assertEquals(out_results, [{
            'db_interest': 0, 'db_other': 1, 'db_seqs_counts': {
                'a': {'NZ_ABEH01000005_641736102': 1,
                      'RESULT-A': 1,
                      'NZ_ABEH01000018_641736102': 1},
                'b': {None: 0, 'RESULT-B': 2}},
            'perfect_interest': 2, 'equal': 1,
            'filename': 'p1_0-a1_50_p2_0-a2_30'}])

if __name__ == "__main__":
    main()
