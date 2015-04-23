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
        base = join(dirname(__file__), 'support_files')
        self.db1 = open(join(base, 'first_db.txt'), 'r')
        self.db2 = open(join(base, 'second_db.txt'), 'r')
        self.blasttest = open(join(base, 'blast_test_output.txt'))
        self.smrtest = open(join(base, 'sortmerna_test_output.txt'))
        self.m9_empty = M9(*[None] * 12)

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
               ('4502804.3.fna_9', fna_9_hits)]

        obs = parse_m9(self.blasttest)
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

        out_results = process_results([0.80], [50], [0.30], [30], best_hits)
        self.assertEquals(out_results, [{
            'db_interest': 0, 'db_other': 1, 'db_seqs_counts': {
                'a': {'NZ_ABEH01000005_641736102': 1,
                      'RESULT-A': 1,
                      'NZ_ABEH01000018_641736102': 1},
                'b': {None: 0, 'RESULT-B': 2}},
            'perfect_interest': 2, 'equal': 1, 'summary':
            ['#SeqId\tFirst\tSecond',
             'HABJ36W02EXF44\tNZ_ABEH01000018_641736102\t',
             'OTHER-BETTER\n\t',
             'SAME-VALUES\tRESULT-A\tRESULT-B',
             'HABJ36W02DLDSY\tNZ_ABEH01000005_641736102\t'],
            'filename': 'p1_0-a1_50_p2_0-a2_30'}])

if __name__ == "__main__":
    main()
