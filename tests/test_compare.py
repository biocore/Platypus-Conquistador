#!/usr/bin/env python
# File created on 17 Jul 2013
from __future__ import division

from platypus.compare import (
    sequences_from_query, PlatypusParseError, PlatypusValueError)

from os.path import dirname, join, abspath
from unittest import TestCase, main


class TopLevelTests(TestCase):

    def setUp(self):
        # get the path of the tests, the support files should be in that dir
        test_dir = dirname(abspath(__file__))
        self.taxonomy_fp = join(test_dir, 'support_files/taxonomy_file.txt')
        self.taxonomy_lines = TAXONOMY_LINES
        self.broken_taxonomy_lines = BROKEN_TAXONOMY_LINES

    def test_sequences_from_query_file(self):
        """Check correct parsing when input is a file path"""
        out_tax_dict = sequences_from_query(self.taxonomy_fp, 'Beggiatoa')
        self.assertEquals(out_tax_dict, {
            'NZ_ABBZ01000843|640963011': 'VibrioBeggiatoa sp. PS',
            'NZ_ABBZ01000613|640963011': 'VibrioBeggiatoa sp. PS',
            'NZ_ABBZ01002042|640963011': 'VibrioBeggiatoa sp. PS',
            'NZ_ABBY01000503|640963012': 'VibrioBeggiatoa sp. SS',
            'NZ_ABBZ01000140|640963011': 'VibrioBeggiatoa sp. PS',
            'NZ_ABBY01000221|640963012': 'VibrioBeggiatoa sp. SS',
            'NZ_ABBZ01006278|640963011': 'VibrioBeggiatoa sp. PS',
            'NZ_ABBZ01005870|640963011': 'VibrioBeggiatoa sp. PS'})

        out_tax_dict = sequences_from_query(self.taxonomy_fp, 'Rumba')
        self.assertEquals(out_tax_dict, {})

    def test_sequences_from_query_string(self):
        """Check correct parsing when input is a string"""
        out_tax_dict = sequences_from_query(self.taxonomy_lines, 'Beggiatoa')
        self.assertEquals(out_tax_dict, {
            'NZ_ABBZ01000843|640963011': 'VibrioBeggiatoa sp. PS',
            'NZ_ABBZ01000613|640963011': 'VibrioBeggiatoa sp. PS',
            'NZ_ABBZ01002042|640963011': 'VibrioBeggiatoa sp. PS',
            'NZ_ABBY01000503|640963012': 'VibrioBeggiatoa sp. SS',
            'NZ_ABBZ01000140|640963011': 'VibrioBeggiatoa sp. PS',
            'NZ_ABBY01000221|640963012': 'VibrioBeggiatoa sp. SS',
            'NZ_ABBZ01006278|640963011': 'VibrioBeggiatoa sp. PS',
            'NZ_ABBZ01005870|640963011': 'VibrioBeggiatoa sp. PS'})

        out_tax_dict = sequences_from_query(self.taxonomy_lines, 'Rumba')
        self.assertEquals(out_tax_dict, {})

    def test_sequences_from_query_exceptions(self):
        """Check exceptions are raised accordingly"""
        # re-format the file to make it "BOOM" delimited
        self.assertRaises(
            PlatypusParseError, sequences_from_query,
            self.taxonomy_lines.replace('\t', 'BOOOM'), 'Beggiatoa')

        # repeated sequence identifiers
        self.assertRaises(
            PlatypusValueError, sequences_from_query,
            self.broken_taxonomy_lines, 'parahaemolyticus')


TAXONOMY_LINES = """NZ_AAOS01000254|638341243\tVibrioYersinia pestis bv Orientalis IP275
NZ_AAOS01000267|638341243\tVibrioYersinia pestis bv Orientalis IP275
NZ_AAVT01000001|639857034\tVibrioMarine gamma proteobacterium sp. HTCC2143
NZ_AAWG01000010|640963003\tVibrioVibrio cholerae 623-39
NZ_AAWG01000117|640963003\tVibrioVibrio cholerae 623-39
NZ_AAWQ01000004|640963009\tVibrioVibrio parahaemolyticus AQ3810
NZ_AAWQ01000613|640963009\tVibrioVibrio parahaemolyticus AQ3810
NZ_AAWQ01001040|640963009\tVibrioVibrio parahaemolyticus AQ3810
NZ_ABBZ01000140|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBZ01000613|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBZ01000843|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBZ01002042|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBZ01005870|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBZ01006278|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBY01000221|640963012\tVibrioBeggiatoa sp. SS
NZ_ABBY01000503|640963012\tVibrioBeggiatoa sp. SS"""

BROKEN_TAXONOMY_LINES = """NZ_AAOS01000254|638341243\tVibrioYersinia pestis bv Orientalis IP275
NZ_AAOS01000267|638341243\tVibrioYersinia pestis bv Orientalis IP275
NZ_AAVT01000001|639857034\tVibrioMarine gamma proteobacterium sp. HTCC2143
NZ_AAWG01000010|640963003\tVibrioVibrio cholerae 623-39
NZ_AAWG01000117|640963003\tVibrioVibrio cholerae 623-39
NZ_AAWQ01001040|640963009\tVibrioVibrio parahaemolyticus AQ3810
NZ_AAWQ01001040|640963009\tVibrioVibrio parahaemolyticus AQ3810
NZ_AAWQ01001040|640963009\tVibrioVibrio parahaemolyticus AQ3810
NZ_ABBZ01000140|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBZ01000613|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBZ01000843|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBZ01002042|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBZ01005870|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBZ01006278|640963011\tVibrioBeggiatoa sp. PS
NZ_ABBY01000221|640963012\tVibrioBeggiatoa sp. SS
NZ_ABBY01000503|640963012\tVibrioBeggiatoa sp. SS"""

if __name__ == "__main__":
    main()
