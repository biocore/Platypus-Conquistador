# ----------------------------------------------------------------------------
# Copyright (c) 2015--, platypus development team.
#
# Distributed under the terms of the GPL License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, dirname, abspath
from shutil import rmtree
from tempfile import gettempdir
from unittest import TestCase, main
from click import BadParameter

from platypus.commands import split_db


class TestCommands(TestCase):
    def setUp(self):
        self.base = abspath(join(dirname(__file__), 'support_files'))
        self.tax_fp = join(self.base, 'small-bacteria.contigs.txt')
        self.seqs_fp = join(self.base, 'small-bacteria.contigs.fna')

        self.to_delete = []

    def tearDown(self):
        for d in self.to_delete:
            try:
                rmtree(d)
            except OSError:
                pass

    def test_split_db(self):
        temp_dir = gettempdir()
        self.to_delete.append(temp_dir)

        split_db(self.tax_fp, self.seqs_fp, 'Streptococcus', temp_dir, None)

        exp_interest = join(self.base, 'interest.fna')
        exp_rest = join(self.base, 'rest.fna')

        out_interest = join(temp_dir, 'interest.fna')
        out_rest = join(temp_dir, 'rest.fna')

        with open(exp_interest) as exp, open(out_interest) as out:
            self.assertItemsEqual(exp.readlines(), out.readlines())

        with open(exp_rest) as exp, open(out_rest) as out:
            self.assertItemsEqual(exp.readlines(), out.readlines())

    def test_split_db_no_results(self):
        with self.assertRaises(BadParameter):
            split_db(self.tax_fp, self.seqs_fp, ":L doesn't exist", 'output',
                     None)

    def test_split_db_split_fp(self):
        temp_dir = gettempdir()
        self.to_delete.append(temp_dir)

        split_fp = join(self.base, 'split_fp.txt')

        split_db(self.tax_fp, self.seqs_fp, None, temp_dir, split_fp)

        exp_interest = join(self.base, 'split-interest.fna')
        exp_rest = join(self.base, 'split-rest.fna')

        out_interest = join(temp_dir, 'interest.fna')
        out_rest = join(temp_dir, 'rest.fna')

        with open(exp_interest) as exp, open(out_interest) as out:
            self.assertItemsEqual(exp.readlines(), out.readlines())

        with open(exp_rest) as exp, open(out_rest) as out:
            self.assertItemsEqual(exp.readlines(), out.readlines())

if __name__ == '__main__':
    main()
