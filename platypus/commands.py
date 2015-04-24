# ----------------------------------------------------------------------------
# Copyright (c) 2015--, platypus development team.
#
# Distributed under the terms of the BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import division

from os.path import join, basename
from operator import itemgetter

from click import BadParameter
from skbio.util import create_dir
from skbio import read

from platypus.compare import (
    sequences_from_query, PlatypusParseError, PlatypusValueError)
from platypus.parse import (parse_first_database, parse_second_database,
                            process_results)


def compare(interest_fp, other_fp, output_dir='blast-results-compare',
            interest_pcts=None, interest_alg_lens=None, other_pcts=None,
            other_alg_lens=None, hits_to_first=False, hits_to_second=False):
    """Compare two databases and write the outputs

    Parameters
    ----------
    interest_fp : str
        BLAST results when searching against the database of interest.
    other_fp : str
        BLAST results when searching against the other database.
    output_dir : str, optional
        Name of the output file path.
    interest_pcts : list, optional
        Minimum percentage identity to be considered as a valid result in the
        interest database search results. If None is passed, it defaults to
        `[70]`.
    interest_alg_lens : list, optional
        Minimum alginment length to be considered a valid result in the
        interest database search results. If None is passed, it defaults to
        `[50]`.
    other_pcts : list, optional
        Minimum percentage identity to be considered as a valid result in the
        other database search results. If None is passed, it defaults to
        `[70]`.
    other_lengths : list, optional
        Minimum alginment length to be considered a valid result in the
        other database search results. If None is passed, it defaults to
        `[50]`.
    hits_to_first : bool, optional defaults to False
        Outputs the labels and counts of the sequences being hit in the first
        database.
    hits_to_second : bool, optional defaults to False
        Outputs the labels and counts of the sequences being hit in the second
        database.

    Raises
    ------
    click.BadParameter
        If the `interest_pcts` and the `other_pcts` lists are of different
        length.
        If the `interest_alg_lens` and the `other_alg_lens` lists are of
        different length.
    """

    if interest_pcts is None:
        interest_pcts = [70]
    if interest_alg_lens is None:
        interest_alg_lens = [50]

    db_a = open(interest_fp, 'U')
    db_b = open(other_fp, 'U')

    # try to create the output directory, if it exists, just continue
    create_dir(output_dir, False)

    # run some validations on the input parameters
    if other_pcts:
        if len(interest_pcts) != len(other_pcts):
            raise BadParameter("The percentage values for both databases "
                               "should be the same length: %s - %s" %
                               (interest_pcts, other_pcts))
    else:
        other_pcts = interest_pcts

    if other_alg_lens:
        if len(interest_alg_lens) != len(other_alg_lens):
            raise BadParameter("The alignment length values for both databases"
                               " should be the length : %s - %s" %
                               (interest_alg_lens, other_alg_lens))
    else:
        other_alg_lens = interest_alg_lens

    # process databases
    total_queries, best_hits = parse_first_database(db_a, interest_pcts,
                                                    interest_alg_lens)
    parse_second_database(db_b, best_hits, other_pcts,
                          other_alg_lens)

    # parse results
    results = process_results(interest_pcts, interest_alg_lens,
                              other_pcts, other_alg_lens, best_hits)

    # Collating output and writing full results
    for i, item in enumerate(results):
        filename = join(output_dir, "summary_" + item['filename'] + ".txt")

        with open(filename, 'w') as fd:
            fd.write('\n'.join(item['summary']))

        if i == 0:
            combined_results = []
            combined_results.append(['filename'])
            combined_results.append(['interest db (%s)' %
                                     basename(interest_fp)])
            combined_results.append(['other db (%s)' % basename(other_fp)])
            combined_results.append(['only interest'])
            combined_results.append(['both dbs'])
            combined_results.append(['no hits in interest db'])

        no_hits = total_queries - item['db_interest'] - item['db_other'] - \
            item['perfect_interest'] - item['equal']

        combined_results[0].append(item['filename'])
        combined_results[1].append(str(item['db_interest']))
        combined_results[2].append(str(item['db_other']))
        combined_results[3].append(str(item['perfect_interest']))
        combined_results[4].append(str(item['equal']))
        combined_results[5].append(str(no_hits))

        # tiny helper function to save hits files
        def save_hits(data, name):

            s_hits = sorted(data, key=itemgetter(1), reverse=True)
            filename = join(output_dir, name)
            with open(filename, 'w') as fd:
                fd.write('\n'.join(['%s\t%d' % (k, v)
                                    for k, v in s_hits if v != 0]))

        if hits_to_first:
            save_hits(item['db_seqs_counts']['a'].items(),
                      "hits_to_first_db_%s.txt" % item['filename'])

        if hits_to_second:
            save_hits(item['db_seqs_counts']['b'].items(),
                      "hits_to_second_db_%s.txt" % item['filename'])

    # saving collated results
    with open(join(output_dir, "compile_output.txt"), 'w') as compiled_output:
        compiled_output.write('\n'.join(['\t'.join(item)
                                         for item in combined_results]))

    fn = join(output_dir, "compile_output_no_nohits.txt")
    with open(fn, 'w') as fd:
        fd.write('\n'.join(['\t'.join(item)
                            for item in combined_results[:-1]]))


def split_db(tax_fp, seqs_fp, query, output_fp, split_fp):
    """Split a database in parts that match a query and parts that don't

    Parameters
    ----------
    tax_fp : str
        Tab-delimited file with two columns, name/identifier of the sequence
        and the taxonomy. The sequence identifier is the longest string before
        a space in the header of the sequence.
    seqs_fp : str
        Path to a FASTA formatted file to split in interest and rest. Note:
        sequence identifiers must match the ones in the taxonomy file.
    query : str
        The query used to split the database, for example: salmonella. The
        query should be an exact match, no wild cards, it can have spaces, and
        it is case insensitive
    output_fp : str
        Output folder path where the results are stored.
    split_fp : str
        The tab delimited query file, where each line is a different sequence
        and the first column is the sequence id.

    Raises
    ------
    BadParameter
        If the Taxonomy file is empty.
        If the query you passed retrieved no results.
    """

    if query is not None:
        # query the taxonomy file for the required sequence identifiers
        try:
            interest_taxonomy = sequences_from_query(open(tax_fp, 'U'),
                                                     query)
        except (PlatypusValueError, PlatypusParseError), e:
            raise BadParameter(e.message)

        if len(interest_taxonomy) == 0:
            raise BadParameter('The query could not retrieve any results, try '
                               'a different one.')
    else:
        interest_taxonomy = {l.strip().split('\t')[0].strip(): ''
                             for l in open(split_fp, 'U')}
        if not interest_taxonomy:
            raise BadParameter('The split_fp is empty!')

    create_dir(output_fp, False)

    interest_fp = open(join(output_fp, 'interest.fna'), 'w')
    rest_fp = open(join(output_fp, 'rest.fna'), 'w')

    for record in read(seqs_fp, format='fasta'):
        full_name = record.id
        seq = record.sequence

        name = full_name.strip().split(' ')[0].strip()

        if name in interest_taxonomy:
            interest_fp.write(">%s\n%s\n" % (full_name, seq))
        else:
            rest_fp.write(">%s\n%s\n" % (full_name, seq))

    interest_fp.close()
    rest_fp.close()
