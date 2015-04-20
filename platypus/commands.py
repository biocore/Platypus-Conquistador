# ----------------------------------------------------------------------------
# Copyright (c) 2015--, platypus development team.
#
# Distributed under the terms of the GPL License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
from operator import itemgetter
from click import BadParameter
from skbio.util import create_dir

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
                              other_pcts, other_alg_lens, best_hits,
                              interest_fp, other_fp)

    # Collating output and writing full results
    for i, item in enumerate(results):
        filename = join(output_dir, "summary_" + item['filename'] + ".txt")

        with open(filename, 'w') as fd:
            fd.write('\n'.join(item['summary']))

        if i == 0:
            combined_results = []
            combined_results.append(['filename'])
            combined_results.append(['interestdb (%s)' % interest_fp])
            combined_results.append(['other db (%s)' % other_fp])
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

        # saving count of hits to the db
        if hits_to_first:
            s_hits = sorted(item['db_seqs_counts']['a'].items(),
                            key=itemgetter(1), reverse=True)

            filename = join(output_dir, "hits_to_first_db_" + item['filename']
                            + ".txt")

            with open(filename, 'w') as fd:
                fd.write('\n'.join(['%s\t%d' % (k, v)
                                    for k, v in s_hits if v != 0]))

        if hits_to_second:
            s_hits = sorted(item['db_seqs_counts']['b'].items(),
                            key=itemgetter(1), reverse=True)

            filename = join(output_dir, "hits_to_second_db_" + item['filename']
                            + ".txt")

            with open(filename, 'w') as fd:
                fd.write('\n'.join(['%s: %d' %
                                    (k, v) for k, v in s_hits if v != 0]))

    # saving collated results
    with open(join(output_dir, "compile_output.txt"), 'w') as compiled_output:
        compiled_output.write('\n'.join(['\t'.join(item)
                                         for item in combined_results]))

    fn = join(output_dir, "compile_output_no_nohits.txt")
    with open(fn, 'w') as fd:
        fd.write('\n'.join(['\t'.join(item)
                            for item in combined_results[:-1]]))

