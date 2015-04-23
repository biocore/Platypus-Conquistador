# ----------------------------------------------------------------------------
# Copyright (c) 2015--, platypus development team.
#
# Distributed under the terms of the BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import division
from itertools import product, izip
from collections import namedtuple
from copy import copy


_header = (('query', str),
           ('subject', str),
           ('percent_id', float),
           ('aln_length', int),
           ('mismatches', int),
           ('gapopenings', int),
           ('q_start', int),
           ('q_end', int),
           ('s_start', int),
           ('s_end', int),
           ('evalue', float),
           ('bitscore', float))

M9 = namedtuple('m9', [h[0] for h in _header])
M9_empty = [M9(**{h: None for h, _ in _header})]


def parse_m9(fp):
    """Parse m9 formatted tabular data

    Parameters
    ----------
    fp : file-like object
        A file pointer that contains the lines to parse. It is expected that
        these lines are in BLAST m9 format, or comparable. Any additional
        columns will be ignored

    Returns
    -------
    list of namedtuple
        The namedtuples describe each field of the m9 format per record.
    """
    res = []
    hits = []
    current_query = None
    start_of_record = False

    for line in fp:
        # Using the header detail from BLAST to differentiate records as this
        # allows us to get a correct count of query sequences from BLAST
        # results. We cannot do this from SortMeRNA
        if line.startswith('#'):
            if line.startswith('# Fields'):
                start_of_record = True
                if hits:
                    res.append((hits[0].query, hits))
                    hits = []

            elif line.startswith('# BLASTN') and start_of_record:
                res.append((None, copy(M9_empty)))

            continue

        # BLAST output contains 12 fields, SortMeRNA has 14. The order is the
        # same, and we only care about the 12 fields common with BLAST.
        parts = line.strip().split('\t')[:12]

        if len(parts) not in (12, 14):  # 12 is BLAST, 14 is SortMeRNA
            raise ValueError("Unexpected number of fields found")

        start_of_record = False
        hit = M9(**{h: c(v) for (h, c), v in zip(_header, parts)})

        # SortMeRNA doesn't have the header output to differentiate records
        if hit.query != current_query and hits:
            res.append((hits[0].query, hits))
            hits = []
            current_query = None
        hits.append(hit)
        current_query = hit.query

    if hits:
        res.append((hits[0].query, hits))
        hits = []
    elif start_of_record:
        res.append((None, copy(M9_empty)))

    return res


def parse_first_database(db, percentage_ids, alignment_lengths):
    """Find hits above a given threshold

    Parameters
    ----------
        db : file-like object
            file pointer to 1st database
        percentage_ids : iterable of ints
            Iterable with percentage ids
        alignment_lengths : iterable of ints
            Iterable with alignment length values

    Returns
    -------
        int
            total number of seqs in the db
        dict
            A dictionary of seqs and hits, of the form:
                {'seq_id':
                    [{
                        'a': {'evalue':%f, 'percentageId':%f, 'bitScore':%f,
                              'subjectId':%s, 'algLength':%f},
                              'evalue': float(h[evalue])},
                        'b': { 'subject_id': None, 'bit_score': -1 },
                       # One element for each combination of %id and aln length
                    ]
                }
    """
    # try blast parser object
    results = parse_m9(db)

    options = list(product(percentage_ids, alignment_lengths))

    best_hits = {}
    for total_queries, (query, hits) in enumerate(results, 1):
        if query is None:
            continue

        best_hits[query] = []
        for p, a in options:
            # best bit score
            bbs = 0
            result = None

            for h in hits:
                valid = (h.percent_id >= p and h.aln_length >= a and
                         h.bitscore > bbs)
                if valid:
                    result = {'a': {'subject_id': h.subject,
                                    'percentage_id': h.percent_id,
                                    'bit_score': h.bitscore,
                                    'alg_length': h.aln_length,
                                    'evalue': h.evalue},
                              'b': {'subject_id': None,
                                    'bit_score': -1}}
                    bbs = h.bitscore
            best_hits[query].append(result)

    return total_queries, best_hits


def parse_second_database(db, best_hits, percentage_ids_other,
                          alignment_lengths_other):
    """Parses 2nd database, only looking at successful hits of the 1st db

    Parameters
    ----------
        db : str
            Filename of the blast results against a database without first
        best_hits: dict
            A dictionary with the successful results from parse_first_database
        percentage_ids : iterable
            Iterable with percentage id values
        alignment_lengths : iterable
            Iterable with with alignment length values

    Notes
    -----
        There are no return values, the command modifies best_hits, mainly the
        'b' key.
    """
    results = parse_m9(db)

    # create function to return results
    for query, hits in results:
        if query is None:
            continue

        if query in best_hits:
            values = product(percentage_ids_other, alignment_lengths_other)
            for i, (p, a) in enumerate(values):
                if not best_hits[query][i]:
                    continue
                # best bit score
                bbs = 0
                result = None
                for h in hits:
                    valid = (h.percent_id >= p and h.aln_length >= a and
                             h.bitscore > bbs)
                    if valid:
                        result = {'subject_id': h.subject,
                                  'percentage_id': h.percent_id,
                                  'bit_score': h.bitscore,
                                  'alg_length': h.aln_length,
                                  'evalue': h.evalue}
                        bbs = h.bitscore
                if result:
                    best_hits[query][i]['b'] = result


def process_results(percentage_ids, alignment_lengths, percentage_ids_other,
                    alignment_lengths_other, best_hits):
    """Format the results into a summary dictionary

    Parameters
    ----------
    percentage_ids : iterable of ints
        An iterable of ints with the percentage identities for the reference
        database.
    alignment_lengths : iterable of ints
        An iterable of ints with the alignment lengths for the reference
        database.
    percentage_ids_other : iterable of ints
        An iterable of ints with the percentage identities for the 'other'
        database.
    alignment_lengths_other : iterable of ints
        An iterable of ints with the alignment lengths for the 'other'
        database.
    best_hits : dict
        A dictionary with the best hits found in the databases.

    Returns
    -------
    list of dicts
        List of dictionaries with the summarized results.
    """
    results = []

    iter_a = product(percentage_ids, alignment_lengths)
    iter_b = product(percentage_ids_other, alignment_lengths_other)
    for (perc_id_a, aln_len_a), (perc_id_b, aln_len_b) in izip(iter_a, iter_b):
        filename = "p1_%d-a1_%d_p2_%d-a2_%d" % (perc_id_a, aln_len_a,
                                                perc_id_b, aln_len_b)
        results.append({
            'filename': filename,
            'db_interest': 0,
            'db_other': 0,
            'perfect_interest': 0,
            'equal': 0,
            'summary': ['#SeqId\tFirst\tSecond'],
            'db_seqs_counts': {'a': {}, 'b': {}}})

    for seq_name, values in best_hits.items():
        seq_name = seq_name.split(' ')[0].strip()
        for i, vals in enumerate(values):
            if not vals:
                continue
            subject_id_a = vals['a']['subject_id']
            subject_id_b = vals['b']['subject_id']
            db_seqs_counts_a = results[i]['db_seqs_counts']['a']
            db_seqs_counts_b = results[i]['db_seqs_counts']['b']

            # validating duplicated results in the databases
            # do this step in a different script early in the pipeline
            if subject_id_a not in db_seqs_counts_a:
                db_seqs_counts_a[subject_id_a] = 0
                if subject_id_a == db_seqs_counts_b:
                    raise Warning("%s is in both databases" % subject_id_a)
            if subject_id_b not in db_seqs_counts_b:
                db_seqs_counts_b[subject_id_b] = 0
                if subject_id_b == db_seqs_counts_a:
                    raise Warning("%s is in both databases" % subject_id_b)

            # Comparing bit_scores to create outputs
            if vals['a']['bit_score'] == vals['b']['bit_score']:
                results[i]['equal'] += 1
                results[i]['summary'].append('%s\t%s\t%s' % (seq_name,
                                                             subject_id_a,
                                                             subject_id_b))
                db_seqs_counts_a[subject_id_a] += 1
                db_seqs_counts_b[subject_id_b] += 1
            elif vals['a']['bit_score'] > vals['b']['bit_score']:
                if not subject_id_b:
                    results[i]['perfect_interest'] += 1
                    results[i]['summary'].append('%s\t%s\t' % (seq_name,
                                                               subject_id_a))
                db_seqs_counts_a[subject_id_a] += 1
            else:
                results[i]['db_other'] += 1
                results[i]['summary'].append('%s\n\t%s' % (seq_name, ''))

                db_seqs_counts_b[subject_id_b] += 1

    return results
