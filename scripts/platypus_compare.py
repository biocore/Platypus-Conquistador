#!/usr/bin/env python
# File created on 04 Feb 2012
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011-2013, The Platypus Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "0.0.8-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

from os.path import join
from operator import itemgetter

from platypus.parse import (parse_first_database, parse_second_database,
    process_results)

from qiime.util import parse_command_line_parameters, make_option, create_dir

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
    make_option( '-i', '--input_path_interest', type="existing_filepath", help=
    'the blast results file of searching against the database of interest'),
    make_option( '-j', '--input_path_other', type="existing_filepath", help=
    'the blast results file of searching against the other database.')
]
script_info['optional_options'] = [
    make_option( '-o', '--output_dir', type="new_dirpath", help='the output '
    'directory [default: %default]', default='blast_results_compare'),
    make_option( '-p', '--percentage_ids', help='Min percentage id to be '
    'considered a valid result. This can be a list using commas. Default: '
    '%default",default="70"'),
    make_option( '-a', '--alignment_lengths', help='Min alignment length to be '
    'considered a valid result. Can be a comma-sepparated list. Default: '
    '%default', default="50"),
    make_option( '-P', '--percentage_ids_other', help='Min percentage id to be '
    'considered a valid result. This can be a list using commas. These values '
    'are used when you want to specify a different threshold for the other db. '
    'Default: %default', default=None),
    make_option( '-A', '--alignment_lengths_other', help='Min alignment length '
    'to be considered a valid result. This can be a list using commas. These '
    'values are  used when you want to specify a different threshold for the '
    'other db. Default: %default',default=None),
    make_option( '--hits_to_1st', action="store_true", help='Outputs the labels'
    ' and counts of the sequences being hit in the first database. Default: '
    '%default',default=None),
    make_option( '--hits_to_2nd', action='store_true', help='Outputs the labels'
    ' and counts of the sequences being hit in the second database. Default: '
    '%default',default=None)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # turn a comma-separated list of numbers into a 1-D list of integers
    list_of_ints = lambda string: map(int, string.split(','))

    input_path_interest = opts.input_path_interest
    input_path_other = opts.input_path_other
    percentage_ids = list_of_ints(opts.percentage_ids)
    alignment_lengths = list_of_ints(opts.alignment_lengths)
    hits_to_1st = opts.hits_to_1st
    hits_to_2nd = opts.hits_to_2nd
    percentage_ids_other = opts.percentage_ids_other
    alignment_lengths_other = opts.alignment_lengths_other
    output_dir = opts.output_dir

    db_a = open(input_path_interest,'U')
    db_b = open(input_path_other,'U')

    # try to create the output directory, if it exists, just continue
    create_dir(output_dir, False)

    # run some validations on the input parameters
    if percentage_ids_other:
        percentage_ids_other = list_of_ints(percentage_ids_other)
        if len(percentage_ids)!=len(percentage_ids_other):
            option_parser.error("The percentage values for both databases "
                "should be the same length: %s - %s" % (percentage_ids,
                percentage_ids_other))
    else:
        percentage_ids_other = percentage_ids

    if alignment_lengths_other:
        alignment_lengths_other = list_of_ints(alignment_lengths_other)

        if len(alignment_lengths)!=len(alignment_lengths_other):
            option_parser.error("The alignment length values for both databases"
                " should be the length : %s - %s" % (alignment_lengths,
                alignment_lengths_other))
    else:
       alignment_lengths_other = alignment_lengths

    # Process databases
    total_queries, best_hits = parse_first_database(db_a, percentage_ids,
        alignment_lengths)
    parse_second_database(db_b, best_hits, percentage_ids_other,
        alignment_lengths_other)

    # Parse results
    results = process_results(percentage_ids, alignment_lengths,
        percentage_ids_other, alignment_lengths_other, best_hits,
        input_path_interest, input_path_other)

    # Collating output and writing full results
    for i, item in enumerate(results):
        filename = join(output_dir, "summary_" + item['filename'] + ".txt")

        fd = open(filename,'w')
        fd.write('\n'.join(item['summary']))
        fd.close()

        if i==0:
            combined_results = []
            combined_results.append(['filename'])
            combined_results.append(['interestdb (%s)' % input_path_interest])
            combined_results.append(['other db (%s)' % input_path_other])
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

        # Printing count of hits to the db
        if hits_to_1st:
            s_hits = sorted(item['db_seqs_counts']['a'].items(),
                key=itemgetter(1), reverse=True)

            filename = join(output_dir, "hits_to_1st_db_" + item['filename'] +\
                ".txt")

            fd = open(filename,'w')
            fd.write('\n'.join(['%s\t%d' % (k,v) for k,v in s_hits if v!=0]))
            fd.close()

        if hits_to_2nd:
            s_hits = sorted(item['db_seqs_counts']['b'].items(),
                key=itemgetter(1), reverse=True)

            filename = join(output_dir, "hits_to_2nd_db_" + item['filename'] +\
                ".txt")

            fd = open(filename,'w')
            fd.write('\n'.join(['%s: %d' % (k,v) for k,v in s_hits if v!=0]))
            fd.close()

    #Printing collated results
    compiled_output_fd = open(join(output_dir,"compile_output.txt"), 'w')
    compiled_output_fd.write('\n'.join(['\t'.join(item)
        for item in combined_results]))
    compiled_output_fd.close()

    compiled_output_no_hits_fd = open(join(output_dir, 
        "compile_output_no_nohits.txt"), 'w')
    compiled_output_no_hits_fd.write('\n'.join(['\t'.join(item)
        for item in combined_results[:-1]]))
    compiled_output_no_hits_fd.close()


if __name__ == "__main__":
    main()
