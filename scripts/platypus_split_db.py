#!/usr/bin/env python
# File created on 25 Apr 2013
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011-2013, The Platypus Project"
__credits__ = ["Antonio Gonzalez Pena", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "0.0.8-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

from platypus.compare import (sequences_from_query, PlatypusParseError,
    PlatypusValueError)

from os import makedirs
from os.path import join
from cogent.parse.fasta import MinimalFastaParser

from qiime.util import parse_command_line_parameters, make_option

script_info = {}
script_info['brief_description'] = "Given a query, split a database in two: "+\
    "the part that matches the query and the part that doesn't"
script_info['script_description'] = "Split a tab-delimited taxonomy file and "+\
    "a FASTA-formatted file pair using a query to search for the matches with"+\
    " the taxonomy file. The end result are two files one compounded by all "+\
    "the matches to the query and the other with the things that didn't meet "+\
    "the described criteria. If a wildcard is used in the query string, it "+\
    "will not have a 'wildcard-like' meaning during the search process."
script_info['script_usage'] = [("Search for the string 'salmonella'","Passing a"
    "tab-delimited taxonomy file and a FASTA-formatted file with the "
    "corresponding sequences; search and split the data in sequences that match"
    " the string 'salmonella' for their taxonomy assignment and the ones that "
    "doesn't. The output is stored in a directory named databases.","%prog -t "
    "bacteria.contigs.txt -f bacteria.contigs.fna -q 'salmonella' -o databases")]
script_info['script_usage'] += [("Search for the ids in the file: seqs_id.txt",
    "Passing a tab-delimited taxonomy file and a FASTA-formatted file with the "
    "corresponding sequences; search and split the data in sequences that match "
    "the ids in the text file and the ones that doesn't. The output is stored in "
    "a directory named databases.","%prog -t bacteria.contigs.txt -f "
    "bacteria.contigs.fna -q seqs_id.txt -o databases")]
script_info['output_description']= "An output directory with two FASTA "+\
    "formatted files, one with the matches to the query (interest.fna) and "+\
    "one with the rest of the contents (rest.fna)."
script_info['required_options'] = [
    make_option( '-t', '--taxonomy_fp', type="existing_filepath", help='tab '
    'separated file with two columns: name/identifier of the sequence and then '
    'the taxonomy. Note: the sequence identifier is the longest string before a'
    ' space in the header of the sequence.'),
    make_option( '-f', '--fasta_fp', type="existing_filepath", help='path to a '
    'FASTA formatted file to split in interest and rest. Note: sequence '
    'identifiers must match the ones in the taxonomy file.'),
    make_option( '-o', '--output_fp', type="new_dirpath", help='output folder '
    'path where the results are stored')
]
script_info['optional_options'] = [
    make_option( '-q', '--query', type="string", help='the query used to split '
    'the database, for example: salmonella. The query should be an exact match, '
    'no wildcards, it can have spaces, and it is case insensitive.'),
    make_option( '-s', '--split_file', type="existing_filepath", help='the tab '
    'separated query file, where the first column is the id of the sequeces that '
    'you want to keep in your interest db.'),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    taxonomy_fp = opts.taxonomy_fp
    fasta_fp = opts.fasta_fp
    query = opts.query
    output_fp = opts.output_fp
    split_file = opts.split_file
    
    if (query==None and split_file==None) or (query!=None and split_file!=None):
        raise IOError, "You must specify one and only one between query: " +\
            "'%s' and split_file: '%s'" % (query, split_file)
    
    if query!=None:
        # query the taxonomy file for the required sequence identifiers
        try:
            interest_taxonomy = sequences_from_query(open(taxonomy_fp, 'U'), query)
        except (PlatypusValueError, PlatypusParseError), e:
            option_parser.error(e.message)

        if len(interest_taxonomy)==0:
            option_parser.error('The query could no results, try a different one')
    else:
        try:
            interest_taxonomy = { l.strip().split('\t')[0].strip(): '' for l in open(split_file, 'U') }
        except (PlatypusValueError, PlatypusParseError), e:
            option_parser.error(e.message)
    try:
       makedirs(opts.output_fp)
    except OSError:
        pass # hopefully dir already exists

    interest_fp = open(join(output_fp, 'interest.fna'), 'w')
    rest_fp = open(join(output_fp, 'rest.fna'), 'w')

    for full_name, seq in MinimalFastaParser(open(opts.fasta_fp,'U')):
         name = full_name.strip().split(' ')[0].strip()

         if name in interest_taxonomy:
            interest_fp.write(">%s\n%s\n" % (full_name, seq))
         else:
            rest_fp.write(">%s\n%s\n" % (full_name, seq))

    interest_fp.close()
    rest_fp.close()


if __name__ == "__main__":
    main()