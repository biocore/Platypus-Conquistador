#!/usr/bin/env python
# File created on 13 Jul 2013
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011-2013, The Platypus Project"
__credits__ = ["Antonio Gonzalez Pena", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "0.0.8-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"


class PlatypusError(Exception):
    """Base exception for errors in the Platypus module"""
    pass


class PlatypusParseError(PlatypusError):
    """Base exception for incorrectly formatted files in the Platypus module"""
    pass


class PlatypusValueError(PlatypusError):
    """Inconsistency of inputs in the Platypus module"""
    pass


def sequences_from_query(taxonomy, query):
    """Parses and searches for a query in the contents of taxonomy

    Inputs:
    taxonomy: file descriptor, StringIO object or lines that contain
    tab-delimited values of sequence identifier to taxonomy assignment
    query: string to search for in the taxonomy assignments

    Output:
    interest_taxonomy: 1-D dictionary where the keys are the sequence
    identifiers and the values are the taxonomy assignments

    Raises:
    PlatypusParseError if the input is not a two column tab-delimited file
    PlatypusValueError if the input contains repeated sequence identifiers
    """
    # Note this function could be greatly benefited from a C extension

    interest_taxonomy = {}
    try:    # file path
        fd = open(taxonomy, 'U')
    except IOError:  # string with lines, split on new lines
        fd = taxonomy.split('\n')
    except TypeError:  # open file descriptor or StringIO object
        fd = taxonomy

    # read each line searching for matches to the query
    for line in fd:
        try:
            sequence_identifier, taxa_name = line.strip().split('\t')
        except ValueError:
            raise PlatypusParseError(
                "Taxonomy file/string is not tab delimited")

        sequence_identifier = sequence_identifier.strip()

        if sequence_identifier in interest_taxonomy:
            # if possible close the file descriptor before leaving the function
            try:
                fd.close()
            except AttributeError:
                pass

            raise PlatypusValueError(
                "There are duplicated entries in the taxonomy file (%s)." %
                sequence_identifier)

        elif query.lower() in taxa_name.lower():
            interest_taxonomy[sequence_identifier] = taxa_name.strip()

    # not all input types are file descriptors
    try:
        fd.close()
    except AttributeError:
        pass

    return interest_taxonomy
