#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2017-04-03: made available by the Cambridge Crystallographic Data Centre
#

import argparse

from csd.pipeline_pilot_interface import StructureSearcher


# Perform a structure search of the CSD. Structure(s) are supplied as a query file and are assumed
# to be format which the the MoleculeReader can handle (e.g. MOL2, SD, MOL, etc.). The type of query
# is supplied by command line, and should be 'substructure', 'exact', or 'similarity'. If type is
# 'similarity', then the similarity threshold must also be supplied, and hits will be found that
# exceed this level of similarity.
#
# Finally, the query itself produces a file containing the refcodes and the index of the query to
# which the refcode related, should multiple queries have been supplied. This can lead to the same
# refcode appearing multiple times, once per query to which it related. If the query type is
# similarity, then the similarity score will also be added to the output.
def parse_command_line_args():
    """Return the command line arguments.
    Parse the command line arguments and return them."""
    parser = argparse.ArgumentParser(description='Perform a structural search against the CSD.')

    parser.add_argument('-q', '--query_file', type=str,
                        help='File containing structures to search for.')
    parser.add_argument('-o', '--output_file', type=str, default='',
                        help='File to output')
    parser.add_argument('-p', '--preferences_file', type=str, default='',
                        help='File containing the settings to use for the search')
    parser.add_argument('-t', '--search_type', type=str, default='substructure',
                        help='Search type; valid values are "substructure", "exact", and "similarity"')
    parser.add_argument('-s', '--similarity_search', dest='similarity_threshold', type=float,
                        help='Threshold similarity score')
    parser.add_argument('-d', '--databases', type=str, default='',
                        help='The database(s) to perform the search across')

    args = parser.parse_args()

    return args


def main():
    # Parse arguments
    args = parse_command_line_args()

    structure_search = StructureSearcher(args.query_file, args.output_file, args.search_type)
    if args.search_type == 'similarity':
        structure_search.similarity_threshold = args.similarity_threshold

    if args.preferences_file != '':
        structure_search.read_settings(args.preferences_file)

    if args.databases != '':
        structure_search.set_databases(args.databases)

    structure_search.search()

if __name__ == '__main__':
    main()
