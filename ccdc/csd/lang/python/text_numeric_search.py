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

from csd.pipeline_pilot_interface import TextNumericSearcher


def parse_command_line_args():
    """Return the command line arguments.
    Parse the command line arguments and return them.
    :return: Arguments for the command."""
    parser = argparse.ArgumentParser(description='Perform a text-numeric search against the CSD.')

    parser.add_argument('-q', '--query_file', type=str,
                        help='The file containing the JSON which describes the query.')
    parser.add_argument('-p', '--preferences_file', type=str, default='',
                        help='File containing the settings to use for the search')
    parser.add_argument('-o', '--output_file', type=str, default='',
                        help='File to output')
    parser.add_argument('-d', '--databases', type=str, default='',
                        help='The database(s) to perform the search across')
    # Used for debug, this allows this script to be run in a debug mode where the query is supplied in the code
    # and the hits are printed to the screen rather than written to disk.
    parser.add_argument('-t', dest='is_test', action="store_true",
                        help='Is in test mode; means standard query will be run and no output file produced.')
    parser.set_defaults(is_test=False)

    args = parser.parse_args()

    return args


def main():
    """
    Execute the text numeric search specified by the parameters passed.
    """
    # Parse arguments
    args = parse_command_line_args()

    searcher = TextNumericSearcher(args.output_file, args.query_file)

    if args.preferences_file != '':
        searcher.read_settings(args.preferences_file)

    if args.databases != '':
        searcher.set_databases(args.databases)

    searcher.search()

if __name__ == '__main__':
    main()
