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

from csd.pipeline_pilot_interface import JoinAttributes


# Join attributes for the CSD crystals, entries, or molecules found for a list of refcodes. Refcodes as supplied a
# GCD file to the entry reader, and the list of entries which this file refers to are iterated over. The attributes
# to be retrieved as supplied as a comma separated list.
#
# The data read is written to a CSV file, unless an attribute is asked for which will not
def parse_command_line_args():
    """Return the command line arguments.
    Parse the command line arguments and return them.
    :rtype: ArgumentParseObject"""

    parser = argparse.ArgumentParser(
        description='Gather attributes for the refcodes found in the list file. Attributes are supplied as a list '
                    'argument.')

    parser.add_argument('-i', '--input_file', type=str,
                        help='<Required> List file containing CSD identities to look up', required=True)
    parser.add_argument('-a', '--attributes', type=str,
                        help='<Required> Set of attributes to recall for each refcode, expressed as a comma separated '
                             'list.',
                        required=True)
    parser.add_argument('-o', '--output_file', type=str,
                        help='<Required> File to output', required=True)
    parser.add_argument('-t', '--join_type', type=str,
                        help='<Required> Type of attributes to be joined; should be Entry, Molecule, or Crystal',
                        choices=['Entry', 'Molecule', 'Crystal'],
                        required=True)
    parser.add_argument('-d', '--database', type=str, default='CSD',
                        help='The database(s) to perform the search across')

    args = parser.parse_args()

    return args


def main():

    args = parse_command_line_args()

    attributes = args.attributes.split(',')
    joiner = JoinAttributes(args.input_file, args.database, args.join_type, attributes)
    joiner.write_attribute_file(args.output_file)


if __name__ == '__main__':

    main()
