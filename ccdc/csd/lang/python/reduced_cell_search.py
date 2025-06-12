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

from ccdc.crystal import CellAngles
from ccdc.crystal import CellLengths

from csd.pipeline_pilot_interface import ReducedCellSearcher


def parse_command_line_args():
    """Return the command line arguments.
    Parse the command line arguments and return them."""
    parser = argparse.ArgumentParser(description='Perform a reduced cell search against the CSD.')

    parser.add_argument('-a', '--a', type=float, default=-1,
                        help='The A length.')
    parser.add_argument('-b', '--b', type=float, default=-1,
                        help='The B length.')
    parser.add_argument('-c', '--c', type=float, default=-1,
                        help='The C length.')
    parser.add_argument('-alpha', '--alpha', type=float, default=-1,
                        help='The alpha angle.')
    parser.add_argument('-beta', '--beta', type=float, default=-1,
                        help='The beta angle.')
    parser.add_argument('-gamma', '--gamma', type=float, default=-1,
                        help='The gamma angle.')
    parser.add_argument('-l', '--lattice_centring', type=str, default='',
                        help='The lattice centring.')
    parser.add_argument('-p', '--preferences_file', type=str, default='',
                        help='File containing the settings to use for the search')
    parser.add_argument('-o', '--output_file', type=str, default='',
                        help='File to output')
    parser.add_argument('-d', '--databases', type=str, default='',
                        help='The database(s) to perform the search across')

    args = parser.parse_args()

    return args


def main():
    args = parse_command_line_args()
    # Construct cell angles, and cell lengths
    if args.a > -1 and args.b > -1 and args.c > -1:
        cell_lengths = CellLengths(args.a, args.b, args.c)
    else:
        cell_lengths = None

    if args.alpha > -1 and args.beta > -1 and args.gamma > -1:
        cell_angles = CellAngles(args.alpha, args.beta, args.gamma)
    else:
        cell_angles = None

    if args.lattice_centring != "":
        lattice_centring = args.lattice_centring
    else:
        lattice_centring = None

    if cell_lengths is not None and cell_angles is not None and lattice_centring is not None:
        reduced_cell_search = ReducedCellSearcher(args.output_file, cell_angles, cell_lengths, lattice_centring)

        if args.preferences_file != '':
            reduced_cell_search.read_settings(args.preferences_file)

        if args.databases != '':
            reduced_cell_search.set_databases(args.databases)

        reduced_cell_search.search()

if __name__ == '__main__':
    main()
