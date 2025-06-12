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

import ccdc.search
import ccdc.io


def parse_command_line_args():
    """Return the command line arguments.
    Parse the command line arguments and return them.
    :return: Arguments for the command."""
    parser = argparse.ArgumentParser(description='Validate a journal title against the CSD.')

    parser.add_argument('-j', '--journal', type=str,
                        help='The journal title.')

    args = parser.parse_args()

    return args


def main():
    """
    Execute the text numeric search specified by the parameters passed.
    """
    # Parse arguments
    args = parse_command_line_args()

    query = ccdc.search.TextNumericSearch()
    print(query.is_journal_valid(args.journal))

if __name__ == '__main__':
    main()
