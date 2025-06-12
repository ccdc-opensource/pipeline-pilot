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

import csd.pipeline_pilot_interface


def parse_command_line_args():
    """Return the command line arguments.
    Parse the command line arguments and return them."""

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-q', '--query', type=str, default='',
                        help='Query file')
    parser.add_argument('-a', '--actives', type=str, default='',
                        help='Actives set')
    parser.add_argument('-d', '--decoys', type=str, default='',
                        help='Decoys set')
    parser.add_argument('-n', '--num_confs', type=int, default=200,
                        help='Maximum number of conformers [200]')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads [1]')
    parser.add_argument('-o', '--output_dir', type=str, default='Screen_data',
                        help='Output directory')

    args = parser.parse_args()

    return args

###############################################################################

def main():
    # read arguments
    args = parse_command_line_args()

    # Construct the validator.
    validator = csd.pipeline_pilot_interface.ValidateScreening(
        args.query, args.actives, args.decoys, args.output_dir, args.threads, args.num_confs)

    actives_scores = validator.screen_actives()
    decoys_scores = validator.screen_decoys()

    all_data = actives_scores
    all_data.extend(decoys_scores)
    results = sorted(all_data)

    validator.write_scores(results)

if __name__ == '__main__':
    main()
