#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2017-04-03: made available by the Cambridge Crystallographic Data Centre
#

import os
import argparse

import csd.pipeline_pilot_interface


#################################################################################################
def parse_command_line_args():
    """Return the command line arguments.

    Parse the command line arguments and return them."""

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-i', '--input', type=str, default='',
                        help='Input file name')
    parser.add_argument('-n', '--num_confs', type=int, default=25,
                        help='Maximum number of conformers [25]')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads [1]')
    parser.add_argument('-u', '--unusual_torsion', type=int, default=2,
                        help='Maximum unusual torsions [2]')
    parser.add_argument('-d', '--destination', type=str, default='',
                        help='Destination location; may be file path or folder')
    parser.add_argument('-s', '--superimpose', default=False, action='store_true',
                        help='Superimpose conformers onto reference? [Yes]')
    parser.add_argument('-o', '--output_directory', default=False, action='store_true',
                        help='Write each conformational ensemble for each input molecule to a separate file? [Yes]')

    args = parser.parse_args()

    return args


#################################################################################################


def main():
    """Run the main program."""

    args = parse_command_line_args()
    output_dir = None
    if args.destination != '':
        output_dir = args.destination
        # Make directory, in case it does not already exist
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

    engine = csd.pipeline_pilot_interface.CompoundConformerLibraryStorage(
        args.input, args.num_confs, args.unusual_torsion, args.threads, args.superimpose)
    conformers_hit_list = engine.generate_conformers()
    engine.write_report(conformers_hit_list, output_dir)

    if not args.output_directory:
        engine.merge_conformers_output(conformers_hit_list, output_dir)
    else:
        engine.split_conformers_output(conformers_hit_list, output_dir)


if __name__ == '__main__':
    main()
