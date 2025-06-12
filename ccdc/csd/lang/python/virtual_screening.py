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
import shutil
import argparse

from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.conformer import ConformerGenerator, ConformerSettings
from ccdc.screening import Screener

import csd.pipeline_pilot_interface

def setup_screener(output_dir):
    """Return the settings object for the ligand screener."""
    settings = Screener.Settings()
    param_dir = settings.parameter_directory
    new_dir = os.path.join(os.getcwd(), 'parameter_files')
    try:
        # copy the parameter_files directory from the default location into the working directory.
        shutil.copytree(param_dir, new_dir)
    except:
        # if there is already a parameter_files directory, then this will be used for the calculation.
        pass
    settings.parameter_directory = new_dir
    settings.output_directory = os.path.join(os.getcwd(), output_dir)
    if os.path.exists(settings.output_directory):
        shutil.rmtree(settings.output_directory)
    os.makedirs(settings.output_directory)
    return settings


def standardise(mol):
    """Return standardised molecule."""
    #    mol.assign_unknown_bond_types(mol)
    mol.standardise_aromatic_bonds()
    mol.standardise_delocalised_bonds()
    return mol


def generate_confs(mol, nc, pt, nt):
    """Return the conformer hit lists."""
    settings = ConformerSettings()
    settings.max_conformers = nc
    gen = ConformerGenerator(settings, nthreads=nt)
    standardised_mol = standardise(mol)
    conf_hit_lists = gen.generate(standardised_mol)
    confs = [[c.molecule for c in conf_hit_lists]]
    return confs


def screen_molecules(screener, mols_to_screen, nc, pt, nt, output_name):
    """Run the ligand screener and write out the screened conformations.
    Return the list of ranked scores"""

    screen_set = [m for m in MoleculeReader(mols_to_screen)]  # Read the molecules to screen
    scores = []
    with MoleculeWriter(output_name) as molwriter:
        for mol in screen_set:

            # If nconformers=0 the input conformation is used, otherwise the CSD-driven conformer
            # generator is called and a number of conformations equal to nconformers is generated.
            if nc == 0:
                confs = [[standardise(m)]]
            else:
                try:
                    confs = generate_confs(mol, nc, pt, nt)
                except RuntimeError:
                    confs = [[standardise(m)]]
            res = screener.screen(confs)  # Screening step
            scores.extend([(r.score, r.identifier) for r in res])
            # Write results
            for r in res:
                molwriter.write(r.molecule)

    return sorted(scores)


def write_scores(results, output):
    """Write results to an output file."""
    with open(output, 'w') as f:
        for r in results:
            score, identifier = r
            f.write('%.3f, %s\n' % (score, identifier))


def parse_command_line_args():
    """Return the command line arguments.

    Parse the command line arguments and return them."""

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-q', '--query', type=str, default='',
                        help='Query file')
    parser.add_argument('-s', '--screen_set', type=str, default='',
                        help='Library of molecules to screen')
    parser.add_argument('-n', '--num_confs', type=int, default=200,
                        help='Maximum number of conformers [200]')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads [1]')
    parser.add_argument('-o', '--output_dir', type=str, default='Screen_data',
                        help='Output directory')

    args = parser.parse_args()

    return args


def main():
    # read arguments
    args = parse_command_line_args()

    virtual_screening = csd.pipeline_pilot_interface.VirtualScreening(
        args.query, args.screen_set, args.threads, args.num_confs, args.output_dir)

    scores = virtual_screening.screen_molecules()

    results = sorted(scores)

    virtual_screening.write_scores(results)


if __name__ == '__main__':
    main()
