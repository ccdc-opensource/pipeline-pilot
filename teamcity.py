"""
This script is intended to be run by TeamCity
"""

from __future__ import division, absolute_import, print_function

import sys
import os
import subprocess

this_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    if os.environ.get('DISTRIBUTION_DIRECTORY') is None:
        os.environ['DISTRIBUTION_DIRECTORY'] = "\\\\artefact\\windows\\cppbuilds"

    for script in ['create_packages', 'test']:
        try:
            print("##teamcity[blockOpened name='" + script +
                  "' description='" + script + "']")
            subprocess.check_call([
                sys.executable,
                os.path.join(this_dir, script + '.py')
            ], cwd=this_dir)
        finally:
            print("##teamcity[blockClosed name='" + script + "']")
