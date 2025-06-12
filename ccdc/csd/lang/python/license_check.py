#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2017-04-03: made available by the Cambridge Crystallographic Data Centre
#
"""Report on which features are licensed in the current environment."""

import sys
print("Python version:" + sys.version)

# Check for CSD (This uses an internal feature)
from ccdc.utilities import _private_importer
with _private_importer() as pi:
    pi.import_ccdc_module('UtilitiesLib')
if UtilitiesLib.licence_check():
    print('CSD license confirmed.')
else:
    print('Failure: CSD license not found.')

# Check for Discovery license
try:
    from ccdc.docking import Docker

    docker = Docker()
    print('CSD-Discovery license confirmed.')
except RuntimeError:
    print('Failure: CSD-Discovery not found.')


# Check for Materials license
try:
    from ccdc.crystal import PackingSimilarity

    PackingSimilarity()
    print('CSD-Materials license confirmed.')
except RuntimeError:
    print('Failure: CSD-Materials not found.')

# Check for cavity license; NOT NEEDED
#try:
#    from ccdc.cavity import Cavity, CavityDB
#    print('You can use Cavity features')
#except ImportError:
#    print('Cavity is only available to Collaborators (RPs)')
