#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2017-04-03: made available by the Cambridge Crystallographic Data Centre
#

from pathlib import Path
import pytest
import os
import sys
import unittest
from unittest.mock import MagicMock

csd_package_dir = str(Path(__file__).resolve().parent.parent)
print(csd_package_dir)
sys.path.insert(0, csd_package_dir)
from csd.pipeline_pilot_interface import (
        CompoundConformerLibraryStorage, ReducedCellSearcher,
        TextNumericSearcher, PipelinePilotBase, StructureSearcher,
        ValidateScreening, VirtualScreening)

from ccdc.io import MoleculeReader
from ccdc.utilities import _find_test_file, _test_output_dir
from ccdc.io import CrystalReader
from ccdc.utilities import _private_importer
with _private_importer() as pi:
    pi.import_ccdc_module('UtilitiesLib')

csd = MoleculeReader('CSD')

class TestsStructureSearcher(unittest.TestCase):

    def testTelemetry(self):
        UtilitiesLib.ccdc_pipeline_pilot_telemetry = MagicMock()
        aspirin_file = _find_test_file(__file__, 'aspirin.mol2')
        StructureSearcher(aspirin_file, "junk_output_dir", 'similarity')
        self.assertEqual(1, UtilitiesLib.ccdc_pipeline_pilot_telemetry.call_count)


class TestCompoundConformerLibraryStorage(unittest.TestCase):
    def setUp(self):
        self.aspirin_file = _find_test_file(__file__, 'aspirin.mol2')

    def testStandardiseAddHydrogens(self):
        aabhtz = csd.molecule('AABHTZ')
        aabhtz.remove_hydrogens()
        CompoundConformerLibraryStorage.standardise(aabhtz, False)
        self.assertEqual(0, len([a for a in aabhtz.atoms if a.atomic_symbol == 'H']))
        CompoundConformerLibraryStorage.standardise(aabhtz, True)
        self.assertEqual(12, len([a for a in aabhtz.atoms if a.atomic_symbol == 'H']))

    @unittest.skip('Requires Mogul data')
    def testGenerateConformers(self):
        engine = CompoundConformerLibraryStorage(self.aspirin_file, 5, 10, 1, True)
        confs = engine.generate_conformers()
        self.assertEqual(5, len(confs[0]))

    def testTelemetry(self):
        UtilitiesLib.ccdc_pipeline_pilot_telemetry = MagicMock()
        CompoundConformerLibraryStorage(self.aspirin_file, 5, 10, 1, True)
        self.assertEqual(1, UtilitiesLib.ccdc_pipeline_pilot_telemetry.call_count)

class TestReducedCellSearcher(unittest.TestCase):
    def setUp(self):
        self.out_directory = _test_output_dir()
        self.out_filename = os.path.join(self.out_directory, 'out.dat')
        reader = CrystalReader('CSD')
        self.aacmhx10 = reader.crystal('AACMHX10')
        self.cell_angles = self.aacmhx10.reduced_cell.cell_angles
        self.cell_lengths = self.aacmhx10.reduced_cell.cell_lengths
        self.lattice_centring = self.aacmhx10.lattice_centring

    def testSearch(self):
        self.assertFalse(os.path.exists(self.out_filename))
        reduced_cell_search = ReducedCellSearcher(
            self.out_filename, self.cell_angles, self.cell_lengths, self.lattice_centring)
        reduced_cell_search.search()
        self.assertTrue(os.path.exists(self.out_filename))

    def testTelemetry(self):
        UtilitiesLib.ccdc_pipeline_pilot_telemetry = MagicMock()
        ReducedCellSearcher(
            self.out_filename, self.cell_angles, self.cell_lengths, self.lattice_centring)
        self.assertEqual(1, UtilitiesLib.ccdc_pipeline_pilot_telemetry.call_count)

class TestTextNumericSearcher(unittest.TestCase):
    def setUp(self):
        self.out_directory = _test_output_dir()
        self.out_filename = os.path.join(self.out_directory, 'out.dat')
        self.query_filename = os.path.join(self.out_directory, 'query.json')
        with open(self.query_filename, 'w') as query_file:
            query_file.write('{"key": "value"}')
        self.searcher = TextNumericSearcher(self.out_filename, self.query_filename)

    def testGetModeMissingReturnsAnywhere(self):
        field = {'none': ''}
        self.assertEqual('anywhere', self.searcher.get_mode(field))

    def testGetModeEmptyReturnsAnywhere(self):
        field = {'mode': ''}
        self.assertEqual('anywhere', self.searcher.get_mode(field))

    def testGetModeInvalidReturnsAnywhere(self):
        field = {'mode': 'invalid'}
        self.assertEqual('anywhere', self.searcher.get_mode(field))

    def testGetModeValidReturnsValue(self):
        field = {'mode': 'start_of_word'}
        self.assertEqual('start_of_word', self.searcher.get_mode(field))

    def testTelemetry(self):
        UtilitiesLib.ccdc_pipeline_pilot_telemetry = MagicMock()
        TextNumericSearcher(self.out_filename, self.query_filename)
        self.assertEqual(1, UtilitiesLib.ccdc_pipeline_pilot_telemetry.call_count)

class TestValidateScreening(unittest.TestCase):

    @pytest.mark.slow
    def testTelemetry(self):
        UtilitiesLib.ccdc_pipeline_pilot_telemetry = MagicMock()
        query = _find_test_file(__file__, 'aspirin.mol2')
        ValidateScreening(query, [], [], "junk_output_dir", 1, 1)
        self.assertEqual(1, UtilitiesLib.ccdc_pipeline_pilot_telemetry.call_count)

class TestVirtualScreening(unittest.TestCase):

    def testTelemetry(self):
        UtilitiesLib.ccdc_pipeline_pilot_telemetry = MagicMock()
        input_file = _find_test_file(__file__, 'aspirin.mol2')
        VirtualScreening(input_file, [], 1, 1, "junk_output_dir")
        self.assertEqual(1, UtilitiesLib.ccdc_pipeline_pilot_telemetry.call_count)


class TestTelemetry(unittest.TestCase):

    def testBaseClass(self):
        UtilitiesLib.ccdc_pipeline_pilot_telemetry = MagicMock()
        PipelinePilotBase()
        self.assertEqual(1, UtilitiesLib.ccdc_pipeline_pilot_telemetry.call_count)
