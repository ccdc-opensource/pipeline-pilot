#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2017-04-03: made available by the Cambridge Crystallographic Data Centre
#

from abc import ABCMeta, abstractmethod

import ccdc.conformer
import ccdc.io
import ccdc.screening
import ccdc.search
import ccdc.entry

import json
import os
import re
import shutil
import types
import datetime

########################################################################################################################
# 
# Monkey-patching to fix JSON-encoding issues introduced by migrating from Python2 to Python3.
# NB. These fixes are only temporary and, ultimately, the Python API code with need modifying.
# https://stackoverflow.com/questions/18478287/making-object-json-serializable-with-regular-encoder/18561055#18561055


def _deposition_date(self):

    d = self._entry.editors_info().accession_date()
    
    return str(datetime.date(d.year(), d.month(), d.day())) if d.year() > 0 else None


ccdc.entry.Entry.deposition_date = property(fget=_deposition_date)


def _journal_to_json(self):

    return str(self.name if self.full_name == '' else self.full_name)


ccdc.entry.Journal._to_json = _journal_to_json


def _new_jsonencoder_default(_, obj):  # NB. Will be used as method, so first arg will be self

    return getattr(obj.__class__, '_to_json', _old_jsonencoder_default)(obj)


_old_jsonencoder_default = json.JSONEncoder.default  # Save default...
json.JSONEncoder.default = _new_jsonencoder_default  # ...then replace it.

# NB. An alternative would be to subclass JSONEncoder; this is perhaps cleaner, but would mean code using the
# encoder would need altering (i.e. 'json_encoder = json.JSONEncoder()' -> 'json_encoder = myJSONEncoder()').
#
# class myJSONEncoder(json.JSONEncoder):
#     def default(self, obj):
#         if isinstance(obj, Journal):
#             return str(obj.name if obj.full_name == '' else obj.full_name)
#         return super().default(self, obj)

########################################################################################################################


class PipelinePilotBase:
    def __init__(self):
        """Base class constructor that calls Telemetry."""
        from ccdc.utilities import _private_importer
        with _private_importer() as pi:
            pi.import_ccdc_module('UtilitiesLib')
            UtilitiesLib.ccdc_pipeline_pilot_telemetry()


class _AbstractSearch(PipelinePilotBase, metaclass=ABCMeta):

    def __init__(self, output_file, is_similarty=False):
        super().__init__()
        self.db_list = []
        self.file_list = []
        self._is_similarity = is_similarty
        self._output_file = output_file
        self._output_file_h = None
        self._databases = ''
        self._reader = None

    def add_hits_to_file(self, hits, i_query_index, source=None):
        if self._output_file_h is None:
            self._output_file_h = open(self._output_file, 'w')
            if self._is_similarity:
                # _index        - the query index number (used to join data to the delimited text file later)
                # csd_refcode   - the refcode for the "hit" structure
                # similarity    - the similarity between this structure and the query.
                self._output_file_h.write('_index,csd_refcode,database,similarity\n')
            else:
                # _index        - the query index number (used to join data to the delimited text file later)
                # csd_refcode   - the refcode for the "hit" structure
                self._output_file_h.write('_index,csd_refcode,database\n')

        for hit in hits:
            try:
                # Gather refcode.
                csd_refcode = str(hit.identifier)
                # Write to file.
                if source is not None:
                    database_entry = source
                else:
                    database_entry = hit.entry.database_name
                if self._is_similarity:
                    self._output_file_h.write(','.join(map(str, [i_query_index, csd_refcode, database_entry,
                                                                 hit.similarity])) + '\n')
                else:
                    self._output_file_h.write(','.join(map(str, [i_query_index, csd_refcode, database_entry])) + '\n')
            except IOError:
                continue

        self._output_file_h.flush()  # It helps debugging to flush when all hits are written (as file is not closed until self is destroyed)

    @abstractmethod
    def search(self):
        pass

    def set_databases(self, databases):
        self._databases = databases

        databases = self._databases.split(',')
        # Need to test to see if any MOL2, SDF, SD, CIF or MOL files are listed.
        #   These get added to files list
        #   CSD databases get added to the db list.
        # The DB list can be tackled in one go, the file list must be tackled, one record at a time
        for idx, db in enumerate(databases):
            if re.search('.(mol2|sdf|sd|cif|mol)$', db):
                print('Adding ' + db + ' to file list')
                self.file_list.append(ccdc.io.MoleculeReader(db))
            else:
                print('Adding ' + db + ' to DB list')
                self.db_list.append(db)

    def read_settings(self, preference_file):
        settings_file = open(preference_file)
        self._reader = SettingsReader(settings_file)


class CompoundConformerLibraryStorage(PipelinePilotBase):
    """Manage the generation and storage of conformers."""

    def __init__(self, input_file, num_confs, max_unusual_torsions, threads, superimpose):
        super().__init__()
        self._input = input_file
        self._num_confs = num_confs
        self._max_unusual_torsions = max_unusual_torsions
        self._threads = threads
        self._superimpose = superimpose

    @staticmethod
    def standardise(mol, add_hydrogens=True):
        """Return standardised molecule."""
        mol.standardise_aromatic_bonds()
        mol.standardise_delocalised_bonds()
        if add_hydrogens:
            mol.add_hydrogens()
        return mol

    def _set_up_conformer_generator(self):
        """Fish out arguments."""
        settings = ccdc.conformer.ConformerSettings()
        settings.max_conformers = self._num_confs
        settings.max_unusual_torsions = self._max_unusual_torsions
        settings.superimpose_conformers_onto_reference = self._superimpose
        return settings

    def generate_conformers(self):
        """Generate the conformers list."""
        molecules = ccdc.io.MoleculeReader(self._input)
        settings = self._set_up_conformer_generator()
        gen = ccdc.conformer.ConformerGenerator(settings, nthreads=self._threads)
        standardised_molecules = [self.standardise(m) for m in molecules]
        conformers_hit_list = gen.generate(standardised_molecules)
        return conformers_hit_list

    def split_conformers_output(self, conformers_hit_list, output_dir=None):
        """Write out each conformational ensemble for each input molecule to a separate file."""
        print("Split conformers")
        in_path = os.path.split(self._input)[0]
        if not in_path:
            in_path = '.'
        if output_dir is None:
            output_dir = os.path.join(in_path, 'conformers')
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)  # Remove all the subdirectories!
            else:
                os.mkdir(output_dir)

        for conformers in conformers_hit_list:
            output_file_name = '%s/%s_confs.sdf' % (output_dir, conformers.identifier)
            print('Writing to ' + output_file_name)
            with ccdc.io.MoleculeWriter(output_file_name) as mol_writer:
                for conformer in conformers:
                    mol_writer.write(conformer.molecule)

    def merge_conformers_output(self, conformers_hit_list, output_dir=None):
        """Write out the conformer libraries to a single file."""
        print("Merge conformers")
        if output_dir is None:
            output_file_name = '%s_confs.sdf' % os.path.splitext(os.path.basename(self._input))[0]
        else:
            output_file_name = output_dir + ('/%s_confs.sdf' % self._input)

        print('Writing to ' + output_file_name)
        with ccdc.io.MoleculeWriter(output_file_name) as mol_writer:
            for conformers in conformers_hit_list:
                for conformer in conformers:
                    mol_writer.write(conformer.molecule)

    def write_report(self, conformers_hit_list, output_dir=None):
        if output_dir is None:
            output_file_name = '%s_confs.csv' % os.path.splitext(os.path.basename(self._input))[0]
        else:
            output_file_name = output_dir + ('/%s_confs.csv' % os.path.splitext(os.path.basename(self._input))[0])
        print("writing report to " + output_file_name)
        with open(output_file_name, 'w') as writer:
            writer.write('Molecule name,max_log_probability,min_log_probability,n_conf_gen_clust,'
                         'n_rotamers_with_no_observations,rotamers_with_no_observations\n')
            for mol in conformers_hit_list:
                if hasattr(mol, 'unusual_rotamers'):
                    writer.write('%s,%.3f,%.3f,%d,%d,%s\n' % (mol.identifier,
                                                              mol.max_log_probability,
                                                              mol.min_log_probability,
                                                              len(mol.hits),
                                                              mol.n_rotamers_with_no_observations,
                                                              len(mol.unusual_rotamers)))
                else:
                    writer.write('%s,%.3f,%.3f,%d,%d,%s\n' % (mol.identifier,
                                                              mol.max_log_probability,
                                                              mol.min_log_probability,
                                                              len(mol.hits),
                                                              mol.n_rotamers_with_no_observations,
                                                              mol.rotamers_with_no_observations))


class JoinAttributes(object):
    def __init__(self, csd_refcode_list_file, database, join_type, attributes):
        self._input_file = csd_refcode_list_file
        self._join_type = join_type
        self._attributes = attributes
        self._database = database

    def write_attribute_file(self, output_file):
        entry_list = open(self._input_file, 'r')

        # Open entry reader
        entry_reader = ccdc.io.EntryReader(self._database)

        # Open output file; will be CSV file with at least the refcode and other attributes
        fh = open(output_file, 'w')

        # Start file with open square bracket
        fh.write('[')

        # Create JSON encoder for records.
        json_encoder = json.JSONEncoder()

        # Initialize the number written.
        i_num_written = 0

        # For each entry...
        for csd_refcode in entry_list:
            # Strip white space from refcode
            csd_refcode = csd_refcode.strip()

            try:
                record_text = ''
                if self._join_type == 'Molecule':
                    mol = entry_reader.molecule(csd_refcode)
                    record_text = self.write_mol_record(mol)
                elif self._join_type == 'Crystal':
                    crystal = entry_reader.crystal(csd_refcode)
                    record_text = self.write_crystal_record(crystal)
                elif self._join_type == 'Entry':
                    entry = entry_reader.entry(csd_refcode)
                    record_text = self.write_entry_record(entry)

                if i_num_written > 0:
                    fh.write(',\n')
                # Write line to file.
                fh.write(json_encoder.encode(record_text))

                # Increment counter
                i_num_written += 1
            except RuntimeError as err:
                if csd_refcode + ' is not in the database' in format(err):
                    pass
                else:
                    raise

        # Close file with closing square bracket
        fh.write(']')
        fh.close()

        print('=======')
        print('num_written=' + str(i_num_written))

    def write_crystal_record(self, crystal):
        # Construct dictionary
        record = {'__csd_refcode': crystal.identifier}

        # Gather attributes from list
        for attr in self._attributes:
            if attr == 'void_volume':
                value = crystal.void_volume()
            elif attr == 'disordered_molecule_ctab':
                value = crystal.disordered_molecule.to_string('mol')
            elif attr == 'disordered_molecule_mol2':
                value = crystal.disordered_molecule.to_string('mol2')
            elif attr == 'molecule_ctab':
                value = crystal.molecule.to_string('mol')
            elif attr == 'molecule_mol2':
                value = crystal.molecule.to_string('mol2')
            else:
                try:
                    value = getattr(crystal, attr, attr + ' not found')
                except RuntimeError:
                    value = attr + ' not found'
            record[attr] = value

        return record

    def write_entry_record(self, entry):
        # Construct dictionary
        record = {'__csd_refcode': entry.identifier}

        # Gather attributes from list
        for attr in self._attributes:
            try:
                value = getattr(entry, attr, '')
                if attr == 'publications':
                    value = str(value[0])
            except RuntimeError:
                value = ''

            record[attr] = value

        return record

    def write_mol_record(self, mol):
        record = {'__csd_refcode': mol.identifier}

        # Gather attributes from list
        for attr in self._attributes:
            if attr == 'ctab':
                # This is a special attribute. Gather from to_string
                value = mol.to_string('mol')
            elif attr == 'mol2':
                # This is a special attribute. Gather from to_string
                value = mol.to_string('mol2')
            elif attr == 'heavy_atoms': 
                # This is a special attribute, as a list of Atom objects is returned.
                # value = ', '.join(str(x) for x in getattr(mol, attr, []))
                value = len(getattr(mol, attr, []))
            else:
                try:
                    value = getattr(mol, attr, '')

                except RuntimeError:
                    value = ''

            record[attr] = value
        return record


class ReducedCellSearcher(_AbstractSearch):
    def __init__(self, output_file, cell_angles, cell_lengths, lattice_centring):
        super(ReducedCellSearcher, self).__init__(output_file)

        self._cell_angles = cell_angles
        self._cell_lengths = cell_lengths
        self._lattice_centring = lattice_centring

    def search(self):
        query = ccdc.search.ReducedCellSearch.Query(self._cell_lengths, self._cell_angles, self._lattice_centring)
        if self._reader is not None:
            settings = ccdc.search.ReducedCellSearch.Settings()
            self._reader.update_settings(settings)
        else:
            settings = None

        searcher = ccdc.search.ReducedCellSearch(query, settings)
        if self._databases != '':
            # Gather database(s) to use.
            # If DB search requested...
            if len(self.db_list) > 0:
                hits = searcher.search(self.db_list)
                if len(hits) > 0:
                    print(str(len(hits)) + " hits found for db search")
                    self.add_hits_to_file(hits, 1)
                else:
                    print('No database hits')
            # If file search requested...
            if len(self.file_list) > 0:
                # Deal with each file in turn...
                for file_source in self.file_list:
                    hits = searcher.search(file_source)
                    if len(hits) > 0:
                        self.add_hits_to_file(hits, 1, file_source.file_name)
                    else:
                        print('No hits found in file source ' + file_source.file_name)
        else:
            # Defaults with no databases listed.
            hits = searcher.search()
            if len(hits) > 0:
                # Open output file; CSV with just ref code.
                self.add_hits_to_file(hits, 1)
                print(str(len(hits)) + " hits")
            else:
                print('No hits')


class SettingsReader(object):
    """A class to parse settings passed from the component to Python. Settings are passed using a JSON file which is
    created by Pipeline Pilot"""
    max_hit_structures = None
    has_3d_coordinates = None
    max_r_factor = None
    no_disorder = None
    no_errors = None
    no_ions = None
    no_power = None
    not_polymeric = None
    only_organic = None
    only_organometallic = None

    def __init__(self, json_file):
        # Load JSON data file
        data = json.load(json_file)

        # Iterate over the attributes of this class
        for a in dir(self):
            # If the attribute doesn't start __ and isn't a function...
            if not a.startswith('__') and not callable(getattr(self, a)):
                # Does JSON mention this attribute?
                if a in data:
                    # It does; set corresponding class attribute
                    setattr(self, a, data[a])

    def settings(self):
        settings = ccdc.search.Search.Settings()
        self.update_settings(settings)
        return settings

    def update_settings(self, settings):
        """
        Get the settings object which this class represents. If the attributes are all set to None, then no settings
        object is required and None is returned.
        :rtype: ccdc.search.Search.Settings
        """
        is_anything_set = False

        # Anything set to anything other than None?
        for a in dir(self):
            # If the attribute doesn't start __ and isn't a function...
            if not a.startswith('__'):
                if getattr(self, a) is not None:
                    if not isinstance(getattr(self, a), types.FunctionType):  # not callable(getattr(self, a)):
                        is_anything_set = True

        if is_anything_set:
            for a in dir(self):
                # If the attribute doesn't start __ and isn't a function...
                if not a.startswith('__') and not callable(getattr(self, a)):
                    # Is this attribute set to something?
                    if getattr(self, a) is not None:
                        print("setting " + a + " on the settings object to " + str(getattr(self, a)))
                        if a == 'max_hit_structures':
                            setattr(settings, a, int(getattr(self, a)))
                        else:
                            setattr(settings, a, getattr(self, a))
        else:
            settings = None

        return settings


class StructureSearcher(_AbstractSearch):
    def __init__(self, query_file, output_file, search_type):
        super(StructureSearcher, self).__init__(output_file, search_type == 'similarity')
        self._search_type = search_type
        # self._is_similarity = search_type == 'similarity'  # Redundant, as _is_similarity is set in _AbstractSearch.__init__()
        self._query_file = ccdc.io.MoleculeReader(query_file)
        self._similarity_threshold = 0.7

    @property
    def similarity_threshold(self):
        """Similarity threshold"""
        return self._similarity_threshold

    @similarity_threshold.setter
    def similarity_threshold(self, threshold):
        self._similarity_threshold = threshold

    def search(self):
        if self._query_file is None:
            return

        # Initialize the number written.
        i_query_index = 0

        print('Search type is ' + self._search_type)

        # For each query in the file...
        for query in self._query_file:
            # Increment the query index
            i_query_index += 1

            # If exact or similarity, add hydrogens.
            if (self._search_type == 'exact') or (self._search_type == 'similarity'):
                # Fix unknown bond types
                query.assign_bond_types(which='unknown')
                # Add missing hydrogens (NB. 'all' would remove all Hs then add all Hs back).
                query.add_hydrogens(mode='missing')

            # Standardise the aromatic bonds...
            query.standardise_aromatic_bonds()

            if (self._search_type == 'exact') or (self._search_type == 'substructure'):
                if self._reader is not None:
                    settings = ccdc.search.SubstructureSearch.Settings()
                    self._reader.update_settings(settings)
                else:
                    settings = None

                # Initialize the substructure - done here as there's no option to remove a substructure
                search_class_instance = ccdc.search.SubstructureSearch(settings)

                for component in query.components:
                    search_class_instance.add_substructure(ccdc.search.MoleculeSubstructure(component))
            else:
                settings = ccdc.search.SimilaritySearch.Settings()
                settings.threshold = self._similarity_threshold
                if self._reader is not None:
                    self._reader.update_settings(settings)
                search_class_instance = ccdc.search.SimilaritySearch(
                    query, self._similarity_threshold, 'tanimoto', settings
                )

            if self._databases != '':
                # If DB search requested...
                if len(self.db_list) > 0:
                    hits = search_class_instance.search(self.db_list)
                    if len(hits) > 0:
                        self.add_hits_to_file(hits, i_query_index)
                    else:
                        print('No database hits')
                # If file search requested...
                if len(self.file_list) > 0:
                    # Deal with each file in turn...
                    for file_source in self.file_list:
                        hits = search_class_instance.search(file_source)
                        if len(hits) > 0:
                            self.add_hits_to_file(hits, i_query_index, file_source.file_name)
                        else:
                            print('No hits found in file source ' + file_source.file_name)
            else:
                # Defaults with no databases listed.
                hits = search_class_instance.search()
                self.add_hits_to_file(hits, i_query_index)


class TextNumericSearcher(_AbstractSearch):
    MODES = ('not_null',
             'anywhere',
             'separate',
             'is_null',
             'exact',
             'start',
             'start_of_word')

    default_mode = 'anywhere'

    def __init__(self, output_file, query_file):
        super(TextNumericSearcher, self).__init__(output_file)
        self._query_file = query_file

        query_json = open(query_file, 'r')
        self._query = json.load(query_json)

    def get_mode(self, field):
        """
        Derive the mode specified for the search. The field supplied is a JSON object which may or may not have a "mode"
        attribute, therefore this is the first thing which is tested. If not mode attribute is found, then an empty
        string is returned. Assuming a mode attribute is found, it is checked to ensure that it is one of the values
        found in the MODE list. If it is not found in that list and is found to be valid then it is returned. Otherwise,
        the default mode of "anywhere" is returned.

        :param field: The JSON object field.
        :return: The mode specified in the field object, or 'anywhere', if the mode is missing or invalid.
        """
        if 'mode' in field:
            mode = field['mode']

            # Check mode is valid.
            if not any(m in mode for m in self.MODES):
                # unknown mode; clear it to ignore it
                mode = self.default_mode
        else:
            mode = self.default_mode
        return mode

    def search(self):
        text_numeric_query = ccdc.search.TextNumericSearch()
        if self._reader is not None:
            self._reader.update_settings(text_numeric_query.settings)

        if 'AllText' in self._query:
            field = self._query['AllText']
            mode = self.get_mode(field)

            value = field['value']

            if isinstance(value, list):
                for v in value:
                    text_numeric_query.add_all_text(v, mode=mode)
            else:
                text_numeric_query.add_all_text(value, mode=mode)

        if 'Bioactivity' in self._query:
            field = self._query['Bioactivity']
            mode = self.get_mode(field)
            if 'ignore_non_alpha_num' in field:
                ignore_non_alpha_num = bool(field['ignore_non_alpha_num'])
            else:
                ignore_non_alpha_num = False

            value = field['value']

            if isinstance(value, list):
                for v in value:
                    text_numeric_query.add_bioactivity(v, mode=mode, ignore_non_alpha_num=ignore_non_alpha_num)
            else:
                text_numeric_query.add_bioactivity(value, mode=mode, ignore_non_alpha_num=ignore_non_alpha_num)

        if 'Citation' in self._query:
            field = self._query['Citation']

            # Create defaults:
            if 'author' in field:
                author = field['author']
            else:
                author = ''

            if 'journal' in field:
                journal = field['journal']
            else:
                journal = ''

            if 'volume' in field:
                volume = int(field['volume'])
            else:
                volume = None

            if 'year' in field:
                year = field['year']
                if isinstance(year, list):
                    year = [int(y) for y in year]
                else:
                    year = int(year)
            else:
                year = None

            if 'first_page' in field:
                first_page = int(field['first_page'])
            else:
                first_page = None

            if 'ignore_non_alpha_num' in field:
                ignore_non_alpha_num = bool(field['ignore_non_alpha_num'])
            else:
                ignore_non_alpha_num = False

            text_numeric_query.add_citation(author=author,
                                            journal=journal,
                                            year=year,
                                            volume=volume,
                                            first_page=first_page,
                                            ignore_non_alpha_num=ignore_non_alpha_num)

        if 'Color' in self._query:
            field = self._query['Color']
            mode = self.get_mode(field)

            text_numeric_query.add_color(field['value'], mode=mode)

        if 'CompoundName' in self._query:
            field = self._query['CompoundName']
            mode = self.get_mode(field)

            text_numeric_query.add_compound_name(field['value'], mode=mode)

        if 'disorder' in self._query:
            field = self._query['disorder']
            mode = self.get_mode(field)

            value = field['value']

            if isinstance(value, list):
                for v in value:
                    text_numeric_query.add_disorder(v, mode=mode)
            else:
                text_numeric_query.add_bioactivity(value, mode=mode)

        if 'Identifier' in self._query:
            field = self._query['Identifier']
            mode = self.get_mode(field)
            if 'ignore_non_alpha_num' in field:
                ignore_non_alpha_num = bool(field['ignore_non_alpha_num'])
            else:
                ignore_non_alpha_num = False

            value = field['value']

            if isinstance(value, list):
                for v in value:
                    text_numeric_query.add_identifier(v, mode=mode, ignore_non_alpha_num=ignore_non_alpha_num)
            else:
                text_numeric_query.add_identifier(value, mode=mode, ignore_non_alpha_num=ignore_non_alpha_num)

        if 'Synonym' in self._query:
            field = self._query['Synonym']
            mode = self.get_mode(field)

            text_numeric_query.add_synonym(field['value'], mode=mode)

        print(text_numeric_query.queries)

        # Perform query
        if self._databases != '':
            # Gather database(s) to use.
            # If DB search requested...
            if len(self.db_list) > 0:
                hits = text_numeric_query.search(self.db_list)
                if len(hits) > 0:
                    print(str(len(hits)) + " hits found for db search")
                    self.add_hits_to_file(hits, 1)
                else:
                    print('No database hits')
            # If file search requested...
            if len(self.file_list) > 0:
                # Deal with each file in turn...
                for file_source in self.file_list:
                    hits = text_numeric_query.search(file_source)
                    if len(hits) > 0:
                        self.add_hits_to_file(hits, 1, file_source.file_name)
                    else:
                        print('No hits found in file source ' + file_source.file_name)
        else:
            # Defaults with no databases listed.
            hits = text_numeric_query.search()
            if len(hits) > 0:
                print(str(len(hits)) + " hits found for default db search")
                self.add_hits_to_file(hits, 1)
            else:
                print('No database hits')


class ValidateScreening(PipelinePilotBase):
    def __init__(self, query, actives_to_screen, decoys_to_screen, output_dir, threads, num_confs):
        super().__init__()
        query_mols = [m for m in ccdc.io.MoleculeReader(query)]  # Read the query mol or overlay of mols
        self._actives_to_screen = actives_to_screen
        self._decoys_to_screen = decoys_to_screen
        self._settings = ccdc.screening.Screener.Settings()
        self.setup_screener_settings(output_dir)
        self._output_name_actives = os.path.join(self._settings.output_directory, "actives_screened.mol2")
        self._output_name_decoys = os.path.join(self._settings.output_directory, "decoys_screened.mol2")
        self._nt = threads
        self._nc = num_confs
        # Create screener from query molecules plus settings and threads.
        self._screener = ccdc.screening.Screener(query_mols, settings=self._settings, nthreads=threads)

    @staticmethod
    def generate_confs(mol, gen):
        """Return the conformer hit lists."""
        standardised_mol = CompoundConformerLibraryStorage.standardise(mol, False)
        conf_hit_lists = gen.generate(standardised_mol)
        confs = [[c.molecule for c in conf_hit_lists]]
        return confs

    def screen_actives(self):
        return self.screen_molecules(self._screener, self._actives_to_screen, 1,
                                     self._output_name_actives)  # Screen set of actives

    def screen_decoys(self):
        return self.screen_molecules(self._screener, self._decoys_to_screen, 0,
                                     self._output_name_decoys)  # Screen set of decoys

    def screen_molecules(self, screener, mols_to_screen, activity, output_name):
        """Run the ligand screener and write out the screened conformations.
        Return the list of ranked scores"""

        if self._nc != 0:
            settings = ccdc.conformer.ConformerSettings()
            settings.max_conformers = self._nc
            gen = ccdc.conformer.ConformerGenerator(settings, nthreads=self._nt)
        else:
            gen = None

        # Read the molecules to screen
        screen_set = [m for m in ccdc.io.MoleculeReader(mols_to_screen)]
        scores = []
        with ccdc.io.MoleculeWriter(output_name) as mol_writer:
            for mol in screen_set:

                # If nconformers=0 the input conformation is used, otherwise the CSD-driven conformer
                # generator is called and a number of conformations equal to nconformers is generated.
                if self._nc == 0:
                    confs = [[CompoundConformerLibraryStorage.standardise(mol, False)]]
                else:
                    try:
                        confs = ValidateScreening.generate_confs(mol, gen)
                    except Exception:  # TODO: Identify exception types that should be caught
                        confs = [[CompoundConformerLibraryStorage.standardise(mol, False)]]
                res = screener.screen(confs)  # Screening step
                scores.extend([(r.score, activity, r.identifier) for r in res])
                # Write results
                for r in res:
                    mol_writer.write(r.molecule)

        return sorted(scores)

    def setup_screener_settings(self, output_dir):
        """Return the settings object for the ligand screener."""
        # self._settings = ccdc.screening.Screener.Settings()
        param_dir = self._settings.parameter_directory
        new_dir = os.path.join(os.getcwd(), 'parameter_files')
        try:
            # copy the parameter_files directory from the default location into the working directory.
            shutil.copytree(param_dir, new_dir)
        except Exception:  # TODO: Identify exception types that should be caught
            # if there is already a parameter_files directory, then this will be used for the calculation.
            pass
        self._settings.parameter_directory = new_dir
        self._settings.output_directory = os.path.join(os.getcwd(), output_dir)
        if os.path.exists(self._settings.output_directory):
            shutil.rmtree(self._settings.output_directory)
        os.makedirs(self._settings.output_directory)

    def write_scores(self, results):
        """Write results to an output file."""
        output_name_scores = os.path.join(self._settings.output_directory, "screening_scores.csv")

        with open(output_name_scores, 'w') as f:
            for r in results:
                score, a, identifier = r
                f.write('%.3f, %d, %s\n' % (score, a, identifier))


class VirtualScreening(PipelinePilotBase):
    def __init__(self, input_file, mols_to_screen, threads, num_confs, output_dir):
        super().__init__()
        self._query = [m for m in ccdc.io.MoleculeReader(input_file)]  # Read the query mol or overlay of mols
        self._mols_to_screen = mols_to_screen
        self._nt = threads
        self._nc = num_confs
        self._output_dir = output_dir
        self._settings = None
        self.setup_screening_settings(output_dir)

    def setup_screening_settings(self, output_dir):
        """Return the settings object for the ligand screener."""
        self._settings = ccdc.screening.Screener.Settings()
        param_dir = self._settings.parameter_directory
        new_dir = os.path.join(os.getcwd(), 'parameter_files')
        try:
            # copy the parameter_files directory from the default location into the working directory.
            shutil.copytree(param_dir, new_dir)
        except Exception:  # TODO: Identify exception types that should be caught
            # if there is already a parameter_files directory, then this will be used for the calculation.
            pass
        self._settings.parameter_directory = new_dir
        self._settings.output_directory = os.path.join(os.getcwd(), output_dir)
        if os.path.exists(self._settings.output_directory):
            shutil.rmtree(self._settings.output_directory)
        os.makedirs(self._settings.output_directory)

    def screen_molecules(self):
        """Run the ligand screener and write out the screened conformations.
        Return the list of ranked scores"""
        if self._nc > 0:
            settings = ccdc.conformer.ConformerSettings()
            settings.max_conformers = self._nc
            gen = ccdc.conformer.ConformerGenerator(settings, nthreads=self._nt)
        else:
            gen = None

        screener = ccdc.screening.Screener(self._query, settings=self._settings, nthreads=self._nt)
        output_name = os.path.join(self._settings.output_directory, "screened_molecules.mol2")

        screen_set = [m for m in ccdc.io.MoleculeReader(self._mols_to_screen)]  # Read the molecules to screen
        scores = []
        mol_writer = ccdc.io.MoleculeWriter(output_name)
        for mol in screen_set:

            # If nconformers=0 the input conformation is used, otherwise the CSD-driven conformer
            # generator is called and a number of conformations equal to nconformers is generated.
            if self._nc == 0:
                confs = [[CompoundConformerLibraryStorage.standardise(mol, False)]]
            else:
                try:
                    confs = ValidateScreening.generate_confs(mol, gen)
                except RuntimeError:
                    confs = [CompoundConformerLibraryStorage.standardise(mol, False)]
            res = screener.screen(confs)  # Screening step
            scores.extend([(r.score, r.identifier) for r in res])
            # Write results
            for r in res:
                mol_writer.write(r.molecule)

        return sorted(scores)

    def write_scores(self, results):
        """Write results to an output file."""
        output_name_scores = os.path.join(self._settings.output_directory, "screening_scores.csv")
        with open(output_name_scores, 'w') as f:
            for r in results:
                score, _id = r
                f.write('%.3f, %s\n' % (score, _id))
