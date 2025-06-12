#! /usr/bin/env python
########################################################################################################################
#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2025-05-02: shared by the Cambridge Crystallographic Data Centre
#
########################################################################################################################

"""
create_packages.py - Create installation package for the CSD Pipeline Pilot Component Collection.

This creates an archive containing two PP collections (note the standard PP vendor/collection hierarchy)...

ccdc/csd    - The main collection containing CSD components and workflows.
ccdc/csddev - An optional package containing material ancilliary to the main collection.
              It currently only contains the regression suite for the main package, but, in future,
              may contain materials useful in developing further components based on the CSD Software
              Portfolio.

Note that the Pipeline Pilot docs refer to what are here called 'collections' as 'packages'; we're avoiding
that to avoid confusion, as 'package' is used more generally to mean the archive that is distributed to customers.

Warning: If this is used in a developer's buildspace, additional files will be included in the zip
archives, such as '*.pyc' files and '.idea' folders.
"""

########################################################################################################################

import os
import shutil
import sys

########################################################################################################################

# Configuration...

vendor_source_dir = 'ccdc'

package_version = '2023_3'

archive_stem = 'csd-pp-cc'

archive_format = 'zip'

distrib_dir = 'distrib'

build_dir = '_build'

########################################################################################################################

def create_dir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)


def remove_dir(dir_name):
    if os.path.exists(dir_name):
        shutil.rmtree(dir_name)


def remove_file(file_name):
    if os.path.exists(file_name):
        os.remove(file_name)


def create_config(collection_dir):

    conf_file = os.path.join(collection_dir, 'package.conf')
        
    with open(conf_file, 'r') as file:
    
        text = file.read()
    
    text = text.replace('#VER_NUM', package_version)
    
    with open(conf_file, 'w') as file:
    
        file.write(text)

########################################################################################################################

def main():
    
    print("Creating package...", file=sys.stdout)
    
    # Set up the temporary build directory and mirror the sources...
    
    remove_dir(build_dir)
 
    create_dir(build_dir)

    vendor_build_dir = os.path.join(build_dir, os.path.basename(vendor_source_dir))

    shutil.copytree(vendor_source_dir, vendor_build_dir)

    # Configure the collections (currently just sets VER_NUM in conf file)...
    
    for name in os.listdir(vendor_build_dir):  # Collections (AKA PP packages) for vendor

        collection_dir = os.path.join(vendor_build_dir, name)
        
        create_config(collection_dir)
        
    # Create distribution package...
        
    base_name = os.path.join(distrib_dir, f"{archive_stem}-{package_version}")  # NB. No file extension
        
    remove_file(f"{base_name}.{archive_format}")

    archive_file = shutil.make_archive(base_name, archive_format, build_dir)  # NB. make_archive adds file extension 

    remove_dir(build_dir)
    
    print(f"Package created: {archive_file}", file=sys.stdout)
    
    ############

    # Copy the package to the artefacts
    
    if not 'DISTRIBUTION_DIRECTORY' in os.environ:
    
        print("Warning: environment variable 'DISTRIBUTION_DIRECTORY' not set, so will not copy distribution to artefacts.")
    
        sys.exit(0)
    
    distribution_directory = os.environ['DISTRIBUTION_DIRECTORY']
    
    installer_dest_dir = os.path.join(distribution_directory, "pipeline_pilot")
    
    if not os.path.isdir(installer_dest_dir):
    
        os.makedirs(installer_dest_dir)
        
    print(f"Copying package to {installer_dest_dir}...", file=sys.stdout)

    for x in os.listdir(installer_dest_dir):
    
        dest = os.path.join(installer_dest_dir, x)
    
        if os.path.isfile(dest):
            os.remove(dest)
        elif os.path.isdir(dest):
            shutil.rmtree(dest)
    
    if not os.path.isfile(archive_file):
        print("No zip file found", file=sys.stdout)
    else:
        shutil.move(archive_file, installer_dest_dir)

    print("Package moved to artefact", file=sys.stdout)

########################################################################################################################

if __name__ == '__main__':

    main()

########################################################################################################################
# End
########################################################################################################################