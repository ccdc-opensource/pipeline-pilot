#! /bin/bash
#######################################################################################################################
#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2025-05-02: shared by the Cambridge Crystallographic Data Centre
#
#######################################################################################################################
# 
# Install enviroment variable file for the CSD Pipeline Pilot Component Collection on Linux
# 
#######################################################################################################################

[ -z "${PACKAGE_DIR}" ] && {

    echo "ERROR! Required environment variable PACKAGE_DIR is not set or is empty." 1>&2

    exit 1
}

[ -z "${CSDHOME}" ] && {

    echo "ERROR! Required environment variable CSDHOME is not set or is empty." 1>&2

    exit 1
}

############

# Files...

SCIROOT=$(realpath ${PACKAGE_DIR}/../../..)  # Absolute path of PP Server installation

env_file=${SCIROOT}/linux_bin/env.d/ccdc_env.sh  # NB. Any file in this dir is run, so name is arbitrary

installed_files=${PACKAGE_DIR}/logs/installed_files.txt

############

# Install enviroment variable file...

cat > ${env_file} <<EOF
# Ensure the CSD can be found...

export CSDHOME=${CSDHOME}
EOF

# Record installed file for uninstall procedure...

echo ${env_file} >> ${installed_files}

########################################################################################################################
# End
########################################################################################################################