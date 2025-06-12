Scripts placed here will be found by the Derive Script Path component. 

The Derive Script Path component works in the following way:

o If a script file exists which matches exactly the name specified in the Script Filename parameter, this file is used.
o Assuming no exact match is found, the appropriate version is selected for the version of CSD installed:
    - If an exact version match is found, then that version of the script is used.
        For example: CSD 0.8.0 is installed, and 0.8.0 script exists.
    - If no exact version match is found then the highest previous version of the file will be used.
        For example: CSD 0.8.0 is installed, but 0.7.0 is the highest version of the script
    - If no previous version is found (e.g. version 1.0.0 is the lowest script version but 0.8.0 CSD is installed) then an error will be thrown.

Whichever script file path is found to be most appropriate is stored in the global property defined in the Script Global Name.

NOTE:

The version is specified in the file name before the .py file extension, and is in the same form as produced by the following Python script:

import ccdc
ccdc.__version__

For example:

script.1.2.3.py