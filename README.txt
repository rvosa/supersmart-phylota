PRE-REQUISITES

All dependencies are now managed by puppet, see conf/manifests/*.pp, see ./INSTALL

CONFIGURATION

The pipeline API has a single point of configuration: the conf/phylota.ini file.
Locations of 3rd party executables, run-time parameters, input and output file
names etc. are defined there (and only there).

EXAMPLES

The examples/* directories contain instructional examples to demonstrate the input
file formats, the run.sh shell script shows how the steps are chained together.