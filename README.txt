PRE-REQUISITES

All dependencies are now managed by puppet, see conf/manifests/*.pp

CONFIGURATION

The pipeline API has a single point of configuration: the conf/phylota.ini file.
Locations of 3rd party executables, run-time parameters, input and output file
names etc. are defined there (and only there).

PARALLELIZATION

At present, the pipeline is designed to be parallelized on an MPI architecture,
specifically OpenMPI.
