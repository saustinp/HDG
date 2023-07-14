## Readme for directory ionicwind2/swarm

- swarm.txt: text file containing output from the BOLSIG+ solver on LXCAT
- bolsig_data_processing.py: first point of entry, processes the data in swarm.txt
    1. Computed polynomial fit coefficients for the curves in swarm.txt
    2. Validates the curve fits by overlaying the fitted curves on the original data set
- bolsig_data.mat: contains the relevant swarm parameters in a matlab-readable file
- swarm_params.m: function, outputs the swarm parameters given an input E, N -> to be used by a driver script -> DEPRECATED: USE THE INDIVIDUAL FUNCTIONS INSTEAD
- get_swarm_validate.m: script to validate the output of the function to calculate the swarm paramters
