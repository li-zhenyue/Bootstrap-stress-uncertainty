# Bootstrap-stress-uncertainty
MATLAB source code for the manuscript "Reliability of Bootstrap-derived stress uncertainty bounds in focal mechanism inversion"

Run the M-files named "figureN.m" (e.g., "figure1.m", "figure2.m") to directly reproduce the corresponding figures in the manuscript.

For figures based on synthetic data, the corresponding MATLAB scripts generate the synthetic focal mechanisms as part of the plotting workflow. The generated synthetic data can therefore be obtained by running these scripts, and the data generation procedures can be examined in the source code.

Because the seismogenic faults in the synthetic experiments are generated randomly, repeated runs of the program will not yield exactly identical results.

The stress inversion results required to reproduce Figures 7 and 8 are provided in AVB_results.mat and TYM_results.mat, respectively. These results can also be reproduced by running real_data.m.

The code was written in MATLAB R2024a.

The fault-slip data for the AVB and TYM regions were extracted from Table A1 of J. Angelier (1990):

Angelier, J., 1990. Inversion of field data in fault tectonics to obtain the regional stress-III. A new rapid direct inversion method by analytical means. Geophys. J. Int., 103(2), 363-376.
