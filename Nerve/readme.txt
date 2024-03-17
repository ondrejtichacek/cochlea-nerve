%%% readme.txt for bruce_stochHH %%%

This is Version 1.0 of the public distribution of the code for the
stochastic Hodgkin-Huxley models of:

    Bruce, I. C. (2007). "Implementation issues in approximate methods for
    stochastic Hodgkin–Huxley models," Annals of Biomedical Engineering
    35(2):315–318.

Please cite this paper if you publish any research results obtained with
this code or any modified versions of this code.

The Matlab scripts

    bruceABME_figs1and2.m

    bruceABME_fig3.m
and
    bruceABME_fig4.m

can be run to generate figures like those of Bruce (2007).
Note that the stochastic models will give slightly different results every
time they are run, so the figures generated with this code will not match
the published figures exactly.

A single model simulation trial can be run with the function
stochasticHH.m, and action potential statistics for multiple trials can be
obtained with the function APstats.m.

Type

    help stochasticHH
or

    help APstats

for instructions on how to run the respective functions.

%%% Ian C. Bruce (ibruce@ieee.org), Faheem Dinath & Melissa T. Perri © 2007 %%%
