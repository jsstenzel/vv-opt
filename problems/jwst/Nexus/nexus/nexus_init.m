% nexus_init.m
% Initialize NEXUS SIMULINK SIMULATION

% flags and options
diagnostics=1;       % print diagnostic messages on command line
verification=1;      % verification steps executed in subroutines (slow)
plotflag=0;          % plots generated in subroutines


%nexus_params
nexus_parfinals

% execute assembly of NEXUS dynamics model

nexus_assy

Azdo=Azd; Bzdo=Bzd; Czdo=Czd; Dzdo=Dzd;


