FOCAL PLANES:
Experiment models are based off of historically-performed experiments, as described in Stenzel 2022.

Historical valdiation data is located in a plotting script. This is not the latest data; after I performed my testing, a number of the focal planes needed to be refitted and retested, but it will do for validation data.

According to the huaso repository llamasConfig, CCDs are on the following benches:
1A
	red - 012
	green - 014
	blue - 019
1B
	red - 026
	green - 023
	blue - 007
2A
	red - 009
	green - 008
	blue - 010
2B
	red - 018
	green - 011
	blue - 002
3A
	red - 020
	green - 016
	blue - 024
3B
	red - 024
	green - 017
	blue - 022
4A
	red - 021
	green - 005
	blue - 013
4B
	red - 006
	green - 025
	blue - 004

QUANTUM EFFICIENCY:
QE was never historically measured in-house for the focal planes.

VPH:
All VPH for all 8 spectrographs were tested in-house. The data for the 24 VPH are in this spreadsheet.

At one point, I used a parabolic interpolation algorithm to make all 24 coating curves. TODO grab those and put them here, for completeness.

TODO think about how to represent testing of 1 spectrograph vs 8 spectrographs in my experiment model. The most straightforward answer right now, considering the scope of a single spectrograph, is that it doesn't make a difference and any difference is out-of-scope for the verification problem.

TODO the inhouse data shows a standard deviation of about 0.0267 for throughput measurements. I should use that for the experiment model.

SL DICHROIC:
There were 2 lots of the dichroics provided, A and B. I'm not sure how many of each we received, but a dichroic from each lot was tested, at 801 measurement points. I threw them together in one excel sheet, SL_run_comparison, to see how they differed.

BG DICHROIC:
There was 1 lot of dichroics provided, 1 dichroic was tested at 801 measurement points.

COLLIMATOR:
In-house testing of a collimator was done, with 1501 measurement points.

LENSES:
It does not seem that camera throughput testing was done for LLAMAS. Plenty of optical tolerancing work was done on LLAMAS, but that is a separate set of parameters that serve other quantities of interest.

Also located here are test results for detailed camera testing of the protoLLAMAS lenses. Here, the camera lens assemblies were tested all together, at 13 points each for the red and blue cameras. This constitutes the full historical camera throughput testing that was done.

FRD:
From evaluating_cleavage...
The average of the sample of 10 stripped is 92.3%, and the scatter is 2.4%.
Previously I used this as the prior, but this wouldn't have been known before testing.
This is more appropriately used as validation data.

GENERAL CONSIDERATIONS:
"Measurement Considerations When Specifying Optical Coatings," Pete Kupinski and Angus Macleod. This paper indicates a best case +- 0.1% T for commercial measurements of highly transmissive coatings.
Therefore, I should assume that any historical data has that uncertainty, unless otherwise stated.
