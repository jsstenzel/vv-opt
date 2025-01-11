FOCAL PLANES:
Experiment models are based off of historically-performed experiments, as described in Stenzel 2022.

Historical valdiation data is TODO

QUANTUM EFFICIENCY:
QE was never historically measured in-house for the focal planes.

VPH:
All VPH for all 8 spectrographs were tested in-house. The data for the 24 VPH are in this spreadsheet.

At one point, I used a parabolic interpolation algorithm to make all 24 coating curves. TODO grab those and put them here, for completeness.

TODO think about how to represent testing of 1 spectrograph vs 8 spectrographs in my experiment model. The most straightforward answer right now, considering the scope of a single spectrograph, is that it doesn't make a difference and any difference is out-of-scope for the verification problem.

SL DICHROIC:
There were 2 lots of the dichroics provided, A and B. I'm not sure how many of each we received, but a dichroic from each lot was tested, at 801 measurement points. I threw them together in one excel sheet, SL_run_comparison, to see how they differed.

TODO update historical test plan.

BG DICHROIC:
There was 1 lot of dichroics provided, 1 dichroic was tested at 801 measurement points.

TODO update historical test plan.

COLLIMATOR:
In-house testing of a collimator was done, with 1501 measurement points.

LENSES:

TODO I am not sure how camera testing was done for LLAMAS, check Box for that.

Also located here are test results for detailed camera testing of the protoLLAMAS lenses. Here, the camera lens assemblies were tested all together, at 11 points each for the red and blue cameras.

FRD:
From evaluating_cleavage...
The average of the sample of 10 stripped is 92.3%, and the scatter is 2.4%.
Previously I used this as the prior, but this wouldn't have been known before testing.
This is more appropriately used as validation data.