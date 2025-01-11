FOCAL PLANES:
For the purposes of priors, the SN20001 focal plane is the prior for red and green, and the SN20003 focal plane is the prior for blue. These were the  "qualification" units used in protoLLAMAS.

TODO check to make sure the values in this spreadsheet are used exactly in my code.

QUANTUM EFFICIENCY:
According to the qualification report, all focal planes are ccd42-40, different coating processes which result in different QE:
Red cameras use "NIR"
Green cameras use "Basic Midband"
Blue cameras use "Enhanced Broadband"

These correspond to different QE curves in ccd42-40, which correspond to the different throughput files.

TODO check to make sure those indeed match closely.

VPH:
Wasach is the vendor. I assume these coating files, part of the ETC, came from the vendor.

SL DICHROIC:
Kind of a design prior. The ETC included an ECI FusedSilica coating curve to represent the dichroics. I modified this one to build in the dichroic wavelength in the design.

BG DICHROIC:
Kind of a design prior. The ETC included an ECI FusedSilica coating curve to represent the dichroics. I modified this one to build in the dichroic wavelength in the design.

COLLIMATOR:
dielectric_mirror.txt was the coating file originally included in ETC. Its provenance is unknown, and there's no strong reason to think it represents the best prior knowledge on the LLAMAS collimator.

The ECI coating curve file, LLAMAS_collimator_coatings_ECI_20210504, is from the vendor. This is the better prior.

TODO get a coating curve out of this?

LENSES:
We have Optimax throughput curves from the vendor for each lens in a spectrograph.

TODO sadly, I have to turn all of these into throughput curves and use em

FRD:
Previously I used evaluating_cleaving_through_bugger.xlsx as the prior, but this wouldn't have been known before testing.
So instead, I'll use a design prior. Requirements say throughput loss has to be less than 5%. So, I'll make the prior a gamma distribution that has a CDF near 1.0 at x=5, and is somewhat linear below that.
This gives:
gamma(alpha=1.575, beta=1.25)

TODO implement that change