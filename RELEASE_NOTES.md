# CubeLFM-IGRF v0.3.0

The 128 × 48 × 48 × 6 native-seam Earth event release adds the corrected CT induction reference, the YYG global PARDISO potential solver and its branch-local electrostatic wall return. Event-driven EUV and precipitation conductance, the dated IGRF orientation and compact MHD/MI output are included.

The Sugon cold start on 192 physical cores completed **04:00–16:00 UT on 23 October 2002**. It advanced 43,200 s in 495,712 steps. The complete MHD process took 87,718.14 s (**0.49249× real time**).

Use `python3 event128/verify.py` from a fresh checkout to reproduce and check the 14 compiled inputs, 47 inherited source mappings and both compact-output writers. The [event guide](event128/README.md) lists cluster prerequisites, submission and output fields. Full runs currently require the compatible native Sugon installation. The earlier 64-grid release is available at its original tag.

# CubeLFM-IGRF v0.2.0

Independent public export of the accepted Earth64 IGRF/tilt/12-hour rotating-MI source snapshot. All49 selected C/Fortran translation units and12 headers match the accepted source hashes.

The full12-hour case passed its frozen numerical audit and complete MHD+MI word comparison, at4.003764681173281× real time on96 cores. MI remains native LFM100/1e-3 with honest warning/continue for finite nonconvergence. Rotation period is12h, Bz is+5nT for2h then−5nT for10h, SigmaP is5S.

The6–7× optimization is ongoing. A portable clean rebuild is not yet qualified; see README and BUILD_DEPENDENCIES.json. This repository contains the accepted source snapshot with independent history.
