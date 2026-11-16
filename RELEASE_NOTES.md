# CubeLFM-IGRF v0.2.0

Independent public export of the accepted Earth64 IGRF/tilt/12-hour rotating-MI source snapshot. All49 selected C/Fortran translation units and12 headers match the accepted source hashes.

The full12-hour case passed its frozen numerical audit and complete MHD+MI word comparison, at4.003764681173281× real time on96 cores. MI remains native LFM100/1e-3 with honest warning/continue for finite nonconvergence. Rotation period is12h, Bz is+5nT for2h then−5nT for10h, SigmaP is5S.

The6–7× optimization is ongoing. A portable clean rebuild is not yet qualified; see README and BUILD_DEPENDENCIES.json. This repository contains the accepted source snapshot with independent history.
