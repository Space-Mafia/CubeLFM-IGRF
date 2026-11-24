"""GSE/GSM to the fixed-spin-axis GAMERA event frame.

This is a proper, right-handed rotation for each vector, not an X reflection.
GSE velocity and GSM IMF must use their respective bases.  The native frame
has +Z along Earth's spin axis and -X along the initial Sun projection.
"""

from __future__ import annotations

from datetime import datetime
import math
from typing import Mapping

from prepare_event_model import solar_and_gmst, utc


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def unit(a):
    n = math.sqrt(dot(a, a))
    if not n > 1e-12:
        raise ValueError("degenerate event coordinate basis")
    return tuple(x / n for x in a)


def model_axes(start: datetime):
    ra, _, _, _ = solar_and_gmst(start)
    return ((-math.cos(ra), -math.sin(ra), 0.0),
            (math.sin(ra), -math.cos(ra), 0.0),
            (0.0, 0.0, 1.0))


def inertial_to_model(v, axes):
    return tuple(dot(v, axis) for axis in axes)


def basis_at(value: datetime, manifest: Mapping):
    start = utc(manifest["start_utc"])
    ra, dec, _, eps = solar_and_gmst(value)
    axes = model_axes(start)
    sun_eq = (math.cos(dec) * math.cos(ra),
              math.cos(dec) * math.sin(ra), math.sin(dec))
    x = unit(inertial_to_model(sun_eq, axes))
    ecliptic_north_eq = (0.0, -math.sin(eps), math.cos(eps))
    ecliptic_north = unit(inertial_to_model(ecliptic_north_eq, axes))
    gse_y = unit(cross(ecliptic_north, x))
    gse_z = unit(cross(x, gse_y))

    angle = (manifest["initial_rotation_phase_rad"] +
             7.2921150e-5 * (value - start).total_seconds())
    d = manifest["dipole_axis_body"]
    dipole = (math.cos(angle) * d[0] - math.sin(angle) * d[1],
              math.sin(angle) * d[0] + math.cos(angle) * d[1], d[2])
    projection = dot(dipole, x)
    gsm_z = unit(tuple(dipole[i] - projection * x[i] for i in range(3)))
    gsm_y = unit(cross(gsm_z, x))
    gsm_z = unit(cross(x, gsm_y))
    return (x, gse_y, gse_z), (x, gsm_y, gsm_z)


def vector_in_model(components, basis):
    return tuple(sum(components[j] * basis[j][i] for j in range(3))
                 for i in range(3))
