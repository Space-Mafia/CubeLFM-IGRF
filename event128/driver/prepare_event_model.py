#!/usr/bin/env python3
"""Prepare a date-specific degree-5 IGRF and phase for one fresh event run.

The model frame has +Z along Earth's spin axis and -X along the *initial*
equatorial projection of the Sun.  No mesh tilt is introduced.  This script
does not silently turn an existing run into an event run: it creates an
immutable model directory that must be named at CMake configure time.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys


COEFF_SHA256 = "8f8d88403028fc4ee92c4f38d97b46e0a87e2cfc496045b43c9e26c1d6b0903c"
COEFF_FILE = Path(__file__).resolve().parents[2] / "examples/earth_magnetosphere/common/igrf14coeffs.txt"
DEGREE = 5


def utc(text: str) -> datetime:
    try:
        value = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError("start_utc must be YYYY-MM-DDTHH:MM:SSZ") from exc
    if value.tzinfo is None or value.utcoffset().total_seconds() != 0:
        raise ValueError("start_utc must explicitly specify UTC")
    if value.microsecond:
        raise ValueError("start_utc must be an integral UTC second")
    if not datetime(1900, 1, 1, tzinfo=timezone.utc) <= value <= datetime(2030, 1, 1, tzinfo=timezone.utc):
        raise ValueError("IGRF-14 date must be within 1900-01-01 to 2030-01-01")
    return value


def julian_day(value: datetime) -> float:
    return value.timestamp() / 86400.0 + 2440587.5


def solar_and_gmst(value: datetime) -> tuple[float, float, float, float]:
    """Approximate geocentric solar RA/declination and Greenwich mean sidereal angle.

    Meeus low-precision solar position and IAU-1982 GMST are adequate for an
    initial orientation.  UTC is used as an approximation to UT1; the
    resulting sub-degree uncertainty is documented in event_model.json.
    """
    jd = julian_day(value)
    t = (jd - 2451545.0) / 36525.0
    m = math.radians((357.52911 + t * (35999.05029 - 0.0001537 * t)) % 360.0)
    l0 = (280.46646 + t * (36000.76983 + 0.0003032 * t)) % 360.0
    c = ((1.914602 - t * (0.004817 + 0.000014 * t)) * math.sin(m)
         + (0.019993 - 0.000101 * t) * math.sin(2 * m)
         + 0.000289 * math.sin(3 * m))
    omega = math.radians(125.04 - 1934.136 * t)
    lam = math.radians(l0 + c - 0.00569 - 0.00478 * math.sin(omega))
    eps0 = 23.0 + (26.0 + (21.448 - t * (46.8150 + t * (0.00059 - 0.001813 * t))) / 60.0) / 60.0
    eps = math.radians(eps0 + 0.00256 * math.cos(omega))
    ra = math.atan2(math.cos(eps) * math.sin(lam), math.cos(lam)) % (2 * math.pi)
    dec = math.asin(math.sin(eps) * math.sin(lam))
    gmst = math.radians((280.46061837 + 360.98564736629 * (jd - 2451545.0)
                         + 0.000387933 * t * t - t**3 / 38710000.0) % 360.0)
    return ra, dec, gmst, eps


def coefficients(value: datetime) -> tuple[list[list[float]], list[list[float]], dict]:
    raw = COEFF_FILE.read_bytes()
    if hashlib.sha256(raw).hexdigest() != COEFF_SHA256:
        raise ValueError("official IGRF-14 coefficient asset checksum mismatch")
    lines = raw.decode("ascii").splitlines()
    year_line = next(line for line in lines if line.startswith("g/h n m "))
    epochs = [float(x) for x in year_line.split()[3:-1]]
    if epochs != [float(y) for y in range(1900, 2030, 5)]:
        raise ValueError("unexpected IGRF-14 epoch columns")
    g = [[0.0] * 15 for _ in range(15)]
    h = [[0.0] * 15 for _ in range(15)]
    year = value.year
    if year >= 2025:
        lo = 2025
        hi = 2030
        fraction = (value - datetime(2025, 1, 1, tzinfo=timezone.utc)).total_seconds() / (
            datetime(2030, 1, 1, tzinfo=timezone.utc) - datetime(2025, 1, 1, tzinfo=timezone.utc)).total_seconds()
    else:
        lo = 1900 + 5 * ((year - 1900) // 5)
        hi = lo + 5
        fraction = (value - datetime(lo, 1, 1, tzinfo=timezone.utc)).total_seconds() / (
            datetime(hi, 1, 1, tzinfo=timezone.utc) - datetime(lo, 1, 1, tzinfo=timezone.utc)).total_seconds()
    seen = set()
    for line in lines:
        fields = line.split()
        if len(fields) < 5 or fields[0] not in ("g", "h"):
            continue
        n, m = int(fields[1]), int(fields[2])
        if n > DEGREE:
            continue
        if (fields[0], n, m) in seen:
            raise ValueError("duplicate IGRF coefficient")
        seen.add((fields[0], n, m))
        data = [float(v) for v in fields[3:]]
        if len(data) != len(epochs) + 1:
            raise ValueError("unexpected IGRF coefficient row length")
        lower = data[epochs.index(lo)]
        upper = lower + 5.0 * data[-1] if lo == 2025 else data[epochs.index(hi)]
        (g if fields[0] == "g" else h)[n][m] = lower + fraction * (upper - lower)
    if len(seen) != 35:
        raise ValueError(f"expected 35 degree-1..5 g/h coefficients, got {len(seen)}")
    return g, h, {"lower_epoch": lo, "upper_epoch": hi, "fraction": fraction,
                   "method": "linear IGRF interpolation" if lo < 2025 else "IGRF 2025 plus 2025-30 SV"}


def header(g: list[list[float]], h: list[list[float]]) -> str:
    result = ["/* Generated from official IGRF-14, degrees 1-5 only. */"]
    for name, values in (("g", g), ("h", h)):
        result.append(f"static const double sh_coeffs_{name}[SH_N_MAX + 1][SH_N_MAX + 1] = {{")
        for n in range(15):
            row = ", ".join(format(values[n][m], ".17g") for m in range(n + 1)) if n <= DEGREE else "0.0"
            result.append("    {" + row + "},")
        result.append("};")
    return "\n".join(result) + "\n"


def prepare(start: datetime, output: Path, table: bool = True) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite prepared model: {output}")
    output.mkdir(parents=True)
    try:
        g, h, interpolation = coefficients(start)
        ra, dec, gmst, eps = solar_and_gmst(start)
        phase = (gmst - ra + math.pi) % (2 * math.pi)
        body_sun = [-math.cos(dec), 0.0, math.sin(dec)]
        canonical_utc = start.strftime("%Y-%m-%dT%H:%M:%SZ")
        (output / "event_model_config.h").write_text(
            "/* Generated event orientation. Do not edit by hand. */\n"
            f'#define GAMERA_EVENT_START_UTC "{canonical_utc}"\n'
            f"#define GAMERA_EVENT_PHASE0_RAD {phase:.17g}\n"
            f"#define GAMERA_EVENT_SUN_X {body_sun[0]:.17g}\n"
            f"#define GAMERA_EVENT_SUN_Y {body_sun[1]:.17g}\n"
            f"#define GAMERA_EVENT_SUN_Z {body_sun[2]:.17g}\n")
        coeff_header = output / "event_igrf_coeffs.h"
        coeff_header.write_text(header(g, h))
        if table:
            subprocess.run([sys.executable, str(Path(__file__).with_name("make_cgm_table.py")),
                            "--coeff-header", str(coeff_header), "--output-dir", str(output)], check=True)
        manifest = {"start_utc": canonical_utc,
                    "igrf14_asset_sha256": COEFF_SHA256, "degree": DEGREE,
                    "coeff_header_sha256": hashlib.sha256(coeff_header.read_bytes()).hexdigest(),
                    "orientation_header_sha256": hashlib.sha256((output / "event_model_config.h").read_bytes()).hexdigest(),
                    "interpolation": interpolation, "initial_rotation_phase_rad": phase,
                    "sunward_direction_model": body_sun,
                    "dipole_axis_body": [v / math.sqrt(sum(u*u for u in (-g[1][1], -h[1][1], -g[1][0])))
                                         for v in (-g[1][1], -h[1][1], -g[1][0])],
                    "solar_right_ascension_rad": ra, "solar_declination_rad": dec,
                    "gmst_rad": gmst, "obliquity_rad": eps,
                    "frame": "+Z spin; -X initial projected Sun; right-handed Y=Z cross X",
                    "geometry_limit": "Sun fixed at start; solar RA/declination drift during run omitted; UTC approximates UT1",
                    "cgm_table_sha256": hashlib.sha256((output / "event_cgm_table.h").read_bytes()).hexdigest() if table else None}
        (output / "event_model.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        return manifest
    except Exception:
        # Leave the failed candidate and its generated evidence for inspection.
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-utc", required=True, help="YYYY-MM-DDTHH:MM:SSZ")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--skip-cgm-for-unit-test", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    result = prepare(utc(args.start_utc), args.output.resolve(), not args.skip_cgm_for_unit_test)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
