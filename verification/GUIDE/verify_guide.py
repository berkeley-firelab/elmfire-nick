#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated verification harness for the ELMFIRE "GUIDE" template cases.

This script:
  1. Copies every case in ./TEMPLATE into a fresh run directory.
  2. Generates the input rasters each case needs (126 x 126 @ 30 m), built
     from the parameters in Table 3.1 of the ELMFIRE guide.
  3. Runs ELMFIRE for every case (and sub-case).
  4. Extracts the relevant quantity from the outputs (head ROS, fireline
     intensity, crown-fire class, stop time, ember statistics, ...).
  5. Compares against the theoretical / BEHAVE-derived solutions reported in
     Chapter 3 of the guide.
  6. Writes a CSV summary, a markdown report, and prints a console table.

The theoretical targets are taken verbatim from the guide's verification
chapter. Where Table 3.1 and the narrative text disagree on an input (e.g.
the "Windy" wind speed), the value that reproduces the published target is
used and the discrepancy is noted in CASES below.

Requires (conda env "elmfire"): rasterio, numpy, pandas, and an `elmfire`
binary on PATH (or set ELMFIRE_BIN).

Usage:
    python verify_guide.py                  # generate inputs, run, compare
    python verify_guide.py --cases Point Windy
    python verify_guide.py --skip-run       # only re-compare existing outputs
    python verify_guide.py --no-generate    # reuse inputs already present
    python verify_guide.py --ia-ensemble 2000   # shrink the heavy suppression run
"""

import argparse
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

try:
    import rasterio
    from rasterio.transform import from_origin
    import pandas as pd
except ImportError as exc:  # pragma: no cover
    sys.exit(f"Missing python dependency ({exc}). Activate the elmfire conda env first.")


# ---------------------------------------------------------------------------
# Static configuration
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.join(HERE, "TEMPLATE")

NCOLS = NROWS = 126               # All guide rasters are 126 x 126 ...
CELLSIZE = 30.0                   # ... at 30 m resolution.
XLL = YLL = 0.0                   # Lower-left domain corner.
# Embed the EPSG code directly (UTM zone 35N) rather than a proj4 string, so
# ELMFIRE reads a clean EPSG from the input rasters without a proj4->EPSG
# database lookup (which is fragile across GDAL/PROJ installations).
CRS_EPSG = "EPSG:32635"
NODATA = -9999.0

# Cell-index of the domain centre (used for quadrant masks / point ignition).
HALF = NROWS // 2                 # 63

# Default pass/fail tolerance on a percentage error (matches FBP harness).
DEFAULT_TOL_PCT = 10.0


# ---------------------------------------------------------------------------
# Raster grid helpers
# ---------------------------------------------------------------------------
def _transform():
    # Origin is the upper-left corner; rows increase downward (decreasing y).
    top = YLL + NROWS * CELLSIZE
    return from_origin(XLL, top, CELLSIZE, CELLSIZE)


def _quadrant_masks():
    """Return NW, NE, SW, SE boolean masks on the (row, col) grid.

    Row 0 is the north (top) edge because the raster origin is upper-left.
    """
    rows = np.arange(NROWS)[:, None] * np.ones((1, NCOLS))
    cols = np.ones((NROWS, 1)) * np.arange(NCOLS)[None, :]
    north = rows < HALF
    west = cols < HALF
    return {
        "NW": north & west,
        "NE": north & ~west,
        "SW": ~north & west,
        "SE": ~north & ~west,
    }


def _half_masks():
    """Return east/west half masks (split on the central column)."""
    cols = np.ones((NROWS, 1)) * np.arange(NCOLS)[None, :]
    west = cols < HALF
    return {"W": west, "E": ~west}


def write_raster(path, array, dtype, nbands=1):
    """Write a (single- or multi-band) GeoTIFF. `array` is 2-D (row, col)."""
    arr = np.asarray(array)
    profile = dict(
        driver="GTiff",
        height=NROWS,
        width=NCOLS,
        count=nbands,
        dtype=dtype,
        crs=CRS_EPSG,
        transform=_transform(),
        nodata=NODATA,
        compress="deflate",
    )
    with rasterio.open(path, "w", **profile) as dst:
        for b in range(1, nbands + 1):
            dst.write(arr.astype(dtype), b)


# ---------------------------------------------------------------------------
# Case specification
# ---------------------------------------------------------------------------
@dataclass
class Layer:
    """A constant, quadrant-varying, or half-varying input field.

    `value` is one of:
        scalar                       -> uniform field
        {"NW":.., "NE":.., ...}      -> per-quadrant field
        {"E":.., "W":..}             -> per-half field
    """
    value: object

    def to_array(self):
        if isinstance(self.value, dict):
            arr = np.zeros((NROWS, NCOLS), dtype=float)
            if set(self.value) <= {"NW", "NE", "SW", "SE"}:
                masks = _quadrant_masks()
            elif set(self.value) <= {"E", "W"}:
                masks = _half_masks()
            else:
                raise ValueError(f"Unrecognized region keys: {set(self.value)}")
            for k, v in self.value.items():
                arr[masks[k]] = v
            return arr
        return np.full((NROWS, NCOLS), float(self.value))


@dataclass
class Target:
    name: str            # human label, e.g. "head ROS (15 deg)"
    value: float         # theoretical value
    units: str
    region: Optional[str] = None   # quadrant/half key to restrict the metric
    tol_pct: float = DEFAULT_TOL_PCT
    abs_tol: Optional[float] = None  # absolute tolerance (overrides pct if set)


@dataclass
class Case:
    name: str
    data_rel: str                       # path to the .data file relative to TEMPLATE
    layers: dict = field(default_factory=dict)   # raster var -> Layer
    targets: list = field(default_factory=list)
    metric: str = "max_ros"             # extractor key (see METRICS)
    notes: str = ""
    quantitative: bool = True


# Integer-valued rasters (everything else is written Float32).
INT_LAYERS = {"slp", "asp", "dem", "fbfm", "cc", "ch", "cbh", "cbd"}
# Rasters that live in the weather directory and are read band-by-band.
WEATHER_LAYERS = {"ws", "wd", "m1", "m10", "m100", "lh", "lw"}

# Default constant fields shared by most cases (overridden per case as needed).
BASE_LAYERS = dict(
    slp=0, asp=0, dem=0, cc=0, ch=0, cbh=0, cbd=0,
    adj=1.0, phi=1.0, ws=0.0, wd=0.0,
    m1=6.0, m10=7.0, m100=8.0, lh=30.0, lw=60.0,
    fbfm=3,
)


def L(**overrides):
    """Build a {var: Layer} dict from BASE_LAYERS plus overrides."""
    merged = dict(BASE_LAYERS)
    merged.update(overrides)
    return {k: Layer(v) for k, v in merged.items()}


# ---------------------------------------------------------------------------
# The verification cases (targets from guide Chapter 3)
# ---------------------------------------------------------------------------
CASES = [
    Case(
        name="Point",
        data_rel="Point/point.data",
        layers=L(fbfm=3, ws=0.0, m1=6.0),
        targets=[Target("head ROS", 1.51, "m/min")],
        metric="max_ros",
        notes="Zero wind/slope -> circular front; BEHAVE max ROS 1.51 m/min.",
    ),
    Case(
        name="Windy",
        data_rel="Windy/windy.data",
        # Table 3.1 (actual sim parameter) = 5 mph 20-ft wind. The guide's BEHAVE
        # target used 6 mph + WAF 0.418; ELMFIRE computes WAF internally.
        layers=L(fbfm=5, ws=6.0, m1=6.0, m10=7.0),
        targets=[Target("head ROS", 3.06, "m/min")],
        metric="max_ros",
        notes=("Wind-driven ellipse (5 mph 20-ft per Table 3.1). Residual error is "
               "a WAF convention gap: ELMFIRE computes WAF internally, the guide's "
               "BEHAVE target used WAF=0.418. Tune wind/WS_AT_10M to match the guide run."),
    ),
    Case(
        name="Valley",
        data_rel="Valley/valley.data",
        # Aspect = downslope azimuth. East half descends west (asp 270) so its
        # upslope head fire develops eastward (outward, away from the central
        # ignition); west half descends east (asp 90) -> head fire westward.
        layers=L(fbfm=3, ws=0.0, m1=6.0,
                 slp={"E": 5, "W": 15}, asp={"E": 270, "W": 90}),
        targets=[
            Target("head ROS (5 deg, east)", 1.92, "m/min", region="E"),
            Target("head ROS (15 deg, west)", 5.37, "m/min", region="W"),
        ],
        metric="max_ros_region",
        notes="Slope-driven spread; SLP raster in degrees; two slopes vs BEHAVE.",
    ),
    Case(
        name="Quadrant",
        data_rel="Quadrant/quadrant.data",
        layers=L(fbfm={"NW": 8, "NE": 7, "SW": 4, "SE": 2},
                 ws=0.0, m1=6.0, m10=7.0, m100=8.0, lh=60.0, lw=90.0),
        targets=[
            Target("FB8 head ROS", 0.080, "m/min", region="NW"),
            Target("FB7 head ROS", 0.472, "m/min", region="NE"),
            Target("FB4 head ROS", 1.492, "m/min", region="SW"),
            Target("FB2 head ROS", 0.871, "m/min", region="SE"),
        ],
        metric="max_ros_region",
        notes="Four fuel models; per-quadrant BEHAVE ROS.",
    ),
    Case(
        name="MoistureQuad",
        data_rel="Moisture Quad/moisture_quad.data",
        layers=L(fbfm=3, ws=0.0,
                 m1={"NW": 6, "NE": 12, "SW": 18, "SE": 25}),
        targets=[
            Target("ROS @ 6%", 1.511, "m/min", region="NW"),
            Target("ROS @ 12%", 1.089, "m/min", region="NE"),
            Target("ROS @ 18%", 0.801, "m/min", region="SW"),
            Target("ROS @ 25%", 0.000, "m/min", region="SE", abs_tol=0.05),
        ],
        metric="max_ros_region",
        notes="Uniform fuel, varying 1-h FMC; 25% = moisture of extinction (no spread).",
    ),
    Case(
        name="Canopy-2mph",
        data_rel="Canopy/2/canopy_lowwind.data",
        layers=L(fbfm=5, ws=2.0, m1=6.0, m10=7.0,
                 cc=80, ch=16, cbh=16, cbd=55),
        targets=[Target("head ROS (no crown)", 1.96, "m/min")],
        metric="max_ros_crown",
        notes="Surface fire only; insufficient intensity for crown fire.",
    ),
    Case(
        name="Canopy-3mph",
        data_rel="Canopy/3/canopy_mediumwind.data",
        layers=L(fbfm=5, ws=3.0, m1=6.0, m10=7.0,
                 cc=80, ch=16, cbh=16, cbd=55),
        targets=[Target("head ROS (passive crown)", 6.1, "m/min")],
        metric="max_ros_crown",
        notes="Passive crown fire expected toward the head.",
    ),
    Case(
        name="Canopy-6mph",
        data_rel="Canopy/6/canopy_highwind.data",
        layers=L(fbfm=5, ws=6.0, m1=6.0, m10=7.0,
                 cc=80, ch=16, cbh=16, cbd=55),
        targets=[Target("head ROS (active crown)", 21.3, "m/min")],
        metric="max_ros_crown",
        notes="Active crown fire expected; theoretical RC 21.5 m/min.",
    ),
    Case(
        name="Complex",
        data_rel="Complex/complex.data",
        layers=L(fbfm={"NW": 8, "NE": 7, "SW": 4, "SE": 2},
                 ws=5.0, m1=6.0, m10=7.0, m100=8.0, lh=60.0, lw=90.0,
                 slp={"E": 5, "W": 15}, asp={"E": 270, "W": 90}),
        targets=[],
        metric="qualitative_ros",
        quantitative=False,
        notes="No analytical solution; reported for qualitative review only.",
    ),
    Case(
        name="Firebrands",
        data_rel="Firebrands/firebrands.data",
        layers=L(fbfm=3, ws=5.0, m1=6.0),
        targets=[Target("max spotting distance", 182.5, "m", tol_pct=40.0)],
        metric="spotting_stats",
        notes="Lognormal ember transport; top-1% downwind distance ~182.5 m (best-effort).",
    ),
    Case(
        name="Overnight",
        data_rel="Overnight/overnight.data",
        layers=L(fbfm=3, ws=0.0, m1=6.0),
        targets=[Target("night/day growth ratio", 0.1, "-", abs_tol=0.15)],
        metric="overnight",
        notes="Diurnal slowdown (x0.1) between sunset/sunrise (best-effort).",
    ),
    Case(
        name="Suppression-extended",
        data_rel="Suppression/extended_attack/extended_attack.data",
        layers=L(fbfm=3, ws=0.0, m1=6.0),
        targets=[Target("containment stop time", 10800.0, "s", tol_pct=15.0)],
        metric="stop_time",
        notes="Analytical containment at ~10800 s (ELMFIRE reported ~10000 s).",
    ),
    Case(
        name="Suppression-initial",
        data_rel="Suppression/initial_attack/initial_attack.data",
        layers=L(fbfm=3, ws=0.0, m1=6.0),
        targets=[Target("fraction contained", 0.54, "-", abs_tol=0.05)],
        metric="initial_attack",
        notes="20000-member ensemble; ~52% contained vs 54% predicted (HEAVY).",
    ),
    Case(
        name="BANDTHICKNESS-1",
        data_rel="BandThickness/1/bandthickness1.data",
        layers=L(fbfm=3, ws=5.0, m1=6.0),
        targets=[],
        metric="consistency_ros",
        quantitative=False,
        notes="Band-thickness invariance check (compared across 1/2/3).",
    ),
    Case(
        name="BANDTHICKNESS-2",
        data_rel="BandThickness/2/bandthickness2.data",
        layers=L(fbfm=3, ws=5.0, m1=6.0),
        targets=[],
        metric="consistency_ros",
        quantitative=False,
        notes="Band-thickness invariance check (compared across 1/2/3).",
    ),
    Case(
        name="BANDTHICKNESS-3",
        data_rel="BandThickness/3/bandthickness3.data",
        layers=L(fbfm=3, ws=5.0, m1=6.0),
        targets=[],
        metric="consistency_ros",
        quantitative=False,
        notes="Band-thickness invariance check (compared across 1/2/3).",
    ),
]


# ---------------------------------------------------------------------------
# .data namelist parsing / patching
# ---------------------------------------------------------------------------
def parse_namelist(path):
    """Very small namelist reader: returns {KEY: raw_string_value}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("!")[0].strip()
            if "=" not in line or line.startswith("&") or line.startswith("/"):
                continue
            key, _, val = line.partition("=")
            out[key.strip().upper()] = val.strip()
    return out


def _unquote(s):
    return s.strip().strip("'").strip('"').strip()


def case_input_filenames(nml):
    """Map each input raster variable -> the base filename the .data expects."""
    mapping = {
        "ws": "WS_FILENAME", "wd": "WD_FILENAME",
        "m1": "M1_FILENAME", "m10": "M10_FILENAME", "m100": "M100_FILENAME",
        "lh": "MLH_FILENAME", "lw": "MLW_FILENAME",
        "asp": "ASP_FILENAME", "slp": "SLP_FILENAME", "dem": "DEM_FILENAME",
        "fbfm": "FBFM_FILENAME", "cc": "CC_FILENAME", "ch": "CH_FILENAME",
        "cbh": "CBH_FILENAME", "cbd": "CBD_FILENAME",
        "adj": "ADJ_FILENAME", "phi": "PHI_FILENAME",
    }
    files = {}
    for var, key in mapping.items():
        if key in nml:
            files[var] = _unquote(nml[key])
    return files


def num_weather_bands(nml):
    try:
        return max(1, int(float(nml.get("NUM_METEOROLOGY_TIMES", "1"))))
    except ValueError:
        return 1


def override_namelist_keys(path, overrides):
    """Override given KEY=value lines in the run copy, in place.

    The committed templates are never touched. Applied to each run copy:
      * PATH_TO_GDAL        -> a GDAL that actually works in this environment
                               (auto-detected by find_working_gdal(), or set via
                               --gdal-path / ELMFIRE_GDAL). The templates ship
                               '/usr/bin/', which is correct on a healthy system
                               but may be broken on a given box (system/conda PROJ
                               db mismatch), so we resolve a working one at run time.
      * NUM_ENSEMBLE_MEMBERS -> shrink the heavy initial-attack run (--ia-ensemble)
    """
    if not overrides:
        return
    out = []
    with open(path) as fh:
        for line in fh:
            key = line.split("=", 1)[0].strip().upper()
            if key in overrides:
                out.append(f"{key} = {overrides[key]}\n")
            else:
                out.append(line)
    with open(path, "w") as fh:
        fh.writelines(out)


_WORKING_GDAL = None  # cache: None=not probed, ""=none found, str=working bin dir


def _gdal_dir_works(bindir, test_tif, tmpdir):
    """True if this GDAL bin dir can BOTH identify and apply EPSG:32635 in the
    current environment -- the two operations ELMFIRE relies on for georeferenced
    output (gdalsrsinfo detection, then gdal_translate -a_srs)."""
    translate = os.path.join(bindir, "gdal_translate")
    srsinfo = os.path.join(bindir, "gdalsrsinfo")
    if not (os.path.exists(translate) and os.path.exists(srsinfo)):
        return False
    try:
        r = subprocess.run([srsinfo, "-o", "epsg", test_tif],
                           capture_output=True, text=True, timeout=30)
        if r.returncode != 0 or "32635" not in r.stdout or "PROJ" in r.stderr:
            return False
        out = os.path.join(tmpdir, "probe_out.tif")
        r2 = subprocess.run([translate, "-q", "-a_srs", "EPSG:32635", test_tif, out],
                            capture_output=True, text=True, timeout=30)
        if r2.returncode != 0 or "PROJ" in r2.stderr or "Failed to process SRS" in r2.stderr:
            return False
    except (OSError, subprocess.SubprocessError):
        return False
    return True


def find_working_gdal():
    """Probe candidate GDAL bin dirs and return the first one (with trailing
    separator) that works in this environment, or None if none do.

    Tries, in order: the gdal_translate on PATH (typically the active conda env),
    then /usr/bin and /usr/local/bin. Result is cached for the process."""
    global _WORKING_GDAL
    if _WORKING_GDAL is not None:
        return _WORKING_GDAL or None

    import tempfile
    tmpdir = tempfile.mkdtemp(prefix="gdalprobe_")
    test_tif = os.path.join(tmpdir, "probe.tif")
    profile = dict(driver="GTiff", height=2, width=2, count=1, dtype="float32",
                   crs=CRS_EPSG, transform=from_origin(0.0, 60.0, 30.0, 30.0),
                   nodata=NODATA)
    with rasterio.open(test_tif, "w", **profile) as dst:
        dst.write(np.ones((2, 2), "float32"), 1)

    candidates = []
    onpath = shutil.which("gdal_translate")
    if onpath:
        candidates.append(os.path.dirname(onpath))
    candidates += ["/usr/bin", "/usr/local/bin"]

    seen, ordered = set(), []
    for c in candidates:
        c = os.path.normpath(c)
        if c not in seen:
            seen.add(c)
            ordered.append(c)

    _WORKING_GDAL = ""
    for d in ordered:
        if _gdal_dir_works(d, test_tif, tmpdir):
            _WORKING_GDAL = d + os.sep
            break
    shutil.rmtree(tmpdir, ignore_errors=True)
    return _WORKING_GDAL or None


# ---------------------------------------------------------------------------
# Input generation
# ---------------------------------------------------------------------------
def generate_case_inputs(case_dir, data_path, case):
    nml = parse_namelist(data_path)
    files = case_input_filenames(nml)
    nbands = num_weather_bands(nml)

    fueltop = _unquote(nml.get("FUELS_AND_TOPOGRAPHY_DIRECTORY", "./inputs"))
    weather = _unquote(nml.get("WEATHER_DIRECTORY", "./inputs"))
    fueltop_dir = os.path.normpath(os.path.join(case_dir, fueltop))
    weather_dir = os.path.normpath(os.path.join(case_dir, weather))
    os.makedirs(fueltop_dir, exist_ok=True)
    os.makedirs(weather_dir, exist_ok=True)

    for var, fname in files.items():
        layer = case.layers.get(var)
        if layer is None:
            # Field not configured for this case; fall back to base default.
            layer = Layer(BASE_LAYERS.get(var, 0))
        arr = layer.to_array()
        dtype = "int16" if var in INT_LAYERS else "float32"
        bands = nbands if var in WEATHER_LAYERS else 1
        target_dir = weather_dir if var in WEATHER_LAYERS else fueltop_dir
        write_raster(os.path.join(target_dir, f"{fname}.tif"), arr, dtype, bands)


# ---------------------------------------------------------------------------
# Running ELMFIRE
# ---------------------------------------------------------------------------
def run_elmfire(case_dir, data_path, elmfire_bin):
    rel_data = os.path.basename(data_path)
    for sub in ("outputs", "scratch"):
        d = os.path.join(case_dir, sub)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d, exist_ok=True)
    log = os.path.join(case_dir, "outputs", "elmfire_run.log")
    with open(log, "w") as lf:
        proc = subprocess.run([elmfire_bin, rel_data], cwd=case_dir,
                              stdout=lf, stderr=subprocess.STDOUT)
    ok = False
    with open(log) as lf:
        ok = "End of simulation reached successfully" in lf.read()
    return ok, log, proc.returncode


# ---------------------------------------------------------------------------
# Output metric extraction
# ---------------------------------------------------------------------------
def _read_stack(paths):
    """Read one or more single-band rasters into a list of masked 2-D arrays."""
    arrays = []
    for p in sorted(paths):
        with rasterio.open(p) as src:
            a = src.read(1, masked=True).astype(float)
            arrays.append(a)
    return arrays


def _glob(out_dir, prefix):
    import glob
    return glob.glob(os.path.join(out_dir, f"{prefix}*.tif"))


def _max_over(arrays, mask=None):
    best = None
    for a in arrays:
        data = a.filled(np.nan)
        if mask is not None:
            data = np.where(mask, data, np.nan)
        data = data[np.isfinite(data)]
        if data.size:
            m = float(np.nanmax(data))
            best = m if best is None else max(best, m)
    return best


def metric_max_ros(out_dir, case):
    arrs = _read_stack(_glob(out_dir, "vs_"))
    return {"head ROS": _max_over(arrs)}


def metric_max_ros_region(out_dir, case):
    arrs = _read_stack(_glob(out_dir, "vs_"))
    masks = {**_quadrant_masks(), **_half_masks()}
    res = {}
    for t in case.targets:
        res[t.name] = _max_over(arrs, masks.get(t.region))
    return res


def metric_max_ros_crown(out_dir, case):
    arrs = _read_stack(_glob(out_dir, "vs_"))
    crown = _read_stack(_glob(out_dir, "crown_fire_"))
    crown_max = _max_over(crown) if crown else None
    label = case.targets[0].name
    return {label: _max_over(arrs), "_crown_class_max": crown_max}


def metric_qualitative_ros(out_dir, case):
    arrs = _read_stack(_glob(out_dir, "vs_"))
    return {"max ROS (qualitative)": _max_over(arrs)}


def metric_consistency_ros(out_dir, case):
    arrs = _read_stack(_glob(out_dir, "vs_"))
    return {"max ROS": _max_over(arrs)}


def metric_stop_time(out_dir, case):
    toa = _read_stack(_glob(out_dir, "time_of_arrival_"))
    return {"containment stop time": _max_over(toa)}


def metric_overnight(out_dir, case):
    """Estimate the ratio of mean night-time to day-time fire-area growth.

    Uses the time-of-arrival raster: counts newly burned cells per hour, then
    compares the average hourly growth during the night window (per guide,
    ~19:00-11:00 local) against the daytime window.
    """
    toa = _glob(out_dir, "time_of_arrival_")
    if not toa:
        return {"night/day growth ratio": None}
    with rasterio.open(sorted(toa)[-1]) as src:
        arr = src.read(1, masked=True).astype(float).filled(np.nan)
    valid = arr[np.isfinite(arr)]
    if valid.size == 0:
        return {"night/day growth ratio": None}
    # Hourly bins of newly arrived cells.
    tmax = np.nanmax(valid)
    hours = int(np.ceil(tmax / 3600.0))
    growth = np.array([
        np.sum((valid >= h * 3600.0) & (valid < (h + 1) * 3600.0))
        for h in range(hours)
    ], dtype=float)
    # Simulation starts at 10:00 (guide). Night window 19:00-11:00 -> sim hours
    # [9, 25) modulo 24 from start; approximate with hours-from-start.
    sim_hour_of_day = (10 + np.arange(hours)) % 24
    night = (sim_hour_of_day >= 19) | (sim_hour_of_day < 11)
    day_growth = growth[~night]
    night_growth = growth[night]
    day_mean = np.mean(day_growth[day_growth > 0]) if np.any(day_growth > 0) else np.nan
    night_mean = np.mean(night_growth[night_growth > 0]) if np.any(night_growth > 0) else 0.0
    ratio = (night_mean / day_mean) if day_mean and np.isfinite(day_mean) else None
    return {"night/day growth ratio": ratio}


def metric_spotting_stats(out_dir, case):
    import glob
    csvs = glob.glob(os.path.join(out_dir, "spotting_stats_*.csv"))
    if not csvs:
        return {"max spotting distance": None, "_ember_count": 0}
    frames = []
    for c in csvs:
        try:
            df = pd.read_csv(c)
            df.columns = df.columns.str.strip()
            frames.append(df)
        except Exception:
            pass
    if not frames:
        return {"max spotting distance": None, "_ember_count": 0}
    df = pd.concat(frames, ignore_index=True)
    dist = pd.to_numeric(df.get("DIST"), errors="coerce").dropna()
    if dist.empty:
        return {"max spotting distance": None, "_ember_count": int(len(df))}
    # Top-1% distance is the published comparison point.
    return {
        "max spotting distance": float(dist.quantile(0.99)),
        "_ember_count": int(len(df)),
        "_mean_distance": float(dist.mean()),
    }


def metric_initial_attack(out_dir, case):
    """Fraction of ensemble members contained at initial attack.

    A contained member burns far less area than an uncontained one; we read
    the fire-size-stats CSV and classify members below 25% of the max final
    size as "contained".
    """
    import glob
    csvs = glob.glob(os.path.join(out_dir, "fire_size_stats*.csv"))
    if not csvs:
        return {"fraction contained": None, "_n_members": 0}
    df = pd.read_csv(sorted(csvs)[0])
    df.columns = df.columns.str.strip()
    size_col = next((c for c in df.columns
                     if "acres" in c.lower() or "area" in c.lower()), None)
    if size_col is None or df.empty:
        return {"fraction contained": None, "_n_members": int(len(df))}
    sizes = pd.to_numeric(df[size_col], errors="coerce").dropna()
    if sizes.empty:
        return {"fraction contained": None, "_n_members": 0}
    threshold = 0.25 * sizes.max()
    contained = float((sizes <= threshold).mean())
    return {"fraction contained": contained, "_n_members": int(len(sizes))}


METRICS = {
    "max_ros": metric_max_ros,
    "max_ros_region": metric_max_ros_region,
    "max_ros_crown": metric_max_ros_crown,
    "qualitative_ros": metric_qualitative_ros,
    "consistency_ros": metric_consistency_ros,
    "stop_time": metric_stop_time,
    "overnight": metric_overnight,
    "spotting_stats": metric_spotting_stats,
    "initial_attack": metric_initial_attack,
}


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------
def pct_error(actual, target):
    if actual is None or target is None or target == 0:
        return None
    return 100.0 * (actual - target) / target


def evaluate_case(case, measured):
    """Return a list of result rows comparing measured vs target."""
    rows = []
    if not case.quantitative or not case.targets:
        # Report measured values, no pass/fail.
        for k, v in measured.items():
            if k.startswith("_"):
                continue
            rows.append(dict(case=case.name, metric=k, target=None,
                             measured=v, pct_error=None, match=None,
                             units="", notes=case.notes))
        return rows

    for t in case.targets:
        actual = measured.get(t.name)
        err = pct_error(actual, t.value)
        if actual is None:
            match = None
        elif t.abs_tol is not None:
            match = abs(actual - t.value) <= t.abs_tol
        else:
            match = (err is not None) and (abs(err) <= t.tol_pct)
        rows.append(dict(case=case.name, metric=t.name, target=t.value,
                         measured=actual, pct_error=err, match=match,
                         units=t.units, notes=case.notes))
    return rows


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def fmt(v, nd=3):
    if v is None:
        return "-"
    if isinstance(v, bool):
        return "PASS" if v else "FAIL"
    try:
        return f"{v:.{nd}f}"
    except (TypeError, ValueError):
        return str(v)


def write_reports(rows, run_dir):
    df = pd.DataFrame(rows)
    csv_path = os.path.join(run_dir, "verification_summary.csv")
    df.to_csv(csv_path, index=False)

    md_path = os.path.join(run_dir, "verification_report.md")
    with open(md_path, "w") as fh:
        fh.write("# ELMFIRE GUIDE verification report\n\n")
        n_eval = df["match"].notna().sum()
        n_pass = (df["match"] == True).sum()  # noqa: E712
        fh.write(f"- Cases run: {df['case'].nunique()}\n")
        fh.write(f"- Quantitative checks: {n_eval} ({n_pass} passed, "
                 f"{n_eval - n_pass} failed)\n\n")
        fh.write("| Case | Metric | Target | ELMFIRE | % err | Result | Units |\n")
        fh.write("|------|--------|-------:|--------:|------:|:------:|-------|\n")
        for r in rows:
            res = "n/a" if r["match"] is None else ("PASS" if r["match"] else "FAIL")
            fh.write(f"| {r['case']} | {r['metric']} | {fmt(r['target'])} | "
                     f"{fmt(r['measured'])} | {fmt(r['pct_error'], 1)} | {res} | "
                     f"{r['units']} |\n")
        fh.write("\n## Notes per case\n\n")
        for case in CASES:
            note = next((r["notes"] for r in rows if r["case"] == case.name), "")
            if note:
                fh.write(f"- **{case.name}**: {note}\n")
    return csv_path, md_path


def print_console(rows):
    hdr = f"{'Case':<20}{'Metric':<28}{'Target':>10}{'ELMFIRE':>10}{'%err':>8}  Result"
    print("\n" + hdr)
    print("-" * len(hdr))
    for r in rows:
        res = "n/a" if r["match"] is None else ("PASS" if r["match"] else "FAIL")
        print(f"{r['case']:<20}{r['metric']:<28}{fmt(r['target']):>10}"
              f"{fmt(r['measured']):>10}{fmt(r['pct_error'], 1):>8}  {res}")


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Run + verify the GUIDE template cases.")
    ap.add_argument("--run-dir", default=os.path.join(HERE, "verification_run"),
                    help="Working directory (TEMPLATE is copied here).")
    ap.add_argument("--cases", nargs="*", help="Subset of case names to run.")
    ap.add_argument("--elmfire-bin", default=os.environ.get("ELMFIRE_BIN", "elmfire"))
    ap.add_argument("--gdal-path", default=os.environ.get("ELMFIRE_GDAL"),
                    help="Force PATH_TO_GDAL in the run copies (e.g. the conda env's "
                         "bin/). Default: auto-detect a GDAL that works in this env.")
    ap.add_argument("--skip-run", action="store_true",
                    help="Do not run ELMFIRE; only re-compare existing outputs.")
    ap.add_argument("--no-generate", action="store_true",
                    help="Do not regenerate input rasters.")
    ap.add_argument("--ia-ensemble", type=int, default=None,
                    help="Override NUM_ENSEMBLE_MEMBERS for the initial-attack case.")
    args = ap.parse_args()

    selected = [c for c in CASES if (not args.cases or c.name in args.cases)]
    if not selected:
        sys.exit(f"No matching cases. Available: {[c.name for c in CASES]}")

    print(f"ELMFIRE GUIDE verification: {len(selected)} case(s) -> "
          f"{', '.join(c.name for c in selected)}", flush=True)
    if not shutil.which(args.elmfire_bin) and not os.path.isfile(args.elmfire_bin):
        print(f"[WARN] elmfire binary '{args.elmfire_bin}' not found on PATH; "
              f"runs will fail. Set ELMFIRE_BIN or use --skip-run.", flush=True)

    # Resolve a GDAL path to write into each run copy's PATH_TO_GDAL: an explicit
    # --gdal-path/ELMFIRE_GDAL wins, otherwise auto-detect one that works here.
    gdal_override = None
    if not args.skip_run:
        gdal_override = args.gdal_path or find_working_gdal()
        if gdal_override:
            print(f"PATH_TO_GDAL -> {gdal_override} (run copies only; templates unchanged)",
                  flush=True)
        else:
            print("[WARN] no working GDAL found; using each template's PATH_TO_GDAL "
                  "as-is (georeferenced .tif outputs may fail).", flush=True)

    if args.skip_run:
        run_dir = args.run_dir
    else:
        run_dir = args.run_dir
        os.makedirs(run_dir, exist_ok=True)

    all_rows = []
    consistency_ros = {}   # for BANDTHICKNESS invariance check

    for case in selected:
        src_data = os.path.join(TEMPLATE_DIR, case.data_rel)
        if not os.path.isfile(src_data):
            print(f"[WARN] {case.name}: data file missing ({case.data_rel}), skipping.")
            continue
        # Each case runs inside run_dir/<case.data_rel parent>.
        case_dir = os.path.join(run_dir, os.path.dirname(case.data_rel))
        os.makedirs(case_dir, exist_ok=True)
        dst_data = os.path.join(case_dir, os.path.basename(case.data_rel))

        if not args.skip_run:
            # Use the template verbatim; only apply opt-in CLI overrides.
            shutil.copy(src_data, dst_data)
            overrides = {}
            if gdal_override:
                overrides["PATH_TO_GDAL"] = f"'{gdal_override}'"
            if args.ia_ensemble is not None and "initial" in case.name.lower():
                overrides["NUM_ENSEMBLE_MEMBERS"] = str(args.ia_ensemble)
            override_namelist_keys(dst_data, overrides)
            if not args.no_generate:
                print(f"[{case.name}] generating inputs ...")
                generate_case_inputs(case_dir, dst_data, case)
            print(f"[{case.name}] running ELMFIRE ...")
            ok, log, rc = run_elmfire(case_dir, dst_data, args.elmfire_bin)
            if not ok:
                print(f"[{case.name}] RUN FAILED (rc={rc}); see {log}")

        out_dir = os.path.join(case_dir, "outputs")
        if not os.path.isdir(out_dir):
            print(f"[WARN] {case.name}: no outputs dir, skipping comparison.")
            continue

        extractor = METRICS[case.metric]
        try:
            measured = extractor(out_dir, case)
        except Exception as exc:  # keep going across cases
            print(f"[{case.name}] metric extraction error: {exc}")
            measured = {}

        if case.metric == "consistency_ros":
            consistency_ros[case.name] = measured.get("max ROS")
        all_rows.extend(evaluate_case(case, measured))

    # Band-thickness invariance: max-ROS should agree across 1/2/3.
    bt = {k: v for k, v in consistency_ros.items() if v is not None}
    if len(bt) >= 2:
        ref = list(bt.values())[0]
        spread = max(bt.values()) - min(bt.values())
        rel = 100.0 * spread / ref if ref else None
        match = (rel is not None) and (rel <= 5.0)
        all_rows.append(dict(case="BANDTHICKNESS", metric="ROS invariance (max-min)",
                             target=0.0, measured=spread, pct_error=rel,
                             match=match, units="m/min",
                             notes="Max ROS should be invariant to band thickness."))

    if not all_rows:
        sys.exit("No results produced.")

    print_console(all_rows)
    csv_path, md_path = write_reports(all_rows, run_dir)
    print(f"\nWrote: {csv_path}\nWrote: {md_path}")


if __name__ == "__main__":
    main()
