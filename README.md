# vesselTracker

Minimal, standalone example of **sub-voxel vessel tracking in fMRI timeseries** with
`fitVesselPatchTimeSeries.m`: a 2D Gaussian (+ flat background) is fit to every frame of a small patch
around a vessel, with the granularity of every parameter chosen independently -- one value per vessel,
per run, per frame, or a low-order polynomial in time per run.

The code is a scoped-down port of the `vesselFit/` folder of the (private) `huMoMain2` project; this
repository exists to share the fitting core with a working example and to start packaging it as a
tool of its own.

## Requirements

MATLAB (developed on R2025a) with the Image Processing Toolbox (`niftiread`) and the Optimization
Toolbox (`lsqnonlin`). No other dependency: everything runs from this repository's own folders.
Optional: `ffmpeg` on the PATH for `.mp4` QA movies (the example writes `.gif`, which needs nothing).

## Running the example

Open `vesselTracker.code-workspace` (or just the folder) and run `doIt.m`, top to bottom (F5) or
section by section (Ctrl+Enter in the MATLAB editor; the sections fold).

```
doIt.m
  Set up environment            addpath of the two code folders
  Load data                     loadNiftiTs   -- one 4D NIfTI per run -> tsIm.im{r} [ny x nx x 1 x T]
  Select vessel centers         drawVesselCenter -- click vessel centers on the time-averaged image
  Build the vessel variable     makeVessel    -- one 11x11 patch per center -> vessel{v}.tsIm
  Fit -- time-averaged anatomy  fitVesselPatchTimeSeries, opts.timeAvg=true          -> .gaussAnat
  Fit -- perVessel/perRun/perFrame/perRunPoly   seeded from .gaussAnat  -> .gaussPerVessel/.gaussPerRun/
                                                                           .gaussPerFrame/.gaussPerRunPoly
  Visualize the fit             showVessel (patch + fitted contour), fitted-position timecourses,
                                fitVesselTimeSeriesDiag (per-frame QA movie of the perFrame fit)
```

The example data is one run of a preprocessed single-slice vfMRI timeseries
(`data/example_preproc_volTs.nii.gz`: 400 x 400 pixels at 0.4 mm, 350 frames, TR 0.84 s). The vessel
centers clicked for it are cached in `drawVesselCenterCache/example.mat`, so the example runs without
any interaction; set `forceThis = 1` at the top of `doIt.m` to re-open the click figure (seeded from the
cache), or `2` to start blank. Outputs land in `figures/` and `qaMovies/` (both gitignored).

## Layout

- `doIt.m` -- the example script (entry point).
- `vesselPatch/` -- the three functions written for this tool: `loadNiftiTs.m`, `drawVesselCenter.m`,
  `makeVessel.m`. Each supports a no-argument call that prints its options and returns the defaults.
- `vesselFit/` -- `fitVesselPatchTimeSeries.m`, its dependencies, and two visualizers
  (`showVessel.m`, `fitVesselTimeSeriesDiag.m`), copied **unmodified** from `huMoMain2/vesselFit/`
  (commit `a443faa`, 2026-09-24), plus `loadNiftiVol.m` from `huMoMain2/humanVessel/`. Only the files
  the four fit granularities and the two visualizers actually execute were brought over; see below.
- `data/` -- the example timeseries.
- `drawVesselCenterCache/` -- the cached manual clicks (git-tracked on purpose: a click is not
  regenerable).

## What was deliberately left out of the port

`fitVesselPatchTimeSeries.m` is copied whole, so its header documents features whose helper files are
**not** here. Calling them fails with an "Undefined function" error rather than silently doing
something else:

- bootstrap / null replicates (`opts.boot`, `opts.null`, the replicate data cache),
- the IRF-image path (`vessel.(fld).irfIm`, `opts.addBase0`, `mode='perKnot'`),
- motion-seeded fits on raw data (`opts.seedMotion`),
- `showVessel.m`'s on-the-fly ROI derivation and stat-map overlays (hand-drawn ROIs, `drawVessel.m`),
  and the `colormap_divergingHue` colormap of `util` (the visualizers fall back to their own maps).

Their code lives in `huMoMain2` and can be brought over the same way if needed.

## Conventions worth knowing

- Images are stored `[rows x cols x 1 x frames]` with rows = NIfTI dim 2 and cols = NIfTI dim 1
  (FreeSurfer `MRIread` orientation, see `loadNiftiVol.m`). A pixel position is `(x, y) = (col, row)`.
- Fitted positions `x0`/`y0` are in mm relative to the patch center pixel.
- A fitted parameter that differs between runs is a `{1 x nRun}` cell, a parameter that is the same in
  every run is stored once, bare; read run `r` with `perRunVal(x, r)`.
- Every user-facing function called with no arguments prints its options and returns the defaults
  (`opts = fitVesselPatchTimeSeries;`).
