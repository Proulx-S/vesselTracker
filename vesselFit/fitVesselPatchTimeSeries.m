function vessel = fitVesselPatchTimeSeries(vessel, fld, opts)
% FITVESSELPATCHTIMESERIES  Jointly fit EVERY peak in a vessel patch, plus one shared background,
% frame by frame -- with per-parameter granularity (perVessel/perFrame/perRun/perKnot) on every
% parameter of every peak independently.
%
% Where fitVesselTimeSeries.m fits ONE Gaussian per frame and has to hide the other peaks behind a
% tight includeMask, this fits them all at once. The practical consequence (see doIt_human.m's own
% "Estimate individual vessel motion" comment, which today restricts the per-frame fit to a narrow ROI
% "to limit the impact of spurious peaks, at the expense of non-validity of other parameters"):
% modelling the secondary peaks explicitly removes the reason for that restriction.
%
% OPTS.MODEL IS THE MODEL EQUATION. Every field of opts.model names a registered component model (see
% componentModelSpec.m); every ELEMENT of that field is one additive component. The prediction is the
% sum of all of them:
%
%   opts = fitVesselPatchTimeSeries();
%   opts.outFld = 'gaussPatchPerFrame';
%   sf = 'tsIm.patchVesselsSimultaneous';
%   opts.model.gaussian(1).seed    = [sf '(main)'];   % the main peak
%   opts.model.gaussian(1).mode.x0 = 'perFrame';      % ... tracked per frame
%   opts.model.gaussian(1).mode.y0 = 'perFrame';
%   opts.model.gaussian(2).seed    = sf;              % every REMAINING peak, one component each
%   opts.model.background.seed     = sf;              % the one shared background
%   vessel = fitVesselPatchTimeSeries(vessel, 'tsIm', opts);
%
% Adding a future model type is a new field (opts.model.parab(k)) plus a componentModelSpec.m entry --
% no change here and none in the solver.
%
% SEED PATHS CARRY THEIR OWN PEAK SELECTOR. A component entry's .seed is a dot-path from the VESSEL
% ROOT to a prior fit -- a fitVesselPatchTimeSeries.m result (read from its .model.gaussian, at
% whatever granularity it was fit, see SEED GRANULARITY below) or a fitPatchVessels.m detector result
% (its .fits) -- optionally suffixed with a selector:
%   'tsIm.gaussAnat(main)'        -- the main peak alone; '(1)' is the same thing
%   'tsIm.gaussAnat(secondary)'   -- every peak BUT the main one, by ROLE; empty (and not an error) on a
%                                    single-peak vessel
%   'tsIm.gaussAnat'              -- every peak NOT already claimed by an earlier entry
%   'tsIm.gaussAnat(2 3)'         -- those peaks, by ARRAY POSITION
%   'tsIm.gaussAnat.main'         -- accepted, but the paren form is canonical (see below)
% '(main)' + '(secondary)' is the pair the 2-arg call emits (see 2-ARG CALL): one entry per ROLE, so ONE
% opts struct works across a batch whose vessels have different peak counts, without the caller ever
% writing K, and without depending on which entry is processed first. The unindexed form still works
% by CLAIMING: entries are processed in order and each takes the peaks earlier ones did not, so
% '(main)' then unindexed means the same thing as '(main)' then '(secondary)' -- but only because entry 1
% claimed peak 1 first (see 2-ARG CALL's THE BUG THIS SIDESTEPS).
%
%   - An UNINDEXED or '(secondary)' entry matching zero peaks contributes NOTHING and is not an error:
%     that is exactly the single-peak vessel under the two-entry template.
%   - An EXPLICIT index past the end IS an error (fitVesselPatchTimeSeries:seedPeakOutOfRange) -- asking
%     for peak 4 of a 3-peak vessel is a real mistake, unlike asking for "the rest" and finding none.
%   - Indices are ARRAY POSITION, not .peakNum. Array position 1 == the main peak is an invariant read
%     unconditionally across this codebase (alignVessel2.m, makeRoiFromVesselFit.m, plotVessels.m), and
%     fitPatchVessels.m already compacts positions when it drops excluded peaks, so positions are dense
%     by the time a result is seeded from. Each peak's own .peakNum still travels into the output.
%   - The paren form is canonical because getNestedField.m walks dot segments blindly, so a trailing
%     '.main' is ambiguous with a genuine field called main; '(main)' never is. The dot form is resolved
%     only after a real field lookup fails.
%   - Peaks flagged .excluded=true (fitPatchVesselsDiag.m's own curation) are DROPPED before any of this
%     -- same contract fitPatchVessels.m's own resolveSeedFits honours -- so selectors and "remaining"
%     both count INCLUDED peaks only.
%
% 'path.heuristic' / 'path(heuristic)' -- a THIRD selector spelling, for when there is no prior
% fitPatchVessels.m result at all yet (a caller's very first "get me a starting anatomical fit" call):
%   'tsIm.heuristic'      -- dot form
%   'tsIm(heuristic)'     -- paren form, mirroring '(main)' exactly; the two are equally accepted and
%                            trigger the IDENTICAL pathway
% Bypasses fitPatchVessels.m entirely: 'tsIm' must instead name a real ts-track (a struct with its own
% .im/.vSize, e.g. vessel.tsIm), and ONE ad hoc seed peak is computed directly from that path's own
% time-averaged image via the legacy single-peak start/bound heuristic (gaussianFitDefaults.m,
% generalized to this file's radius/aspectRatio parameterization by gaussianFitDefaultsRadiusAspect.m),
% seeded at the patch center. Always yields exactly one synthetic peak -- there is no compound/indexed
% heuristic selector (no 'tsIm.heuristic(2)'). RESPECTS opts.includeMask/.excludeMask (2026-09-09, a
% real bug found investigating a reported off-FOV convergence): excluded/non-included pixels are
% replaced with the valid region's own median before any heuristic statistic (amplitude, background) is
% computed from the image, so a bright excluded feature (e.g. a nearby sinus) can no longer bias the
% starting amplitude/background estimate the way it could before this fix -- the SAME mask the real fit
% that follows will use, resolved identically. See resolveHeuristicSeedPeak's own note.
%
% BARE 'heuristic' (2026-09-08, Seb's own ask, and this is the DEFAULT .seed -- see defaultOpts's own
% emptyEntryFor) -- a seed of exactly 'heuristic', with no path at all, expands to THIS CALL's own
% '<fld>.heuristic' (the same fld passed to fitVesselPatchTimeSeries itself), so a caller does not have
% to repeat fld in every entry's seed just to ask for "the heuristic on my own data." A caller fitting
% a DIFFERENT track's heuristic still spells the full path out (e.g. 'tsImMotionCorrected.heuristic').
%
% START VALUES AND BOUNDS ARE LITERAL NUMBERS. There is no .source dot-path mechanism here. Each
% parameter's start comes from its component's own seed peak; each bound comes from .seedBounds applied
% to that same seed (additive for x0/y0/theta, multiplicative for a/radius/aspectRatio/b -- the exact
% vocabulary and formulas fitPatchVessels.m's 'simultaneous' mode already uses, via the shared
% additiveSeedBound.m/multiplicativeSeedBound.m). Either can be overridden with a literal number via
% .start/.lower/.upper. Two things follow that were not previously possible:
%   - .fixed=true works with NO seed at all, by pairing it with a literal .start.
%   - .start/.lower/.upper mean the same thing on input as on output, so a paramSpec ROUND-TRIPS: feed a
%     previous fit's own reported values straight back in to pin a new fit at them.
%
% SEED GRANULARITY (2026-09-10, Seb's own rule). A seed peak read from a fitVesselPatchTimeSeries.m
% result carries each parameter at that fit's own granularity (per vessel, per run, or per frame), and
% the parameter's requested .mode HERE decides what becomes of it (collapseSeedToMode):
%   - target COARSER than the seed -> AVERAGE: a per-frame seed into a 'perRun' fit starts each run at
%     that run's mean; anything into a 'perVessel' (or 'perRunPoly'/irf-path) fit starts at one mean.
%   - target AS FINE OR FINER -> keep every value the seed has: a per-run seed into a 'perRun' fit
%     starts each run at its own value, into a 'perFrame' fit each run's value is repeated over that
%     run's frames; a per-frame seed into a 'perFrame' fit is used frame for frame.
% A seed that differs between runs can only feed a fit of the SAME runs (run/frame counts must match;
% fitVesselPatchTimeSeries:seedRunCountMismatch / :seedFrameCountMismatch otherwise). .seedBounds is
% applied elementwise to whatever shape the start ends up with. So a chain gaussAnat (time-averaged)
% -> gaussAnatPerRun -> a per-frame fit loses nothing at any step.
%
% GRANULARITY DEFAULTS TO 'perVessel', NOT 'perFrame'. This is a deliberate difference from
% fitVessel.m/fitComponentFit.m, whose default is 'perFrame'. The seed here is always an
% already-converged time-averaged fit, so the natural default action is "refine that jointly", with
% per-frame tracking opted INTO for the specific parameters of interest. The other way round, a
% default-template call on a 4-peak, 6-run, 350-frame vessel would silently launch a fit with tens of
% thousands of free parameters. Cheap and conservative by default; expensive only on request.
%
% RUN HANDLING. Every run is concatenated into one series and fit JOINTLY (fitVessel.m's own runCat
% semantics, always on -- there is no opts.runCat here). 'perVessel' therefore means one value for the
% whole vessel across every run, and 'perRun' means one value per original run. Output is then SPLIT
% BACK per run by the same three-way length convention fitVesselTimeSeries.m's own PER-RUN OUTPUT
% established (1 = perVessel/fixed scalar, nRun = perRun, T_total = perFrame) -- and then IDENTICAL
% VALUES ARE STORED ONCE (2026-09-10, Seb's own ask -- see collapseIdentical): a per-parameter field is
% a {1 x nRun} cell only when its value actually DIFFERS between runs, and a bare value when it is the
% same in every run (every 'perVessel'/'fixed' parameter, and a per-frame row that is constant within a
% run is likewise one number). A cell therefore always means "varies by run". Readers index runs with
% perRunVal.m (bare for any run, x{r} for a cell) rather than assuming a cell.
%
% opts.timeAvg (default FALSE) collapses time BEFORE fitting -- how much of it, depends on the FINEST
% mode requested anywhere in opts.model (2026-09-08, Seb's own design), not a blanket full collapse:
%   - every parameter 'perVessel' -> every frame of every run pooled into ONE time-averaged frame --
%     the same reduction fitMultiVessel.m performs, so a timeAvg call with everything 'perVessel' is
%     the time-averaged multi-peak refit, expressed in this framework.
%   - at least one parameter 'perRun' (and none finer) -> average WITHIN each run only, keeping runs
%     distinct as nRun pseudo-frames, so a 'perRun' parameter still gets one value per run instead of
%     being silently pooled into a single global value it never asked for.
%   - 'perFrame'/'perKnot'/'perRunPoly' anywhere in opts.model -> REJECTED outright
%     (fitVesselPatchTimeSeries:timeAvgFinerThanPerRun): once ANY time-averaging happens there is no
%     frame/knot axis left to vary over, and silently collapsing them too would invent values nobody
%     asked for.
% Note the default (false) is the OPPOSITE of fitVessel.m's (true), because a genuine per-frame fit is
% this function's whole point.
%
% .fits COMPATIBILITY VIEW (see assembleResult/attachFitsView) is attached whenever every parameter
% ends up scalar -- timeAvg=true (any mode), OR timeAvg=false with every parameter 'perVessel' (a
% joint fit across raw per-frame data, no pre-averaging, still exactly one value per parameter).
% timeAvg=false with any coarser-than-scalar mode (e.g. 'perRun') has no .fits at all.
%
% WHICH TIME SERIES -- TWO kinds, selected STRUCTURALLY by what sits at vessel.(fld), not by fld's own
% spelling (the same dispatch fitVesselTimeSeries.m uses):
%   - vessel.(fld) carries an .irfIm.im (a fitIRF.m per-track IRF container) -> fits that image, whose
%     nKnot delay-slices stand in for time frames. e.g. fld='tsIm.irf'. The result stores one level
%     deeper, at vessel.(fld).irfIm.(outFld). .vSize/.rois are borrowed from fld's own PARENT path --
%     the real ts track the IRF was fit from, which shares its exact crop and voxel grid.
%   - otherwise vessel.(fld) carries .im directly -> fits its real acquired frames.
% On the irf-path a "frame" is a knot DELAY, so per-run concepts do not apply: the output is NOT split
% into per-run cells (a knot delay does not belong to a scan run), and opts.timeAvg is refused outright
% (averaging across knot delays destroys the response shape that image exists to represent).
%
% opts.addBase0 (default true, IRF-PATH ONLY) -- fitIRF.m bakes the degree-0 (DC) baseline into
% irfIm.im at fit time, so true is a no-op: the fit sees the absolute signal level, and .a/.b mean what
% you would expect. false SUBTRACTS it back out (irfAddBase0.m, sign=-1) to fit the pure
% stimulus-evoked deviation instead. Not a skip -- it has to UNDO fitIRF.m's own addition -- and it
% errors if there was never a base0 to remove, which is correct: asking to subtract a baseline that was
% never fit is a real caller mistake.
%
% mode='perKnot' (TS-PATH ONLY) -- one value per knot-group, with each real frame assigned to the knot
% its IRF design gives it. fitComponentFit.m already implements the granularity; the only thing this
% level adds is the frame->knot map, read from vessel.(fld).irf.irfMat.X (the solver has no concept of
% a knot). Requires a one-hot design. KNOT COVERAGE: a frame the design does not cover at all (run-edge
% dead time before the first / after the last modelled knot) is EXCLUDED from the fit rather than
% forced into a group it is not in, then restored as NaN in the output so every reported per-frame
% array still spans the full raw frame count in the original order. Pair with .fixed=true and a
% literal nKnot-long .start to PIN a parameter per knot instead of estimating it.
%
% mode='perRunPoly([<degrees>])' (2026-09-07) -- an INDEPENDENT low-order Legendre polynomial in
% within-run frame time, with its own coefficients PER ORIGINAL RUN, filling the gap between 'perRun'
% (one constant per run) and 'perFrame' (fully free per frame) -- e.g. a vessel position that drifts
% smoothly within a run, with its own drift trajectory per run. fitComponentFit.m owns the entire
% mechanism (parsing, basis construction, seeding, packing, output expansion) -- see that file's own
% mode='perRunPoly' doc for the full spec; only the parts relevant AT THIS LEVEL are repeated here:
%   - The degree list is embedded in the mode string itself, e.g. .mode.x0 = 'perRunPoly([0 1])'
%     (intercept + linear drift); 'perRunPoly(1,3)', 'perRunPoly[1 2]', 'perRunPoly(0:4)' are all
%     accepted spellings of the same comma/space-separated, optionally-ranged degree list.
%   - SEEDING follows this file's own ordinary seeding exactly: the SAME scalar seed value/bound every
%     other mode would get (from the seed peak, .seedBounds, or a literal .start/.lower/.upper override)
%     is what fitComponentFit.m broadcasts onto the degree-0 coefficient of every run (every other
%     coefficient starts at 0, symmetric default bound) -- so switching a parameter's mode to
%     'perRunPoly([...])' needs NO other opts.model.<type>(k) field touched.
%   - ANGLE WRAPPING: a periodic parameter (currently just 'gaussian'.theta, period pi -- see
%     componentModelSpec.m's own .wrapParamPeriod) has its fitted per-frame trajectory wrapped into its
%     principal range automatically -- no caller action needed.
%   - Requires OPTS.TIMEAVG=FALSE (rejected exactly like 'perFrame'/'perKnot' -- see opts.timeAvg above)
%     and, like every other mode, does NOT support .fixed=true.
%   - KNOWN LIMITATION, not fixed here: fitVesselTimeSeriesDiag.m's own eligibility check for its
%     animated-panel figure tests literally for a mode='perFrame' parameter
%     (~any(strcmp(struct2cell(r.paramSpec.mode),'perFrame'))), so a fit whose only per-frame variation
%     comes from a perRunPoly parameter is currently reported as ineligible ("nothing to animate") even
%     though it genuinely varies per frame -- a real gap, deliberately left for a separate change rather
%     than reached for here.
%
% A PARALLEL POOL IS STARTED ONLY WHEN REPLICATES ARE REQUESTED. The point-estimate fit has no parfor
% anywhere -- the per-vessel loop is a plain for and the solver is one lsqnonlin call -- so a plain
% call never touches ensureParPool.m at all. An earlier version called it unconditionally, copied from
% fitVesselTimeSeries.m, where it did nothing here but spin up a 32-worker pool no code used: ~20 s of
% startup per call, and a hard failure of the WHOLE call when the pool could not open a port
% (observed, which is how it was noticed). The replicate loops below are the one genuine parfor, so
% the pool call now sits with them, gated on actually having replicates to run.
%
% BOOTSTRAP AND NULL REPLICATES (opts.bootRep / opts.nullRep). Two distinct STATISTICAL PURPOSES,
% named for the purpose rather than the method:
%   .bootRep -- resampling the REAL signal, so the spread across replicates is this fit's own SAMPLING
%               UNCERTAINTY. Sourced from fitIRF.m's own trial-cluster bootstrap (irfIm.imBoot).
%               IRF-PATH ONLY: there is no raw-data bootstrap here (resampling frames would change each
%               replicate's own frame count, leaving a 'perFrame' parameter with no fixed-shape output
%               to collect -- a genuinely harder problem, deliberately not attempted).
%   .nullRep -- SURROGATE data carrying no real temporal structure, so the spread describes what the
%               fitted parameters look like under the NULL HYPOTHESIS. Compare the point estimate
%               against it for an empirical p-value. Two sources, same quantity:
%                 irf-path: fitIRF.m's own pre-built surrogates (irfIm.imPerm). opts.nullRep=true.
%                 ts-path:  generated HERE from the raw data. opts.nullRep=N, the surrogate count.
%
% ON THE NAMING (2026-09-04, Seb's own call, and it fixes a real misnomer). This vocabulary replaces
% fitVesselTimeSeries.m's own .boot/.perm pair. "Permutation" there does not mean a permutation test:
% permutationNullIRF.m calls phaseScrambleRun.m, i.e. FFT phase randomization (Theiler et al. 1992),
% and there is no randperm anywhere in vesselFit/ -- so "perm" and "phase-randomize" were two names for
% ONE mechanism, which is what made them read as two concepts. boot/null names the PURPOSE and leaves
% the method an implementation detail, so a genuine permutation-based null could be added later as
% another source of .nullRep with nothing renamed. It also turns a caveat into a fact: the irf-path and
% ts-path nulls write the same fields because they ARE the same quantity, not because they collide.
% fitIRF.m's own imBoot/imPerm keep their names -- they are inputs read here, and renaming them is a
% separate change touching a live function with its own consumers. permutationNullIRF.m is likewise
% still misnamed; a separable cleanup, noted rather than done here.
%
% SHAPE -- .bootRep.<param> / .nullRep.<param> are [nRep x pLen] (raw and UNSUMMARISED, which is the
% whole contract -- see REPLICATES ONLY below -- so any std, custom interval or empirical p-value can
% be computed from them downstream). They live per component (opts.model.<type>(k).bootRep) AND on the
% flat main-peak alias. Deliberately NOT split per run: the leading replicate dimension is not a run
% axis. There are no .bootSE/.nullSE fields; replicateSE.m is the one-liner for that.
%
% SCOPE. Deliberately NOT ported from fitVesselTimeSeries.m, because nothing in this project exercises
% any of it and each carries its own substantial contract (this project's convention is scoped-down
% ports -- port what a real call path uses, not the whole original):
%   - the RUN-EDGE BASELINE and 0BASE sub-fits (opts.outFld_runEdgeBase/opts.outFld_0base with their
%     own seeds and paramSpecs) -- two additional, separately-seeded fits with their own output shapes
%     and .colLabels design-matching, not a variation on this one;
%   - SHAPE-B seeding (seeding from a prior per-knot fit and mapping it onto real frames);
%   - opts.removeBaselineDrift;
% fitVesselTimeSeries.m remains the function to use for any of those.
%
% SEED MOTION (opts.seedMotion, ts-path only) -- lets a gaussian be fit on RAW, un-motion-corrected
% patches while still being seeded from a fit made on the PREPROCESSED ones. A char dot-path or a
% cellstr of them, each naming a getPreprocMotion.m-shaped struct; several are SUMMED
% (resolveMotionSeed.m, shared with the predecessor so both resolve identically).
%
%   opts.seedMotion = 'tsIm.preprocMotion';
%   opts.model.gaussian(1).seed    = 'tsIm.patchVesselsSimultaneous(main)';
%   opts.model.gaussian(1).mode.x0 = 'perFrame';   % REQUIRED for any position it corrects
%   opts.model.gaussian(1).mode.y0 = 'perFrame';
%
% TOP-LEVEL, NOT PER COMPONENT (Seb's call, 2026-09-05): motion is a property of the TRACK, not of any
% one peak, so one setting applies to every component and opts.model stays purely the model equation.
% The alternative -- a per-component .seedMotion -- was available and declined; nothing needs per-peak
% motion, and it would invite two peaks of one patch being corrected inconsistently.
%
% SIGN: getPreprocMotion.m reports the correction preprocessing APPLIED, so preprocessed-position
% MINUS that correction is where the same feature sits in the raw data. Subtracting it turns a scalar
% position seed into a [1 x T] row, which is the whole point -- and the corresponding .lower/.upper
% then become per-frame too, each bounded around its own frame's corrected seed rather than around one
% global value.
%
% TWO THINGS IT REFUSES TO DO, both deliberately loud rather than silent:
%   - opts.timeAvg=true is REJECTED. A motion trace is inherently per-frame and timeAvg collapses
%     every frame before fitting, so there is no frame axis left to correct.
%   - a corrected position whose mode is NOT 'perFrame' is REJECTED. Proceeding would mean either
%     averaging the correction away or keeping only frame 1 of it, and both silently discard the very
%     thing the option was set for.
%
% REPLICATE DATA CACHE -- TWO ORTHOGONAL AXES, per kind (opts.boot.* / opts.null.*), plus one shared
% opts.cacheDir. The general convention is written up once in
% .bass/pkm/patterns/data-cache-offload-recompute.md and is meant for ANY heavy data-processing step,
% not just this one:
%   .offloadToCache (default TRUE)  -- MEMORY axis. true leaves the raw replicates on disk and puts a
%                    dataCacheDescriptor.m POINTER in the result; false materialises them in it.
%   .reComputeCache (default FALSE) -- CPU axis. false reuses whatever already exists (in the vessel,
%                    else in the cache file); true recomputes regardless, overwriting the file when
%                    offloading too.
% So: WRITE iff offloadToCache. READ iff ~reComputeCache and something exists. POINTER iff
% offloadToCache. The full table, including the already-in-vessel case:
%
%   offload  reCompute  exists          result
%   -------  ---------  --------------  ------------------------------------------------
%      1         0      yes             POINTER ONLY -- the file is never even read
%      1         0      no              compute -> write -> pointer
%      1         1      either          compute -> OVERWRITE -> pointer
%      0         0      yes             load into the vessel struct
%      0         0      no              compute -> vessel (nothing written)
%      0         1      either          compute -> vessel (nothing written)
%   already in the vessel: reused unless reComputeCache, offloaded if offloadToCache.
%
% The default pair (offload on, recompute off) is the case that should almost always be wanted: never
% recompute, never hold in RAM.
%
% REPLICATES ONLY -- NOTHING DERIVED IS COMPUTED HERE (Seb's rule, 2026-09-05, applying to every
% replicate producer including fitIRF.m). This function computes and stores the raw replicate fits and
% nothing else: no .se, no confidence interval, no empirical p-value, on EITHER offload path.
% Downstream code that needs something derived computes it, at the point of use:
%
%   se = replicateSE(fit.bootRep);              % accepts a pointer OR raw values
%   ci = prctile(dataCacheLoad(fit.nullRep).val{1}.a, [2.5 97.5]);
%
% WHY, because it is worth not relitigating: a derived value stored next to (or instead of) its parent
% is a second representation of one thing, and it is the derived one that goes stale silently. .se is
% the spread over whatever nRep held when it was written, so resuming to a larger nRep leaves it
% describing a subset with nothing marking it wrong. Keeping the producer to raw replicates makes that
% class of bug unrepresentable, and it keeps the choice of summary (std? percentile interval? one- or
% two-sided p?) with the analysis that actually knows which one it wants. fitIRF.m already reached the
% same conclusion for its own replicate mode -- opts.stats.enable=false + opts.nPerm>0 accumulates
% raw replicate images "with NO derived statistics/p-values at all, for a caller who intends to derive
% its own significance testing from the replicates at a later stage" (its header, 2026-08-24). This
% rule generalises that from an opt-in combination to the only behaviour.
%
% The only things travelling alongside the replicates are cache METADATA -- .nBootRep/.nNullRep,
% .bootSource/.nullSource, .bootKey/.nullKey -- which describe the replicate SET, not its values, and
% which the cache policy needs in order to decide reuse. Nothing derived is written to the cache file
% either, so a pointer-only hit reads key/n and nothing more.
%
% ON THE NAMING, because the predecessor's is actively misleading: fitVesselTimeSeries.m and fitIRF.m
% both use a single opts.forceCache, default TRUE, which means "force RECOMPUTATION" -- it ignores the
% cache and overwrites it. That conflated the two axes above behind one flag whose name reads as its
% own opposite, and its consequence was that caching was implemented but its read half never ran
% unless a caller explicitly passed false, which doIt_human.m never does: written every run, never
% read. Splitting the axes also reunites the two PURPOSES the predecessor had split across separate
% mechanisms -- fitIRF.m's replicateImCacheDir for memory offload, fitVesselTimeSeries.m's
% replicateTsCacheDir for recompute-avoidance -- into one policy, leaving the choice to the caller.
%
% DATA CACHE vs ANALYSIS CACHE (Seb's own distinction, worth keeping): this is a DATA cache -- a
% per-process file holding one heavy result. The ANALYSIS cache is the whole-workspace checkpoint
% checkCache.m/saveCache.m/loadCache.m manage in a doIt_*.m script. They compose in a useful
% direction: the more the data cache is used, the less actual data the analysis cache contains, which
% makes the analysis cache a true cache of the analysis PROCESS rather than a dump of its data.
%
% RESUME is by replicate INDEX, and only ever a contiguous TAIL -- structurally, not by convention:
% the reuse count is always a PREFIX count (never a list of indices), the refit is
% `parfor r = nHave+1:nRep`, and the merge concatenates base-then-new, so row order stays 1..nRep and
% the replicateIm(:,:,:,:,r) slice always sees the bare loop variable (which is what keeps parfor's
% slicing classification provable -- see the OOM note at that loop).
% Replicate r is index-stable -- slice r of imBoot/imPerm is fixed, and the ts-path null's own
% substream r is fully determined by opts.dataSeed. VERIFIED, not assumed: stream r out of
% RandStream.create(...,'NumStreams',N,...) is byte-identical for N=3 and N=6, so a resumed run's new
% tail belongs to the same stream family as the cached prefix rather than to a re-partitioned set
% (asserted in the stage-6 test suite, since it is a property of MATLAB rather than of this code).
% Asking for more replicates later therefore refits only the new ones and merges. NOTE this also
% makes the TS-PATH null cacheable, which the
% predecessor explicitly could not do: it drew fresh randomness per replicate with no substream
% bookkeeping, so "replicate r" was not a stable thing to cache. One substream per replicate fixes it.
%
% THE KEY IS THE WHOLE RESOLVED SPEC, and it guards in-vessel reuse as well as cache reuse. This
% matters far more with reuse ON by default: anything that changes the fitted values must invalidate,
% or a re-run after tightening a bound or flipping a mode silently hands back the OLD replicates. The
% key covers every component's model/label/mode/fixed/start/lower/upper, the exact mask (by linear
% index -- size and count alone would let two different equal-size masks collide),
% nIsochromatPerVoxDim, the runIdx/knotIdx granularity vectors, addBase0, dataSeed, and the replicate
% SOURCE's own identity (a file-backed stack by path+size+mtime, an in-memory one by shape+checksum).
% Stored verbatim rather than hashed, so a miss can be diagnosed by reading the file. The predecessor
% compared only the stored parameter NAME LIST -- which would not have caught any of those.
%
% NO-ARG CALL -- fitVesselPatchTimeSeries() prints every opts field with its allowed values and returns
% opts fully populated with defaults, including the DEFAULT MODEL TEMPLATE: two gaussian entries (main,
% then all remaining) plus one background, with empty .seed placeholders for the caller to fill. Same
% convention as fitVessel.m/fitPatchVessels.m.
%
% NAMED-MODEL CALL -- fitVesselPatchTimeSeries(modelName), with modelName a char/string and NO other
% argument, is the same no-real-work opts-help-and-return contract as the no-arg call above, but with
% opts.model.gaussian built from the named template instead of the default two-entry one:
%   fitVesselPatchTimeSeries('N gaussian + background')  -- opts.model.gaussian is an N-ELEMENT array
%     (N a positive integer), every entry an IDENTICAL, independent, fully-populated placeholder
%     (2026-09-08, Seb's own ask: "The first one will be the main and the other N-1 will be
%     secondary"). Array position 1 == main peak is the same invariant read unconditionally
%     throughout this codebase (SEED PATHS above); positions 2..N are each their own SEPARATE slot,
%     not one shared "everything else" entry -- assign each an explicit peak selector once there is
%     real multi-peak fitPatchVessels.m data to seed from (e.g. 'tsIm.patchVesselsSimultaneous(2)').
%     Left untouched, every entry still defaults to the bare 'heuristic' seed (see SEED PATHS'
%     default-seed note) -- since the heuristic only ever produces ONE synthetic peak, position 1
%     claims it and positions 2..N each find nothing left, contributing no component (the same
%     "unindexed entry matching zero remaining peaks is not an error" convention SEED PATHS already
%     documents), not an error.
%   fitVesselPatchTimeSeries('1 gaussian + background')  -- N=1: opts.model.gaussian is a SINGLE
%     entry (no secondary slots at all). For a caller who knows up front there is only one peak.
%   fitVesselPatchTimeSeries('main gaussian + background')  -- alias for N=1, same as
%     '1 gaussian + background' (2026-09-08, Seb's own ask): reads better standalone, or as the base
%     for the secondary-gaussian suffix below.
%   fitVesselPatchTimeSeries('gaussian + background')  -- another alias for N=1, same as
%     '1 gaussian + background' (2026-09-09, Seb's own ask -- the name he actually reaches for day to
%     day). A LITERAL alias, not a suffix base -- does not accept the secondary-gaussian suffix below.
%   fitVesselPatchTimeSeries('gaussians + background')  -- alias for N=2, same as
%     '2 gaussian + background' (2026-09-09, Seb's own ask, same day-to-day motivation as the singular
%     form above). Also a LITERAL alias, not a suffix base.
%   Either NUMBERED/MAIN base ('N gaussian + background' or 'main gaussian + background' -- NOT the two
%   literal aliases just above) may carry an optional trailing suffix that adds exactly ONE secondary
%   slot (N = N+1), i.e. populates opts.model.gaussian(2) in preparation for a fit with one (or, once
%   assigned real per-peak seeds, more) secondary peaks -- the three spellings below are
%   interchangeable, not different peak counts:
%     ' + sec gaussian', ' + secondary gaussian', ' + gaussians'
%   e.g. 'main gaussian + background + sec gaussian' and '1 gaussian + background + gaussians' both
%   produce the same 2-entry opts.model.gaussian as '2 gaussian + background' (and as the literal
%   'gaussians + background' alias above -- all three are equivalent end states, just spelled
%   differently depending on how a caller arrived there).
% Any name not matching one of the forms above errors (fitVesselPatchTimeSeries:unknownModelName)
% rather than guessing -- exact pattern match only. The no-arg call's own two-gaussian template ('')
% is unaffected by any of this (though it is now exactly what 'N gaussian + background' with N=2, and
% 'gaussians + background', also produce).
%
% 2-ARG CALL (2026-09-08, Seb's own ask, generalizing drawVesselPatchPeaks.m's own established 2-arg
% call to this function) -- fitVesselPatchTimeSeries(vessel, fld), with NO opts and no output beyond
% opts itself, is the SAME no-real-work contract as the no-arg/NAMED-MODEL calls above, but derives
% opts.model.gaussian/.background from vessel(1).(fld)'s ALREADY-EXISTING starting fit instead of a
% generic template -- at most two gaussian entries, one per ROLE: .seed='<fld>(main)' and, when that fit
% has any secondary peak, .seed='<fld>(secondary)' (see SEED PATHS), each with the corresponding
% component's own already-converged .mode/.fixed/.seedBounds copied verbatim (the result carries the
% spec it was built from, see OUTPUT), plus opts.includeMask/.excludeMask copied from that fit too.
% .start/.lower/.upper are deliberately left EMPTY (2026-09-09, Seb's own correction -- setting them at
% this "just give me a template" stage would bake in a frozen snapshot instead of a genuine USER
% override; the .seed already resolves to the SAME values at actual fit-call time, freshly, from
% whatever vessel.(fld) holds then, at whatever granularity it has -- see SEED GRANULARITY). See
% defaultOptsFromFit.
%   fld HERE MEANS SOMETHING DIFFERENT than in the 3-arg real-fit call below: a dot-path to the FIT
%   RESULT itself (e.g. 'tsIm.gaussAnat'), one level DEEPER than the 3-arg call's own fld (the raw
%   image TRACK, e.g. 'tsIm') -- the same fld-meaning divergence drawVesselPatchPeaks.m's own file
%   header already documents relative to THIS function, now also true across this one function's own
%   call forms. Typical use:
%     opts = fitVesselPatchTimeSeries(vessel, 'tsIm.gaussAnat');   % derive from an existing fit
%     opts.model.gaussian(2).mode.x0 = 'perRun';                   % override just what differs
%     opts.model.gaussian(2).mode.y0 = 'perRun';
%     vessel = fitVesselPatchTimeSeries(vessel, 'tsIm', opts);     % then actually run it (3-arg,
%                                                                   % fld back to the TRACK)
%   THE BUG THIS SIDESTEPS (2026-09-08, real incident): cloning one entry onto a second by hand
%   (`opts.model.gaussian(2) = opts.model.gaussian(1);`) copies its UNINDEXED seed along with it, so
%   both entries end up pointing at the exact same unclaimed-peak pool. Peaks are claimed in array
%   order (see SEED PATHS above), so entry 1 -- processed first -- silently claims EVERY peak the path
%   has to offer, and entry 2 then finds nothing left and contributes ZERO components: any override on
%   it (like a different .mode) never manifests anywhere, with no error at all (an unindexed entry
%   matching zero remaining peaks is ordinary, valid behaviour). This call avoids the trap
%   structurally by naming each entry's ROLE ('(main)' / '(secondary)'), neither of which depends on
%   claiming order, so two derived entries can never collide. Deriving the template this way, THEN
%   overriding only what should differ, is the recommended way to seed a multi-entry model from a
%   real prior multi-peak fit; do not give two entries of the same type an identical unindexed seed.
%
% INPUT
%   vessel : vessel struct (scalar or array), or a cell array of vessel structs.
%   fld    : dot-path (any depth) to the image track, e.g. 'tsIm'. Needs vessel.(fld).im/.vSize, and
%            .rois if any mask label is used.
%   opts   : struct --
%       .outFld  REQUIRED (default ''). Result stored at vessel.(fld).(outFld); overwritten with a
%                warning if present.
%       .model   the model equation (see above). Default: gaussian(1)+gaussian(2)+background, all with
%                empty .seed. Each entry:
%                  .seed    dot-path + optional selector '(main)'|'(secondary)'|'(k ...)'|none (see
%                           SEED PATHS). REQUIRED per entry, but never actually empty in practice:
%                           defaults to 'heuristic' (see SEED PATHS' own bare-'heuristic' note), which
%                           needs no prior fit at all.
%                  .mode    struct, per parameter: 'perVessel'|'perFrame'|'perRun'|'perKnot'|
%                           'perRunPoly([<degrees>])' (see mode='perKnot'/mode='perRunPoly' above;
%                           'perKnot' is ts-path only). Default 'perVessel'.
%                  .fixed   struct, per parameter, logical. Default false.
%                  .start/.lower/.upper  struct, per parameter, LITERAL numeric overrides. Default: from
%                           the seed peak / from .seedBounds.
%                  .seedBounds  PER MODEL TYPE (2026-09-08) -- only the bound names that type's own
%                           componentModelSpec.m parameters need; setting any other name errors
%                           (fitVesselPatchTimeSeries:badSeedBoundsField). gaussian: .x0y0/.theta
%                           (additive) and .a/.radius/.aspectRatio (multiplicative). background: .b
%                           (multiplicative) only. Same additive/multiplicative semantics as
%                           fitPatchVessels.m's own opts.seedBounds. Defaults: x0y0 [] (= 2 voxels, per
%                           axis), theta Inf (unbounded), the rest (a/radius/aspectRatio/b) 2.
%       .timeAvg logical, default FALSE -- see RUN HANDLING/opts.timeAvg above: collapses fully
%                (perVessel) or per-run only (perRun), by the finest mode requested;
%                perFrame/perKnot/perRunPoly are rejected. Refused on the irf-path.
%       .addBase0 logical, default TRUE -- IRF-PATH ONLY, see opts.addBase0 above.
%       .includeMask/.excludeMask  char | cellstr of vessel.(fld).rois label(s). Default {}.
%       .nIsochromatPerVoxDim  positive ODD integer. Default 7.
%       .progress logical, default true.
%
% OUTPUT
%   vessel : input with vessel.(fld).(outFld) added, per vessel --
%       .model.<type>(k)  mirrors opts.model. Per entry: .params (one field per model parameter, each a
%                 {1 x nRun} cell when it varies between runs, bare otherwise -- see RUN HANDLING),
%                 .start/.lower/.upper (same shape), .mode/.fixed (echoed), .seedBounds (the spec the
%                 bounds were built from, defaults filled in -- what the 2-arg call copies back),
%                 .peakNum/.label/.seed (provenance), .priorFit -- the seed peak verbatim -- and, for any
%                 parameter whose mode is 'perRunPoly([<degrees>])': .runPolyCoef (one field per such
%                 parameter, a {1 x nRun} cell of [1 x nCoef] raw fitted coefficient vectors -- the ONLY
%                 place the actual free fit variables of that parameter are visible directly; .params
%                 itself is the per-frame EXPANSION, same shape/contract as an ordinary 'perFrame'
%                 parameter) and .runPolyDeg (the resolved degree list). See fitComponentFit.m's own
%                 mode='perRunPoly' doc for the full seeding/expansion/wrapping contract.
%       .a/.x0/.y0/.radius/.aspectRatio/.theta/.b  FLAT ALIAS of the MAIN gaussian's own parameters
%                 (plus b from the background), so numeric consumers written against
%                 fitVesselTimeSeries.m's own flat output keep working unchanged. .runPolyCoef/
%                 .runPolyDeg are mirrored onto this flat alias too, for any of the MAIN gaussian's (or
%                 background's) own parameters using mode='perRunPoly' -- same contract as the
%                 per-component .model.<type>(k) fields above.
%                 CAVEAT, stated rather than buried: the alias describes the main peak ONLY. A renderer
%                 that REBUILDS a predicted image from it (buildGaussianFitDiag.m, and so
%                 fitVesselTimeSeriesDiag.m) will therefore draw a main-peak-only model and show the
%                 secondary peaks as residual. The numbers are right; the reconstructed picture is
%                 incomplete. A genuinely multi-component diagnostic is future work, not this stage.
%       .fits     TIMEAVG ONLY -- a fitPatchVessels.m-shaped struct array (one entry per gaussian
%                 component, .a/.x0/.y0/.radius/.aspectRatio/.theta/.b/.start/.lower/.upper/.peakNum/
%                 .includeMask/.excludeMask/.nIso/.priorFit), so every existing consumer of a
%                 'simultaneous'-mode result -- fitPatchVesselsDiag.m, makeRoiFromVesselFit.m,
%                 alignVessel2.m, fitVesselTimeSeries.m's Shape-A seed -- reads this output unchanged.
%                 Absent for a per-frame fit: a .fits entry holds SCALAR parameter values, which is
%                 only well defined when there is one frame. The shared background's b is REPLICATED
%                 onto each entry (a derived alias of the one background, not K baselines) because
%                 drawGaussianPeakContour.m reads peak.b directly -- the same thing fitMultiVessel.m
%                 did with its own shared b.
%       .resnorm/.resnormPerFrame/.exitFlag/.iterations/.maxIterations/.functionTolerance
%       .fitMode/.timeAvg/.nIso/.includeMask/.excludeMask/.method  provenance. With timeAvg, .K and
%                 .bShared too, matching the old 'simultaneous' container shape.
%       .validMask  the ACTUAL resolved logical array fitComponentFit.m was fit against (included &
%                 ~excluded) -- NOT the same thing as .includeMask/.excludeMask above, which are only
%                 the caller's own input label list(s). .nValidVox is nnz(this).
%
% See also fitComponentFit, componentModelSpec, fitPatchVessels, fitVesselTimeSeries, fitVessel.

    if nargin==0
        printOptsHelp();
        vessel = defaultOpts();
        return
    end
    if nargin==1 && (ischar(vessel) || isstring(vessel))
        % NAMED-MODEL CALL (see file header) -- a single char/string argument names a model template
        % rather than supplying a vessel. Same no-real-work contract as the no-arg call, just with
        % opts.model.gaussian built from the requested template.
        modelName = strtrim(char(vessel));
        namedOpts = defaultOpts(modelName);   % validates modelName first -- throws before any printing
        fprintf('fitVesselPatchTimeSeries: returning default opts for named model template ''%s''.\n', modelName);
        printOptsHelp();
        vessel = namedOpts;
        return
    end
    if nargin==2
        % 2-ARG CALL (see file header) -- vessel + fld together point at an EXISTING starting fit to
        % derive opts.model from, rather than supplying a real vessel + track to fit. Same "just give
        % me opts, no real work" contract as the no-arg/NAMED-MODEL calls, mirroring
        % drawVesselPatchPeaks.m's own established 2-arg call. NOTE fld means something DIFFERENT here
        % than in the 3-arg real-fit call below (a fit-result path, e.g. 'tsIm.gaussAnat', one level
        % DEEPER than the track path, e.g. 'tsIm', the 3-arg call's fld) -- same divergence
        % drawVesselPatchPeaks.m's own file header already documents for its own fld, now also true
        % across this one function's own call forms; see defaultOptsFromFit.
        vessel = defaultOptsFromFit(vessel, fld);
        return
    end
    if nargin<2 || isempty(fld); error('fitVesselPatchTimeSeries:noFld', 'fld is required.'); end
    if nargin<3; opts = struct(); end
    opts = fillOptsDefaults(opts);
    assert(~isempty(opts.outFld), 'fitVesselPatchTimeSeries:noOutFld', ...
        'opts.outFld is required (no default) -- name the result explicitly.');

    wasCell = iscell(vessel);
    if wasCell; cellSz = size(vessel); vessel = reshape([vessel{:}], cellSz); end
    nV = numel(vessel);

    for v = 1:nV
        vessel(v) = fitOneVessel(vessel(v), fld, opts, v, nV);
    end

    if wasCell; vessel = reshape(num2cell(vessel), cellSz); end
end

% ---------------------------------------------------------------------------
function vesselV = fitOneVessel(vesselV, fld, opts, v, nV)
    idStr = strtrim([char(string(vesselV.sId)) ' ' char(string(vesselV.label))]);
    % Printed IMMEDIATELY, before any setup work below (mask resolution, seed/motion resolution,
    % buildComponents) -- a caller staring at a silent prompt during a long fit has no way to tell
    % the call even started until now, since the more informative pre-solve message further down
    % (component/frame counts) has to wait on that setup finishing first.
    if opts.progress
        fprintf('fitVesselPatchTimeSeries: starting vessel %d/%d (%s) ...\n', v, nV, idStr);
    end
    base = getNestedField(vesselV, fld);
    assert(isstruct(base), 'fitVesselPatchTimeSeries:noFld', ...
        'vessel %s: vessel.%s is missing or not a struct.', idStr, fld);

    % WHICH TIME SERIES -- structural dispatch, the same test fitVesselTimeSeries.m uses: an irf-path
    % container is one carrying its own .irfIm.im; anything else must carry .im directly. Not a naming
    % convention on fld.
    isIrf = isfield(base,'irfIm') && isstruct(base.irfIm) && isfield(base.irfIm,'im') && ~isempty(base.irfIm.im);
    if isIrf
        % The IRF image's "frames" are KNOT DELAYS, not acquired frames, and its own container has no
        % .vSize/.rois -- those come from fld's PARENT path, the real ts<Name> track the IRF was fit
        % from, which shares the identical crop and voxel grid.
        segs = strsplit(fld, '.');
        if numel(segs) > 1; srcFldPath = strjoin(segs(1:end-1), '.'); else; srcFldPath = fld; end
        srcBase = getNestedField(vesselV, srcFldPath);
        assert(isstruct(srcBase) && isfield(srcBase,'vSize'), 'fitVesselPatchTimeSeries:noSrcBase', ...
            ['vessel %s: vessel.%s (fld''s own parent path, expected to be the ts track the IRF was fit ' ...
             'from) has no .vSize -- needed to fit vessel.%s.irfIm.'], idStr, srcFldPath, fld);
        imIrf = base.irfIm.im;
        if ~opts.addBase0
            % fitIRF.m BAKES the degree-0 (DC) baseline into irfIm.im at fit time, so false must
            % SUBTRACT it back out (sign=-1) to recover the pure stimulus-locked deviation -- it is not
            % a no-op skip. irfAddBase0.m hard-errors if there was never a base0 to remove, which is
            % correct: asking to subtract a baseline that was never fit is a real caller mistake.
            imIrf = irfAddBase0(base.irfIm, imIrf, [], -1);
        end
        imC = {imIrf};              % ONE pseudo-run: an IRF image is already pooled across runs
        metaBase = srcBase;         % .vSize/.rois borrowed from the source track
        maskFld = srcFldPath;
    else
        looksLikeFitResult = isstruct(base) && isfield(base,'model') && isstruct(base.model) ...
            && isfield(base.model,'gaussian') && isfield(base,'fits');
        if looksLikeFitResult
            % fld NAMED A FIT RESULT, NOT A TRACK (2026-09-10, real user mistake caught at the prompt) --
            % the 3-arg call's fld is the TRACK path (e.g. 'tsIm'), one level SHALLOWER than a fit-result
            % path (e.g. 'tsIm.gaussAnat'); only the 2-arg "give me opts" call and an
            % opts.model.<comp>.seed value take the deeper fit-result form (see file header 2-ARG CALL,
            % defaultOptsFromFit). A fit result never carries its own .im -- the image lives on fld's
            % PARENT, the track that was actually fit -- so that mismatch is caught explicitly here,
            % instead of letting the generic "neither ... is present" assert below fire with no hint of
            % the actual mistake.
            segs = strsplit(fld, '.');
            if numel(segs) > 1; parentFld = strjoin(segs(1:end-1), '.'); else; parentFld = fld; end
            error('fitVesselPatchTimeSeries:fldIsFitResultNotTrack', ...
                ['vessel %s: vessel.%s looks like a PRIOR fitVesselPatchTimeSeries.m result (it has ' ...
                 '.model.gaussian/.fits), not a ts-track -- the 3-arg call''s fld must name the TRACK ' ...
                 'itself (e.g. ''%s''), not the fit result. ''%s'' is only valid as the 2-arg call''s fld ' ...
                 '(opts = fitVesselPatchTimeSeries(vessel, ''%s'')) or as an opts.model.<comp>.seed ' ...
                 'value. Did you mean fitVesselPatchTimeSeries(vessel, ''%s'', opts)?'], ...
                idStr, fld, parentFld, fld, fld, parentFld);
        end
        assert(isfield(base,'im') && ~isempty(base.im), 'fitVesselPatchTimeSeries:noIm', ...
            'vessel %s: neither vessel.%s.irfIm.im nor vessel.%s.im is present.', idStr, fld, fld);
        imC = base.im; if ~iscell(imC); imC = {imC}; end
        metaBase = base;
        maskFld = fld;
    end
    runLens = cellfun(@(x) size(x,4), imC);
    nRun = numel(imC);
    voxSz = metaBase.vSize(:);
    voxSz2 = [voxSz(1) voxSz(min(2,numel(voxSz)))];

    % Runs concatenated into one series; runIdx records which original run each frame came from, which
    % is what makes mode='perRun' meaningful downstream (see file header RUN HANDLING).
    if opts.timeAvg
        assert(~isIrf, 'fitVesselPatchTimeSeries:timeAvgIrfPathNotAllowed', ...
            ['vessel %s: opts.timeAvg is ts-path only -- averaging an IRF image across its KNOT DELAYS ' ...
             'collapses the response shape that image exists to represent.'], idStr);
        % WHICH AXIS timeAvg is allowed to collapse depends on the FINEST mode requested anywhere in
        % opts.model (2026-09-08, Seb's own design): every parameter 'perVessel' -> collapse everything,
        % runs included, into ONE image (today's original behaviour, and the same reduction
        % fitMultiVessel.m performs: cat then mean over dim 4). At least one parameter 'perRun' -> only
        % average WITHIN each run, keeping runs separate as nRun pseudo-frames, so a 'perRun' parameter
        % still ends up with one value per run instead of being silently pooled into a single global
        % value it never asked for. 'perFrame'/'perKnot'/'perRunPoly' have no time axis left once ANY
        % time-averaging happens, so they are rejected outright rather than silently collapsed too --
        % same deliberately loud convention opts.seedMotion+timeAvg already uses just below.
        finestMode = finestModeAcrossModel(opts.model);
        assert(~any(strcmp(finestMode, {'perFrame','perKnot'})) && ~startsWith(finestMode, 'perRunPoly'), ...
            'fitVesselPatchTimeSeries:timeAvgFinerThanPerRun', ...
            ['vessel %s: opts.timeAvg=true collapses time, but opts.model requests mode=''%s'' on at ' ...
             'least one parameter -- there is no frame/knot axis left for it to vary over once time is ' ...
             'averaged. Set every parameter to ''perVessel'' or ''perRun'', or drop opts.timeAvg.'], ...
            idStr, finestMode);
        if strcmp(finestMode, 'perRun')
            % Average WITHIN each run only -- one pseudo-frame per run, runs kept distinct.
            im4 = zeros(size(imC{1},1), size(imC{1},2), 1, nRun);
            for r = 1:nRun; im4(:,:,1,r) = mean(imC{r}(:,:,1,:), 4, 'omitnan'); end
            runIdx = 1:nRun;
        else
            % Every parameter 'perVessel' -- pool every run's every frame into ONE image.
            acc = cat(4, imC{:});
            im4 = mean(acc(:,:,1,:), 4, 'omitnan');
            runIdx = 1;
        end
    else
        im4 = cat(4, imC{:});
        im4 = im4(:,:,1,:);
        runIdx = repelem(1:nRun, runLens);
    end
    [ny,nx,~,T] = size(im4);

    largeSz = [ny nx];
    validMask = resolveIncludeMaskLocal(vesselV, maskFld, opts.includeMask, largeSz) ...
                & ~unionMasksLocal(vesselV, maskFld, opts.excludeMask, largeSz);
    assert(any(validMask(:)), 'fitVesselPatchTimeSeries:emptyValidRegion', ...
        'vessel %s: opts.includeMask/excludeMask leave no valid voxel.', idStr);

    % SEED MOTION, resolved once per vessel and applied to EVERY component's position seed. Resolved
    % against the PRE-timeAvg run cell, because a motion trace is inherently per frame -- see the
    % assert in buildComponents for why timeAvg and seedMotion are mutually exclusive.
    motionSeed = struct();
    if ~isempty(opts.seedMotion)
        motionSeed = resolveMotionSeed(vesselV, opts.seedMotion, imC, idStr);
    end
    [components, meta] = buildComponents(vesselV, opts, idStr, voxSz2, motionSeed, T, fld, runLens, isIrf);

    % Replicate requests, resolved BEFORE the (potentially slow) point-estimate fit so a skip note
    % prints before the wait rather than after. Same silently-skip-with-a-note convention
    % fitVesselTimeSeries.m uses: one opts struct is routinely reused across a mix of ts-path and
    % irf-path calls, and erroring on the path that cannot serve a request would make that impossible.
    useBoot = opts.boot.nRep > 0 && isIrf && isfield(base.irfIm,'imBoot') && ~isempty(base.irfIm.imBoot);
    if isIrf
        useNull = opts.null.nRep > 0 && isfield(base.irfIm,'imPerm') && ~isempty(base.irfIm.imPerm);
    else
        useNull = opts.null.nRep > 0;   % ts-path: a COUNT of surrogates to generate here
    end
    if opts.progress
        if opts.boot.nRep > 0 && ~useBoot
            fprintf(['  fitVesselPatchTimeSeries: vessel %s: opts.boot.nRep>0 but this is not an ' ...
                     'irf-path fld with a populated .irfIm.imBoot -- skipping.\n'], idStr);
        end
        if opts.null.nRep > 0 && ~useNull
            fprintf(['  fitVesselPatchTimeSeries: vessel %s: opts.null.nRep>0 but this irf-path fld ' ...
                     'has no populated .irfIm.imPerm -- skipping.\n'], idStr);
        end
    end
    if useBoot || useNull; ensureParPool(opts.nPool); end
    % Anything already sitting at the output field from a PREVIOUS call -- see the file header on the
    % already-in-vessel case: reused unless reComputeCache says otherwise, and validated against the
    % same spec key the cache uses so a stale reuse is impossible.
    priorRes = getNestedField(vesselV, [fld '.' opts.outFld]);
    if isIrf; priorRes = getNestedField(vesselV, [fld '.irfIm.' opts.outFld]); end
    priorBoot = priorReplicates(priorRes, 'boot');
    priorNull = priorReplicates(priorRes, 'null');

    % mode='perKnot' -- resolve the frame->knot map, and drop any frame the design does not cover.
    % fitComponentFit.m already implements the granularity itself (verified in stage 1); the only thing
    % missing at this level is the mapping, which needs the IRF design matrix this function can see and
    % the solver cannot.
    knotIdx = [];
    keepFrame = [];
    if anyPerKnot(components)
        assert(~isIrf, 'fitVesselPatchTimeSeries:perKnotIrfPathNotAllowed', ...
            ['vessel %s: mode=''perKnot'' is ts-path only -- an irf-path fit''s own "frames" ALREADY ' ...
             'are the knots, so there is nothing left to group.'], idStr);
        assert(~opts.timeAvg, 'fitVesselPatchTimeSeries:perKnotWithTimeAvg', ...
            'vessel %s: mode=''perKnot'' is meaningless with opts.timeAvg (one frame, no knots).', idStr);
        [knotIdx, keepFrame] = resolveKnotIdx(vesselV, fld, nRun, runLens, T, idStr);
        if ~all(keepFrame)
            % KNOT COVERAGE -- run-edge dead time before the first / after the last modelled knot has
            % no knot to belong to (getIRFmat.m's own unmodelled gaps). Those frames are EXCLUDED from
            % the fit rather than forced into a group they are not in, and their rows are restored as
            % NaN afterwards so the reported output still spans the full raw frame count in the
            % original order.
            if opts.progress
                fprintf(['  fitVesselPatchTimeSeries: vessel %s: %d of %d frame(s) have no knot ' ...
                         'coverage -- excluded from the fit, restored as NaN in the output.\n'], ...
                        idStr, nnz(~keepFrame), T);
            end
            im4 = im4(:,:,1,keepFrame);
            runIdx = runIdx(keepFrame);
            T = size(im4,4);
        end
    end

    if opts.progress
        fprintf('fitVesselPatchTimeSeries: vessel %d/%d (%s), %d component(s), %d frame(s) ...\n', ...
            v, nV, idStr, numel(components), T);
    end
    t0 = tic;
    fcOpts = struct('nIsochromatPerVoxDim',opts.nIsochromatPerVoxDim, 'runIdx',runIdx, ...
                    'knotIdx',knotIdx, 'cacheFrameInvariant',true);
    fit = fitComponentFit(im4, validMask, voxSz2, components, fcOpts);
    if opts.progress
        fprintf('fitVesselPatchTimeSeries: vessel %d/%d (%s) done in %.1f s -- resnorm=%.6g, %d free\n', ...
            v, nV, idStr, toc(t0), fit.resnorm, fit.nFree);
    end

    % Restore the frames KNOT COVERAGE excluded, as NaN, so every reported per-frame array spans the
    % full raw frame count in the original order -- what the per-run split below, and every consumer
    % after it, assumes.
    if ~isempty(keepFrame) && ~all(keepFrame)
        fit = reinflateFit(fit, keepFrame);
        T = numel(keepFrame);
    end

    % With timeAvg's all-perVessel reduction the runs are POOLED into one frame, so there are no
    % per-run values to report -- one run, so every value comes out bare (collapseIdentical), never
    % {1xnRun} broadcasting one number as if it were nRun of them. With timeAvg's perRun reduction,
    % T==nRun (one pseudo-frame per run) and
    % outRunLens is a vector of nRun ones, so splitPerRun's own length==nRun branch gives each run
    % exactly its own value, and a 'perVessel' parameter (length 1) still broadcasts across all of them
    % -- the SAME three-way split contract as the non-timeAvg path, just with 1 real frame per run
    % instead of runLens(r). The IRF path is NOT split at all: its "frames" are knot delays, and a knot
    % delay does not belong to a scan run, so a per-run cell would be a category error. Same call
    % fitVesselTimeSeries.m makes (its own PER-RUN OUTPUT is ts-path only, irf-path fields stay flat).
    outRunLens = runLens;
    % ones(1,T), not T itself: when T==1 (all-perVessel) this IS just the scalar 1, identical to
    % today's behaviour; when T==nRun (perRun) it is a proper nRun-long vector of 1s, which is what
    % splitPerRun's length==nRun branch needs to hand each run its own single value.
    if opts.timeAvg; outRunLens = ones(1,T); end
    % REPLICATES -- see the file header's own BOOTSTRAP AND NULL REPLICATES section for the vocabulary.
    % Both kinds refit the EXACT SAME component spec (same components, same fcOpts) against a different
    % image, so neither can drift from the point estimate in what it is actually fitting.
    repOut = struct();
    if useBoot
        [repOut.boot, repOut.bootStore] = runReplicates(base.irfIm.imBoot, [], components, validMask, ...
            voxSz2, fcOpts, opts, base.irfIm, isIrf, idStr, 'boot', v, nV, [fld '_' opts.outFld], priorBoot);
    end
    if useNull
        if isIrf
            [repOut.null, repOut.nullStore] = runReplicates(base.irfIm.imPerm, [], components, validMask, ...
                voxSz2, fcOpts, opts, base.irfIm, isIrf, idStr, 'null', v, nV, [fld '_' opts.outFld], priorNull);
        else
            [repOut.null, repOut.nullStore] = runReplicates([], imC, components, validMask, ...
                voxSz2, fcOpts, opts, [], isIrf, idStr, 'null', v, nV, [fld '_' opts.outFld], priorNull);
        end
    end

    nKnotGroups = 0; if ~isempty(knotIdx); nKnotGroups = numel(unique(knotIdx)); end
    res = assembleResult(fit, meta, opts, outRunLens, T, validMask, ~isIrf, nKnotGroups, repOut);
    if isIrf
        outPath = [fld '.irfIm.' opts.outFld];
    else
        outPath = [fld '.' opts.outFld];
    end
    if ~isempty(getNestedField(vesselV, outPath))
        warning('fitVesselPatchTimeSeries:outFldExists', ...
            'vessel %d/%d (%s): vessel.%s already exists -- overwriting.', v, nV, idStr, outPath);
    end
    vesselV = setNestedField(vesselV, outPath, res);
end

% ---------------------------------------------------------------------------
function [components, meta] = buildComponents(vesselV, opts, idStr, voxSz2, motionSeed, T, fld, runLens, isIrf)
    % Expand opts.model into the flat component list fitComponentFit.m consumes, resolving each entry's
    % own seed path + selector into concrete peaks. `claimed` is keyed by seed PATH so that two entries
    % pointing at DIFFERENT fits don't interfere with each other's "remaining" accounting.
    % runLens/isIrf describe the TARGET of this fit -- needed by collapseSeedToMode, which maps a seed
    % peak's own granularity (a prior per-run/per-frame fit) onto each parameter's requested mode.
    modelTypes = fieldnames(opts.model);
    components = struct('model',{}, 'label',{}, 'mode',{}, 'fixed',{}, ...
                        'start',{}, 'lower',{}, 'upper',{});
    meta = struct('type',{}, 'entry',{}, 'peakNum',{}, 'seed',{}, 'priorFit',{}, 'isMainGaussian',{}, ...
                  'seedBounds',{});
    claimed = struct();   % sanitized path -> already-claimed array positions

    for mt = 1:numel(modelTypes)
        typeName = modelTypes{mt};
        entries = opts.model.(typeName);
        spec = componentModelSpec(typeName);
        for e = 1:numel(entries)
            en = entries(e);
            assert(isfield(en,'seed') && ~isempty(en.seed), 'fitVesselPatchTimeSeries:noSeed', ...
                ['vessel %s: opts.model.%s(%d).seed is empty -- every component entry needs a dot-path ' ...
                 'to a fitPatchVessels.m result (optionally with a (main)/(k) peak selector), or the ' ...
                 'bare word ''heuristic'' (defaults to THIS fld''s own ''.heuristic'').'], ...
                idStr, typeName, e);
            % A BARE 'heuristic' (no path at all -- defaultOpts()'s own default seed, 2026-09-08,
            % Seb's own ask) means "run the heuristic on THIS CALL's own fld", i.e. exactly
            % [fld '.heuristic'] -- sparing a caller from repeating fld in every entry's seed just to
            % ask for the same thing parseSeedPath already recognizes with it spelled out.
            seedStr = en.seed;
            if strcmpi(strtrim(seedStr), 'heuristic'); seedStr = [fld '.heuristic']; end
            [pathStr, sel] = parseSeedPath(seedStr);
            if sel.heuristic
                peaks = resolveHeuristicSeedPeak(vesselV, pathStr, idStr, opts.includeMask, opts.excludeMask);
            else
                peaks = resolveSeedPeaks(vesselV, pathStr, idStr);
            end
            key = matlab.lang.makeValidName(pathStr);
            if ~isfield(claimed, key); claimed.(key) = []; end

            % CLAIMING applies only to models that expand PER PEAK (spec.perPeak -- see
            % componentModelSpec.m). A SINGLETON model like 'background' is not selecting peaks to
            % model at all; its seed only says where to read a starting value from. Letting it claim
            % would be actively wrong in both directions: run last (as the default template does) and
            % an unindexed background finds every peak already claimed, resolving to nothing and
            % seeding b from an empty mean (NaN) -- caught on real data, not by inspection; run first
            % and it would starve the gaussian entries of the peaks they were meant to model.
            switch sel.kind
                case 'main';  idx = 1;
                case 'secondary'
                    % Every peak but the main one, by ROLE rather than by number (2026-09-10, Seb's own
                    % ask): EMPTY on a single-peak vessel, which is not an error -- the same "matches
                    % nothing, contributes nothing" contract an unindexed entry has, so ONE opts struct
                    % carrying a '(secondary)' entry runs unchanged across a batch whose vessels have
                    % different peak counts. Unlike the unindexed form it does not depend on an earlier
                    % entry having claimed the main peak first.
                    idx = 2:numel(peaks);
                case 'index'; idx = sel.idx;
                case 'all'
                    if spec.perPeak
                        idx = setdiff(1:numel(peaks), claimed.(key), 'stable');
                    else
                        idx = 1:numel(peaks);   % singleton: "all peaks" means all, never "the remainder"
                    end
            end
            if ~any(strcmp(sel.kind,{'all','secondary'}))
                assert(all(idx>=1 & idx<=numel(peaks)), 'fitVesselPatchTimeSeries:seedPeakOutOfRange', ...
                    ['vessel %s: opts.model.%s(%d).seed selects peak(s) %s, but ''%s'' has only %d ' ...
                     'INCLUDED peak(s). An explicit selector past the end is a real mistake; use an ' ...
                     'unindexed seed if you meant "whatever remains".'], ...
                    idStr, typeName, e, mat2str(idx), pathStr, numel(peaks));
            end
            if spec.perPeak; claimed.(key) = union(claimed.(key), idx); end

            if ~spec.perPeak
                % Structurally singleton: ONE additive constant for the patch, regardless of peak count.
                % A second one would be exactly degenerate with the first.
                assert(isscalar(entries), 'fitVesselPatchTimeSeries:multipleBackgrounds', ...
                    ['vessel %s: opts.model.%s must have exactly one entry (got %d) -- it is a ' ...
                     'SINGLETON model (componentModelSpec.m .perPeak=false); a second one would be ' ...
                     'degenerate with the first.'], idStr, typeName, numel(entries));
                assert(~isempty(idx), 'fitVesselPatchTimeSeries:backgroundNoSeed', ...
                    ['vessel %s: opts.model.%s.seed resolved to no peak to take a starting value ' ...
                     'from.'], idStr, typeName);
                % b seeds from the MEAN of the resolved peaks' own b. For a 'simultaneous' seed every
                % peak carries the identical shared b, so the mean is exactly that value; for a
                % 'sequential' seed (a running per-peak background estimate) it is the sensible summary.
                % Each peak's b is first mapped onto THIS entry's own mode (see collapseSeedToMode), so
                % the mean is taken per run/frame where the target keeps that granularity.
                lbl = sprintf('%s', typeName);
                bSeeds = cell(1, numel(idx));
                for ii = 1:numel(idx)
                    assert(isfield(peaks(idx(ii)),'b') && ~isempty(peaks(idx(ii)).b), ...
                        'fitVesselPatchTimeSeries:backgroundSeedNoB', ...
                        ['vessel %s: seed ''%s'' carries no background value to seed opts.model.%s from ' ...
                         '(a fitVesselPatchTimeSeries.m result fit without a background component?).'], ...
                        idStr, pathStr, typeName);
                    bSeeds{ii} = collapseSeedToMode(peaks(idx(ii)).b, subOr(en,'mode','b','perVessel'), ...
                                                    runLens, isIrf, idStr, lbl, 'b');
                end
                seedVals = struct('b', mean(cat(1, bSeeds{:}), 1));
                [components(end+1), sbSpec] = mkComponent(typeName, spec, en, seedVals, voxSz2, lbl); %#ok<AGROW>
                meta(end+1) = struct('type',typeName, 'entry',e, 'peakNum',[], 'seed',en.seed, ...
                                     'priorFit',[], 'isMainGaussian',false, 'seedBounds',sbSpec); %#ok<AGROW>
            else
                for ii = 1:numel(idx)
                    pk = peaks(idx(ii));
                    pkNum = idx(ii);
                    if isfield(pk,'peakNum') && ~isempty(pk.peakNum); pkNum = pk.peakNum; end
                    lbl = sprintf('%s#peak%d', typeName, pkNum);
                    seedVals = struct();
                    for pn = 1:numel(spec.paramNames)
                        nm = spec.paramNames{pn};
                        assert(isfield(pk,nm), 'fitVesselPatchTimeSeries:seedMissingParam', ...
                            ['vessel %s: seed peak %d of ''%s'' has no .%s -- it does not look like a ' ...
                             '''%s''-shaped fit result.'], idStr, idx(ii), pathStr, nm, typeName);
                        seedVals.(nm) = collapseSeedToMode(pk.(nm), subOr(en,'mode',nm,'perVessel'), ...
                                                           runLens, isIrf, idStr, lbl, nm);
                    end
                    seedVals = applyMotionToSeed(seedVals, motionSeed, en, T, idStr, lbl, opts);
                    [components(end+1), sbSpec] = mkComponent(typeName, spec, en, seedVals, voxSz2, lbl); %#ok<AGROW>
                    meta(end+1) = struct('type',typeName, 'entry',e, 'peakNum',pkNum, 'seed',en.seed, ...
                                         'priorFit',pk, 'isMainGaussian', ...
                                         strcmp(typeName,'gaussian') && idx(ii)==1, ...
                                         'seedBounds',sbSpec); %#ok<AGROW>
                end
            end
        end
    end
    assert(~isempty(components), 'fitVesselPatchTimeSeries:noComponents', ...
        'vessel %s: opts.model expanded to zero components -- nothing to fit.', idStr);
end

% ---------------------------------------------------------------------------
function seedVals = applyMotionToSeed(seedVals, motionSeed, en, T, idStr, lbl, opts)
% Subtract the per-frame motion from this component's POSITION seeds, turning a scalar seed into a
% [1 x T] row. SIGN: getPreprocMotion.m reports the correction preprocessing APPLIED, so
% preprocessed-position minus correction is where the feature sits in the RAW data -- see
% resolveMotionSeed.m. Only x0/y0/z0 are touched; amplitude/width/baseline seeds are left alone,
% because a motion trace has none of those.
    if isempty(fieldnames(motionSeed)); return; end
    assert(~opts.timeAvg, 'fitVesselPatchTimeSeries:seedMotionWithTimeAvg', ...
        ['vessel %s: opts.seedMotion and opts.timeAvg=true are mutually exclusive. A motion trace is ' ...
         'inherently per-frame and timeAvg collapses every frame into one before fitting, so there is ' ...
         'no frame axis left for it to correct. Fit the time series, or drop opts.seedMotion.'], idStr);
    for pc = {'x0','y0','z0'}
        p = pc{1};
        if ~isfield(seedVals, p) || ~isfield(motionSeed, p) || isempty(motionSeed.(p)); continue; end
        m = motionSeed.(p)(:).';
        assert(numel(m) == T, 'fitVesselPatchTimeSeries:seedMotionFrameCount', ...
            ['vessel %s: motion seed .%s covers %d frame(s) but the fit target has %d -- the motion ' ...
             'source and the fit target disagree on the frame axis.'], idStr, p, numel(m), T);
        b = seedVals.(p); b = b(:).';
        assert(isscalar(b) || numel(b) == T, 'fitVesselPatchTimeSeries:seedMotionSeedShape', ...
            'vessel %s: base seed .%s has %d value(s), expected 1 or %d.', idStr, p, numel(b), T);
        % A scalar base broadcasts against the [1 x T] motion row, so the seed becomes genuinely
        % PER-FRAME -- which is the entire point of the feature.
        seedVals.(p) = b - m;
        % THE MODE HAS TO AGREE, and this must ERROR rather than silently collapse. A per-frame seed
        % under a perVessel/perRun parameter has more values than that parameter has groups, and the
        % only ways to proceed would be to average it (inventing a number the caller never asked for)
        % or to keep frame 1 (silently discarding the correction that was the whole reason for setting
        % opts.seedMotion). Both are the silent-fallback shape this project rejects -- fail loudly.
        md = 'perVessel';
        if isfield(en,'mode') && isstruct(en.mode) && isfield(en.mode, p) && ~isempty(en.mode.(p))
            md = en.mode.(p);
        end
        assert(strcmp(md,'perFrame'), 'fitVesselPatchTimeSeries:seedMotionModeMismatch', ...
            ['vessel %s: opts.seedMotion makes the .%s seed of component ''%s'' per-frame, but its ' ...
             'mode is ''%s''. Set opts.model.<type>(k).mode.%s = ''perFrame'' for any position a ' ...
             'motion seed corrects, or drop opts.seedMotion -- averaging the correction away, or ' ...
             'keeping only frame 1 of it, would silently discard the very thing it was set for.'], ...
            idStr, p, lbl, md, p);
    end
end

% ---------------------------------------------------------------------------
function [c, sbSpec] = mkComponent(typeName, spec, en, seedVals, voxSz2, label)
    % One fitComponentFit.m component: granularity from the entry, start from the seed peak (or a
    % literal override), bounds from .seedBounds applied to that same seed (or a literal override).
    % A seed value may be a VECTOR here (a per-run/per-frame start mapped from a prior fit of the same
    % granularity -- see collapseSeedToMode); additiveSeedBound.m/multiplicativeSeedBound.m are
    % elementwise, so the bound comes out with the same shape and fitComponentFit.m's own fitToGroups
    % accepts it directly. sbSpec is the caller-facing .seedBounds this component was built from
    % (defaults filled in), echoed into the output so a later 2-arg call can copy it verbatim.
    [sb, sbSpec] = resolveSeedBounds(en, voxSz2, typeName);
    c = struct('model',typeName, 'label',label, 'mode',struct(), 'fixed',struct(), ...
               'start',struct(), 'lower',struct(), 'upper',struct());
    for pn = 1:numel(spec.paramNames)
        nm = spec.paramNames{pn};
        % 'perVessel' default -- deliberately NOT fitComponentFit.m's own 'perFrame' (see file header
        % GRANULARITY DEFAULTS).
        c.mode.(nm)  = subOr(en, 'mode',  nm, 'perVessel');
        c.fixed.(nm) = logical(subOr(en, 'fixed', nm, false));
        sVal = subOr(en, 'start', nm, seedVals.(nm));
        c.start.(nm) = sVal;
        [lo, hi] = seedBoundFor(nm, sVal, sb);
        c.lower.(nm) = subOr(en, 'lower', nm, lo);
        c.upper.(nm) = subOr(en, 'upper', nm, hi);
    end
end

% ---------------------------------------------------------------------------
function [lo, hi] = seedBoundFor(name, seedVal, sb)
    % Same additive/multiplicative split, and the same per-parameter Inf special-casing, that
    % fitPatchVessels.m's 'simultaneous' mode uses -- via the shared additiveSeedBound.m/
    % multiplicativeSeedBound.m rather than a local copy of either formula.
    switch name
        case 'x0';    [lo,hi] = additiveSeedBound(seedVal, sb.x0);
        case 'y0';    [lo,hi] = additiveSeedBound(seedVal, sb.y0);
        case 'theta'; [lo,hi] = additiveSeedBound(seedVal, sb.theta);
        case 'a'
            % 'a' is genuinely bipolar (an amplitude can be a negative dip), so scalar Inf must mean the
            % whole real line rather than whatever the general formula makes of dividing by Inf.
            if isscalar(sb.a) && isinf(sb.a); lo = -Inf; hi = Inf;
            else; [lo,hi] = multiplicativeSeedBound(seedVal, sb.a); end
        case 'b'
            % 'b' is a background level, non-negative by physical meaning -- scalar Inf means [0,Inf).
            if isscalar(sb.b) && isinf(sb.b); lo = 0; hi = Inf;
            else; [lo,hi] = multiplicativeSeedBound(seedVal, sb.b); end
        case 'radius';      [lo,hi] = multiplicativeSeedBound(seedVal, sb.radius);
        case 'aspectRatio'
            [lo,hi] = multiplicativeSeedBound(seedVal, sb.aspectRatio);
            % Symmetrized around 1 for the same reason fitPatchVessels.m/fitMultiVessel.m symmetrize
            % theirs: theta is free here too, so canonicalizeRadiusAspectShape's reciprocal can push a
            % converged aspectRatio outside a non-symmetric bound.
            [lo,hi] = symmetrizeAspectRatioBound(lo, hi);
        otherwise;    lo = -Inf; hi = Inf;   % an unknown//future parameter is left unbounded
    end
end

% ---------------------------------------------------------------------------
function [sb, spec] = resolveSeedBounds(en, voxSz2, typeName)
    % PER-TYPE (2026-09-08, Seb's own ask): gaussian's and background's .seedBounds previously shared
    % ONE struct carrying every field of both (x0y0/theta/a/radius/aspectRatio AND b), so
    % opts.model.background.seedBounds showed five fields it can never use (background has no x0/y0/
    % theta/a/radius/aspectRatio parameter at all) -- confusing, and it silently accepted a caller
    % setting one of them with no effect. Now scoped to exactly the bound names THIS type's own
    % componentModelSpec.m paramNames need (defaultSeedBounds(typeName)), and any field name outside
    % that set is a loud error rather than a silent no-op.
    d = defaultSeedBounds(typeName);
    if isfield(en,'seedBounds') && isstruct(en.seedBounds)
        fn = fieldnames(en.seedBounds);
        for i = 1:numel(fn)
            assert(isfield(d, fn{i}), 'fitVesselPatchTimeSeries:badSeedBoundsField', ...
                ['opts.model.%s.seedBounds.%s is not a valid bound for a ''%s'' component -- only {%s} ' ...
                 'apply here (see that model''s own componentModelSpec.m paramNames).'], ...
                typeName, fn{i}, typeName, strjoin(fieldnames(d), ', '));
            if ~isempty(en.seedBounds.(fn{i})); d.(fn{i}) = en.seedBounds.(fn{i}); end
        end
    end
    spec = d;   % the merged caller-facing spec, BEFORE x0y0 is expanded into x0/y0 below
    sb = d;
    % x0y0's own default is "2 voxels", which is PER AXIS and so not knowable until the voxel size is --
    % same asymmetric-default/symmetric-override convention fitPatchVessels.m uses (a caller-set spec
    % applies identically to both axes). GAUSSIAN ONLY -- background has no x0y0 field to resolve.
    if isfield(sb, 'x0y0')
        if isempty(sb.x0y0); sb.x0 = 2*voxSz2(2); sb.y0 = 2*voxSz2(1);
        else;                sb.x0 = sb.x0y0;     sb.y0 = sb.x0y0;
        end
    end
end

% ---------------------------------------------------------------------------
function d = defaultSeedBounds(typeName)
    % One bound spec per REAL parameter of typeName (componentModelSpec.m's own paramNames) -- x0/y0
    % share the single grouped 'x0y0' spec (see resolveSeedBounds), never split into x0/y0 separately.
    switch typeName
        case 'gaussian';   d = struct('x0y0',[], 'theta',Inf, 'a',2, 'radius',2, 'aspectRatio',2);
        case 'background'; d = struct('b',2);
        otherwise
            error('fitVesselPatchTimeSeries:noSeedBoundsForType', ...
                ['no default .seedBounds defined for model type ''%s'' -- add a case here alongside a ' ...
                 'new componentModelSpec.m entry.'], typeName);
    end
end

% ---------------------------------------------------------------------------
function s = collapseSeedToMode(val, mode, runLens, isIrf, idStr, lbl, nm)
    % Map ONE seed parameter's value, at whatever granularity the prior fit stored it, onto the
    % granularity THIS fit requests for it (2026-09-10, Seb's own rule): AVERAGE a finer seed into a
    % coarser target (per-frame -> per-run: mean within each run; anything -> per-vessel: one mean),
    % but NEVER average when the target is as fine as or finer than the seed -- keep every value the
    % seed has (a per-run seed into a per-run fit stays per run; into a per-frame fit each run's value
    % is repeated over that run's frames).
    %
    % val is either a {1 x nRunSeed} cell (one entry per run, each a scalar or a [1 x Tr] per-frame
    % row -- how every fitVesselPatchTimeSeries.m result stores a value that differs between runs) or a
    % BARE numeric (a value identical across runs, stored once -- see splitPerRun; also every .fits
    % entry, every heuristic seed, and every result of the pre-collapse era's {1x1} cell, which is the
    % same thing spelled longer). Output is what fitComponentFit.m's own fitToGroups accepts for the
    % target mode: a scalar (broadcast), or exactly one value per group.
    nRun = numel(runLens);
    if iscell(val); runs = cellfun(@(x) double(x(:).'), val, 'UniformOutput', false);
    else;           runs = {double(val(:).')};
    end
    allV = [runs{:}];
    assert(~isempty(allV), 'fitVesselPatchTimeSeries:emptySeedValue', ...
        'vessel %s: %s seed .%s is empty.', idStr, lbl, nm);
    if isIrf
        % An irf-path fit's "frames" are knot delays and its one pseudo-run is not a scan run, so a
        % ts-path seed's run/frame structure has no counterpart there: one value.
        s = mean(allV, 'omitnan'); return
    end
    switch mode
        case 'perRun'
            if isscalar(runs)
                s = mean(runs{1}, 'omitnan');   % identical across runs (or a single-run seed): broadcast
            else
                assert(numel(runs) == nRun, 'fitVesselPatchTimeSeries:seedRunCountMismatch', ...
                    ['vessel %s: %s seed .%s has %d run(s) but this fit has %d -- a per-run seed can ' ...
                     'only feed a fit of the SAME runs.'], idStr, lbl, nm, numel(runs), nRun);
                s = cellfun(@(x) mean(x, 'omitnan'), runs);
            end
        case 'perFrame'
            if isscalar(runs) && isscalar(runs{1})
                s = runs{1};
            elseif isscalar(runs)
                % One per-frame row identical across runs: every run must have that many frames.
                assert(all(runLens == numel(runs{1})), 'fitVesselPatchTimeSeries:seedFrameCountMismatch', ...
                    ['vessel %s: %s seed .%s has %d value(s) per run but this fit''s runs have %s ' ...
                     'frame(s).'], idStr, lbl, nm, numel(runs{1}), mat2str(runLens(:).'));
                s = repmat(runs{1}, 1, nRun);
            else
                assert(numel(runs) == nRun, 'fitVesselPatchTimeSeries:seedRunCountMismatch', ...
                    ['vessel %s: %s seed .%s has %d run(s) but this fit has %d -- a per-run seed can ' ...
                     'only feed a fit of the SAME runs.'], idStr, lbl, nm, numel(runs), nRun);
                s = cell(1, nRun);
                for r = 1:nRun
                    if isscalar(runs{r})
                        s{r} = repmat(runs{r}, 1, runLens(r));
                    else
                        assert(numel(runs{r}) == runLens(r), 'fitVesselPatchTimeSeries:seedFrameCountMismatch', ...
                            ['vessel %s: %s seed .%s run %d has %d value(s) but this fit''s run %d has ' ...
                             '%d frame(s).'], idStr, lbl, nm, r, numel(runs{r}), r, runLens(r));
                        s{r} = runs{r};
                    end
                end
                s = [s{:}];
            end
        case 'perKnot'
            % A knot-shaped seed (one row, identical across runs) passes through; fitComponentFit.m
            % checks its length against the resolved knot count. Anything else: one value.
            if isscalar(runs) && ~isscalar(runs{1}); s = runs{1}; else; s = mean(allV, 'omitnan'); end
        otherwise
            % 'perVessel', and 'perRunPoly([...])' (fitComponentFit.m requires a scalar start there,
            % broadcast onto every run's degree-0 coefficient): one value.
            s = mean(allV, 'omitnan');
    end
end

% ---------------------------------------------------------------------------
function v = subOr(en, group, name, dflt)
    v = dflt;
    if isfield(en,group) && isstruct(en.(group)) && isfield(en.(group),name) && ~isempty(en.(group).(name))
        v = en.(group).(name);
    end
end

% ---------------------------------------------------------------------------
function [pathStr, sel] = parseSeedPath(seedStr)
    % Split 'a.b.c(main)' / 'a.b.c(2 3)' / 'a.b.c' into a plain dot-path plus a selector. The trailing
    % '.main' spelling is also accepted, resolved only where a paren form would have been -- see file
    % header on why the paren form is canonical. 'a.b.c(heuristic)' / 'a.b.c.heuristic' are accepted the
    % same way -- see file header "'path.heuristic' / 'path(heuristic)'" note; sel.kind stays 'all' for
    % this case (there is no compound/indexed heuristic selector), with sel.heuristic=true marking it.
    seedStr = strtrim(char(string(seedStr)));
    sel = struct('kind','all', 'idx',[], 'heuristic',false);
    tok = regexp(seedStr, '^(.*?)\s*\(\s*([^)]*)\s*\)$', 'tokens', 'once');
    if ~isempty(tok)
        pathStr = strtrim(tok{1});
        inner = strtrim(tok{2});
        if strcmpi(inner,'main')
            sel.kind = 'main';
        elseif strcmpi(inner,'secondary')
            sel.kind = 'secondary';
        elseif strcmpi(inner,'heuristic')
            sel.heuristic = true;
        else
            idx = str2num(inner); %#ok<ST2NM>  -- accepts '2', '2 3', '2:4'
            assert(~isempty(idx) && all(idx==round(idx)) && all(idx>=1), ...
                'fitVesselPatchTimeSeries:badSeedSelector', ...
                ['seed selector ''(%s)'' must be ''main'', ''secondary'', ''heuristic'', or a vector of ' ...
                 'positive integer array positions (e.g. ''(2)'', ''(2 3)'', ''(2:4)'').'], inner);
            sel.kind = 'index'; sel.idx = idx(:).';
        end
        return
    end
    segs = strsplit(seedStr, '.');
    if numel(segs)>1 && strcmpi(segs{end},'main')
        pathStr = strjoin(segs(1:end-1), '.');
        sel.kind = 'main';
    elseif numel(segs)>1 && strcmpi(segs{end},'secondary')
        pathStr = strjoin(segs(1:end-1), '.');
        sel.kind = 'secondary';
    elseif numel(segs)>1 && strcmpi(segs{end},'heuristic')
        pathStr = strjoin(segs(1:end-1), '.');
        sel.heuristic = true;
    else
        pathStr = seedStr;
    end
end

% ---------------------------------------------------------------------------
function peaks = resolveSeedPeaks(vesselV, pathStr, idStr)
    % Resolve a seed path to a fitPatchVessels.m .fits array, dropping curated exclusions. The '.main'
    % fallback is applied HERE (not in parseSeedPath) precisely so a genuine field called 'main' wins:
    % parseSeedPath only strips it when the full path fails to resolve.
    S = getNestedField(vesselV, pathStr);
    if isempty(S) || ~isstruct(S)
        error('fitVesselPatchTimeSeries:seedNotFound', ...
            ['vessel %s: seed path ''%s'' resolves to nothing on this vessel -- run fitPatchVessels.m ' ...
             'first, or check the path (it is relative to the VESSEL ROOT, not to fld).'], idStr, pathStr);
    end
    % TWO SOURCE SHAPES, told apart structurally (2026-09-10):
    %   - a fitVesselPatchTimeSeries.m result carries .model.gaussian(k).params, at whatever granularity
    %     that fit used -- read from THERE, so a per-run/per-frame result (which has no .fits at all,
    %     see file header .fits COMPATIBILITY VIEW) can seed the next fit, keeping its granularity
    %     (collapseSeedToMode does the mapping onto the target's own modes);
    %   - a fitPatchVessels.m detector result carries only .fits, every value scalar.
    if isfield(S,'model') && isstruct(S.model) && isfield(S.model,'gaussian') && ~isempty(S.model.gaussian)
        peaks = peaksFromModel(S, pathStr, idStr);
    else
        assert(isfield(S,'fits'), 'fitVesselPatchTimeSeries:seedNotAFit', ...
            ['vessel %s: seed path ''%s'' has neither .model.gaussian nor .fits -- it must name a ' ...
             'fitVesselPatchTimeSeries.m or fitPatchVessels.m result (e.g. ''tsIm.gaussAnat'', ' ...
             '''tsIm.patchVesselsSequential'').'], idStr, pathStr);
        peaks = S.fits;
    end
    % Same contract fitPatchVessels.m's own resolveSeedFits honours: a peak the user (or an auto
    % criterion) already judged spurious is not refit as if it were a real vessel. Curation lives on
    % .fits (fitPatchVesselsDiag.m), by array position, so it is read from there whichever shape the
    % values came from. isfield-guarded -- a result never run through fitPatchVesselsDiag.m has no
    % .excluded field at all.
    if isfield(S,'fits') && isfield(S.fits,'excluded')
        assert(numel(S.fits) == numel(peaks), 'fitVesselPatchTimeSeries:seedFitsModelMismatch', ...
            'vessel %s: seed path ''%s'' has %d .fits entries but %d .model.gaussian entries.', ...
            idStr, pathStr, numel(S.fits), numel(peaks));
        drop = arrayfun(@(p) ~isempty(p.excluded) && p.excluded, S.fits);
        peaks = peaks(~drop);
    end
    assert(~isempty(peaks), 'fitVesselPatchTimeSeries:seedAllExcluded', ...
        'vessel %s: seed path ''%s'' has no INCLUDED peak left.', idStr, pathStr);
end

% ---------------------------------------------------------------------------
function peaks = peaksFromModel(S, pathStr, idStr)
    % One seed peak per .model.gaussian(k), values taken from .params VERBATIM (per-run cell or bare,
    % see collapseSeedToMode for what each means), b from the one background component, .peakNum from
    % the entry. Array position 1 == the main peak, the same invariant .fits carries.
    g = S.model.gaussian;
    bVal = [];
    if isfield(S.model,'background') && ~isempty(S.model.background)
        bVal = S.model.background(1).params.b;
    end
    peaks = struct('a',{}, 'x0',{}, 'y0',{}, 'radius',{}, 'aspectRatio',{}, 'theta',{}, 'b',{}, 'peakNum',{});
    for k = 1:numel(g)
        assert(isfield(g(k),'params') && isstruct(g(k).params), 'fitVesselPatchTimeSeries:seedNotAFit', ...
            'vessel %s: seed path ''%s'' .model.gaussian(%d) has no .params.', idStr, pathStr, k);
        p = g(k).params;
        pkNum = k;
        if isfield(g(k),'peakNum') && ~isempty(g(k).peakNum); pkNum = g(k).peakNum; end
        % Every value braced: struct() would otherwise expand a {1 x nRun} cell into nRun structs.
        peaks(k) = struct('a',{p.a}, 'x0',{p.x0}, 'y0',{p.y0}, 'radius',{p.radius}, ...
                          'aspectRatio',{p.aspectRatio}, 'theta',{p.theta}, 'b',{bVal}, 'peakNum',{pkNum});
    end
end

% ---------------------------------------------------------------------------
function peaks = resolveHeuristicSeedPeak(vesselV, pathStr, idStr, includeMaskArg, excludeMaskArg)
    % Sibling to resolveSeedPeaks, for the 'path.heuristic'/'path(heuristic)' spelling: pathStr must
    % name a real ts-track (e.g. 'tsIm') carrying its own .im/.vSize DIRECTLY -- NOT a fitPatchVessels.m
    % .fits-shaped result (resolveSeedPeaks' own job). Builds ONE ad hoc seed peak straight from that
    % path's own time-averaged image, via the legacy single-peak start/bound heuristic
    % (gaussianFitDefaultsRadiusAspect.m), seeded at the patch center -- for a caller with no prior
    % fitPatchVessels.m result to seed from at all (the very first "get me a starting anatomical fit"
    % call in the pipeline).
    %
    % includeMaskArg/excludeMaskArg -- the CALLER's own opts.includeMask/.excludeMask (2026-09-09, a
    % real bug found investigating a reported off-FOV convergence, drawVesselPatchPeaks patch idx=27:
    % "a clear peak present at the patch center... but it should have been properly excluded by
    % opts.excludeMask='excludeManualSinus'"). BEFORE this fix, this function computed its own
    % amplitude/background estimate (b0 = mean of the WHOLE smoothed image; a0 = local value at the
    % seed pixel MINUS b0) from the image UNCONDITIONALLY, completely ignoring opts.excludeMask -- even
    % though the REAL fit that follows correctly excludes those voxels from its own cost function
    % (confirmed by direct code read of the caller's own validMask/resolveIncludeMaskLocal/
    % unionMasksLocal). A bright excluded region (a sinus) inflates b0, which can make a0 badly
    % underestimated (even clamped near zero) -- a heuristic that no longer looks like "there is a real
    % peak here". Combined with a caller-widened seedBounds.x0y0 (e.g. Inf, removing
    % gaussianFitDefaults.m's own built-in patch-extent x0/y0 bound -- see that file's own d.lower.x0/
    % d.upper.x0), a badly-biased starting amplitude/shape can let the optimizer converge somewhere far
    % from the real peak, including completely outside the displayed patch, since nothing bounds
    % position once that safety net is removed. Excluded (and non-included) pixels are now replaced
    % with the INCLUDED region's own median before ANY heuristic statistic is computed from the image --
    % a neutral, locally-representative fill value that does not bias b0/a0 toward content the real fit
    % will never see either. Position (seedXY, always the patch center) and the shape/radius starting
    % guess (a FIXED voxel-size-derived constant, not derived from image content at all) are unaffected
    % either way.
    base = getNestedField(vesselV, pathStr);
    % 2026-09-10, Seb's own ask: mixing '.heuristic' with a dot-path to an EXISTING fit must error
    % loudly, not silently misbehave. Verified this is ALREADY the case, with no code change needed --
    % fitVesselPatchTimeSeries.m's own output never stores an .im field at all (only a genuine raw
    % ts-track does, e.g. vessel.tsIm), so a path that is actually an existing fit result fails this
    % SAME assert on .im alone, already with a clear, specific message ("...not a fitPatchVessels.m
    % result"). A candidate ~isfield(base,'fits') addition here would be unreachable dead code: if base
    % lacked .im it already threw above, and a genuine raw ts-track never ALSO carries .fits (that
    % lives one level deeper, e.g. vessel.tsIm.gaussAnat.fits, never on vessel.tsIm itself).
    assert(isstruct(base) && isfield(base,'im') && ~isempty(base.im) && isfield(base,'vSize'), ...
        'fitVesselPatchTimeSeries:heuristicSeedNoIm', ...
        ['vessel %s: seed path ''%s.heuristic''/''%s(heuristic)'' requires ''%s'' to be a real ts-track ' ...
         '(a struct with a non-empty .im and a .vSize), not a fitPatchVessels.m result -- it is relative ' ...
         'to the VESSEL ROOT, not to fld.'], idStr, pathStr, pathStr, pathStr);

    % Same time-averaging reduction fitOneVessel itself performs for opts.timeAvg (see that function,
    % ~line 403-417) -- matched exactly, then squeezed to a plain 2-D array before the heuristic call.
    im = base.im; if ~iscell(im); im = {im}; end
    acc = cat(4, im{:});
    im2d = squeeze(mean(acc(:,:,1,:), 4, 'omitnan'));

    % Resolved FRESH from base.vSize rather than assumed to match the caller's own possibly-different
    % fld -- see fitOneVessel's identical derivation from metaBase.vSize (~line 398-399).
    voxSz2 = [base.vSize(1) base.vSize(min(2,numel(base.vSize)))];

    % Mask out excluded/non-included pixels (see includeMaskArg/excludeMaskArg note above) -- resolved
    % against pathStr itself (where these labels actually live), the SAME resolveIncludeMaskLocal/
    % unionMasksLocal this file's own real fit-region mask uses, so "excluded from the heuristic" and
    % "excluded from the fit" can never silently disagree.
    validMask = resolveIncludeMaskLocal(vesselV, pathStr, includeMaskArg, size(im2d)) ...
                & ~unionMasksLocal(vesselV, pathStr, excludeMaskArg, size(im2d));
    if any(~validMask(:)) && any(validMask(:))
        im2d(~validMask) = median(im2d(validMask), 'omitnan');
    end

    % No known peak location yet -- these vessel patches are already COM-aligned upstream, so dead
    % center is the reasonable default. gaussianFitDefaults.m's own seedXY convention is [xPix,yPix]
    % ([col,row]); this choice makes the returned x0/y0 come out as exactly 0 ("centered") by
    % construction, matching gaussianFitDefaults.m's own ctrCol/ctrRow formula exactly.
    [ny,nx] = size(im2d);   % matches gaussianFitDefaults.m's own [ny,nx]=size(im2d) exactly
    seedXY = [(nx+1)/2, (ny+1)/2];

    d = gaussianFitDefaultsRadiusAspect(im2d, voxSz2, seedXY);
    % Only d.start feeds the peak's values -- bounds for this component still go through mkComponent's
    % own generic .seedBounds machinery unchanged, exactly like every other seed source.
    peaks = struct('a',d.start.a, 'x0',d.start.x0, 'y0',d.start.y0, 'radius',d.start.radius, ...
                   'aspectRatio',d.start.aspectRatio, 'theta',d.start.theta, 'b',d.start.b, 'peakNum',1);
end

% ---------------------------------------------------------------------------
function res = assembleResult(fit, meta, opts, runLens, T, validMask, doSplit, nKnotGroups, repOut)
    % Reshape the solver's flat component list back into opts.model's own shape, split every
    % per-parameter field per run, and add the flat main-peak alias.
    res = struct();
    res.model = struct();
    mainIdx = [];
    bgIdx = [];
    for c = 1:numel(meta)
        tn = meta(c).type;
        entry = struct();
        entry.params = splitStructPerRun(fit.components(c).params, runLens, T, doSplit, nKnotGroups);
        entry.start  = splitStructPerRun(fit.components(c).start,  runLens, T, doSplit, nKnotGroups);
        entry.lower  = splitStructPerRun(fit.components(c).lower,  runLens, T, doSplit, nKnotGroups);
        entry.upper  = splitStructPerRun(fit.components(c).upper,  runLens, T, doSplit, nKnotGroups);
        entry.mode   = fit.components(c).mode;
        entry.fixed  = fit.components(c).fixed;
        entry.label  = fit.components(c).label;
        entry.peakNum = meta(c).peakNum;
        entry.seed    = meta(c).seed;
        entry.priorFit = meta(c).priorFit;
        % The caller-facing bound SPEC this component's .lower/.upper were derived from (defaults
        % filled in), so the 2-arg call can hand it back verbatim rather than inverting the bounds.
        entry.seedBounds = meta(c).seedBounds;
        % .runPolyCoef/.runPolyDeg (mode='perRunPoly' parameters only) -- ALREADY a {1 x nRunGroups}
        % cell / a resolved degree list respectively (fitComponentFit.m's own output), so passed through
        % verbatim rather than through splitStructPerRun above: that function's three/four-way length
        % convention is for a PER-FRAME-shaped value, and this is neither -- it is already shaped per
        % run (or, on the irf-path, one degenerate "run" spanning the whole knot axis). Present (as an
        % empty struct, no sub-fields) even when no parameter of this component used the mode.
        if ~isempty(fieldnames(fit.components(c).runPolyCoef))
            entry.runPolyCoef = fit.components(c).runPolyCoef;
            entry.runPolyDeg  = fit.components(c).runPolyDeg;
        end
        % Replicate arrays carry an extra LEADING replicate dimension, so they are deliberately left
        % UNSPLIT by run -- [nRep x pLen] means the same thing on either path, and splitting a
        % replicate axis by run would conflate two different dimensions. Same call
        % fitVesselTimeSeries.m made for its own .boot/.perm.
        % .bootRep/.nullRep hold either the raw [nRep x pLen] values or a dataCacheDescriptor.m
        % POINTER to them, per opts.<kind>.offloadToCache -- read either form with dataCacheLoad.m,
        % which passes materialised data straight through. PER COMPONENT the offloaded form is one
        % pointer to the whole file (the file holds every component), so resolving it once gives
        % .val{c} for all.
        % NO .bootSE/.nullSE, on either path -- see the file header's REPLICATES ONLY. Call
        % replicateSE.m (or your own summary) on .bootRep/.nullRep at the point of use.
        if isfield(repOut,'boot')
            if isDataCacheDescriptor(repOut.bootStore)
                entry.bootRep = repOut.bootStore;
            else
                entry.bootRep = repOut.bootStore{c};
            end
        end
        if isfield(repOut,'null')
            if isDataCacheDescriptor(repOut.nullStore)
                entry.nullRep = repOut.nullStore;
            else
                entry.nullRep = repOut.nullStore{c};
            end
        end
        if isfield(res.model, tn); res.model.(tn)(end+1) = entry; else; res.model.(tn) = entry; end
        if meta(c).isMainGaussian; mainIdx = c; end
        if strcmp(tn,'background'); bgIdx = c; end
    end

    % FLAT ALIAS of the main gaussian (+ background b) -- see file header OUTPUT for the caveat that it
    % describes the main peak only.
    if ~isempty(mainIdx)
        fn = fieldnames(fit.components(mainIdx).params);
        for i = 1:numel(fn)
            res.(fn{i}) = splitPerRun(fit.components(mainIdx).params.(fn{i}), runLens, T, doSplit, nKnotGroups);
        end
        res.start = splitStructPerRun(fit.components(mainIdx).start, runLens, T, doSplit, nKnotGroups);
        res.lower = splitStructPerRun(fit.components(mainIdx).lower, runLens, T, doSplit, nKnotGroups);
        res.upper = splitStructPerRun(fit.components(mainIdx).upper, runLens, T, doSplit, nKnotGroups);
    end
    if ~isempty(bgIdx)
        res.b = splitPerRun(fit.components(bgIdx).params.b, runLens, T, doSplit, nKnotGroups);
        if isfield(res,'start')
            res.start.b = splitPerRun(fit.components(bgIdx).start.b, runLens, T, doSplit, nKnotGroups);
            res.lower.b = splitPerRun(fit.components(bgIdx).lower.b, runLens, T, doSplit, nKnotGroups);
            res.upper.b = splitPerRun(fit.components(bgIdx).upper.b, runLens, T, doSplit, nKnotGroups);
        end
    end
    % .runPolyCoef/.runPolyDeg flat alias -- mirrors the main gaussian's own (plus background's b),
    % same convention as every other flat-alias field above; passed through verbatim (see the identical
    % per-component note above for why this skips splitStructPerRun).
    if ~isempty(mainIdx) && ~isempty(fieldnames(fit.components(mainIdx).runPolyCoef))
        res.runPolyCoef = fit.components(mainIdx).runPolyCoef;
        res.runPolyDeg  = fit.components(mainIdx).runPolyDeg;
    end
    if ~isempty(bgIdx) && ~isempty(fieldnames(fit.components(bgIdx).runPolyCoef))
        if ~isfield(res,'runPolyCoef'); res.runPolyCoef = struct(); res.runPolyDeg = struct(); end
        res.runPolyCoef.b = fit.components(bgIdx).runPolyCoef.b;
        res.runPolyDeg.b  = fit.components(bgIdx).runPolyDeg.b;
    end
    % .paramSpec completes the flat alias. fitVessel.m-shaped results carry one, and consumers use it
    % to reason about the fit WITHOUT reading the values -- fitVesselTimeSeriesDiag.m's own eligibility
    % check reads .paramSpec.mode to decide whether anything varies per frame and is therefore worth
    % animating. Without it that check rejects this function's output outright, which is how a movie of
    % a genuinely per-frame multi-peak fit would otherwise be refused. Mirrors the MAIN gaussian's own
    % modes/fixed flags plus b from the background, i.e. exactly the parameters the flat alias reports.
    if ~isempty(mainIdx)
        res.paramSpec = struct('mode', fit.components(mainIdx).mode, ...
                               'fixed', fit.components(mainIdx).fixed);
        if ~isempty(bgIdx)
            res.paramSpec.mode.b  = fit.components(bgIdx).mode.b;
            res.paramSpec.fixed.b = fit.components(bgIdx).fixed.b;
        end
    end

    % Flat alias for the replicates too, mirroring the main gaussian (+ background b), same convention
    % as the point-estimate alias above. RAW REPLICATES ONLY -- no summaries, see the file header.
    if ~isempty(mainIdx) && isfield(repOut,'boot')
        if isDataCacheDescriptor(repOut.bootStore)
            res.bootRep = repOut.bootStore;
        else
            res.bootRep = repOut.bootStore{mainIdx};
            if ~isempty(bgIdx); res.bootRep.b = repOut.bootStore{bgIdx}.b; end
        end
        res.nBootRep = repOut.boot.n;
        res.bootSource = repOut.boot.source;   % 'cache' | 'vessel' | 'computed' | 'cache+computed'
        res.bootKey = repOut.boot.key;         % validates an already-in-vessel reuse on a later call
    end
    if ~isempty(mainIdx) && isfield(repOut,'null')
        if isDataCacheDescriptor(repOut.nullStore)
            res.nullRep = repOut.nullStore;
        else
            res.nullRep = repOut.nullStore{mainIdx};
            if ~isempty(bgIdx); res.nullRep.b = repOut.nullStore{bgIdx}.b; end
        end
        res.nNullRep = repOut.null.n;
        res.nullSource = repOut.null.source;
        res.nullKey = repOut.null.key;
    end

    res.resnorm           = fit.resnorm;
    res.resnormPerFrame   = splitPerRun(fit.resnormPerFrame, runLens, T, doSplit, nKnotGroups);
    res.exitFlag          = fit.exitFlag;
    res.iterations        = fit.iterations;
    res.maxIterations     = fit.maxIterations;
    res.functionTolerance = fit.functionTolerance;
    res.nFree             = fit.nFree;
    res.method            = 'gaussian';
    res.fitMode           = 'componentJoint';

    % IDENTIFIABILITY INPUTS, passed through RAW -- see the file header's REPLICATES ONLY, which is
    % the same rule: the jacobian and bound multipliers are solver output, so they travel; the
    % diagnosis is fitIdentifiability.m's job, at the point of use. This is what makes
    % "is this fitted number a measurement?" answerable from the fit itself rather than from a size
    % heuristic -- and it is the only route that can see CROSS-COMPONENT degeneracy, since the
    % jacobian spans every component's free parameters at once.
    %   d = fitIdentifiability(vessel.tsIm.myFit.solve);
    % Kept in a .solve sub-struct rather than flattened, because .jacobian is [m x nFree] and would
    % otherwise sit in the same namespace as the fitted values, inviting exactly the confusion
    % between a parameter and a diagnostic that this codebase has already paid for once.
    res.solve = struct('jacobian', fit.jacobian, 'lambda', fit.lambda, ...
                       'freeLabel', {fit.freeLabel}, 'freeComp', fit.freeComp, ...
                       'freeParam', {fit.freeParam}, 'freeGroup', fit.freeGroup, ...
                       'freeStart', fit.freeStart, 'freeLower', fit.freeLower, ...
                       'freeUpper', fit.freeUpper, 'freeFit', fit.freeFit, ...
                       'resnorm', fit.resnorm, 'nFree', fit.nFree, 'components', fit.components);
    res.timeAvg           = opts.timeAvg;
    res.nIso              = opts.nIsochromatPerVoxDim;
    res.nIsochromatPerVoxDim = opts.nIsochromatPerVoxDim;
    res.includeMask       = opts.includeMask;
    res.excludeMask       = opts.excludeMask;
    res.nValidVox         = nnz(validMask);
    % .validMask (2026-09-08, Seb's own ask): the ACTUAL resolved logical array fitComponentFit.m was
    % given -- .includeMask/.excludeMask above are only the caller's own INPUT label list(s), not what
    % was actually fit against once they're combined (validMask = included & ~excluded, via this
    % function's own resolveIncludeMaskLocal/unionMasksLocal). Kept alongside, not instead of, the
    % input lists -- both are useful provenance and neither can be derived from the other without the
    % vessel's own .rois still being available.
    res.validMask         = validMask;

    % .fits COMPATIBILITY VIEW -- timeAvg, OR any other combination that still resolves every
    % parameter to a genuine scalar (2026-09-10, Seb's own ask: "even for perVessel fits we want
    % either a fit on the time average or fitting a single parameter on time-resolved data" -- a
    % mode='perVessel' fit against raw per-frame data, opts.timeAvg=false, is jointly fit across
    % every frame and comes out just as scalar as a timeAvg=true one; only the DATA differs
    % (pre-averaged vs. not), never the parameter count). See attachFitsView's own comment for why
    % scalar-ness, not opts.timeAvg itself, is the real precondition. isAllScalarComponents checks
    % fit.components DIRECTLY (the raw, pre-run-split solver output attachFitsView itself reads),
    % since that is the one place "scalar" is unambiguous -- res.model's own fields are already
    % {1 x nRun} cells by this point (splitStructPerRun, above), always, even for a perVessel
    % param, so checking res itself would find every field "non-scalar" regardless of mode.
    if opts.timeAvg || isAllScalarComponents(fit)
        res = attachFitsView(res, fit, meta, opts);
    end
end

% ---------------------------------------------------------------------------
function tf = isAllScalarComponents(fit)
    tf = true;
    for c = 1:numel(fit.components)
        pn = fieldnames(fit.components(c).params);
        for k = 1:numel(pn)
            if numel(fit.components(c).params.(pn{k})) ~= 1
                tf = false; return;
            end
        end
    end
end

% ---------------------------------------------------------------------------
function res = attachFitsView(res, fit, meta, opts)
    % The fitPatchVessels.m-shaped .fits struct array, so every existing consumer of a
    % 'simultaneous'-mode result -- fitPatchVesselsDiag.m, makeRoiFromVesselFit.m, alignVessel2.m,
    % fitVesselTimeSeries.m's Shape-A seed -- reads this function's output with NO change at all.
    %
    % SCALAR RESULTS ONLY, deliberately (2026-09-10: broadened from "timeAvg only" -- see this
    % function's own caller, isAllScalarComponents). A .fits entry is a struct of SCALAR parameter
    % values; that is only well defined when every parameter resolved to exactly one value, which
    % happens both when timeAvg pre-collapses the data AND when every parameter is mode='perVessel'
    % on raw per-frame data (opts.timeAvg=false) -- the solver still estimates one joint value either
    % way, only the DATA differs. A genuinely per-run/per-frame fit would mean either silently
    % reporting one frame as if it were the answer, or putting arrays where every consumer expects
    % scalars. Neither is worth a compatibility shim, so such a result simply has no .fits and the
    % .model view is the only one.
    %
    % b IS REPLICATED ONTO EVERY PEAK. In this framework the background is one separate component, but
    % fitPatchVessels.m's own contract gives each peak its own .b, and drawGaussianPeakContour.m reads
    % peak.b DIRECTLY when it evaluates the model for a contour. So the shared value is copied onto
    % each entry -- a derived alias of the one background, not K independent baselines. Same thing
    % fitMultiVessel.m already did with its own shared b.
    gIdx = find(strcmp({meta.type},'gaussian'));
    bIdx = find(strcmp({meta.type},'background'), 1);
    if isempty(gIdx); return; end
    bVal = []; bStart = []; bLower = []; bUpper = [];
    if ~isempty(bIdx)
        bVal   = fit.components(bIdx).params.b;
        bStart = fit.components(bIdx).start.b;
        bLower = fit.components(bIdx).lower.b;
        bUpper = fit.components(bIdx).upper.b;
    end

    % Built by assigning fields on an indexed element rather than growing a pre-templated struct array:
    % that gives every entry the same field order by construction, with no orderfields() call needed
    % (and no dependence on an empty struct array being a valid orderfields template).
    fits = struct([]);
    for i = 1:numel(gIdx)
        c = gIdx(i);
        pn = fieldnames(fit.components(c).params);
        for k = 1:numel(pn); fits(i).(pn{k}) = fit.components(c).params.(pn{k}); end
        fits(i).b = bVal;
        fits(i).peakNum = meta(c).peakNum;
        s = fit.components(c).start; s.b = bStart; fits(i).start = s;
        l = fit.components(c).lower; l.b = bLower; fits(i).lower = l;
        u = fit.components(c).upper; u.b = bUpper; fits(i).upper = u;
        % Provenance fitPatchVesselsDiag.m ASSERTS on (.includeMask/.excludeMask/.nIso/.lower/.upper) --
        % it recomputes validMask from these rather than trusting a caller to pass the same masks again.
        fits(i).includeMask = opts.includeMask;
        fits(i).excludeMask = opts.excludeMask;
        fits(i).nIso = opts.nIsochromatPerVoxDim;
        fits(i).priorFit = meta(c).priorFit;
    end
    res.fits = fits;
    % Container-level fields matching fitPatchVessels.m's own 'simultaneous' output shape.
    if opts.timeAvg; scalarTag = 'timeAvg'; else; scalarTag = 'scalarPerVessel'; end
    res.fitMode = sprintf('componentJoint (%s, seeded from %s)', scalarTag, meta(gIdx(1)).seed);
    res.K = numel(gIdx);
    if ~isempty(bIdx); res.bShared = bVal; end
end

% ---------------------------------------------------------------------------
function [rep, store] = runReplicates(replicateIm, rawImC, components, validMask, voxSz2, fcOpts, ...
                                      opts, irfCore, isIrf, idStr, kind, v, nV, fldTag, prior)
    % One replicate set, under the DATA-CACHE POLICY (see the file header's own REPLICATE CACHE
    % section, and .bass/pkm/patterns/data-cache-offload-recompute.md for the general convention).
    %
    % TWO ORTHOGONAL AXES, both per-kind (opts.boot.* / opts.null.*):
    %   .offloadToCache (default TRUE)  -- MEMORY. true leaves the raw replicates on disk and puts a
    %                    dataCacheDescriptor.m POINTER in the result; false materialises them.
    %   .reComputeCache (default FALSE) -- CPU. false reuses what already exists; true recomputes
    %                    regardless, overwriting the cache file when offloading too.
    % So: WRITE iff offloadToCache; READ iff ~reComputeCache and something exists; POINTER iff
    % offloadToCache. The default pair -- offload on, recompute off -- is the case that should almost
    % always be wanted: never recompute, never hold in RAM.
    %
    % REPLICATES ONLY -- see the file header. This returns raw replicate fits and cache metadata, and
    % nothing derived from the values, on either offload path. replicateSE.m is the downstream helper.
    %
    % `store` (2nd output) is what the caller writes into the result: either the raw values or a
    % pointer to them, per offloadToCache.
    spec = opts.(kind);
    nReq = spec.nRep;

    % --- resolve the replicate SOURCE and its identity ----------------------
    % srcTag identifies the source for the cache key. A file-backed stack by path+size+mtime (cheap,
    % and it invalidates if fitIRF.m rewrites it); an in-memory one by shape+checksum; the ts-path null
    % has no stack at all, so the raw image's shape/checksum stands in (dataSeed is keyed separately).
    srcTagIn = '';
    if ~isempty(replicateIm)
        if isReplicateImDescriptor(replicateIm)
            d = dir(replicateIm.file);
            if isempty(d)
                srcTagIn = sprintf('file:%s:missing', replicateIm.file);
            else
                srcTagIn = sprintf('file:%s:%d:%.6f', replicateIm.file, d.bytes, d.datenum);
            end
        end
    end

    cacheFile = fullfile(opts.cacheDir, sprintf('%s_%s_%s.mat', ...
        matlab.lang.makeValidName(idStr), matlab.lang.makeValidName(fldTag), kind));

    % --- POLICY, resolved BEFORE any expensive work -------------------------
    % Three possible sources, in cost order: already in the vessel, in the cache file, or computed.
    % ~reComputeCache honours the first two; reComputeCache=true skips straight to computing.
    haveInVessel = false; inVesselVal = {}; nInVessel = 0;

    % The key needs the resolved stack for its shape/checksum, so a file-backed stack that is going to
    % be REUSED from cache never has to be read at all -- which is the "don't even load it" case.
    if isempty(replicateIm)
        nRep = nReq;
        srcTag = sprintf('ts:%s:%.12g', mat2str(cellfun(@numel, rawImC)), ...
                         sum(cellfun(@(x) sum(double(x(:))), rawImC)));
    elseif ~isempty(srcTagIn)
        srcTag = srcTagIn;              % file-backed: identity WITHOUT reading the file
        nRep = nReq;
    else
        srcTag = sprintf('mem:%s:%.12g', mat2str(size(replicateIm)), sum(double(replicateIm(:))));
        nRep = nReq;
    end
    key = cacheKeyString(components, validMask, voxSz2, fcOpts, opts, kind, srcTag);
    % In-vessel reuse is validated against the SAME key as the cache -- otherwise re-running with a
    % changed bound would happily reuse the previous call's replicates out of the struct, which is the
    % identical stale-reuse hazard the cache key exists to prevent.
    if ~spec.reComputeCache && ~isempty(prior) && isfield(prior,'key') && strcmp(prior.key, key) ...
            && numel(prior.val) == numel(components)
        inVesselVal = prior.val; nInVessel = prior.n; haveInVessel = true;
    end

    nCached = 0; cacheKeyOk = false;
    if ~spec.reComputeCache && ~haveInVessel
        [nCached, cacheKeyOk] = peekReplicateCache(cacheFile, key);
    end

    % --- CASE 1: reuse, nothing to compute ----------------------------------
    if haveInVessel && nInVessel >= nRep
        if opts.progress
            fprintf('  fitVesselPatchTimeSeries: vessel %d/%d (%s): %s -- %d replicate(s) already in vessel\n', ...
                    v, nV, idStr, kind, nRep);
        end
        [rep, store] = finaliseReplicates(inVesselVal, nRep, numel(components), key, ...
                                          spec, cacheFile, isIrf, 'vessel');
        return
    end
    if cacheKeyOk && nCached == nRep && spec.offloadToCache
        % THE DEFAULT PATH, and the whole point of the memory axis: the file already holds exactly what
        % was asked for, so it is never read, never materialised, never summarised -- only pointed at.
        % NOTHING is loaded here. Error bars come from replicateSE.m on the pointer, when and if
        % someone actually wants them.
        if opts.progress
            fprintf('  fitVesselPatchTimeSeries: vessel %d/%d (%s): %s -- %d replicate(s) in cache, NOT loaded (pointer only)\n', ...
                    v, nV, idStr, kind, nRep);
        end
        rep = struct('n', nRep, 'isIrf', isIrf, 'source', 'cache', 'key', key);
        store = dataCacheDescriptor(cacheFile);
        return
    end
    if cacheKeyOk && nCached >= nRep
        % Either not offloading (so it has to be materialised anyway), or more is stored than asked
        % for, which needs a trim -- and a spread over a superset is not the spread over the subset,
        % so any summary has to come from the trimmed values.
        if opts.progress
            fprintf('  fitVesselPatchTimeSeries: vessel %d/%d (%s): %s -- %d replicate(s) loaded from cache (%d stored)\n', ...
                    v, nV, idStr, kind, nRep, nCached);
        end
        full = dataCacheLoad(dataCacheDescriptor(cacheFile));
        [rep, store] = finaliseReplicates(full.val, nRep, numel(components), key, ...
                                          spec, cacheFile, isIrf, 'cache');
        return
    end

    % --- CASE 2: compute (possibly only a tail) -----------------------------
    % PARTIAL RESUME. Not in the specified truth table, added because the binary reading would return
    % 512 replicates to a caller asking for 1024. Replicate r is index-stable -- slice r of
    % imBoot/imPerm is fixed, and the ts-path null's substream r is determined by opts.dataSeed -- so
    % the already-computed ones never need redoing. Steady state (a full cache) still never loads.
    baseVal = {}; nHave = 0; fromWhere = 'computed';
    if haveInVessel
        baseVal = inVesselVal; nHave = nInVessel; fromWhere = 'vessel+computed';
    elseif cacheKeyOk && nCached > 0
        full = dataCacheLoad(dataCacheDescriptor(cacheFile));
        baseVal = full.val; nHave = nCached; fromWhere = 'cache+computed';
    end

    % Now the stack genuinely has to be in memory. A file-backed one is resolved into a LOCAL and
    % dropped when this function returns -- the same discipline replicateImData.m documents.
    if ~isempty(replicateIm)
        replicateIm = replicateImData(replicateIm);
        nAvail = size(replicateIm, 5);
        assert(nRep <= nAvail, 'fitVesselPatchTimeSeries:tooManyReplicatesRequested', ...
            ['vessel %s: opts.%s.nRep=%d but the replicate stack only has %d -- rebuild it with a ' ...
             'larger opts.nBoot/opts.nPerm in fitIRF.m, or ask for fewer.'], idStr, kind, nRep, nAvail);
    end
    assert(nRep >= 2, 'fitVesselPatchTimeSeries:tooFewReplicates', ...
        'vessel %s: %s asked for %d replicate(s) -- a spread across replicates needs at least 2.', ...
        idStr, kind, nRep);

    streams = {};
    if isempty(replicateIm)
        streams = RandStream.create('mrg32k3a', 'NumStreams', nRep, 'Seed', opts.dataSeed, 'CellOutput', true);
    end
    addB0 = opts.addBase0;
    r0 = nHave + 1;
    fits = cell(1, nRep);
    if opts.progress
        fprintf('  fitVesselPatchTimeSeries: vessel %d/%d (%s): %s -- %d replicate(s) in parallel (%d reused) ...\n', ...
                v, nV, idStr, kind, nRep-nHave, nHave);
    end
    t0 = tic;
    % TWO SEPARATE LOOPS, not one with a branch inside. A parfor variable gets ONE role: writing
    % `if isempty(replicateIm)` inside the body uses it WHOLE while `replicateIm(:,:,:,:,r)` uses it
    % SLICED, and MATLAB cannot classify it as both -- a hard parfor error, not a silent fallback (hit
    % while writing this). Splitting the branch out also keeps the slicing PROVABLE, which is the
    % point: index that array with anything but the loop variable itself and parfor reverts to
    % broadcasting a full copy to every worker -- the OOM this file family has already root-caused
    % (819MB vs 205MB of worker traffic, killing a real nRep=1024 x 32-worker batch).
    if isempty(replicateIm)
        parfor r = r0:nRep
            imR = phaseRandomizeRuns(rawImC, streams{r});
            fits{r} = fitComponentFit(imR, validMask, voxSz2, components, fcOpts);
        end
    else
        parfor r = r0:nRep
            imR = replicateIm(:,:,:,:,r);      % DIRECT slice -- see above
            imR = imR(:,:,1,:);
            % A replicate arrives with base0 already baked in by fitIRF.m, exactly like the point
            % estimate, so addBase0=false must subtract it back out here too -- reusing the POINT
            % ESTIMATE's own base0 (fitIRF.m carries no per-replicate baseline coefficients).
            if ~addB0; imR = irfAddBase0(irfCore, imR, [], -1); end
            fits{r} = fitComponentFit(imR, validMask, voxSz2, components, fcOpts);
        end
    end
    if opts.progress
        fprintf('  fitVesselPatchTimeSeries: vessel %d/%d (%s): %s done in %.1f s\n', ...
                v, nV, idStr, kind, toc(t0));
    end

    newVal = stackReplicates(fits(r0:nRep), numel(components));
    val = mergeReplicateVal(baseVal, newVal, numel(components));
    [rep, store] = finaliseReplicates(val, nRep, numel(components), key, ...
                                      spec, cacheFile, isIrf, fromWhere);
end

% ---------------------------------------------------------------------------
function prior = priorReplicates(priorRes, kind)
    % Replicates ALREADY sitting at the output field from an earlier call -- the "odd chance" case.
    % Returned only when they are genuinely materialised AND carry the spec key that produced them, so
    % the caller can verify them exactly as it verifies a cache file. A POINTER left by a previous
    % offloading call is deliberately NOT treated as in-vessel data: the cache path handles that, and
    % it handles it without loading anything.
    prior = [];
    if isempty(priorRes) || ~isstruct(priorRes); return; end
    fRep = [kind 'Rep']; fKey = [kind 'Key']; fN = ['n' upper(kind(1)) kind(2:end) 'Rep'];
    if ~isfield(priorRes, fRep) || ~isfield(priorRes, fKey) || ~isfield(priorRes, fN); return; end
    if isDataCacheDescriptor(priorRes.(fRep)); return; end
    % The flat alias only carries the MAIN component, so it cannot reconstruct every component's
    % values -- the per-component .model view is the authoritative copy. Recovered from there.
    if ~isfield(priorRes,'model') || ~isstruct(priorRes.model); return; end
    val = {};
    tn = fieldnames(priorRes.model);
    for i = 1:numel(tn)
        entries = priorRes.model.(tn{i});
        for e = 1:numel(entries)
            if ~isfield(entries(e), fRep) || isDataCacheDescriptor(entries(e).(fRep)); return; end
            val{end+1} = entries(e).(fRep); %#ok<AGROW>
        end
    end
    prior = struct('val', {val}, 'n', priorRes.(fN), 'key', priorRes.(fKey));
end

% ---------------------------------------------------------------------------
function [rep, store] = finaliseReplicates(val, nRep, nComp, key, spec, cacheFile, isIrf, source)
    % Apply the MEMORY axis: write-and-point (offloadToCache=true) or keep resident (false).
    %
    % NOTHING DERIVED FROM THE REPLICATES IS COMPUTED HERE -- see the file header's REPLICATES ONLY.
    % No .se in either case. `n`/`key`/`source` are cache METADATA (they describe the replicate SET,
    % not its values) and are the only things that travel alongside.
    rep = finishReplicates(val, nRep, nComp);
    rep.n = nRep; rep.isIrf = isIrf; rep.source = source; rep.key = key;
    if spec.offloadToCache
        saveDataCache(cacheFile, key, struct('val',{val}), nRep);
        store = dataCacheDescriptor(cacheFile);
    else
        store = rep.val;
    end
end

% ---------------------------------------------------------------------------
function key = cacheKeyString(components, validMask, voxSz2, fcOpts, opts, kind, srcTag)
    % Canonical text description of EVERYTHING that determines the replicate values. Stored verbatim in
    % the cache file and compared as a string, so a miss is readable rather than an opaque hash
    % difference. Anything omitted here is something a caller could change while silently getting the
    % previous run's replicates back -- the failure mode reuse-by-default introduces -- so the list is
    % deliberately generous. The predecessor's cache compared only the stored parameter NAME LIST,
    % which would not have caught a changed bound, mode, or mask at all.
    parts = {sprintf('kind=%s', kind), sprintf('src=%s', srcTag), ...
             sprintf('nIso=%d', fcOpts.nIsochromatPerVoxDim), ...
             sprintf('vox=%.12g,%.12g', voxSz2(1), voxSz2(2)), ...
             sprintf('addBase0=%d', opts.addBase0), ...
             sprintf('dataSeed=%d', opts.dataSeed), ...
             sprintf('runIdx=%s', mat2str(fcOpts.runIdx)), ...
             sprintf('knotIdx=%s', mat2str(fcOpts.knotIdx)), ...
             sprintf('maskIdx=%s', mat2str(find(validMask(:)).'))};
    for c = 1:numel(components)
        cc = components(c);
        mdl = cc.model; if ~ischar(mdl); mdl = mdl.name; end
        parts{end+1} = sprintf('c%d:model=%s,label=%s', c, mdl, cc.label); %#ok<AGROW>
        pn = fieldnames(cc.mode);
        for i = 1:numel(pn)
            n = pn{i};
            parts{end+1} = sprintf('c%d:%s:mode=%s,fixed=%d,start=%s,lower=%s,upper=%s', ...
                c, n, cc.mode.(n), cc.fixed.(n), mat2str(cc.start.(n),12), ...
                mat2str(cc.lower.(n),12), mat2str(cc.upper.(n),12)); %#ok<AGROW>
        end
    end
    key = strjoin(parts, '|');
end

% ---------------------------------------------------------------------------
function [nStored, keyOk] = peekReplicateCache(cacheFile, key)
    % Read ONLY the key and the stored count -- never the raw replicates, and (since 2026-09-05) no
    % stored summary either, because none is written any more. That is what makes offloadToCache's
    % "don't even load it" case real: enough comes back to decide the policy without touching the
    % bulk. Saved as separate top-level variables for exactly this reason (see saveDataCache) --
    % load() with an explicit variable list does not read the rest of the file.
    %
    % Reports only what is IN the file; the caller decides what to do with it, because the three cases
    % differ (exact count -> pointer; more than asked -> load and trim; fewer -> resume from it).
    %
    % A miss -- absent, unreadable, or key mismatch -- is never an error, just keyOk=false. Same
    % "stale/incompatible cache is a miss, not a failure" convention fitIRF.m's own cache uses.
    nStored = 0; keyOk = false;
    if ~exist(cacheFile, 'file'); return; end
    try
        S = load(cacheFile, 'key', 'n');
    catch
        return
    end
    if ~isfield(S,'key') || ~strcmp(S.key, key); return; end
    if ~isfield(S,'n'); return; end
    keyOk = true; nStored = S.n;
end

% ---------------------------------------------------------------------------
function saveDataCache(cacheFile, key, payload, n)
    d = fileparts(cacheFile);
    if ~exist(d, 'dir'); mkdir(d); end   % created lazily, only when there is something to write
    data = payload;      % dataCacheLoad.m expects a variable named 'data'
    % THREE SEPARATE top-level variables, not one struct: that is what lets a pointer-only cache hit
    % load key/n and leave the [nRep x pLen] arrays on disk. -v7.3 for partial loading.
    % No 'se' is stored -- see finaliseReplicates on why the summary is not precomputed. Files written
    % by the earlier version carry an extra 'se' variable; it is simply never read, so they stay
    % valid and no cache is invalidated by this change.
    save(cacheFile, 'data', 'key', 'n', '-v7.3');
end

% ---------------------------------------------------------------------------
function val = mergeReplicateVal(baseVal, newVal, nComp)
    if isempty(baseVal); val = newVal; return; end
    val = baseVal;
    for c = 1:nComp
        fn = fieldnames(newVal{c});
        for i = 1:numel(fn)
            val{c}.(fn{i}) = [val{c}.(fn{i}); newVal{c}.(fn{i})];
        end
    end
end

% ---------------------------------------------------------------------------
function imOut = phaseRandomizeRuns(rawImC, stream)
    % One surrogate dataset: every run phase-randomized independently, then concatenated on the same
    % frame axis the point estimate used. Each run's [ny nx 1 T] block becomes [T x nVox] (phase
    % randomization is a per-CHANNEL transform over time, with one shared phase draw across voxels so
    % cross-voxel covariance survives -- see phaseScrambleRun.m) and back again. An FFT transform, not
    % a resample: frame count and order are preserved exactly, which is why every parameter's replicate
    % output has the same shape as the point estimate's, 'perFrame' included.
    nRun = numel(rawImC);
    out = cell(1, nRun);
    for r = 1:nRun
        blk = double(rawImC{r});
        [ny,nx,~,Tr] = size(blk);
        flat = reshape(blk, [], Tr).';            % [T x nVox]
        flatP = phaseScrambleRun(flat, stream);   % DC left untouched (that file's own default)
        out{r} = reshape(flatP.', ny, nx, 1, Tr);
    end
    imOut = cat(4, out{:});
end

% ---------------------------------------------------------------------------
function val = stackReplicates(fits, nComp)
    % Stack a set of replicate fits into [nRep x pLen] per component per parameter. VALUES ONLY -- the
    % spread is computed once, on the MERGED set, by finishReplicates below, so a resumed run's std is
    % over every replicate rather than only the newly-fitted tail.
    nRep = numel(fits);
    val = cell(1, nComp);
    for c = 1:nComp
        pn = fieldnames(fits{1}.components(c).params);
        vals = struct();
        for i = 1:numel(pn)
            pLen = numel(fits{1}.components(c).params.(pn{i}));
            M = zeros(nRep, pLen);
            for r = 1:nRep; M(r,:) = fits{r}.components(c).params.(pn{i})(:).'; end
            vals.(pn{i}) = M;
        end
        val{c} = vals;
    end
end

% ---------------------------------------------------------------------------
function rep = finishReplicates(val, nRep, nComp)
    % Trim to the requested count, and that is all. Raw replicates are kept UNSUMMARISED so a caller
    % can compute any std / custom interval / empirical p-value it wants -- the same "don't summarise
    % prematurely" convention fitIRF.m's own imBoot mode follows. See the file header's REPLICATES
    % ONLY for why no summary is produced anywhere in this function.
    rep = struct('val', {cell(1,nComp)});
    for c = 1:nComp
        fn = fieldnames(val{c});
        vals = struct();
        for i = 1:numel(fn)
            vals.(fn{i}) = val{c}.(fn{i})(1:nRep, :);
        end
        rep.val{c} = vals;
    end
end

% ---------------------------------------------------------------------------
function tf = anyPerKnot(components)
    tf = false;
    for c = 1:numel(components)
        fn = fieldnames(components(c).mode);
        for i = 1:numel(fn)
            if strcmp(components(c).mode.(fn{i}), 'perKnot'); tf = true; return; end
        end
    end
end

% ---------------------------------------------------------------------------
function finest = finestModeAcrossModel(model)
    % The FINEST granularity requested anywhere in RAW opts.model (before buildComponents resolves
    % seeds/claiming), ranked perVessel < perRun < perFrame == perKnot == perRunPoly -- what
    % opts.timeAvg's own reduction-strategy choice needs (see that call site). Scans opts.model directly
    % rather than the resolved component list on purpose: a caller's un-set '.mode.<param>' (an empty
    % struct field, per defaultOpts()'s emptyEntry) correctly contributes nothing here, exactly like it
    % correctly defaults to 'perVessel' once mkComponent's own subOr() resolves it later -- scanning
    % post-resolution would require building components before knowing how much of the time axis to
    % keep, which timeAvg's own reduction has to decide FIRST.
    %
    % mode='perRunPoly([...])' (2026-09-07) is matched by PREFIX, not exact strcmp against `rank`'s own
    % literal keys -- its mode string carries an embedded, per-parameter degree list (see
    % fitComponentFit.m's own mode='perRunPoly' doc), so no single literal can name it. It ranks with
    % perFrame/perKnot: a within-run polynomial trajectory has no time axis left once ANY time-averaging
    % happens, exactly like they do.
    rank = struct('perVessel',1, 'perRun',2, 'perFrame',3, 'perKnot',3);
    best = 1; finest = 'perVessel';
    types = fieldnames(model);
    for t = 1:numel(types)
        entries = model.(types{t});
        for e = 1:numel(entries)
            if ~isstruct(entries(e).mode); continue; end
            pn = fieldnames(entries(e).mode);
            for p = 1:numel(pn)
                m = entries(e).mode.(pn{p});
                if startsWith(m, 'perRunPoly')
                    r = 3;
                elseif isfield(rank, m)
                    r = rank.(m);
                else
                    continue
                end
                if r > best; best = r; finest = m; end
            end
        end
    end
end

% ---------------------------------------------------------------------------
function [knotIdx, keepFrame] = resolveKnotIdx(vesselV, fld, nRun, runLens, T, idStr)
    % Frame -> knot-group, from the IRF's own FIR design matrix. This is the ONE thing
    % fitComponentFit.m cannot work out for itself: it has no concept of a "knot" at all, only of a
    % caller-supplied grouping vector (see that file's own granularity note).
    %
    % vessel.(fld).irf.irfMat.X is ONE RUN's design ([nDes x nKnot], getIRFmat.m); the real frame axis
    % is runs concatenated in ascending order, which is exactly what fitIRF.m itself tiles
    % (Xfir = repmat(X, nRun, 1)), so row t here lines up with frame t of this function's own
    % concatenation.
    irfBase = getNestedField(vesselV, fld);
    ok = isstruct(irfBase) && isfield(irfBase,'irf') && isstruct(irfBase.irf) && ...
         isfield(irfBase.irf,'irfMat') && isstruct(irfBase.irf.irfMat) && isfield(irfBase.irf.irfMat,'X');
    assert(ok, 'fitVesselPatchTimeSeries:noKnotDesign', ...
        ['vessel %s: mode=''perKnot'' requires vessel.%s.irf.irfMat.X (the design matrix mapping real ' ...
         'frames onto knots) -- run fitIRF.m on this track first.'], idStr, fld);
    X = irfBase.irf.irfMat.X;
    assert(size(X,1)*nRun == T || all(runLens == size(X,1)), 'fitVesselPatchTimeSeries:knotDesignLength', ...
        ['vessel %s: vessel.%s.irf.irfMat.X has %d row(s) per run but this fld''s runs are %s frame(s) ' ...
         '-- the design and the data disagree about run length.'], idStr, fld, size(X,1), mat2str(runLens));
    Xfir = repmat(X, nRun, 1);
    keepFrame = sum(Xfir, 2) > 0;   % a frame with NO knot coverage belongs to no group at all
    keepFrame = keepFrame(:).';

    Xkept = Xfir(keepFrame, :);
    % A HARD group assignment needs a one-hot design. With overlapping knot weights a frame genuinely
    % belongs to several knots at once, and there is no well-defined single group to estimate it
    % against -- so this fails loudly rather than picking a winner by argmax and hiding the ambiguity.
    assert(all(sum(Xkept > 0, 2) <= 1), 'fitVesselPatchTimeSeries:perKnotRequiresOneHot', ...
        ['vessel %s: mode=''perKnot'' requires a one-hot (non-overlapping) knot design -- ' ...
         'vessel.%s.irf.irfMat.X has at least one covered frame with nonzero weight on more than one ' ...
         'knot.'], idStr, fld);
    [~, gi] = max(Xkept, [], 2);
    knotIdx = gi(:).';
    % fitComponentFit.m requires contiguous 1:nGroups labels and says so clearly if they are not, so a
    % genuinely unused knot column (zero weight on every frame) surfaces there rather than being
    % silently renumbered here -- the same fail-loud-never-relabel call made for runIdx.
end

% ---------------------------------------------------------------------------
function fit = reinflateFit(fit, keepFrame)
    % Put the KNOT-COVERAGE-excluded frames back as NaN. Only PER-FRAME-length arrays are touched: a
    % coarser granularity (perVessel/perRun/perKnot) has no row per frame to restore, and rewriting it
    % would corrupt a value that was never frame-indexed in the first place.
    nKept = nnz(keepFrame);
    for c = 1:numel(fit.components)
        for grp = {'params','start','lower','upper'}
            g = grp{1};
            fn = fieldnames(fit.components(c).(g));
            for i = 1:numel(fn)
                fit.components(c).(g).(fn{i}) = reinflateVec(fit.components(c).(g).(fn{i}), keepFrame, nKept);
            end
        end
    end
    fit.resnormPerFrame = reinflateVec(fit.resnormPerFrame, keepFrame, nKept);
end

% ---------------------------------------------------------------------------
function v = reinflateVec(v, keepFrame, nKept)
    if numel(v) ~= nKept; return; end   % not a per-frame array -- see reinflateFit
    out = nan(1, numel(keepFrame));
    out(keepFrame) = v;
    v = out;
end

% ---------------------------------------------------------------------------
function s = splitStructPerRun(s, runLens, T, doSplit, nKnotGroups)
    fn = fieldnames(s);
    for i = 1:numel(fn); s.(fn{i}) = splitPerRun(s.(fn{i}), runLens, T, doSplit, nKnotGroups); end
end

% ---------------------------------------------------------------------------
function segs = splitPerRun(val, runLens, T, doSplit, nKnotGroups)
    % doSplit=false (irf-path) returns the value unsplit -- see assembleResult's own note -- but still
    % passes through the identical-values collapse below.
    if nargin>=4 && ~doSplit; segs = collapseIdentical(val(:).'); return; end
    if nargin<5; nKnotGroups = 0; end
    % The SAME three-way length convention fitVesselTimeSeries.m's own splitResultPerRun uses:
    % length 1 broadcasts to every run, length nRun is already per-run, length T splits at the real
    % per-run frame boundaries. THEN (2026-09-10, Seb's own ask) any run of IDENTICAL values is
    % stored once -- see collapseIdentical -- so a value that does not vary between runs is a bare
    % number, not nRun copies of it, and a reader cannot mistake a per-vessel value for a per-run one.
    nRun = numel(runLens);
    val = val(:).';
    n = numel(val);
    segs = cell(1,nRun);
    if n == 1
        for r = 1:nRun; segs{r} = val; end
    elseif n == nRun
        for r = 1:nRun; segs{r} = val(r); end
    elseif n == T
        edges = [0 cumsum(runLens(:).')];
        for r = 1:nRun; segs{r} = val(edges(r)+1:edges(r+1)); end
    elseif nKnotGroups > 0 && n == nKnotGroups
        % mode='perKnot' -- a FOURTH length the original three-way convention did not cover. The knot
        % grouping is tiled IDENTICALLY across runs (one run's FIR design, repeated -- see
        % resolveKnotIdx), so every run genuinely shares the same per-knot values: the whole nKnot-long
        % array is broadcast into each run's cell, exactly as a perVessel scalar is. Not split, because
        % a knot group is not a subset of one run.
        for r = 1:nRun; segs{r} = val; end
    else
        error('fitVesselPatchTimeSeries:badPerRunSplitLength', ...
            ['a field resolved to %d value(s) -- expected 1 (perVessel/fixed), %d (perRun), %d ' ...
             '(perFrame) or %d (perKnot) -- cannot split per run.'], n, nRun, T, nKnotGroups);
    end
    segs = collapseIdentical(segs);
end

% ---------------------------------------------------------------------------
function v = collapseIdentical(v)
    % IDENTICAL VALUES ARE STORED ONCE (2026-09-10, Seb's own ask -- a human cannot see that three
    % printed cells hold the same number, and a per-vessel value printed as nRun copies was mistaken
    % for a per-run result on real data). Applied to every per-parameter output field, at BOTH levels:
    %   - WITHIN a run: a vector whose entries are all equal becomes that scalar;
    %   - ACROSS runs: a {1 x nRun} cell whose entries are all equal becomes that one entry, BARE.
    % So in a result, a CELL always means "this varies between runs" and a bare value always means
    % "one value for the whole vessel" -- by VALUE, not by mode (a 'perRun' parameter that converged to
    % the same number in every run collapses too; a 'fixed' one always does). Readers index runs with
    % perRunVal.m, which returns a bare value for any run. Note for floating-point fits this is only
    % ever structural in practice: two runs' free estimates never coincide exactly. A single-run
    % vessel (or a timeAvg fit's one pooled pseudo-run) collapses the same way: nothing varies by run,
    % so nothing is a cell -- no more {1x1} wrapper.
    if iscell(v)
        for i = 1:numel(v); v{i} = collapseIdentical(v{i}); end
        if ~isempty(v) && all(cellfun(@(x) isequaln(x, v{1}), v(2:end)))
            v = v{1};
        end
    elseif isnumeric(v) && numel(v) > 1 && all(isequaln(v(:), repmat(v(1), numel(v), 1)))
        v = v(1);
    end
end

% ---------------------------------------------------------------------------
function mUnion = unionMasksLocal(vessel, fld, labelList, largeSz)
    % Verbatim copy of fitVessel.m's own unionMasks (this codebase's established
    % small-helper-duplication convention for vessel-struct plumbing -- fitMultiVessel.m and
    % buildGaussianFitDiag.m carry the same pair).
    mUnion = false(largeSz);
    if isempty(labelList); return; end
    base = getNestedField(vessel, fld);
    if ~isfield(base,'rois') || isempty(base.rois); return; end
    labels = {base.rois.label};
    for i = 1:numel(labelList)
        idx = find(strcmp(labels, labelList{i}), 1);
        assert(~isempty(idx), 'fitVesselPatchTimeSeries:maskLabelNotFound', ...
            'mask label ''%s'' is not on vessel.%s.rois.', labelList{i}, fld);
        mUnion = mUnion | base.rois(idx).mask;
    end
end

% ---------------------------------------------------------------------------
function m = resolveIncludeMaskLocal(vessel, fld, includeMaskArg, largeSz)
    if isempty(includeMaskArg); m = true(largeSz); return; end
    m = unionMasksLocal(vessel, fld, includeMaskArg, largeSz);
end

% ---------------------------------------------------------------------------
function opts = defaultOpts(modelName)
    % The DEFAULT MODEL TEMPLATE (see file header NO-ARG CALL): main peak, then every remaining peak,
    % plus one background. Seeds are empty placeholders -- REQUIRED fields are present-but-empty rather
    % than absent, per the self-populating-default-opts convention, so a caller starting from this
    % struct can see every field that exists.
    %
    % modelName (default '') selects the NAMED-MODEL CALL template (see file header). '' is today's
    % unnamed two-gaussian-plus-background template, UNCHANGED. Any other non-empty name must match
    % 'N gaussian + background' for a positive integer N (parseGaussianCount below) or errors --
    % pattern match only, no fuzzy/guessed names.
    %
    % FULLY POPULATED, per parameter (2026-09-08, Seb's own ask): a caller should never have to write
    % out `.mode.a = 'perVessel'; .mode.x0 = 'perVessel'; ...` for every parameter of a model type just
    % to see/confirm the default -- that belongs in the template already. Every entry's .mode/.fixed/
    % .start/.lower/.upper carries one field per componentModelSpec.m paramName for that model type
    % (6 for gaussian, 1 for background): .mode.<param>='perVessel', .fixed.<param>=false,
    % .start/.lower/.upper.<param>=[] (present, not fixed to a literal, matching "no override" exactly
    % as an ABSENT field did before -- mkComponent's own subOr() treats [] and absent identically).
    if nargin<1; modelName = ''; end
    modelName = strtrim(modelName);
    gaussianEntry   = emptyEntryFor('gaussian');
    backgroundEntry = emptyEntryFor('background');
    opts.outFld  = '';
    opts.model   = struct();
    if strcmp(modelName, '')
        opts.model.gaussian = [gaussianEntry, gaussianEntry];   % (1) main, (2) all remaining
    else
        n = parseGaussianCount(modelName);
        opts.model.gaussian = repmat(gaussianEntry, 1, n);      % (1) main, (2..n) each its own slot
    end
    opts.model.background = backgroundEntry;
    opts.timeAvg = false;
    opts.addBase0 = true;
    % TOP-LEVEL, not per component (Seb's call, 2026-09-05): motion is a property of the TRACK, not of
    % any one peak, so one setting applies to every component and opts.model stays purely the model
    % equation. See file header SEED MOTION.
    opts.seedMotion = {};
    % Per-kind replicate groups. .nRep replaces the old opts.nBoot/opts.nPerm spelling; the two cache
    % axes are per-kind too, since a caller may well want the (cheap, pre-built) bootstrap offloaded
    % and pointed at while actively re-running the null, or vice versa.
    opts.boot = struct('nRep',0, 'offloadToCache',true, 'reComputeCache',false);
    opts.null = struct('nRep',0, 'offloadToCache',true, 'reComputeCache',false);
    opts.dataSeed = 0;
    opts.cacheDir = fullfile(pwd, 'fitVesselPatchTimeSeriesCache');
    opts.nPool = '';
    opts.includeMask = {};
    opts.excludeMask = {};
    opts.nIsochromatPerVoxDim = 7;
    opts.progress = true;
end

% ---------------------------------------------------------------------------
function n = parseGaussianCount(modelName)
    % Matches 'N gaussian + background' for a positive integer N (e.g. '1 gaussian + background',
    % '3 gaussian + background'), or 'main gaussian + background' as an alias for N=1 (2026-09-08,
    % Seb's own ask) -- see file header NAMED-MODEL CALL. Either base may carry an optional trailing
    % secondary-gaussian suffix -- ' + sec gaussian', ' + secondary gaussian', or ' + gaussians' --
    % which adds exactly ONE secondary slot (N = N+1), i.e. it populates opts.model.gaussian(2) in
    % preparation for a fit with one (or, once assigned real per-peak seeds, more) secondary peaks;
    % the three spellings are interchangeable, not different peak counts. Pattern match only, no
    % fuzzy/guessed names; N=0 or a non-integer/negative count is rejected the same way an
    % unrecognised name is.
    %
    % TWO LITERAL ALIASES (2026-09-09, Seb's own ask -- reads more naturally than spelling out a count
    % for the two he actually uses day to day): 'gaussian + background' == '1 gaussian + background'
    % (N=1), 'gaussians + background' == '2 gaussian + background' (N=2). Checked as plain exact-string
    % matches BEFORE the general pattern below, deliberately NOT folded into that regex, so they never
    % interact with (or get misread as) that pattern's own ' + gaussians' SUFFIX -- a different concept
    % entirely (SUFFIX adds one secondary slot to a base count; these two ARE the base count itself,
    % just spelled without a leading number).
    if strcmp(modelName, 'gaussian + background'); n = 1; return; end
    if strcmp(modelName, 'gaussians + background'); n = 2; return; end

    tok = regexp(modelName, ...
        '^(main|\d+) gaussian \+ background( \+ sec gaussian| \+ secondary gaussian| \+ gaussians)?$', ...
        'tokens', 'once');
    assert(~isempty(tok), 'fitVesselPatchTimeSeries:unknownModelName', ...
        ['unrecognised model name ''%s'' -- expected ''gaussian + background'', ''gaussians + ' ...
         'background'', ''N gaussian + background'' (N a positive integer), or ''main gaussian + ' ...
         'background'', the last two optionally followed by '' + sec gaussian'', '' + secondary ' ...
         'gaussian'', or '' + gaussians'' (e.g. ''1 gaussian + background'', ''main gaussian + ' ...
         'background + sec gaussian'').'], modelName);
    if strcmp(tok{1}, 'main')
        n = 1;
    else
        n = str2double(tok{1});
        assert(n >= 1, 'fitVesselPatchTimeSeries:unknownModelName', ...
            'unrecognised model name ''%s'' -- N must be a positive integer.', modelName);
    end
    if ~isempty(tok{2})
        n = n + 1;
    end
end

% ---------------------------------------------------------------------------
function e = emptyEntryFor(typeName)
    % One default opts.model.<typeName>(k) entry, fully populated per parameter (see defaultOpts's own
    % note) from componentModelSpec.m's own paramNames for that type -- 6 for gaussian, 1
    % (background's own 'b') for background. .start/.lower/.upper are [] per parameter, not absent:
    % mkComponent's own subOr() treats an empty field exactly like a missing one (falls back to the
    % seed/seedBounds default), so this changes nothing about behaviour, only what a caller SEES when
    % inspecting the template.
    %
    % .seed defaults to 'heuristic' (2026-09-08, Seb's own ask: "the default seed of everything should
    % be 'heuristic'"), not '' -- buildComponents expands a bare 'heuristic' to THIS call's own
    % '<fld>.heuristic', so a caller gets a working, no-prior-fit-required seed for every entry with no
    % setup at all. Two default gaussian entries seeding from the SAME single heuristic peak still
    % behaves correctly unchanged (entry 1 claims it, entry 2's own 'all' selector then finds nothing
    % left to claim and contributes no component -- exactly the "every remaining peak" contract with a
    % heuristic source that only ever has one peak to give).
    spec = componentModelSpec(typeName);
    e = struct('seed','heuristic', 'mode',struct(), 'fixed',struct(), ...
               'start',struct(), 'lower',struct(), 'upper',struct(), ...
               'seedBounds',defaultSeedBounds(typeName));
    for i = 1:numel(spec.paramNames)
        nm = spec.paramNames{i};
        e.mode.(nm)  = 'perVessel';
        e.fixed.(nm) = false;
        e.start.(nm) = [];
        e.lower.(nm) = [];
        e.upper.(nm) = [];
    end
end

% ---------------------------------------------------------------------------
function opts = defaultOptsFromFit(vessel, fld)
    % 2-ARG CALL (2026-09-08, Seb's own ask, generalizing drawVesselPatchPeaks.m's own established
    % 2-arg call to this function) -- mirrors the NAMED-MODEL/no-arg calls' own "just give me opts, no
    % real work" contract, but derives opts.model.gaussian/.background from vessel(1).(fld)'s ALREADY-
    % EXISTING starting fit (fld one level DEEPER than the real 3-arg call's own fld -- e.g.
    % 'tsIm.gaussAnat', not 'tsIm' -- see file header 2-ARG CALL and drawVesselPatchPeaks.m's own FLD
    % note, the same fld-meaning divergence, now also across THIS function's own call forms).
    %
    % THE BUG THIS EXISTS TO SIDESTEP (2026-09-08, real incident): a caller who clones one entry onto a
    % second (`opts.model.gaussian(2) = opts.model.gaussian(1); opts.model.gaussian(2).mode.x0 =
    % 'perRun';`) ends up with BOTH entries pointing at the exact same UNINDEXED seed path. Peaks are
    % claimed in array order (see buildComponents/emptyEntryFor's own note on the heuristic case), so
    % entry 1 -- processed first -- silently claims EVERY peak the path has to offer (there was nothing
    % claimed yet to subtract), and entry 2 then finds nothing left and contributes ZERO components: its
    % distinct .mode never manifests anywhere, with no error at all (an unindexed entry matching zero
    % remaining peaks is valid, ordinary behaviour -- see file header SEED PATHS). This function avoids
    % the trap structurally, by naming each entry's ROLE ('(main)' / '(secondary)', 2026-09-10 -- before
    % that, an explicit index `'<fld>(<k>)'` on entry 1 and an unindexed entry 2, which worked only
    % through claiming order), a pattern in the same spirit as the one drawVesselPatchPeaks.m's own
    % SEED CONSTRUCTION note uses for its own re-fit: role/explicit selectors are never subject to
    % claiming, so N entries into the SAME source can never collide. (Older text follows.) Exactly the
    % pattern drawVesselPatchPeaks.m's own SEED CONSTRUCTION note already uses for the same reason:
    % explicit indices are never subject to claiming, so N entries into the SAME source can never
    % collide.
    %
    % AT MOST 2 ENTRIES (2026-09-10, Seb's own ask, harmonizing with drawVesselPatchPeaks.m's own
    % always-2-entry contract) -- entry 1 is the MAIN peak, EXPLICITLY indexed ('<fld>(1)'); entry 2 is
    % ONE UNINDEXED template ('<fld>', no selector) representing every secondary peak collectively, not
    % one entry per secondary peak actually present. buildComponents' own existing claiming machinery
    % (see file header SEED PATHS) is what expands that single unindexed entry back into one component
    % per remaining peak at the real 3-arg call, deriving each one's own .start/.lower/.upper FRESH from
    % that peak's own seed data then -- exactly the mechanism (1) above already relies on, so this
    % template need not (and structurally cannot, being ONE entry) carry N separately-derived copies.
    %   PRE-2026-09-10 BEHAVIOR (removed): one entry PER PEAK actually present, each with its own
    %   explicit index and its own individually-derived .seedBounds. Correct, but redundant with what
    %   buildComponents already does for an unindexed entry, and out of step with
    %   drawVesselPatchPeaks.m's own single-shared-secondary-template convention.
    %   REPRESENTATIVE SOURCE FOR THE SHARED SECONDARY TEMPLATE: entry 2's own .mode/.fixed/.seedBounds
    %   are derived from the starting fit's OWN peak 2 (the first actual secondary), when one exists --
    %   a real behavior change from copying every secondary's own individual values: if secondaries in
    %   the starting fit do not all already share identical .mode/.fixed/.seedBounds (unusual in
    %   practice, since they typically came from this same shared-template convention in the first
    %   place), only peak 2's own values survive into the derived template; the others are not silently
    %   dropped from the FIT itself (buildComponents still claims and fits every one of them), only from
    %   what this TEMPLATE reports about their individual .mode/.fixed/.seedBounds. "MAXIMUM 2" IS A CAP,
    %   NOT A FLOOR (2026-09-10 correction, caught by testEmptyStartStillProducesAValidRefit's own
    %   failure): a starting fit with NO secondary peak at all (nPeak==1) gets exactly ONE entry, not a
    %   manufactured gaussian(2) -- unlike drawVesselPatchPeaks.m's own always-2 template, whose
    %   heuristic-seeded secondary sits inert until a human actually clicks, THIS function's own 2-arg
    %   call feeds directly into an unconditional 3-arg re-fit, where a bare '.seed=heuristic' secondary
    %   ALWAYS resolves to a real synthetic peak -- including one here for a genuinely single-peak
    %   starting fit would invent a second peak nobody asked for.
    %
    % .start/.lower/.upper are LEFT EMPTY (2026-09-09, Seb's own correction of this function's earlier
    % behavior -- "setting them at this high level should be a user override", not something this
    % template call bakes in on its own) -- each entry's .seed already points at that SAME peak via an
    % explicit index, so mkComponent's own seed/seedBounds resolution (see emptyEntryFor's own note --
    % an empty .start/.lower/.upper is treated exactly like a genuinely missing one) reproduces the
    % SAME literal numeric values a caller would have gotten from a snapshot copy, but resolved FRESH
    % at the actual 3-arg fit call, from whatever vessel.(fld) holds AT THAT TIME -- not frozen to
    % whatever it happened to be when this 2-arg call ran. A caller who genuinely wants a hard literal
    % override (pinning a value regardless of any later change) can still set .start/.lower/.upper by
    % hand afterward, same as building opts.model from scratch -- this call's job is only to hand back
    % a correctly-seeded, collision-free template, not to pre-decide which values are frozen.
    % (Lifts an earlier restriction to a single-run starting fit -- that restriction existed only
    % because .start/.lower/.upper's own per-run cell shape needed collapsing to one literal value;
    % with nothing to collapse anymore, a multi-run starting fit derives just as well.)
    %
    % opts.includeMask/.excludeMask ARE copied from the starting fit (2026-09-08, Seb's own real
    % script always re-set these by hand after the 1-arg call) -- DELIBERATELY UNLIKE
    % drawVesselPatchPeaks.m's own 2-arg call, which leaves them at {} because ITS internal re-fit
    % wrapper always unions with the starting fit's own masks automatically regardless of caller opts.
    % This function's own OUTPUT is a real, standalone opts struct headed for the ordinary 3-arg call,
    % which has no such automatic inheritance -- so copying here is what makes the returned opts usable
    % as a real starting template (a caller can still clear or extend either field afterward, same as
    % every other opts field).
    %
    % Every other opts field (.timeAvg, .nIsochromatPerVoxDim, ...) is left at defaultOpts()'s own plain
    % default -- this call is about "prepopulating MODEL options" (Seb's own words), not the surrounding
    % fit configuration, which usually differs on purpose between the summary starting fit and whatever
    % comes next (e.g. Seb's own real case: a timeAvg=true starting fit feeding a timeAvg=false, finer-
    % grained follow-up).
    if iscell(vessel); vessel = vessel{1}; else; vessel = vessel(1); end
    startFit = getNestedField(vessel, fld);
    idStr = strtrim([char(string(vessel.sId)) ' ' char(string(vessel.label))]);
    assert(isstruct(startFit) && isfield(startFit,'model') && isfield(startFit.model,'gaussian') ...
           && ~isempty(startFit.model.gaussian) && isfield(startFit.model,'background'), ...
        'fitVesselPatchTimeSeries:noStartingFit', ...
        ['vessel %s: vessel.%s has no .model.gaussian/.background -- it does not look like a ' ...
         'fitVesselPatchTimeSeries.m result (run one first, e.g. with opts.timeAvg=true).'], idStr, fld);

    nPeak = numel(startFit.model.gaussian);

    % SEEDS ARE BY ROLE, NOT BY NUMBER (2026-09-10, Seb's own ask): '<fld>(main)' and
    % '<fld>(secondary)' -- the roles the two entries actually carry -- rather than the earlier
    % '<fld>(1)' + unindexed pair, which encoded the same thing through array positions and claiming
    % order (see file header SEED PATHS for both selectors). '(secondary)' resolves to every peak but
    % the main one at the real 3-arg call, on each vessel separately, and to nothing on a single-peak
    % vessel, so one template still runs across a batch with mixed peak counts.
    %
    % .seedBounds IS COPIED, not reconstructed (2026-09-10): every fitVesselPatchTimeSeries.m result
    % now carries the spec each component was built from (.model.<type>(k).seedBounds), so the
    % earlier lower/upper inversion (lossy for aspectRatio, whose bound is symmetrized after the
    % fact) is gone. A result predating that field cannot be used here -- re-run the fit that made it.
    gMain = emptyEntryFor('gaussian');
    mainFit = startFit.model.gaussian(1);
    gMain.seed       = sprintf('%s(main)', fld);
    gMain.mode       = mainFit.mode;
    gMain.fixed      = mainFit.fixed;
    gMain.seedBounds = seedBoundsFromFit(mainFit, sprintf('%s.model.gaussian(1)', fld), idStr);

    opts = defaultOpts();
    if nPeak >= 2
        % "maximum 2" is a CAP, not a floor (2026-09-10 correction, caught by
        % testEmptyStartStillProducesAValidRefit's own failure): a single-peak starting fit must stay a
        % single ENTRY here, not gain a manufactured gaussian(2). Unlike drawVesselPatchPeaks.m's own
        % always-2 template (whose heuristic-seeded secondary sits inert until a human actually clicks),
        % this function's own 2-arg call feeds directly into an unconditional 3-arg re-fit -- a bare
        % '.seed=heuristic' secondary there ALWAYS resolves to a real synthetic peak, so including one
        % here for a genuinely single-peak starting fit would invent a second peak nobody asked for.
        gSecondary = emptyEntryFor('gaussian');
        secFit = startFit.model.gaussian(2);   % representative -- see AT MOST 2 ENTRIES note above
        gSecondary.seed       = sprintf('%s(secondary)', fld);
        gSecondary.mode       = secFit.mode;
        gSecondary.fixed      = secFit.fixed;
        gSecondary.seedBounds = seedBoundsFromFit(secFit, sprintf('%s.model.gaussian(2)', fld), idStr);
        opts.model.gaussian = [gMain, gSecondary];
    else
        opts.model.gaussian = gMain;
    end

    bFit = startFit.model.background;
    bE = emptyEntryFor('background');
    bE.seed  = fld;
    bE.mode  = bFit.mode;
    bE.fixed = bFit.fixed;
    bE.seedBounds = seedBoundsFromFit(bFit, sprintf('%s.model.background', fld), idStr);
    opts.model.background = bE;

    opts.includeMask = startFit.includeMask;
    opts.excludeMask = startFit.excludeMask;
end

% ---------------------------------------------------------------------------
function sb = seedBoundsFromFit(entryFit, what, idStr)
    % The .seedBounds a fit result's component carries (its own spec, defaults filled in -- see
    % mkComponent/assembleResult), copied verbatim. No reconstruction from .lower/.upper: a result
    % without the field predates 2026-09-10 and has to be re-fit, which is the loud answer rather than
    % a silently lossy inversion (aspectRatio's bound is symmetrized after the spec is applied, so
    % .lower/.upper alone cannot give the spec back exactly).
    assert(isfield(entryFit,'seedBounds') && isstruct(entryFit.seedBounds), ...
        'fitVesselPatchTimeSeries:noSeedBoundsInFit', ...
        ['vessel %s: vessel.%s carries no .seedBounds -- fitVesselPatchTimeSeries.m results written ' ...
         'before 2026-09-10 predate that field. Re-run the fit that produced it.'], idStr, what);
    sb = entryFit.seedBounds;
end

% ---------------------------------------------------------------------------
function opts = fillOptsDefaults(opts)
    opts = fillOptsFromDefaults(opts, defaultOpts());
    if ischar(opts.includeMask); opts.includeMask = {opts.includeMask}; end
    if ischar(opts.excludeMask); opts.excludeMask = {opts.excludeMask}; end
    if ischar(opts.seedMotion); opts.seedMotion = {opts.seedMotion}; end
    if isempty(opts.seedMotion); opts.seedMotion = {}; end
    assert(iscellstr(opts.seedMotion), 'fitVesselPatchTimeSeries:badSeedMotion', ...  %#ok<ISCLSTR>
        'opts.seedMotion must be a dot-path char or a cellstr of them (got %s).', class(opts.seedMotion));
    assert(isstruct(opts.model) && ~isempty(fieldnames(opts.model)), ...
        'fitVesselPatchTimeSeries:emptyModel', ...
        'opts.model must have at least one component field (e.g. .gaussian) -- it IS the model equation.');
end

% ---------------------------------------------------------------------------
function printOptsHelp()
    fprintf('fitVesselPatchTimeSeries opts -- allowed value(s) per field (default in parentheses):\n');
    fprintf('  opts.outFld       : any nonempty char (REQUIRED)                                 ('''')\n');
    fprintf('  opts.model.<type>(k) : the MODEL EQUATION -- every field a registered model (componentModelSpec.m),\n');
    fprintf('                      every element one additive component. Default: gaussian(1)+gaussian(2)+background.\n');
    fprintf('                      Named-model call, fitVesselPatchTimeSeries(''N gaussian + background''), returns an\n');
    fprintf('                      N-entry template instead (1=main only; N-1 secondary slots) -- see NAMED-MODEL CALL.\n');
    fprintf('    .seed           : ''path(main)'' | ''path(secondary)'' | ''path(2 3)'' | ''path'' -- dot-path from the VESSEL\n');
    fprintf('                      ROOT to a prior fit (a fitVesselPatchTimeSeries.m result, seeded at its own granularity,\n');
    fprintf('                      or a fitPatchVessels.m .fits); (secondary) = every peak but the main; unindexed = every\n');
    fprintf('                      peak not claimed by an earlier entry                              (REQUIRED, ''heuristic'')\n');
    fprintf('                      | ''path.heuristic'' | ''path(heuristic)'' -- no prior fit needed; ONE ad hoc peak computed\n');
    fprintf('                      from path''s own time-averaged image via the legacy single-peak heuristic, seeded at center\n');
    fprintf('                      | ''heuristic'' (bare, the DEFAULT) -- shorthand for ''<fld>.heuristic'', this call''s own fld\n');
    fprintf('    .mode.<param>   : ''perVessel'' | ''perFrame'' | ''perRun'' | ''perKnot'' |\n');
    fprintf('                      ''perRunPoly([<degrees>])'' (embedded degree list, e.g. ''perRunPoly([0 1])'' -- see file\n');
    fprintf('                      header mode=''perRunPoly'')                                          (''perVessel'')\n');
    fprintf('    .fixed.<param>  : true | false                                                  (false)\n');
    fprintf('    .start/.lower/.upper.<param> : LITERAL numeric override                         (from seed / seedBounds)\n');
    fprintf('    .seedBounds     : PER MODEL TYPE -- only that type''s own parameters, other names error:\n');
    fprintf('                      gaussian:   .x0y0 (mm, additive; [] = 2 voxels per axis), .theta (rad, additive; Inf),\n');
    fprintf('                                  .a/.radius/.aspectRatio (multiplicative; 2)\n');
    fprintf('                      background: .b (multiplicative; 2)\n');
    fprintf('  opts.timeAvg      : true | false -- collapse time by the FINEST mode requested (perVessel=all;\n');
    fprintf('                      perRun=within-run only; perFrame/perKnot/perRunPoly rejected)      (false)\n');
    fprintf('  opts.addBase0     : true | false -- IRF-PATH ONLY; false subtracts fitIRF.m''s baked-in DC baseline (true)\n');
    fprintf('  opts.seedMotion   : '''' | dot-path | cellstr -- getPreprocMotion.m source(s) SUBTRACTED from every\n');
    fprintf('                      component''s x0/y0 seed, per frame; requires mode.x0/.y0=''perFrame''    ({})\n');
    fprintf('  opts.boot.nRep    : 0 = off | N -- IRF-PATH ONLY, refit N of fitIRF.m''s own imBoot replicates (0)\n');
    fprintf('  opts.null.nRep    : 0 = off | N -- irf-path: N of imPerm; ts-path: generate N surrogates here (0)\n');
    fprintf('  opts.<boot|null>.offloadToCache : true | false -- MEMORY axis; true keeps raw replicates on disk, struct holds a pointer (true)\n');
    fprintf('  opts.<boot|null>.reComputeCache : true | false -- CPU axis; false reuses what exists, true recomputes (false)\n');
    fprintf('  opts.dataSeed     : integer, seeds the ts-path null''s own local RandStream            (0)\n');
    fprintf('  opts.cacheDir     : directory for the replicate DATA cache   (<pwd>/fitVesselPatchTimeSeriesCache)\n');
    fprintf('  opts.nPool        : '''' | positive integer -- only when replicates start a pool         ('''')\n');
    fprintf('  opts.includeMask  : char | cellstr of vessel.(fld).rois label(s)                   ({})\n');
    fprintf('  opts.excludeMask  : char | cellstr of vessel.(fld).rois label(s)                   ({})\n');
    fprintf('  opts.nIsochromatPerVoxDim : positive ODD integer                                   (7)\n');
    fprintf('  opts.progress     : true | false                                                   (true)\n');
end
