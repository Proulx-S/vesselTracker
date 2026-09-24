function ensureParPool(nPool)
% ENSUREPARPOOL  Starts a parallel pool if none is running; REUSES any existing pool exactly as-is,
% regardless of its own worker count, never tearing one down to match nPool -- a caller wanting a
% DIFFERENT pool must shut the existing one down itself (delete(gcp)) first, never done automatically
% here. Shared by every parfor-using function in this project (fitVesselTimeSeries.m, fitVessel.m,
% alignVessel2.m) so each independently guarantees SOME pool exists before its own parfor, without
% depending on a caller to have started one -- see those files' own opts.nPool for how this gets called.
%
% Root cause this works around: parfor's own AutoCreate is disabled in this MATLAB installation
% (parallel.Settings().Pool.AutoCreate == 0) -- license/toolbox/cluster profile are otherwise fine, but
% parfor never auto-starts a pool on its own, so every parfor silently runs serially "on the client"
% unless something explicitly starts one first.
%
% INPUT
%   nPool : OPTIONAL positive integer -- worker count to request ONLY when actually STARTING a NEW pool
%           (no effect at all if a pool is already running). Omitted/empty -- parpool()'s own default
%           cluster profile worker count.
%
% No output -- side-effecting only (starts a pool as a side effect, or does nothing if one's already
% there).
%
% WORKER GUARD: gcp('nocreate') called FROM WITHIN a worker (e.g. this file's own callers nested inside
% another already-running parfor -- fitVessel called from fitVesselTimeSeries's own per-vessel parfor,
% or from attachReplicateErrorBars'/attachPhaseRandomizeData's own per-replicate parfor) does NOT see
% the pool that worker is already part of, so the emptiness check below would otherwise wrongly try to
% start a SECOND pool -- which MATLAB refuses outright ("A parallel pool cannot be started from a
% worker, only from a client MATLAB"), erroring instead of just being redundant. getCurrentTask()
% returns non-empty iff running inside a worker right now, regardless of what gcp reports from there --
% checked FIRST, before gcp, so this is a safe no-op whenever called from a nested context: being
% inside a worker at all is, by construction, proof some pool already exists.

    if nargin<1; nPool = []; end
    if ~isempty(getCurrentTask()); return; end
    if ~isempty(gcp('nocreate')); return; end
    if isempty(nPool)
        parpool();
    else
        parpool(nPool);
    end
end
