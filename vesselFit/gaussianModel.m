function v = gaussianModel(X, Y, a, x0, y0, sx, sy, theta, b)
% GAUSSIANMODEL  THE rotated-2D-Gaussian model formula used throughout this pipeline's Gaussian fit
% (fitVesselProfile.m) -- isolated on its own, standalone, mirroring parabolicModel.m (same file
% layout/1D-mode convention, so the two are easy to compare/audit side by side), so it can be
% read/audited on its own and so a diagnostic plot can never silently drift onto a different formula
% than what was actually fit.
%
% ADDITIVE everywhere -- v = a*exp(-(xr^2/(2*sx^2) + yr^2/(2*sy^2))) + b, where (xr,yr) is (X,Y)
% rotated by theta into the fit's own frame around (x0,y0). Unlike parabolicModel.m's two-compartment
% (blood-XOR-tissue) model, `b` here rides UNDER the peak everywhere, including at the vessel center --
% this is DELIBERATE: the Gaussian fit is not meant to be a physically realistic model of blood/tissue
% signal (just a smooth, well-behaved reference/seed fit -- see fitVesselProfile.m), so an additive
% baseline is fine here even though it would be wrong for the parabolic (realistic) model.
%
% NOTE: this file defines the formula for reading/plotting only -- the actual Gaussian FIT itself
% (fitVesselProfile.m) still uses MATLAB's Curve Fitting Toolbox `fittype`/`fit()` mechanism with its
% own literal formula string (a `fittype` string can't call an arbitrary function handle the same way
% `lsqnonlin`'s residual can), so that string and this function must be kept algebraically identical
% by hand -- there are only the two places (the `fittype` string, and this file) where the formula
% needs to match.
%
% 1D MODE (for auditing -- collapse the 2D model to a single EFFECTIVE-DISTANCE profile, without
% building a 2D grid): leave Y empty ([]) to get the profile as a function of X = an EFFECTIVE
% distance from center in MAJOR-axis-equivalent units (whichever of sx/sy is LARGER); leave X empty
% for Y = the same but in MINOR-axis-equivalent units (whichever is SMALLER). X (or Y) is NOT a literal
% on-axis coordinate (the point it came from need not lie anywhere near that axis) -- it is
% `rho*sMajor` (or `rho*sMinor`), where `rho = sqrt((xr/sx)^2+(yr/sy)^2)` is the point's OWN full,
% angle-correct elliptical-normalized radius (see the 2D branch below). This is why the formula
% reduces to exactly `a*exp(-(X/sMajor)^2/2)+b` (etc.): substituting X=rho*sMajor gives
% `a*exp(-(rho*sMajor/sMajor)^2/2)+b = a*exp(-rho^2/2)+b`, algebraically IDENTICAL to the 2D formula
% evaluated at whatever (x,y) produced that rho -- so a caller can pre-compute `rho*sMajor` for an
% arbitrary (possibly off-axis) point and this 1D mode reports the exact same value the full 2D
% evaluation would have. Always non-negative (no separate positive/negative side -- only X^2 enters
% the formula), matching a true "distance from center," not a signed coordinate. x0/y0/theta are
% unused in this mode -- the rotation/offset is already baked into rho by whoever computed it.
% Feeding BOTH X and Y (2D mode, the normal case) evaluates the full rotated-ellipse formula as usual.
%
% INPUT
%   X, Y  : coordinates to evaluate at (any size/shape; typically a meshgrid or a vector). 2D mode:
%           both non-empty, absolute image coordinates. 1D mode: exactly one of X/Y non-empty, an
%           unsigned EFFECTIVE distance from center (`rho*sMajor` for X, `rho*sMinor` for Y -- see
%           above), not a literal on-axis coordinate.
%   a     : peak amplitude above background.
%   x0,y0 : ellipse center. Unused in 1D mode.
%   sx,sy : standard deviations (not radii), along xr/yr respectively (xr = the X-axis rotated by
%           theta, yr = the axis PERPENDICULAR to xr, i.e. rotated by theta+pi/2 -- see the 2D branch
%           below). NOT fixed to "major"/"minor" by this FORMULA's own definition (it evaluates
%           correctly for ANY sx,sy,theta, regardless of which one happens to be larger) -- this file
%           itself is agnostic to convention. This pipeline's OWN FITS (fitVesselProfile.m always;
%           fitVessel.m's 'gaussian' method whenever sx/sy/theta are ALL jointly free with the same
%           mode -- see that file's own THETA CANONICALIZATION note) DO canonicalize their output via
%           canonicalizeEllipseShape.m so that sx>=sy and theta is in (-pi/2,pi/2] always -- i.e. sx
%           IS the major axis, sy IS the minor, in values produced by those fits. This relies on a
%           genuine two-fold degeneracy: the SAME ellipse is equally well described by (theta,sx,sy) or
%           (theta -/+ pi/2, sy,sx) (swap the axis labels, rotate a further quarter turn). Whenever
%           "along the long/short axis" is needed from a VALUE THAT MIGHT NOT BE CANONICALIZED (e.g. a
%           fitVessel.m 'gaussian' result where sx/sy were 'fixed', not jointly fit), use
%           sMajor=max(sx,sy) / sMinor=min(sx,sy) -- exactly what the 1D mode below does -- rather than
%           assuming sx or sy is always one or the other.
%   theta : ellipse rotation (radians). Unused in 1D mode.
%   b     : background level (additive everywhere, jointly fit -- see header note).
%
%   v = gaussianModel(Xf, Yf, a, x0, y0, sx, sy, theta, b)                  % 2D
%   v = gaussianModel(rho*max(sx,sy), [], a, x0, y0, sx, sy, theta, b)      % 1D, major-axis effective distance
%   v = gaussianModel([], rho*min(sx,sy), a, x0, y0, sx, sy, theta, b)      % 1D, minor-axis effective distance

    if ~isempty(X) && ~isempty(Y)
        xr =  (X-x0)*cos(theta) + (Y-y0)*sin(theta);
        yr = -(X-x0)*sin(theta) + (Y-y0)*cos(theta);
        v = a*exp(-(xr.^2/(2*sx^2) + yr.^2/(2*sy^2))) + b;
    elseif ~isempty(X) && isempty(Y)
        sMajor = max(sx,sy);
        v = a*exp(-(X.^2)/(2*sMajor^2)) + b;
    elseif isempty(X) && ~isempty(Y)
        sMinor = min(sx,sy);
        v = a*exp(-(Y.^2)/(2*sMinor^2)) + b;
    else
        error('gaussianModel:emptyXY', 'At least one of X or Y must be non-empty.');
    end
end
