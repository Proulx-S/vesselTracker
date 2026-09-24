function s = contourLevelFracStr(fracs)
% CONTOURLEVELFRACSTR  Human-readable "10%/90%" style string for a vector of peak-amplitude fractions
% (the level(s) a fitted-vessel contour is drawn at, e.g. drawGaussianPeakContour.m's own levelFracs) --
% shared formatting so every panel/legend that reports a contour's own threshold level(s) uses identical
% wording. Promoted once 3 callers (showVessel.m, fitPatchVesselsDiag.m, fitVesselDiag.m) needed the
% identical formatting -- same promote-once-2+-callers-need-it convention as hasIrfField.m.
%
% INPUT
%   fracs : numeric vector, each a fraction of peak amplitude in (0,1) (e.g. [0.1 0.9]).
%
% OUTPUT
%   s : e.g. '10%/90%' for [0.1 0.9], '50%' for a scalar 0.5. Percentages are ROUNDED (round(f*100)) --
%       every caller today uses round-number fractions (0.1/0.3/0.5/0.7/0.9), so this never loses
%       information in practice; a caller with a genuinely fractional-percent level would want its own
%       formatting instead.
%
%   s = contourLevelFracStr([0.1 0.9])   % -> '10%/90%'

    parts = arrayfun(@(f) sprintf('%g%%', round(f*100)), fracs, 'UniformOutput', false);
    s = strjoin(parts, '/');
end
