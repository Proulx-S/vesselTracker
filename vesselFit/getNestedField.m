function val = getNestedField(s, path)
% GETNESTEDFIELD  Dot-path struct GET -- s.(parts{1}).(parts{2})..., walked one segment at a time
% since MATLAB's own dynamic-fieldname syntax (s.(path)) rejects a dotted string outright. Tolerant,
% like isfield/resolveIrfField.m's own style: [] (not an error) the moment any segment is missing or
% the value stops being a struct, rather than only ever being called after a separate existence check.
%
% Promoted 2026-07-31 (previously a local function inside fitVesselTimeSeries.m, added there 2026-07-30
% for that file's own dot-path fld argument -- see that file's own FLD AS DOT-PATH) once a SECOND
% caller (gaussFitToVecField.m, resolving both a per-frame Gaussian fit result AND its own PARENT
% path's metadata) needed the identical resolution -- see this codebase's established
% promote-once-2+-callers-need-it convention (irfContainerCore.m, resolveIrfField.m,
% principalAxisProjection.m, etc.).
%
% INPUT
%   s    : a struct (or anything -- non-struct input is a graceful [], not an error).
%   path : dot-separated field-name string, e.g. 'tsImCentered.irf.irfIm'.
%
% OUTPUT
%   val : s.(path)'s own resolved value, or [] if any segment is missing/not a struct along the way.
%
%   val = getNestedField(vessel, 'tsImCentered.irf.irfIm')
%
% See also fitVesselTimeSeries.

    parts = strsplit(path, '.');
    val = s;
    for i = 1:numel(parts)
        if ~isstruct(val) || ~isfield(val, parts{i})
            val = [];
            return
        end
        val = val.(parts{i});
    end
end
