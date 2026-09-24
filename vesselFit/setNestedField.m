function s = setNestedField(s, path, val)
% SETNESTEDFIELD  Dot-path struct SET -- the STORAGE-side counterpart of getNestedField.m.
% Intermediate segments are created as empty structs if missing (mirrors plain dynamic-fieldname
% assignment's own auto-vivification, e.g. s.a.b = 1 creating s.a on the fly).
%
% Promoted 2026-08-20 (previously a local function inside fitVesselTimeSeries.m, added there
% 2026-07-30 for that file's own dot-path fld argument -- see that file's own FLD AS DOT-PATH) once a
% SECOND caller (irfSvd.m, storing its own SVD result nested under the caller-chosen fldOut path)
% needed the identical dot-path SET -- see this codebase's established promote-once-2+-callers-need-
% it convention (getNestedField.m's own header lists the same pattern for the GET side).
%
% INPUT
%   s    : a struct.
%   path : dot-separated field-name string, e.g. 'tsIm.irf.irfIm.svd'.
%   val  : the value to store at s.(path).
%
% OUTPUT
%   s : the struct with s.(path) set to val (intermediate structs auto-vivified as needed).
%
%   s = setNestedField(vessel, 'tsIm.irf.irfIm.svd', svdResult)
%
% See also getNestedField, hasNestedField, fitVesselTimeSeries, irfSvd.

    parts = strsplit(path, '.');
    s = setNestedFieldParts(s, parts, val);
end

% ---------------------------------------------------------------------------
function s = setNestedFieldParts(s, parts, val)
    if isscalar(parts)
        s.(parts{1}) = val;
        return
    end
    if ~isstruct(s) || ~isfield(s, parts{1}) || ~isstruct(s.(parts{1}))
        s.(parts{1}) = struct();
    end
    s.(parts{1}) = setNestedFieldParts(s.(parts{1}), parts(2:end), val);
end
