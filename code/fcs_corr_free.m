function fcs_corr_free(handle)
% function fcs_corr_free(h)
% frees a correlator handle. If h==0 all handles are destroyed.
%
%       h       The correlator handle

if (nargin == 0)
    handle = -1;
end

omexFcs(4, handle);
end