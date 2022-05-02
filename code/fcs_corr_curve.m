function [time, curve] = fcs_corr_curve(h, varargin)
% function [time, curve] = fcs_corr_curve(h, varargin)
%   returns the estimator for the current correlation curve and
%   calculates the time bin bases.
%
%       time    The centers of the time bins
%       curve   The current correlation curve
%
%       h       The correlator handle
%       n_norm  Normalization method
%               0 : no normalization         : <f(t)g(t+tau)> - <g(t)>*<g(t)>
%               1 : asymmetric normalization : <f(t)g(t+tau)>/(<f(t)><g(t)>) - 1
%                   (the normalization mean is taken over the whole trace)
%               2 : symmetric normalization  : <f(t)g(t+tau)>/(<f(t)><g(t)>) - 1
%                   (the normalization mean is taken only over the overlap, default)
%               Asymmetric normalization approaches the result of symmetric
%               normalization for traces long compared to the maximum lag time.

params.n_norm = 2;
params = omex_read_params(params, varargin);

[time, curve] = omexFcs(3, h, params);
end