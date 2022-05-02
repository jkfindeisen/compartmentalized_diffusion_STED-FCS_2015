function handle = fcs_corr_init(varargin)
% function handle = fcs_corr_init(varargin)
% initializes a correlator object and returns a handle on it. The correlator
% may be of linear type or of multi tau type as provided by ordinary hardware
% correlators. See papers of Schaetzel for details.
% Variable parameters are
%
%   type            The correlator type.
%                   'multi'  : multi tau correlator (default)
%                   'linear' : linear correlator    (via Fast Fourier Transform)
%   f_time          The basic binning time step of the data trace.
%                   The default is 1.
%   n_compression   The compression rate. This describes the factor by
%                   which the original bin is augmented.
%                   The default is also 1.
%
%   only for linear correlator    :
%   n_steps         Number of binning steps     (default : 100)
%
%   only for multi tau correlator :
%   n_octaves       Number of octaves           (default : 20)
%   n_bins_per_oct  Number of bins per octave   (default : 8)

params.type          = 'multi';
params.f_time        = 1;
params.n_compression = 1;
% for linear correlator
params.n_steps           = 100;
% for multi tau
params.n_octaves         = 20;
params.n_bins_per_octave = 8;

params = omex_read_params(params, varargin);

if     strcmp(params.type,'multi')
    params.n_corrtype = 0;
elseif strcmp(params.type,'linear')
    params.n_corrtype = 1;
else
    error('unknown correlator type!')
end

% call the s_mex_simplex_init gateway function. Here we pass the
% parameter array really as a struct
handle = omexFcs(1, params);

end
