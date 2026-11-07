%RNFFTMEX Gateway function to RNFFT module from NFFT3
%   Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
%
%   The NFFT of real-valued functions. The commands are those of nfftmex.
%   f is real, f_hat holds the N_total coefficients of the half spectrum
%   0 <= k(d) <= N(d)/2, -N(t)/2 <= k(t) < N(t)/2 for t < d, the last index
%   running fastest, and 'solver' runs the c2r solver.
%
%   p = rnfftmex('init_guru',{d,N(1),...,N(d),M,n(1),...,n(d),m,flags,fftw_flags})
%   rnfftmex('set_x',p,x)          % x is d x M in [-1/2,1/2)
%   rnfftmex('precompute_psi',p)
%   rnfftmex('set_f_hat',p,fhat)   % fhat is N_total x 1
%   rnfftmex('trafo',p)            % f_j = sum_k c_k Re(fhat_k exp(-2 pi i k x_j))
%   f = rnfftmex('get_f',p)        % c_k = 1 for k(d) in {0, N(d)/2}, else 2
%   rnfftmex('set_f',p,f)
%   rnfftmex('adjoint',p)          % fhat_k = c_k sum_j f_j exp(2 pi i k x_j)
%   fhat = rnfftmex('get_f_hat',p)
%   rnfftmex('finalize',p)
