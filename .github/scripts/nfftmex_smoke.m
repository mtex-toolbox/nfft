% smoke test of a freshly built nfftmex: fast against direct transforms, and
% the half spectrum of a real plan; errors make the CI step fail
fprintf('nfftmex: %s\nthreads: %d\n', which('nfftmex'), nfftmex('get_num_threads'));
% more threads than the runner's default, so the parallel code runs
nfftmex('set_num_threads', 4);
fprintf('threads set to %d\n', nfftmex('get_num_threads'));
rng(1);
for d = 2:3
  N = 8*ones(1,d); M = 60;
  x = rand(d,M) - 0.5;
  fh = complex(randn(prod(N),1), randn(prod(N),1));
  p = nfftmex('init_guru', [{d}, num2cell(N), {M}, num2cell(2*N), {4, 1+2^10+2^12, int8(64)}]);
  nfftmex('set_x', p, x);
  nfftmex('precompute_psi', p);
  nfftmex('set_f_hat', p, fh);
  nfftmex('trafo', p);         f1 = nfftmex('get_f', p);
  nfftmex('trafo_direct', p);  f2 = nfftmex('get_f', p);
  nfftmex('set_f', p, f2);
  nfftmex('adjoint', p);       g1 = nfftmex('get_f_hat', p);
  nfftmex('adjoint_direct', p); g2 = nfftmex('get_f_hat', p);
  nfftmex('finalize', p);
  e = [norm(f1-f2)/norm(f2), norm(g1-g2)/norm(g2)];
  fprintf('%dd trafo error %.1e, adjoint error %.1e\n', d, e);
  assert(all(e < 1e-6), 'nfftmex: %dd transform does not match the direct one', d);
end
% lattices large enough for the tiled parallel adjoint: 4 threads equal 1
for d = 2:3
  N = [64 64 32]; N = N(1:d); M = 20000;
  x = rand(d,M) - 0.5; f = complex(rand(M,1), rand(M,1));
  g = cell(1,2); nt = [1 4];
  for k = 1:2
    nfftmex('set_num_threads', nt(k));
    p = nfftmex('init_guru', [{d}, num2cell(N), {M}, num2cell(2*N), {4, 1+2^10+2^12, int8(64)}]);
    nfftmex('set_x', p, x);
    nfftmex('precompute_psi', p);
    nfftmex('set_f', p, f);
    nfftmex('adjoint', p);
    g{k} = nfftmex('get_f_hat', p);
    nfftmex('finalize', p);
  end
  e = norm(g{1}-g{2})/norm(g{1});
  fprintf('%dd adjoint, 4 against 1 thread %.1e\n', d, e);
  assert(e < 1e-12, 'nfftmex: the parallel %dd adjoint differs from the serial one', d);
end
p = nfftmex('init_guru', {2, 4, 4, 1, 8, 8, 2, 1+2^10+2^14, int8(64)});
n = numel(nfftmex('get_f_hat', p));
nfftmex('finalize', p);
assert(n == 8, 'nfftmex: a real plan holds %d coefficients, not the half spectrum 8', n);
disp('nfftmex smoke test passed');
