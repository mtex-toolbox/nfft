% smoke test of a freshly built rnfftmex: fast against direct transforms, real
% samples, the threaded adjoint against the serial one and the c2r solver;
% errors make the CI step fail
fprintf('rnfftmex: %s\nthreads: %d\n', which('rnfftmex'), rnfftmex('get_num_threads'));
rnfftmex('set_num_threads', 4);
fprintf('threads set to %d\n', rnfftmex('get_num_threads'));
rng(1);
for d = 1:3
  N = [16 12 10]; N = N(1:d); M = 200;
  x = rand(d,M) - 0.5;
  p = rnfftmex('init_guru', [{d}, num2cell(N), {M}, num2cell(2*N), {4, 1+2^10+2^12, int8(64)}]);
  rnfftmex('set_x', p, x);
  rnfftmex('precompute_psi', p);
  fh = rnfftmex('get_f_hat', p);
  fh = complex(randn(size(fh)), randn(size(fh)));
  rnfftmex('set_f_hat', p, fh);
  rnfftmex('trafo', p);         f1 = rnfftmex('get_f', p);
  rnfftmex('trafo_direct', p);  f2 = rnfftmex('get_f', p);
  assert(isreal(f1) && numel(f1) == M, 'rnfftmex: %dd samples are not M real numbers', d);
  rnfftmex('set_f', p, f2);
  rnfftmex('adjoint', p);       g1 = rnfftmex('get_f_hat', p);
  rnfftmex('adjoint_direct', p); g2 = rnfftmex('get_f_hat', p);
  rnfftmex('set_f_hat', p, zeros(size(fh)));
  rnfftmex('solver', p, 20, 0);
  rnfftmex('trafo', p); r = norm(rnfftmex('get_f', p) - f2)/norm(f2);
  rnfftmex('finalize', p);
  e = [norm(f1-f2)/norm(f2), norm(g1-g2)/norm(g2)];
  fprintf('%dd trafo error %.1e, adjoint error %.1e, solver residual %.1e\n', d, e, r);
  assert(all(e < 1e-6), 'rnfftmex: %dd transform does not match the direct one', d);
  assert(r < 0.1, 'rnfftmex: the %dd solver does not reduce the residual', d);
end
% lattices large enough for the parallel adjoint: 4 threads equal 1
for d = 2:3
  N = [64 64 32]; N = N(1:d); M = 20000;
  x = rand(d,M) - 0.5; f = rand(M,1);
  g = cell(1,2); nt = [1 4];
  for k = 1:2
    rnfftmex('set_num_threads', nt(k));
    p = rnfftmex('init_guru', [{d}, num2cell(N), {M}, num2cell(2*N), {4, 1+2^9+2^10+2^12, int8(64)}]);
    rnfftmex('set_x', p, x);
    rnfftmex('precompute_psi', p);
    rnfftmex('set_f', p, f);
    rnfftmex('adjoint', p);
    g{k} = rnfftmex('get_f_hat', p);
    rnfftmex('finalize', p);
  end
  e = norm(g{1}-g{2})/norm(g{1});
  fprintf('%dd adjoint, 4 against 1 thread %.1e\n', d, e);
  assert(e < 1e-12, 'rnfftmex: the parallel %dd adjoint differs from the serial one', d);
end
disp('rnfftmex smoke test passed');
