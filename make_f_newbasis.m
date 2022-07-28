function [u] = make_f_newbasis(f, Nres,  centers, r, N_multi, k)
  fu  = @(u, x) expansion_eval_newbasis(u, Nres, N_multi, centers, x, k)
  err = integral(expansion_eval_newbasis - f(x), )
  
  x0 = ones(Nres*N_multi))
  u  = fminsearch(err, x0)
end
