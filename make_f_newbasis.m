function [u_coeff] = make_f_newbasis(f, Nres,  centers, r, N_multi, k)
  assert(isa(f,'function_handle'));
  
  acc = 200; %nr of random points to estimate difference
  thresh = 10; %optimisation steps until warning
  expansion_fit_error = 0.1; %desired 
  rand_err = 1;
  fu = @(u, x) expansion_eval_newbasis(u, Nres, N_multi, centers, x, k);
  
  steps = 0;
  u_coeff = ones(Nres*N_multi));
  p_test_new = (rand(3, acc) - 0.5.*ones(3, acc)).*2; %random testing points in [-1, 1]^3
  while rand_err > expansion_fit_error
    p_test = p_test_new;
    rand_err = sum_square_diff(u_coeff, fu, f, acc, p_test)
    rand_eval_err = @(u) (sum_square_diff(u, fu, f, acc, p_test));

    u_coeff_new = fminsearch(rand_eval_err, u_coeff);
    
    p_test_new = (rand(3, acc) - 0.5.*ones(3, acc)).*2; %tested on different points than optimised
    rand_err_new = sum_square_diff(u_coeff_new, fu, f, acc, p_test_new);
    if rand_err_new < rand_err
      u_coeff = u_coeff_new;
      rand_err = rand_err_new;
    end
    
    steps = steps + 1;
    if steps = thresh
      printf('expansion fit is taking more than' + str(thresh) + 'steps');
    end
  end
end

function [s] = sum_square_diff(u, fu, f, acc, p_test)
  s = 0;
  for i = 1:acc
    s = s + (fu(u, p_test(i,:))  - f(p_test(i,:)))^2;
  end
end
