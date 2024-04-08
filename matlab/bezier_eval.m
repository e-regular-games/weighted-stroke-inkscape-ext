function F = bezier_eval(C, t)
  M = bezier_creator_matrix(columns(C));
  F = C * M * power(t, 0:columns(C)-1)';
end