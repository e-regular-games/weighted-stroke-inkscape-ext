function s = bezier_slope(C, t)
  n = columns(C) - 1;
  M = bezier_creator_matrix(n);
  Q = (n+1) * (C(:,2:n+1) - C(:,1:n));
  s = Q * M * power(t, [0:n-1]');
end