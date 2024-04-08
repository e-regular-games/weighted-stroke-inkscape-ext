function C = bezier_from_points(P)
  n = columns(P);
  M = bezier_creator_matrix(n);
  T = power([0:1/(n-1):1], [0:n-1]');
  rats(inv(M*T))
  C = P*inv(M*T);
end