function B = bezier_reduce(C)
  n = columns(C);
  M = zeros(n,n-1);
  M(2:n-1,1:n-2) += ([1:n-2]' / (n-1)) .* eye(n-2);
  M(2:n-1,2:n-1) += (1 - [1:n-2]' / (n-1)) .* eye(n-2);
  M(1,1) = 1;
  M(n,n-1) = 1;
  B = (inv(M' * M) * M' * C')';
end