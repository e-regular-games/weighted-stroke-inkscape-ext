function C = bezier_solver_by_derivative(P)
  # find a bezier curve of order 3 to represent the curve between P0 and P1.
  n = columns(P);
  M = bezier_creator_matrix(n);

  t = [0:1/(n-1):1];
  T = power(ones(n) .* t, ones(n) .* [0:n-1]');
  C = P*inv(M*T);

  [Q60, R6] = casteljau(C, t(2));

  s0 = bezier_slope(Q60, 0);
  s1 = bezier_slope(Q60, 0.5);
  s2 = bezier_slope(Q60, 1);

  M3 = bezier_creator_matrix(3);
  s = [s0, s2];
  Q3 = [Q60(:,1), [0 0; 0 0], Q60(:,7)];

  # the following equation is a formulation of the bezier curve of the derivative
  # at t=0.5. It follows the constraint that the C0 and C1 will remain on the same line.
  # ie the unit vector of C1 - C0 will be the same. This applies for C3 - C2 as well.
  # Note: The final control points of the curve in terms of these multipliers, k_a and k_b
  # is as follows: C = subs([q30x s0x s2x q33x; q30y s0y s2y q33y] * [1 1 z 0; 0 k_a 0 0; 0 0 -k_b 0; 0 z 1 1], z, 0)
  # The derivative curve is then D = n * C * [-1 0 0; 1 -1 0; 0 1 -1; 0 0 1]
  # Evaluating at 0.5 to match the slope of Q60 at 0.5
  # s1 = D * M3 * power(0.5, [0:2]')
  # s1 + 2*Q3(:,1) - 2 * Q3(:,3) = [s0, s2] * [-k_a 0; 0 -k_b]
  K = inv(s) * (s1 + 2*Q3(:,1) - 2 * Q3(:,4));
  a = -K(1);
  b = -K(2);

  Q3(:,2) = Q3(:,1) + a * s0;
  Q3(:,3) = Q3(:,4) - b * s2;

  C = Q3;
end