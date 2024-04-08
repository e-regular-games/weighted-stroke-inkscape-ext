function C = bezier3_match_slopes(s0, s3, p0, p1, p2, p3)
  # p1 and p2 must be at times t1=1/3 and t2=2/3  
  # the following equation is a formulation of the bezier curve of the derivative
  # at t=0.5. It follows the constraint that the C0 and C1 will remain on the same line.
  # ie the unit vector of C1 - C0 will be the same. This applies for C3 - C2 as well.
  # Note: The final control points of the curve in terms of these multipliers, k_a and k_b
  # is as follows: C = subs([q30x s0x s2x q33x; q30y s0y s2y q33y] * [1 1 z 0; 0 k_a 0 0; 0 0 -k_b 0; 0 z 1 1], z, 0)
  # Evaluating at t=1/3 and t=2/3 to match the points at Q60
  # S = [s0x s2x; s0y s2y];
  # 27 * [Q60(1/3), Q60(2/3)] = 27 * C * M4 * power([1/3, 2/3], [0:3]')
  # 27 * [Q60(1/3), Q60(2/3)] - [Q3(:,1), Q3(:,4)] * [20 7; 7 20] = 27 * C * M4 * power([1/3, 2/3], [0:3]') - [Q3(:,1), Q3(:,4)] * [20 7; 7 20]
  # 27 * [Q60(1/3), Q60(2/3)] - [Q3(:,1), Q3(:,4)] * [20 7; 7 20] = S * [k_a 0; 0 k_b] * [12 6; -6 -12]
  # inv(S) * (27 * [Q60(1/3), Q60(2/3)] - [Q3(:,1), Q3(:,4)] * [20 7; 7 20]) * inv([12 6; -6 -12]) = [k_a 0; 0 k_b]
  # s1 + 2*Q3(:,1) - 2 * Q3(:,3) = [s0, s2] * [-k_a 0; 0 -k_b]
  
  s = [s0, s3];
  K = inv(s) * (27 * [p1, p2] - [p0, p3] * [20 7; 7 20]) * inv([12 6; -6 -12])
  a = K(1,1);
  b = K(2,2);

  C = [p0, p0 + a * s0, p3 - b * s3, p3];
end