P = [-6 -2 1 1 -1 -2 -1; 0 1 2 3 4 5 6];
n = 7; # 3 * (number of cubic curves) + 1; in this case c = 2

M = bezier_creator_matrix(n);

t = [0:1/(n-1):1];
T = power(ones(n) .* t, ones(n) .* [0:n-1]');

C = P*inv(M*T);

[Q60, R6] = casteljau(C, t(2));

s0 = bezier_slope(Q60, 0)
s1 = bezier_slope(Q60, 0.5)
s2 = bezier_slope(Q60, 1)

M3 = bezier_creator_matrix(3);
k = [k_a 0; 0 k_b].' * [1 0; -1 -1; 0 1]' * M3 * power(0.5, [0:2]');

s = [s0, s2];
Q3 = [Q60(:,1), [0 0; 0 0], Q60(:,7)];

# s1 = C'(0.5)
# s1 / n = [s0, s2] * [k_a 0; 0 k_b].' * [1 0; -1 -1; 0 1]' * M3 * power(0.5, [0:2]');
#K = -4 * inv(s) * (s1./n)

# subs([q30x s0x s2x q33x; q30y s0y s2y q33y] * [1 1 z 0; 0 k_a 0 0; 0 0 -k_b 0; 0 z 1 1], z, 0) * [-1 0 0; 1 -1 0; 0 1 -1; 0 0 1] * M3 * T
# s1/n*4 + 2*Q3(:,1) - 2 * Q3(:,3) = [s0, s2] * [k_a 0; 0 k_b]
K = inv(s) * (s1/n*4 + 2*Q3(:,1) - 2 * Q3(:,4))

a = -K(1);
b = -K(2);


Q3(:,2) = Q3(:,1) + a * s0;
Q3(:,3) = Q3(:,4) - b * s2

bezier_eval(Q3, bezier_creator_matrix(4), 0.5)
bezier_eval(Q60, M, 0.5)

norm(e)




