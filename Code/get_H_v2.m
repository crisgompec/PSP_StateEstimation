
% Function to get matrix H evaluated at state x
% This version 2 takes out the integral for equations 1, 11-16, 20, 21


function H = get_H_v2(x, x_ans, params, vouth, ibmh)

% Prepare parameters
n = params(1);
g_m = params(2);
L_1 = params(2);
L_2 = params(4);
L_3 = params(5);
M_23 = params(6);
g_s1 = params(7);
g_s2 = params(8);
g_s3 = params(9);
r_1 = params(10);
r_2 = params(11);
r_3 = params(12);
R_b = params(13);
g_b = params(14);
l_0 = params(15);
i_0 = params(16);
L_0 = params(17);
h = params(18);

% Initial values required for H
% x = [v1 v2 v3 v4 e lambda y1 y2 y3 y4 ip im iL1 iL2 iL3]
lambda = x(6);
lambda_ans = x_ans(6);

y1 = x(7);
y1_ans = x_ans(7);

y2 = x(8);
y2_ans = x_ans(8);

y3 = x(9);
y3_ans = x_ans(9);

y4 = x(10);
y4_ans = x_ans(10);

%%% Create matrix H
H = zeros(15,21);

% Equation 1 - Affected
H(3,1) = 1;
H(4,1) = -1;


% Equation 2
H(5,2) = -g_m*h/2;
H(11,2) = h/2/n;
H(12,2) = -h/2;
H(13,2) = g_s1*L_1 + h/2;

% Equation 3
H(5,3) = g_m*h/2;
H(11,3) = -h/2/n;
H(12,3) = h/2;
H(14,3) = -g_s2*L_2-h/2;
H(15,3) = g_s2*M_23;

% Equation 4
H(5,4) = -g_m*h/2;
H(11,4) = h/2/n;
H(12,4) = -h/2;
H(14,4) = -M_23*g_s3;
H(15,4) = g_s3*L_3 + h/2;

% Equation 5
H(1,5) = - h/2;
H(2,5) = h/2;
H(5,5) = h/2 + r_1*h/2*g_m;
H(11,5) = - r_1*h/2/n;
H(12,5) = r_1*h/2;
H(13,5) = L_1;

% Equation 6
H(1,6) = h/2;
H(3,6) = - h/2;
H(14,6) = r_2 * (h/2 + g_s2*L_2) + L_2;
H(15,6) = - r_2 * g_s2 * M_23 - M_23;

% Equation 7
H(2,7) = -h/2;
H(4,7) = h/2;
H(14,7) = - r_3 * g_s3 * M_23 - M_23;
H(15,7) = r_3 * (h/2 + g_s3 * L_3) + L_3;

% Equation 8
H(3,8) = g_b*h/2;
H(4,8) = - g_b*h/2;
H(14,8) = h/2 + g_s2 * L_2;
H(15,8) = - g_s2 * M_23;

% Equation 9
H(3,9) = - h/2 * g_b;
H(4,9) = h/2 * g_b;
H(14,9) = g_s3 * M_23;
H(15,9) = - h/2 - g_s3 * L_3;

% Equation 10
H(5,10) = h/2;
H(6,10) = -1;

% Equation 11 - Affected
H(6,11) = -2/(l_0^2) * lambda;
H(7,11) = 1;

% Equation 12 - Affected
H(7,12) = -2*y1;
H(8,12) = 1;

% Equation 13 - Affected
H(8,13) = -2*y2;
H(9,13) = 1;

% Equation 14 - Affected
H(7,14) = -y3;
H(9,14) = -y1;
H(10,14) = 1;

% Equation 15 - Affected
H(6,15) = -i_0/l_0*y4 - 1/L_0;
H(10,15) = - i_0/l_0 * (lambda);
H(12,15) = 1;

% Equation 16 - Affected
H(3,16) = -g_b;
H(4,16) = g_b;

% Equation 17 - Affected by left-hand side integration
H(13,17) = h/2 + g_s1 * L_1;
H(13,17) = H(13,17) * 2/h;

% Equation 18 - Affected by left-hand side integration
H(14,18) = h/2 + g_s2 * L_2;
H(14,18) = H(14,18) * 2/h;
H(15,18) = - g_s2 * M_23;
H(15,18) = H(15,18) * 2/h;


% Equation 19 - Affected by left-hand side integration
H(14,19) = g_s3 * M_23;
H(14,19) = H(14,19) * 2/h;
H(15,19) = - h/2 - g_s3 * L_3;
H(15,19) = H(15,19) * 2/h;


% Equation 20 - Affected
H(5,20) = g_m;
H(11,20) = -1/n;
H(12,20) = 1;

% Equation 21 - Affected
H(4,21) = 1;

% Now we transpose the matrix
H = H';

end