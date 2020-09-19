
%%% AVAILABLE FUNCTIONS
% get_measures_eqs(x, x_ans, params)
% get_H(x, x_ans, params)
% get_phasors(v_int)
% get_measures_values(v,i)

% Function to get measure vector from the state
function m = get_measures_eqs(x, x_ans, params, measures_h)
% Task 3: Develop the equations of all measurements, actual, derived, virtual 
% and pseudo as a function of the state of the system. 

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

vouth = measures_h(1);
ibmh = measures_h(2);

% x = [v1 v2 v3 v4 e lambda y1 y2 y3 y4 ip im iL1 iL2 iL3]
v1 = x(1); v2 = x(2); v3 = x(3); v4 = x(4); e = x(5); lambda = x(6);
y1 = x(7); y2 = x(8); y3 = x(9); y4 = x(10); ip = x(11); im = x(12);
iL1 = x(13); iL2 = x(14); iL3 = x(15);

%%% Create previous state variables vector x 
v1h = x_ans(1); v2h = x_ans(2); v3h = x_ans(3); v4h = x_ans(4);
eh = x_ans(5); lambdah = x_ans(6); y1h = x_ans(7); y2h = x_ans(8);
y3h = x_ans(9); y4h = x_ans(10); iph = x_ans(11); imh = x_ans(12);
iL1h = x_ans(13); iL2h = x_ans(14); iL3h = x_ans(15);

%%% Create Measure estimation from state variables vector - h
m = zeros(21,1);

m(1) = h/2 * (v3h + v3) + h/2 * (v4h + v4);

m(2) = -g_m * h/2 * (eh + e) - h/2 * (imh + im) + h/2/n * (iph + ip) + ...
    g_s1 * L_1 * (iL1 - iL1h) + h/2 * (iL1+iL1h);

m(3) = g_m * h/2 * (eh + e) + h/2 * (imh + im) - h/2/n * (iph + ip) - ...
    g_s2 * (L_2 * (iL2- iL2h) - M_23 * (iL3 - iL3h)) - h/2 * (iL2 + iL2h);

m(4) = - g_m * h/2 * (eh + e) - h/2 * (imh + im) + h/2/n * (iph+ip) + ...
    g_s3 * (L_3 * (iL3 - iL3h) - M_23 * (iL2-iL2h)) + h/2 * (iL3 + iL3h);

m(5) = -h/2 * (v1h + v1) + h/2 * (v2h + v2) + h/2 * (eh + e) + L_1 * (iL1 - iL1h) ...
    + r_1 * h/2 * (g_m * (eh + e) + imh + im - 1/n * (iph + ip));

m(6) = - h/2 * (v3h + v3) + h/2 * (v1h + v1) + r_2 * (h/2 * (iL2h + iL2) + ...
    g_s2 * (L_2 * (iL2 - iL2h) - M_23 * (iL3 - iL3h)) + L_2 * (iL2 - iL2h) - ...
    M_23 * (iL3 - iL3h));

m(7) = -h/2 * (v2h + v2) + h/2 * (v4h + v4) + r_3 * (h/2 * (iL3h + iL3) + ...
    g_s3 * (L_3 * (iL3 - iL3h) - M_23 * (iL2 - iL2h) ) ) + L_3 * (iL3 - iL3h) ...
    - M_23 * (iL2 - iL2h);

m(8) = h/2 * (iL2h + iL2) + g_s2 * (L_2 * (iL2 - iL2h) - M_23 * (iL3 - iL3h)) ...
    + g_b * h/2 * (v3h + v3 - (v4h + v4));

m(9) = -h/2 * (iL3h + iL3) - g_s3 * (L_3 * (iL3 - iL3h) - M_23 * (iL2 - iL2h)) ...
    + g_b * (v4h + v4 - (v3h + v3)) * h/2;

m(10) = h/2 * (eh + e) - (lambda - lambdah);

m(11) = h/2 * (y1h - y1) - h/3/l_0^2 * (lambdah*lambda + lambdah^2 + lambda^2);

m(12) = h/2 * (y2h + y2) - h/3 * (y1h*y1 + y1^2 + y1h^2);

m(13) = h/2 * (y3h - y3) - h/3 * (y2h*y2 + y2h^2 + y2^2);

m(14) = h/2 * (y4h + y4) - h/3 * (y3h*y1h + y3h*y1/2 + y3*y1h/2 + y1*y3);

m(15) = h/2 * (imh + im) - i_0*h/l_0/3 * (lambdah*y4h + lambdah*y4/2 + lambda*y4h/2 + ...
    lambda*y4) - h/2/L_0 * (lambdah + lambda);

m(16) = - g_b * ((v3h - v3) - (v4h + v4)) * h/2;

m(17) = h/2 * (iL1h + iL1) + g_s1 * L_1 * (iL1 - iL1h);

m(18) = h/2 * (iL2h + iL2) + g_s2 * (L_2 * (iL2 - iL2h) - M_23 * (iL3 - iL3h));

m(19) = -h/2 * (iL3h + iL3) - g_s3 * (L_3 * (iL3 - iL3h) - M_23 * (iL2 - iL2h));

m(20) = g_m * (e + eh) * h/2 + h/2 * (imh + im) - h/2/n * (iph + ip);

m(21) = h/2 * (v4h + v4);


% Equations 1 and 16-20 are affected by the integration in the left side.
m(1) = m(1) * 2/h - vouth;

m(16) = m(16) * 2/h - ibmh;
m(17) = m(17) * 2/h - ibmh;
m(18) = m(18) * 2/h - ibmh;
m(19) = m(19) * 2/h - ibmh;
m(20) = m(20) * 2/h - ibmh;


end







