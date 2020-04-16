%% Homework 3 - Power Systems Protections
%close all; clear all; clc;

%% Task 1: Read data
% t_data = readtable('data/ECE6323_PROJECT_02_CTCHANNEL_EVENT_01_MAIN.dat');
% cfg_id = fopen('data/ECE6323_PROJECT_02_CTCHANNEL_EVENT_01_MAIN.cfg');
% cfg = textscan(cfg_id, '%s', 'delimiter', '\n');
% 
% config_data = textscan(cfg{1,1}{3,1}, '%s', 'Delimiter', ',');
% sample_data = textscan(cfg{1,1}{6,1}, '%s', 'Delimiter', ',');
% 
% % Data from cfg file
% v_base = str2double(config_data{1,1}{6,1});%v_base = 0.000038357;
% v_offset = str2double(config_data{1,1}{7,1}); %-0.204407170;
% %c_base = 3.4529695416e-002; There is no current to read
% sampling_rate = str2double(sample_data{1,1}{1,1}); 
% number_samples = str2int(sample_data{1,1}{2,1});
% frequency = str2double(cfg{1,1}{4,1});



load('data/saved_data.mat')
n_samples_cycle = sampling_rate/frequency; % per cycle
time_base = 1/sampling_rate;
period = 1/frequency;
factor_samples2time = period/n_samples_cycle;
omega = 2*pi*frequency;

% Extract data 
time = t_data{:,2}/1e6; %in microseconds
voltage = t_data{:,3}*v_base + v_offset;




%% Task 2: Develop the mathematical model of the instrumentation channel
n = 400;
g_m = 0.001;
L_1 = 26.526e-6;
L_2 = 348.0e-6;
L_3 = 348.0e-6;
M_23 = 287.0e-6;
g_s1 = 1.9635;
g_s2 = 0.1497;
g_s3 = 0.1497;
r_1 = 0.005;
r_2 = 0.4469;
r_3 = 0.4469;
R_b = 0.1;
g_b = 1/R_b;
l_0 = 0.1876;
i_0 = 6.09109;
L_0 = 2.36;
h = 1/sampling_rate;
params = [n g_m L_1 L_2 L_3 M_23 g_s1 g_s2 g_s3 r_1 r_2 r_3 R_b g_b l_0 i_0 L_0 h sampling_rate/frequency omega];

% Get useful phasors
%v_out = get_phasors(voltage,time, params);
v_out = voltage;
imb = -v_out / R_b;
% figure;
% plot(time,voltage);
% xlabel('Time (s)');ylabel('Voltage (V)');title('Voltage data record');
% figure;
% plot(v_out);
% %xlabel('Time (s)');º
% ylabel('Vout data phasor magnitude(V)');title('Voltage data record');

%% Task 3: Develop the equations of all measurements [Done in functions file]
% Though we can obtain matrix W for the algorithm

% Create vector of deviations
%sigma = [s_eq1 s_eq2 s_eq3 ... s_eq21]
sigmas = [0.005 0.005 0.005 0.005 0.0005 ...
          0.0005 0.0005 0.005 0.005 0.0005 ...
          0.00005 0.00005 0.00005 0.00005 0.0003 ...
          0.05 0.05 0.05 0.05 0.05 1];
s = size(sigmas);
W = zeros(s(2),s(2));
for i=1:s(2)
    W(i,i) = 1/(sigmas(i))^2;
end

%% Task 4: Develop the algorithm to solve the dynamic state estimation

%%% Create matrixes of state vectors that will be recorded
s = size(v_out);
X = zeros(15,s(2));
Z = zeros(21,s(2));
p = zeros(1,s(2));
NMED = zeros(21,s(2));

% Get all the parameters, matrices and variables needed for the algorithm
e = 1e-4;
deg_free = 21 - 15;

if s(1)>s(2); s = s(1); else; s = s(2); end;
progressbar('Performing Dynamic State Estimation');
for i=2:s
    % X(:,i)
    %fprintf('Iteracion: %i\n',i)
    
    % First iteration
    x_aux = zeros(15,1);
    x_aux_prev = ones(15,1);
    x_aux_next = ones(15,1);
    
    % Define measures vector for this iteration
    z = zeros(21,1);
    z(1) = v_out(i);
    z(16:20) = imb(i);
    
    % Create measures vector-h needed for get_measures_eqs
    measures_h = [v_out(i-1), imb(i-1)];
    
    while sum(abs(x_aux - x_aux_prev)) > e
        %disp(abs(x_aux_next - x_aux))
        % Evaluate H at current point
        H = get_H_v2(x_aux, X(:,i-1), params);
        m = get_measures_eqs_v2(x_aux, X(:,i-1), params, measures_h);
        matrix = inv(H' * W * H) * H' * W;
        % Apply iteration computations
        x_aux_next = x_aux + matrix * (z - m);
        x_aux_prev = x_aux;
        x_aux = x_aux_next;
        
    end
    
    X(:,i) = x_aux;
    Z(:,i) = get_measures_eqs_v2(x_aux, X(:,i-1), params, measures_h);
    % Perform Chi-test
    nmed = (z - get_measures_eqs_v2(x_aux, X(:,i-1), params, measures_h));
    NMED(:,i) = nmed;
    seta = 0;
    for measure=1:21
        seta = seta + nmed(measure)^2*W(measure,measure);
    end
    
    % Get confidence level
    p(i) = 1 - chi2cdf(seta,deg_free);
    
    % Perform 2nd Chi-Test
    i_est = X(11,:)';
    i_med = -imb*n;
    
    progressbar(i/s);
end

save('data/computed_estimation.mat','X','Z','time','p');


%% FINAL ANALYSIS
% for i=1:15
%     figure;
%     plot(X(i,:))
%     ylabel('Estimation');title('Estimation');
% end
%%
%figure; plot(X(11,:)'); ylabel('Primary current Estimation');title('Primary current Estimation');
%hold on 
%plot(-imb*n); plot(X(11,:)' - (-imb*n)); legend('Estimated current','Measured current', 'Difference')
%figure; plot(p); ylabel('Probability');title('Confidence level');

%%
%figure; plot(((X(11,:)')/n)); hold on; plot((-imb));title('Difference of measured and estimated current'); legend('Estimated','Measured')
%figure; plot(X(6,:));title('Flux linkage');
%%
%figure; plot(((X(11,:)')/n - (-imb)).^2./(.05^2));title('Residue');
% figure; 
% for i=1:15
% plot(X(i,:));title('State'); hold on;
% end;
% legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15')
% %% 
% for i=1:21
%     figure; plot(NMED(i,:)');
% end
% 




