
% Function to get the phasor from sinusoidal 
function v_out = get_phasors(v_int, time_int, params)

n_samples_cycle = params(19);
omega = params(20);
size_v = size(v_int);

for n=1:size_v-n_samples_cycle+1
    % Define vectors
    vec_time = time_int(n:n+n_samples_cycle-1);
    vec_cos = cos(omega*vec_time);
    vec_sin = sin(omega*vec_time);
    vec_v_values = (v_int(n:n+n_samples_cycle-1))';
    %vec_i_values = (c_int(n:n+n_samples_cycle-1))';
    
    V(n) = (vec_v_values*vec_cos+j*vec_v_values*vec_sin)*sqrt(2)/n_samples_cycle;
    %I(n) = (vec_i_values*vec_cos+j*vec_i_values*vec_sin)*sqrt(2)/n_samples_cycle;
end

v_out = abs(V);

end 