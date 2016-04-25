function add_noise_to_W(ESN)
% add noise to W
    ESN.W = ESN.W + abs(ESN.W) .* (ESN.WparameterNL * (rand(ESN.DR.size,ESN.DR.size)-0.5));
end