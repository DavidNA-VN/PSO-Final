%% Simulation
function results = simulation(params, output_filename)

    save_filename = output_filename;
    results = {};

    % Parfor Supported loop
    parfor rep = 1:params.repetitions
        fprintf('\n Repetition %i:', rep);

        % ------------------ Generate geometry ------------------
        [UE_pos, AP_pos, target_pos] = generate_positions(params.T, ...
            params.U, params.M_t, params.geo.line_length, ...
            params.geo.target_y, params.geo.UE_y, ...
            params.geo.min_dist, params.geo.max_dist);

        results{rep}.P_comm_ratio = params.P_comm_ratio;
        results{rep}.AP = AP_pos;
        results{rep}.UE = UE_pos;
        results{rep}.Target = target_pos;

        % ------------------ Channel generation ------------------
        H_comm = LOS_channel(AP_pos, UE_pos, params.N_t);

        % ------------------ Sensing beamsteering ------------------
        [sensing_angle, ~] = compute_angle_dist(AP_pos, target_pos);
        sensing_beamsteering = beamsteering(sensing_angle.', params.N_t); % (Target x M_t x N_t)

        % Conjugate BF - power-normalized beamsteering (norm only)
        F_sensing_CB_norm = sensing_beamsteering * sqrt(1 / params.N_t);

        % Nullspace BF - power-normalized (norm only)
        F_sensing_NS_norm = beam_nulling(H_comm, sensing_beamsteering);

        % ------------------ Sweep rho ------------------
        for p_i = 1:length(params.P_comm_ratio)

            % Power ratio of comm and sensing (per AP)
            P_comm    = params.P * params.P_comm_ratio(p_i);
            P_sensing = params.P * (1 - params.P_comm_ratio(p_i));

            % Scale sensing beams with sensing power
            F_sensing_CB = F_sensing_CB_norm * sqrt(P_sensing);
            F_sensing_NS = F_sensing_NS_norm * sqrt(P_sensing);

            solution_counter = 1;

            % ============================================================
            % 1) Baseline: NS Sensing + RZF Comm
            % ============================================================
            F_star_RZF = beam_regularized_zeroforcing(H_comm, P_comm, params.sigmasq_ue) * sqrt(P_comm);

            results{rep}.power{p_i}{solution_counter} = compute_metrics( ...
                H_comm, F_star_RZF, params.sigmasq_ue, sensing_beamsteering, ...
                F_sensing_NS, params.sigmasq_radar_rcs);

            results{rep}.power{p_i}{solution_counter}.name = 'NS+RZF';
            solution_counter = solution_counter + 1;

            % ============================================================
            % 2) Low-dimension PSO for sensing only (Comm fixed = RZF)
            %    - PSO optimizes complex weights per AP applied to base sensing beam
            %    - Keeps SINR above a relaxed threshold to show trade-off
            % ============================================================
            SINR_RZF = compute_SINR(H_comm, F_star_RZF, F_sensing_NS, params.sigmasq_ue);

            % Relaxed SINR threshold so PSO has room to improve sensing
            % (Try 0.6 ~ 0.9 if you want to tune the curve)
            gamma_ref = 0.7 * min(SINR_RZF);

            sens_streams = 1;

            % PSO returns optimized sensing beam; Comm stays fixed (RZF)
            [F_jsc_sensing, feasible, best_SSNR, best_fit] = opt_sens_PSO_lowdim( ...
                H_comm, F_star_RZF, params.sigmasq_ue, sensing_beamsteering, ...
                F_sensing_NS, params.sigmasq_radar_rcs, gamma_ref, P_sensing);

            F_jsc_comm = F_star_RZF;

            results{rep}.power{p_i}{solution_counter} = compute_metrics( ...
                H_comm, F_jsc_comm, params.sigmasq_ue, sensing_beamsteering, ...
                F_jsc_sensing, params.sigmasq_radar_rcs);

            results{rep}.power{p_i}{solution_counter}.feasible = feasible;
            results{rep}.power{p_i}{solution_counter}.name = strcat('JSC+PSO-LD', num2str(sens_streams));
            results{rep}.power{p_i}{solution_counter}.SSNR_opt = best_SSNR;
            results{rep}.power{p_i}{solution_counter}.best_fit = best_fit;
            results{rep}.power{p_i}{solution_counter}.gamma_ref = gamma_ref;

            solution_counter = solution_counter + 1;

        end % for p_i
    end % parfor rep

    % ------------------ Save Results ------------------
    output_folder = './output/';
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    save(strcat(output_folder, save_filename, '.mat'));

end
