function [F_comm_best, F_sensing_best, feasible_best, SSNR_best, J_best] = opt_jsc_PSO( ...
    H_comm, sigmasq_ue, gamma, sensing_beamsteering, sensing_streams, sigmasq_radar_rcs, P_all)
% opt_jsc_PSO: PSO solver replacing opt_jsc_SDP (SDP/CVX).
% Maximizes sensing SSNR subject to SINR>=gamma and per-AP power<=P_all.
%
% Inputs match opt_jsc_SDP usage in simulation.m:
%   H_comm: U x M x N
%   sigmasq_ue: scalar
%   gamma: scalar (SINR threshold)
%   sensing_beamsteering: (Tgt=1) x M x N  (in your sim, T=1 target)
%   sensing_streams: integer >=0
%   sigmasq_radar_rcs: scalar
%   P_all: scalar (power per AP)
%
% Outputs:
%   F_comm_best: U x M x N
%   F_sensing_best: sensing_streams x M x N (or 1xMxN if streams=0 -> return zeros(1,M,N))
%   feasible_best: bool (all users meet SINR threshold)
%   SSNR_best: sensing SNR achieved
%   J_best: best fitness (minimized)

    [U, M, N] = size(H_comm);
    T = sensing_streams;
    if T <= 0
        T = 0;
    end

    % ---------- PSO hyperparameters (tune if needed) ----------
    nPop  = 40;      % swarm size
    MaxIt = 80;      % iterations

    w_max = 0.9;     % inertia start
    w_min = 0.4;     % inertia end
    c1 = 1.6;        % cognitive
    c2 = 1.6;        % social

    % Decision variables: complex beams for comm and sensing
    D = M*N;
    n_comm = U*D;
    n_sens = T*D;
    dim = 2*(n_comm + n_sens);  % real+imag

    % Search bounds (amplitude). Power is enforced by projection anyway.
    xMin = -1; xMax = 1;

    % Initialize
    x = xMin + (xMax-xMin)*rand(nPop, dim);
    v = zeros(nPop, dim);

    pBest = x;
    pBestVal = inf(nPop,1);

    gBest = x(1,:);
    gBestVal = inf;

    % Evaluate initial
    for i=1:nPop
        val = fitness(x(i,:));
        pBestVal(i) = val;
        if val < gBestVal
            gBestVal = val;
            gBest = x(i,:);
        end
    end

    % Main loop
    for it=1:MaxIt
        w = w_max - (w_max-w_min)*(it-1)/(MaxIt-1);

        for i=1:nPop
            r1 = rand(1,dim);
            r2 = rand(1,dim);
            v(i,:) = w*v(i,:) + c1*r1.*(pBest(i,:)-x(i,:)) + c2*r2.*(gBest-x(i,:));
            x(i,:) = x(i,:) + v(i,:);

            % clamp
            x(i,:) = max(min(x(i,:), xMax), xMin);

            val = fitness(x(i,:));

            if val < pBestVal(i)
                pBest(i,:) = x(i,:);
                pBestVal(i) = val;

                if val < gBestVal
                    gBestVal = val;
                    gBest = x(i,:);
                end
            end
        end
    end

    % Decode best
    [F_comm_best, F_sensing_best] = decode_and_project(gBest);

    % Compute final metrics
    SINR = compute_SINR(H_comm, F_comm_best, F_sensing_best, sigmasq_ue);
    SSNR_best = compute_sensing_SNR(sigmasq_radar_rcs, sensing_beamsteering, F_comm_best, F_sensing_best);

    feasible_best = all(SINR >= (gamma - 1e-6));
    J_best = gBestVal;

    % ----------------- nested helpers -----------------

    function [F_comm, F_sensing] = decode_and_project(xvec)
        % Split real/imag
        idx = 0;
        re_comm = xvec(idx+1:idx+n_comm); idx=idx+n_comm;
        im_comm = xvec(idx+1:idx+n_comm); idx=idx+n_comm;
        re_sens = xvec(idx+1:idx+n_sens); idx=idx+n_sens;
        im_sens = xvec(idx+1:idx+n_sens);

        F_comm_st = (re_comm + 1j*im_comm);
        F_comm = reshape(F_comm_st, [U, M, N]);

        if T > 0
            F_sens_st = (re_sens + 1j*im_sens);
            F_sensing = reshape(F_sens_st, [T, M, N]);
        else
            % keep shape consistent with compute_SINR loop (it iterates size(F_sensing,1))
            F_sensing = zeros(0, M, N);
        end

        % Hard projection to satisfy per-AP power <= P_all
        % Power is computed as sum over streams and antennas per AP
        Pm = compute_power(F_comm) + compute_power(F_sensing); % 1 x M

        for m=1:M
            if Pm(m) > P_all + 1e-12
                scale = sqrt(P_all / Pm(m));
                F_comm(:,m,:) = F_comm(:,m,:) * scale;
                if T > 0
                    F_sensing(:,m,:) = F_sensing(:,m,:) * scale;
                end
            end
        end
    end

    function J = fitness(xvec)
        [F_comm, F_sensing] = decode_and_project(xvec);

        SINR = compute_SINR(H_comm, F_comm, F_sensing, sigmasq_ue);
        SSNR = compute_sensing_SNR(sigmasq_radar_rcs, sensing_beamsteering, F_comm, F_sensing);

        % Objective: maximize SSNR -> minimize -SSNR
        obj = -SSNR;

        % SINR penalty (soft). Power constraint already projected hard, but keep tiny penalty for stability
        pen_sinr = sum(max(0, gamma - SINR).^2);

        % Optional tiny power penalty (after projection, should be ~0)
        Pm = compute_power(F_comm) + compute_power(F_sensing);
        pen_pow = sum(max(0, Pm - P_all).^2);

        % Weights (tune). SINR feasibility is usually the hardest, so weight it high.
        lambda_sinr = 1e3;
        lambda_pow  = 1e1;

        J = obj + lambda_sinr*pen_sinr + lambda_pow*pen_pow;
    end

end