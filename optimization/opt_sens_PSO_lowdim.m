function [F_sens_best, feasible_best, SSNR_best, J_best] = opt_sens_PSO_lowdim( ...
    H_comm, F_comm_fixed, sigmasq_ue, sensing_beamsteering, ...
    F_sens_base, sigmasq_radar_rcs, gamma, P_sensing)
% Low-dimension PSO (MATLAB Online friendly):
%   - Comm beam is fixed (e.g., RZF)
%   - PSO optimizes only per-AP complex weights applied to a base sensing beam
%     F_sens(m,:) = c_m * F_sens_base(m,:)
%
% Inputs:
%   H_comm: U x M x N
%   F_comm_fixed: U x M x N (fixed comm beam, e.g., RZF)
%   F_sens_base: 1 x M x N (base sensing beam, e.g., NS or beamsteering)
%   gamma: SINR threshold
%   P_sensing: sensing power per AP (scalar)
%
% Outputs:
%   F_sens_best: 1 x M x N
%   feasible_best: bool
%   SSNR_best: best sensing SNR
%   J_best: best fitness

    [~, M, ~] = size(H_comm);

    % Decision variables: c_m (complex) for m=1..M
    % Encode as [Re(c1..cM), Im(c1..cM)] => dim = 2M
    dim = 2*M;

    % PSO params (small because problem is tiny)
    nPop  = 25;
    MaxIt = 60;
    w_max = 0.9; w_min = 0.4;
    c1 = 1.6; c2 = 1.6;

    % bounds for weights (amplitude)
    % final power is repaired per AP anyway
    xMin = -2; xMax = 2;

    % init
    x = xMin + (xMax-xMin)*rand(nPop, dim);
    v = zeros(nPop, dim);

    pBest = x;
    pBestVal = inf(nPop,1);

    gBest = x(1,:);
    gBestVal = inf;

    for i=1:nPop
        val = fitness(x(i,:));
        pBestVal(i) = val;
        if val < gBestVal
            gBestVal = val; gBest = x(i,:);
        end
    end

    for it=1:MaxIt
        w = w_max - (w_max-w_min)*(it-1)/(MaxIt-1);
        for i=1:nPop
            r1 = rand(1,dim); r2 = rand(1,dim);
            v(i,:) = w*v(i,:) + c1*r1.*(pBest(i,:)-x(i,:)) + c2*r2.*(gBest-x(i,:));
            x(i,:) = x(i,:) + v(i,:);
            x(i,:) = max(min(x(i,:), xMax), xMin);

            val = fitness(x(i,:));
            if val < pBestVal(i)
                pBest(i,:) = x(i,:);
                pBestVal(i) = val;
                if val < gBestVal
                    gBestVal = val; gBest = x(i,:);
                end
            end
        end
    end

    % decode best
    F_sens_best = build_and_repair(gBest);

    SINR = compute_SINR(H_comm, F_comm_fixed, F_sens_best, sigmasq_ue);
    feasible_best = all(SINR >= (gamma - 1e-6));

    SSNR_best = compute_sensing_SNR(sigmasq_radar_rcs, sensing_beamsteering, F_comm_fixed, F_sens_best);

    J_best = gBestVal;

    % ---------------- helpers ----------------
    function F_sens = build_and_repair(xvec)
        re = xvec(1:M);
        im = xvec(M+1:end);
        c = re + 1j*im; % 1xM

        % Apply weights per AP (broadcast along antennas)
        F_sens = F_sens_base;
        for m=1:M
            F_sens(:,m,:) = F_sens(:,m,:) * c(m);
        end

        % Repair sensing power per AP
        Ps = compute_power(F_sens); % 1xM
        for m=1:M
            if Ps(m) > P_sensing + 1e-12
                s = sqrt(P_sensing / Ps(m));
                F_sens(:,m,:) = F_sens(:,m,:) * s;
            end
        end
    end

    function J = fitness(xvec)
        F_sens = build_and_repair(xvec);

        SINR = compute_SINR(H_comm, F_comm_fixed, F_sens, sigmasq_ue);
        SSNR = compute_sensing_SNR(sigmasq_radar_rcs, sensing_beamsteering, F_comm_fixed, F_sens);

        % Feasibility-first: nặng cho violation, sau đó maximize SSNR
        viol = sum(max(0, gamma - SINR).^2);
        J = 1e4*viol - SSNR;
    end

end
