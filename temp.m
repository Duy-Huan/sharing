%% envtemp MODULE - ambient temperature model with slow random fluctuation
% T(t) = setpoint + soft_limit(x(t)) + HVAC cycle + slow drift
%   x(t) : Ornstein-Uhlenbeck (first-order Gauss-Markov) process, stationary std = sigma,
%          correlation time = tau_corr  (slow, random, mean-reverting around the setpoint)
%   soft_limit(x) = max_dev*tanh(x/max_dev)  keeps |deviation| below max_dev
classdef envtemp < handle
    properties
        % Public parameters
        Setpoint        % Target ambient temperature (oC)
        Sigma           % Std of the random component (oC)
        TauCorr         % Correlation time of the random component (s)
        MaxDev          % Soft limit on the random deviation (oC); 0 = no limit
        HvacAmp         % Amplitude of the periodic (air-conditioning) cycle (oC); 0 = off
        HvacPeriod      % Period of the periodic cycle (s)
        DriftRate       % Slow linear drift (oC/hour); 0 = off

        % Private parameters
        dt              % Time step (s)
        Phi             % AR(1) coefficient exp(-dt/TauCorr)
        x               % Current state of the random component (oC)
        t               % Elapsed time (s)
        HvacPhase       % Random initial phase of the periodic cycle (rad)
    end

    methods
        %% Constructor
        function obj = envtemp(config, dt_val)
            % Defaults (used for any missing field)
            d.setpoint = 25.0; d.sigma = 0.5; d.tau_corr = 600; d.max_dev = 2.0;
            d.hvac_amp = 0.3;  d.hvac_period = 900; d.drift_rate = 0.0;

            if nargin == 0 || isempty(config)
                default_file = 'config.json';
                config = struct();
                if exist(default_file, 'file')
                    raw_config = jsondecode(fileread(default_file));
                    if isfield(raw_config, 'envtemp'), config = raw_config.envtemp; end
                end
            end
            names = fieldnames(d);
            for i = 1:numel(names)
                if ~isfield(config, names{i}), config.(names{i}) = d.(names{i}); end
            end

            obj.Setpoint   = config.setpoint;
            obj.Sigma      = config.sigma;
            obj.TauCorr    = config.tau_corr;
            obj.MaxDev     = config.max_dev;
            obj.HvacAmp    = config.hvac_amp;
            obj.HvacPeriod = config.hvac_period;
            obj.DriftRate  = config.drift_rate;

            if nargin < 2 || isempty(dt_val), obj.dt = 1e-3; else, obj.dt = dt_val; end
            obj.Phi = exp(-obj.dt / obj.TauCorr);
            obj.reset();
        end

        %% Reset to a stationary random state
        function reset(obj)
            obj.x = obj.Sigma * randn();            % start inside the stationary distribution
            obj.t = 0;
            obj.HvacPhase = 2*pi*rand();
        end

        %% One sample (use inside the simulation loop)
        function T = step(obj)
            obj.x = obj.Phi * obj.x + obj.Sigma * sqrt(1 - obj.Phi^2) * randn();
            T = obj.compose(obj.x, obj.t);
            obj.t = obj.t + obj.dt;
        end

        %% N samples at once (fast; use before the loop)
        function T = generate(obj, N)
            w  = randn(N, 1);
            b0 = obj.Sigma * sqrt(1 - obj.Phi^2);
            xv = filter(b0, [1 -obj.Phi], w, obj.Phi * obj.x);   % continues from current state
            tv = obj.t + (0:N-1)' * obj.dt;
            T  = obj.compose(xv, tv);
            obj.x = xv(end);
            obj.t = obj.t + N * obj.dt;
        end
    end

    methods (Access = private)
        function T = compose(obj, xv, tv)
            dev = xv;
            if obj.MaxDev > 0
                dev = obj.MaxDev * tanh(xv / obj.MaxDev);
            end
            hv = 0;
            if obj.HvacAmp ~= 0
                hv = obj.HvacAmp * sin(2*pi*tv / obj.HvacPeriod + obj.HvacPhase);
            end
            T = obj.Setpoint + dev + hv + (obj.DriftRate / 3600) * tv;
        end
    end
end
