function [osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp_multi(Y,fs,MAX_OSC,MAX_VAR,algorithm,grad)
%
% Input:
%    Y:             multivariate time series (J times T)
%    fs:            sampling frequency (scalar)
%    MAX_OSC:       maximum number of oscillation components (scalar)
%    MAX_VAR:       maximum VAR degree (scalar), default = ceil(2*MAX_OSC/J)
%    algorithm:     optimization algorithm (string), default = 'quasi-Newton'
%    grad:          specify gradient or not (boolean), default = false
%
% Output:
%    osc_param:     estimated parameters of the oscillator model (MAX_OSC times (2*J+1)*MAX_OSC+1)
%                   osc_param(K,1:K) is the estimated a_1,...,a_K of the oscillator model with K oscillation components
%                   osc_param(K,K+1:2*K) is the estimated f_1,...f_K of the oscillator model with K oscillation components
%                   osc_param(K,2*K+1:3*K) is the estimated sigma_1^2,...,sigma_K^2 of the oscillator model with K oscillation components
%                   osc_param(K,3*K+1:(2*J+1)*K) is the estimated c_{21,1},c_{21,2},...,c_{J1,1},c_{J1,2},...,c_{JK,1},c_{JK,2} of the oscillator model with K oscillation components
%                   osc_param(K,(2*J+1)*K+1) is the estimated tau^2 of the oscillator model with K oscillation components
%    osc_AIC:       AIC of the oscillator model (1 times MAX_OSC)
%                   AIC_osc(K) is the AIC of the oscillator model with K oscillation components
%    osc_mean:      smoothed coordinate of each oscillator (2*MAX_OSC times T times MAX_OSC)
%                   decomp_mu(2*k-1:2*k,:,K) is the smoothed coordinate of the k-th oscillator in the decomposition into K components
%    osc_cov:       smoothed covariance of each oscillator (2*MAX_OSC times 2*MAX_OSC times T times MAX_OSC)
%                   decomp_cov(2*k-1:2*k,2*k-1:2*k,:,K) is the smoothed covariance of the k-th oscillator in the decomposition into K components
%    osc_phase:     estimated phase of each oscillation component (MAX_OSC times T times MAX_OSC)
%                   phi_prop(k,:,K) is the estimated phase of the k-th oscillator in the decomposition into K components
%
	J = size(Y,1);
	T = size(Y,2);
    if nargin == 3
        MAX_VAR = ceil(2*MAX_OSC/J);
    else
        MAX_VAR = max(MAX_VAR,ceil(2*MAX_OSC/J));
    end
    [VARwithnoise_A,VARwithnoise_E,VARwithnoise_r,VARwithnoise_AIC] = VAR_fit(Y,MAX_VAR);
    osc_param = zeros(MAX_OSC,(2*J+1)*MAX_OSC+1);
    osc_AIC = zeros(1,MAX_OSC);
    osc_mean = zeros(2*MAX_OSC,T,MAX_OSC);
    osc_cov = zeros(2*MAX_OSC,2*MAX_OSC,T,MAX_OSC);
    osc_phase = zeros(MAX_OSC,T,MAX_OSC);
	for K=1:MAX_OSC
    	K
	    minAIC = inf;
    	minAIC2 = inf;
	    minK = inf;
    	for ARdeg=ceil(K/J):MAX_VAR
	        [Vtmp,tmp] = polyeig_VAR(VARwithnoise_A(:,1:J*ARdeg,ARdeg));
    	    if J*ARdeg-nnz(imag(tmp))/2 == K && VARwithnoise_AIC(ARdeg) < minAIC
        	    V0 = Vtmp;
            	z0 = tmp;
            	E0 = VARwithnoise_E(:,:,ARdeg);
            	R0 = VARwithnoise_r(ARdeg);
            	minAIC = VARwithnoise_AIC(ARdeg);
    		    optARdeg = ARdeg;
	        end
	        if J*ARdeg-nnz(imag(tmp))/2 > K && (J*ARdeg-nnz(imag(tmp))/2 < minK || (J*ARdeg-nnz(imag(tmp))/2 == minK && VARwithnoise_AIC(ARdeg) < minAIC2))
    	        V1 = Vtmp;
        	    z1 = tmp;
            	E1 = VARwithnoise_E(:,:,ARdeg);
            	R1 = VARwithnoise_r(ARdeg);
            	minAIC2 = VARwithnoise_AIC(ARdeg);
        	    optARdeg2 = ARdeg;
    	        minK = J*ARdeg-nnz(imag(tmp))/2;
	        end
	    end
    	if minAIC == inf
        	if minAIC2 == inf
            	warning('no VAR model with %d oscillators',K);
	        end
    	    V0 = V1;
        	z0 = z1;
	        E0 = E1;
    	    R0 = R1;
        	optARdeg = optARdeg2;
	    end
    	[~,ARdeg] = min(VARwithnoise_AIC);
	    [Vtmp,tmp] = polyeig_VAR(VARwithnoise_A(:,1:J*ARdeg,ARdeg));
    	if nnz(imag(tmp)>=0) >= K
	        V0 = Vtmp;
    	    z0 = tmp;
        	E0 = VARwithnoise_E(:,:,ARdeg);
        	R0 = VARwithnoise_r(ARdeg);
        	minAIC = VARwithnoise_AIC(ARdeg);
    	    optARdeg = ARdeg;
	    end
	    VV = zeros(J*optARdeg,J*optARdeg);
    	for j=1:J*optARdeg
        	for i=1:optARdeg
            	VV((i-1)*J+1:i*J,j) = z0(j)^(1-i)*V0(:,j);
    	    end
	    end
    	QQ = inv(VV)*[E0 zeros(J,J*(optARdeg-1)); zeros(J*(optARdeg-1),J*optARdeg)]*inv(VV)';
	    [~,I] = sort(diag(real(QQ))./(1-abs(z0).^2),'descend');
    	V0 = V0(:,I);
    	z0 = z0(I);
    	if R0 == 0
    	    R0 = 10^-5;
        end
        
	    init_a = zeros(1,K);
    	init_theta = zeros(1,K);
    	init_c = zeros(1,2*(J-1)*K);
    	kk = 1;
    	for k=1:K
        	init_a(k) = abs(z0(kk));
    	    init_theta(k) = abs(angle(z0(kk)));
    	    for ii=1:J-1
        	    init_c(2*(J-1)*(k-1)+2*ii-1) = real(V0(ii+1,kk)/V0(1,kk));
            	if imag(z0(kk)) < 0
                	init_c(2*(J-1)*(k-1)+2*ii) = imag(V0(ii+1,kk)/V0(1,kk));
	            else
    	            init_c(2*(J-1)*(k-1)+2*ii) = -imag(V0(ii+1,kk)/V0(1,kk));
        	    end
	        end
    	    if imag(z0(kk)) == 0
        	    kk = kk+1;
	        else
    	        kk = kk+2;
            end
        end
        
        if mod(T,2) == 0
        	freq = [2*pi/T*(0:T/2-1) 4];
	    else
    	    freq = [2*pi/T*(0:(T-1)/2)];
        end        
        P = zeros(J*length(freq),K);
	    for k=1:K
    	    a = init_a(k);
        	theta = init_theta(k);
        	A = (1-2*a^2*cos(theta)^2+a^4*cos(2*theta))/a/(a^2-1)/cos(theta);
        	b = (A-2*a*cos(theta)+sign(cos(theta))*sqrt((A-2*a*cos(theta))^2-4))/2;
        	for j=1:length(freq)
            	P(j,k) = -a*cos(theta)/b*abs(1+b*exp(-1i*freq(j)))^2/abs(1-2*a*cos(theta)*exp(-1i*freq(j))+a^2*exp(-2*1i*freq(j))).^2;
	            for ii=1:J-1
    	            P(ii*length(freq)+j,k) = norm(init_c(2*(J-1)*(k-1)+2*ii-1:2*(J-1)*(k-1)+2*ii))^2*P(j,k);
        	    end
    	    end
	    end
	    [~,bestARdeg] = min(VARwithnoise_AIC);
    	A0 = VARwithnoise_A(:,1:J*bestARdeg,bestARdeg);
	    E0 = VARwithnoise_E(:,:,bestARdeg);
    	p = zeros(J*length(freq),1);
	    for j=1:length(freq)
    	    Phi = eye(J);
        	for k=1:bestARdeg
            	Phi = Phi-A0(:,(k-1)*J+1:k*J)*exp(-1i*k*freq(j));
	        end
    	    tmp = inv(Phi)*E0*inv(Phi)';
        	for ii=1:J
            	p((ii-1)*length(freq)+j) = real(tmp(ii,ii));
    	    end
	    end
	    p = zeros(J*length(freq),1);
    	d1 = zeros(J*length(freq),1);
    	for j=1:length(freq)
        	if freq(j) == 0 || freq(j) == 4
            	d1(j:length(freq):J*length(freq)) = 1;
	        end
    	    for ii=1:J
        	    p((ii-1)*length(freq)+j) = abs(Y(ii,:)*exp(-1i*freq(j)*(0:T-1)'))^2/T;
    	    end
        end
        if cond(P'*P) < 10^6
            init_sigma2 = (P'*P)\(P'*p);
        else
            init_sigma2 = expGLMfit(P,p);
        end
        init_sigma2(init_sigma2<0) = R0;
        init_sigma2 = init_sigma2';
        init_tau2 = R0;

        if nargin >= 5
            switch lower(algorithm)
                case 'quasi-newton'
                    if nargin == 6 && grad == true
                        options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000,'SpecifyObjectiveGradient',true);
                        param = fminunc(@(param)osc_multi_prof_ll(Y,param,init_theta,true),[atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2) init_c],options);
                    else
                        options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000,'SpecifyObjectiveGradient',false);
                        param = fminunc(@(param)osc_multi_prof_ll(Y,param,init_theta,false),[atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2) init_c],options);
                    end
                case 'trust-region'
                    options = optimoptions('fminunc','Algorithm','trust-region','MaxFunEvals',10000,'SpecifyObjectiveGradient',true);
                    param = fminunc(@(param)osc_multi_prof_ll(Y,param,init_theta,true),[atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2) init_c],options);
                case 'cg'
                    param = minimize([atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2) init_c]',@(param)osc_multi_prof_ll(Y,param,init_theta,true),-10000)';
                case 'em'
                    param = osc_em([atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2) init_c]);
                otherwise
                    error('osc_decomp_multi: unknown algorithm')
            end
        else
            options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000,'SpecifyObjectiveGradient',false);
            param = fminunc(@(param)osc_multi_prof_ll(Y,param,init_theta,false),[atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2) init_c],options);
        end
    	[mll,~,osc_tau2] = osc_multi_prof_ll(Y,param,init_theta,false);
    	osc_AIC(K) = 2*mll+2*((2*J+1)*K+1);

        param(K+1:2*K) = init_theta+tanh(param(K+1:2*K))*pi;
	    [~,I] = sort(abs(param(K+1:2*K)));
    	osc_a = (tanh(param(I))+1)/2;
    	osc_f = abs(param(K+I))*fs/2/pi;
    	osc_sigma2 = exp(param(2*K+I))*osc_tau2;
        tmp = zeros(1,2*(J-1)*K);
        for k=1:K
            tmp((k-1)*2*(J-1)+1:k*2*(J-1)) = I(k)-1;
        end
        osc_c = param(3*K+tmp*2*(J-1)+repmat(1:2*(J-1),1,K));
		osc_param(K,1:(2*J+1)*K+1) = [osc_a osc_f osc_sigma2 osc_c osc_tau2];
        [osc_mean(1:2*K,:,K),osc_cov(1:2*K,1:2*K,:,K)] = osc_smooth(Y,fs,osc_param(K,1:(2*J+1)*K+1));
	    for k=1:K
        	osc_phase(k,:,K) = atan2(osc_mean(2*k,:,K),osc_mean(2*k-1,:,K));
    	end
	end
end


