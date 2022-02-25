function [osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp(Y,fs,MAX_OSC,MAX_AR,algorithm,grad)
    arguments
        Y
        fs (1,1)
        MAX_OSC (1,1) {mustBeInteger(MAX_OSC),mustBePositive(MAX_OSC)} = 1
        MAX_AR (1,1) {mustBeInteger(MAX_AR),mustBePositive(MAX_AR)} = ceil(2*MAX_OSC/size(Y,1))
        algorithm = 'quasi-newton'
        grad = true
    end
    if size(Y,1) == 1
        [osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp_uni(Y,fs,MAX_OSC,MAX_AR,algorithm,grad);
    else
        [osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp_multi(Y,fs,MAX_OSC,MAX_AR,algorithm,grad);    
    end
end


