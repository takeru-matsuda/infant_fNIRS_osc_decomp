function []=plot_Vdecomp(Y,fs,decomp_mu,decomp_cov,osc_param,num_component)
	T = size(Y,2);
	T1 = 1;
	T2 = T;
	Tstart = 0;
    osc_a = osc_param(num_component,1:num_component);
    osc_sigma = osc_param(num_component,2*num_component+1:3*num_component);
	subplot(num_component+1,1,1);
	plot(Tstart+(T1:T2)/fs,Y(1,:));
	xlim(Tstart+[T1 T2]/fs);
	set(gca,'FontSize',12);
	for i=1:num_component
    	c = sqrt(osc_sigma(i))./sqrt(1-osc_a(i).^2);
    	subplot(num_component+1,1,i+1);
	    hold on
	    plot(Tstart+(T1:T2)/fs,c*decomp_mu(2*i-1,T1:T2,num_component));
    	std_bar = reshape(sqrt(decomp_cov(2*i-1,2*i-1,T1:T2,num_component)),1,T2-T1+1);
	    xx = [Tstart+(T1:T2)/fs Tstart+(T2:-1:T1)/fs]';
    	yy = [c*(decomp_mu(2*i-1,T1:T2,num_component)-chi2inv(0.9,1)*std_bar) c*(decomp_mu(2*i-1,T2:-1:T1,num_component)+chi2inv(0.9,1)*std_bar(T2-T1+1:-1:1))]';
    	fill(xx,yy,[.8,.8,.8],'EdgeColor','none');
    	plot(Tstart+(T1:T2)/fs,c*decomp_mu(2*i-1,T1:T2,num_component));
	    xlim(Tstart+[T1 T2]/fs);
    	set(gca,'FontSize',12);
	end
end
