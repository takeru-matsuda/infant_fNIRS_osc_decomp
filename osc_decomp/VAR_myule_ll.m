function mll = VAR_myule_ll(Y,r,C,c)
    d = size(Y,1);
    ARdeg = size(C,1)/d;
    B = (C-r*eye(d*ARdeg))\c;
    A = zeros(d,d*ARdeg);
    for k=1:ARdeg
        A(:,(k-1)*d+1:k*d) = B((k-1)*d+1:k*d,:)';
    end
    E = C(1:d,1:d)-r*eye(d);
    for k=1:ARdeg
        E = E-A(:,(k-1)*d+1:k*d)*c((k-1)*d+1:k*d,:);
    end
    vecE = zeros(1,d*(d+1)/2);
    k = 1;
    for i=1:d
        for j=1:i
            vecE(k) = E(i,j);
            k = k+1;
        end
    end
    mll = VARwithnoise_ll(Y,A,E,r*eye(d));
end

