function [A,E,r,mll] = VAR_myule(Y,ARdeg)
    d = size(Y,1);
    T = size(Y,2);
    acov = zeros(d,d,ARdeg+1);
    Y = Y-mean(Y,2)*ones(1,T);
    for k=1:ARdeg+1
        acov(:,:,k) = Y(:,1:T-k+1)*Y(:,k:T)'/T;
    end
    C = zeros(d*ARdeg,d*ARdeg);
    for k=1:ARdeg
        for j=1:k
            C((k-1)*d+1:k*d,(j-1)*d+1:j*d) = acov(:,:,k-j+1);
        end
        for j=k+1:ARdeg
            C((k-1)*d+1:k*d,(j-1)*d+1:j*d) = acov(:,:,j-k+1)';
        end
    end
    c = zeros(d*ARdeg,d);
    for k=1:ARdeg
        c((k-1)*d+1:k*d,:) = acov(:,:,k+1);
    end
    D = zeros(d*(ARdeg+1),d*(ARdeg+1));
    D(1:d*ARdeg,1:d*ARdeg) = C;
    for k=1:ARdeg
        D((k-1)*d+1:k*d,ARdeg*d+1:(ARdeg+1)*d) = c((ARdeg-k)*d+1:(ARdeg-k+1)*d,:);
        D(ARdeg*d+1:(ARdeg+1)*d,(k-1)*d+1:k*d) = c((ARdeg-k)*d+1:(ARdeg-k+1)*d,:)';
    end
    D(ARdeg*d+1:(ARdeg+1)*d,ARdeg*d+1:(ARdeg+1)*d) = acov(:,:,1);
    eigs = eig(D);
    if min(eigs) <= 0
        r = 0;
        mll = VAR_myule_ll(Y,r,C,c);
    else
        [r,mll] = fminbnd(@(r)VAR_myule_ll(Y,r,C,c),0,min(eigs));
    end
    R = r*eye(d);
    B = (C-r*eye(d*ARdeg))\c;
    A = zeros(d,d*ARdeg);
    for k=1:ARdeg
        A(:,(k-1)*d+1:k*d) = B((k-1)*d+1:k*d,:)';
    end
    E = acov(:,:,1)-R;
    for k=1:ARdeg
        E = E-A(:,(k-1)*d+1:k*d)*acov(:,:,k+1);
    end
    A = [eye(d) -A];
end

