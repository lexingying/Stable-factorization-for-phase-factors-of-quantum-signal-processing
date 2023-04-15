if(0)
    d = 5; %5 times
    N = max(40*d,1024);
    %L = round(N/2.2);
    L = 2*d+2;
    D = 2*d+1;
    
    ts = [0:(N-1)]'/N * 2*pi;
    zs = exp(i*ts);
    xs = cos(ts);
    
    as = 0.3*cos(1*ts);    %plot(xs,as);
    bs = 0.5*sin(d*ts);
end

Ctime = tic;
if(1)
    vs = (1-as.^2-bs.^2);
    us = 1./(1-as.^2-bs.^2);
    uk = fft(us);
    
    hk = real(uk(N+1-[1:(L+D)])); %Fourier modes -1,-2,...

    switch nullAlg
      case 'svd' %svd
        M = zeros(L,D);    for a=1:D;        M(:,a) = hk(a:(a-1+L));    end
        [tU,tS,tV] = svd(M); p = tV(:,D); p = p/p(end);
      case 'itr' %iterative, but form matrix
        M = zeros(L,D);    for a=1:D;        M(:,a) = hk(a:(a-1+L));    end
        eta = 1e+0;
        Bfunc = @(x) M'*(M*x) + eta*x;
        p0 = zeros(D,1);
        p = randn(D,1); p=p/norm(p);
        cnt= 0;
        while(norm(p-p0)>1e-16)
            p0 = p;
            [p,flag,relres,iter,resvec] = pcg(Bfunc,p,1e-6,1000,[],[],p); 
            fprintf(1,'%d %d %d | ', flag, cnt, iter);
            p=p/norm(p);
            cnt = cnt+1;
        end
        p = p/p(end);
      case 'fft' %iterative, but use FFT for matrix multvec
        amphat = fft(hk([end,1:end-1])); %length D+L
        Dfunc = @(x) x(1:D);
        Lfunc = @(x) x(1:L);
        Mfunc = @(x) real(Lfunc(ifft(fft([zeros(L,1);x(end:-1:1)]).*amphat)));
        Mpfunc= @(x) real(Dfunc(ifft(fft([zeros(D,1);x(end:-1:1)]).*amphat)));
        %z = randn(D,1);        norm(M*z - Mfunc(z))        %z = randn(L,1);        norm(M'*z - Mpfunc(z))
        
        eta = 1e+0;
        Bfunc = @(x) Mpfunc(Mfunc(x)) + eta*x;
        p0 = zeros(D,1);
        p = randn(D,1); p=p/norm(p);
        cnt= 0;
        while(norm(p-p0)>1e-16)
            p0 = p;
            [p,flag,relres,iter,resvec] = pcg(Bfunc,p,1e-6,1000,[],[],p); 
            fprintf(1,'%d %d %d | ', flag, cnt, iter);
            p=p/norm(p);
            cnt = cnt+1;
        end
        p = p/p(end);
    end
    
    pk = zeros(L,1); pk(1:D) = p;
    pk = circshift(pk, -d);
    ez = ifft(pk) * L;
    e1z = fft(pk);
    
    as = interpft(as,L);
    bs = interpft(bs,L);
    vs = interpft(vs,L);
    
    alpha = real( (ez.*e1z) \ vs );    fprintf(1, 'alpha err %1.1d\n', norm( (ez.*e1z)*alpha-vs)/norm(vs));
    cs = real( (ez+e1z)/2     * sqrt(alpha) );
    ds = real( (ez-e1z)/(2*i) * sqrt(alpha) );
    
    %compute all angles
    ps = as + i*cs;
    qs = bs + i*ds;
    
    %ak = fft(as)/N;
    %bk = fft(bs)/N; 
    %ck = fft(cs)/N;
    %dk = fft(ds)/N;
    %pk = ak + i*ck;
    %qk = bk + i*dk;

    ts = [0:(L-1)]'/L * 2*pi;
    coss = cos(ts);
    sins = sin(ts);
    
    angles = zeros(d+1,1);
    for g=0:(d-1)
        
        pk = fft(ps);
        qk = fft(qs);
        
        %get angle
        idx = d+1-g;        %tmp = (ak(idx)+i*ck(idx))/(bk(idx)+i*dk(idx));
        tmp = pk(idx)/qk(idx);
        ang = atan2(imag(tmp),real(tmp)) / 2; 
        angles(d+1-g) = ang; %store angle
        pk = pk * exp(-i*ang);
        qk = qk * exp(+i*ang);
        
        ps = ifft(pk);
        qs = ifft(qk);
        
        %lower degree
        po = ps;
        qo = qs;
        ps = po.*coss + qo.*(-i*sins);
        qs = po.*(-i*sins) + qo.*coss;
        
    end
    %last one
    if(1)
        pk = fft(ps);
        qk = fft(qs);
        
        tmp = pk(1);
        ang = atan2(imag(tmp),real(tmp));
        angles(1) = ang;
    end
    %angles;
    Ctime = toc(Ctime);
end

if(1)
    I = [i 0; 0 -i];
    J = [0 1; -1 0];
    K = [0 i; i 0];
    R = zeros(2,2,L);
    
    angMat = zeros(2,2,d+1);
    for h = 1:(d+1)
        angMat(:,:,h) = expm(angles(h)*I);
    end
    
    [tmpval,tmpidx] = max(abs(as));
    PK = sort([ceil(rand(1,20)*L) tmpidx]);
    
    for g=PK
        t = ts(g);
        eKt = expm(K*t);
        res = angMat(:,:,1);
        for h = 2:(d+1)            %tmp = angMat(:,:,h);
            res = res * eKt * angMat(:,:,h);
        end
        R(:,:,g) = res;
    end

    ps = as + i*cs;
    qs = bs + i*ds;
    
    pa = squeeze(R(1,1,PK));
    qa = squeeze(R(1,2,PK));
    perr = norm(pa-ps(PK))/norm(ps(PK));
    qerr = norm(qa-qs(PK))/norm(qs(PK));
    
    Caerr = norm(real(pa-ps(PK)),Inf)/norm(as,Inf);
    fprintf(1, 'perr %1.1d, qerr %1.1d aerr %1.1d\n', perr, qerr, Caerr);
    
end


