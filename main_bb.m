%weight is not correct

clear all;
N=10;
K=4;
%M=5;
Rk = 1;
rhoP = 1; %primary users' tranmit power
sigma2_dbm= -90;%+10*log10(BW)+Nf; %Thermal noise in dBm
sigma = 10^((sigma2_dbm-30)/10);
Pmax = rhoP;  
xi=100000000;%1.0000e+15;
taux = 0.0000001;tauy = 50; epsilon_cbv=0.01;
al=2; %path loss
c0=3*10^8; fc = 300*10^9; rl = 5*exp(-3);%molecular absorption
ct=500;
eps_BB = 0.1;
max_itr=200;

Mvec = [   1 2 4 6 8    ];
%Kvec = [1 3 5 7 10];
for iM = 1: length(Mvec)
    M = Mvec(iM);
    weight = ones(1,M);
    sum4=0;
    for ict = 1 : ct
        %primary users' channels
        akP = complex(sqrt(0.5)*randn(K,1),sqrt(0.5)*randn(K,1));
        rP = 10; %size of the square
        %rkP = 10;
        coodP = [rP*rand(K,1).*sign(randn(K,1)) rP*rand(K,1).*sign(randn(K,1))]; %locations
        rkP = sqrt(abs(coodP(:,1)).^2+abs(coodP(:,2)).^2);
        PLkP = (c0/4/pi/fc)^2*exp(-rl*rkP)./(1+rkP.^al);
        scale = max(1./PLkP);%need to scale the term inside of the log for optimization

        thetak = -pi/2+pi/K*[1:1:K]';%pi*rand(K,1)-pi/2;%within -pi/2 and pi/2 %sign(randn)*2*pi/3*rand(K,1);
        NN = [0:1:N-1]';
        H = exp(-complex(0,1)*pi*sin(thetak)'.*NN)*(diag(akP).*sqrt(PLkP));%channel matrix, N by K

        %build the codebook
        NQ = 10; %number of codewords/beams in the codebook
        %build the analog beamforming
        theta_vector = pi*[0:1:NQ-1]/NQ-pi/2;
        ABF=[]; beam=[];
        for i = 1 : K
            [temp1,indexk] = min(abs(thetak(i)-theta_vector));
            ABF(:,i) = exp(-complex(0,1)*pi*sin(theta_vector(indexk))*NN)/sqrt(N);
            theta_vector(indexk) = [];%to avoid the use of the same vector
        end 
        %build digital beamforming
        Habf = H'*ABF; %channel matrix after analog beamforming - using H-ABF*sqrt(N) to check
        D_power = diag(1./[diag(inv(Habf'*Habf))]);
        DBF = inv(Habf)*sqrt(D_power);%/sqrt(K); %normalized digital beamforming 

        for k = 1 : K
            beam(:,k) = ABF*DBF(:,k);
        end
        ABF=beam; %beam is the composite precoding

        %secondary users' channels
        akS = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));
        %rkS = 5;
        rS = 10; %size of square
        coodS = [rS*rand(M,1).*sign(randn(M,1)) rS*rand(M,1).*sign(randn(M,1))];%locations
        rkS = sqrt(abs(coodS(:,1)).^2+abs(coodS(:,2)).^2);

        PLkS = (c0/4/pi/fc)^2*exp(-rl*rkS)./(1+rkS.^al);

        thetak_s = pi*rand(M,1)-pi/2;%within -pi/2 and pi/2 %sign(randn)*2*pi/3*rand(K,1);
        G = exp(-complex(0,1)*pi*sin(thetak_s)'.*NN)*(diag(akS).*sqrt(PLkS));%channel matrix, N by M 

        %G'*ABF*DBF
        %H'*ABF*DBF

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %generate h^P_ki
        hPki=[];hSjk=[];HPkk=[];ck=[];HSjk=[];bjk=[];tjk=[];
        for k =1 : K
            hPki(k,:) = abs(H(:,k)'*ABF).^2; %K rows, 
        end
        %generate h^S_jk
        for j =1 : M
            hSjk(j,:) = abs(G(:,j)'*ABF).^2; %M rows,  K column
        end
        %generate H^P_kk, ck, a vector
        for k =1 : K
            HPkk(k) = sum(hPki(k,:)) - hPki(k,k); %remove the kth one, 
            ck(k) = (sum(hPki(k,:)) - hPki(k,k))/hPki(k,k) - rhoP/(2^Rk-1)+sigma/hPki(k,k);
        end
        %generate H^S_jk, bjk, and tjk a matrix
        for j =1 : M
            for k = 1 : K
                HSjk(j,k) = sum(hSjk(j,:)) - hSjk(j,k); %M rows, K column, jk 
                bjk(j,k) = (sum(hSjk(j,:)) - hSjk(j,k))/hSjk(j,k)*rhoP - rhoP/(2^Rk-1)+sigma/hSjk(j,k);
                tjk(j,k) = (sum(hSjk(j,:)) - hSjk(j,k))*rhoP +sigma;
            end
        end

       %  load('data_x.mat')
       %  eps_BB = 0.01;
       % max_itr=300;
        %  Rev_mat=[];   xi=1.0000e+15;
        % taux = 0.00000001;tauy = 100; epsilon_cbv=0.01;
        %for data_wrong3.mat, reduce xi to 10, why?
        %%%%%%%%%%%%% algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        %initilization
        structure = zeros(M,K);
        for j = 1 : M
            for k = 1 : K
                if bjk(j,k) <=0 & ck(k) <=0;
                    structure(j,k) = 1;
                end
            end
        end
        
% % % %         % user scheduling
%         for k = 1 : K
%             index_uk = find(structure(:,k));
%             if length(index_uk)~=0 %non empty empty
%                 [temp1, uk] = max(hSjk(index_uk,k));%uk will be scheduled at beam k
%                 structure(:,k)=0;
%                 structure(index_uk(uk),k)=1; %only one user is scheduled
%             end
%         end
        
        [row,col] = find(structure); %row is j - [1,M], col is k - [1,K]
        S = [reshape(row,[],1) reshape(col,[],1)];
        number_sij = size(S,1);
        if number_sij == 0%no beam available
             continue;
         end


        %generate the matrices to recover each row
        R = []; Rev_mat=[];
        for k =1 : K %  K  submatrices for K users/subcarrier
            Rev_mat(1:M,:,k) = zeros(M,number_sij); %each submatrix is M by number_sij
            for n = 1 : number_sij
                if S(n,2)==k % this secondary user is active on k-th primary uer
                    Rev_mat(S(n,1),n,k) = 1; %this user is m-th user, S(n,2)
                end
            end    
            R = [R;Rev_mat(:,:,k)];
        end

        %generate cp, dp, ep, fp and dpx
        cp=[];dp=[];ep=[];fp=[];dpx=[];
        for p = 1 : number_sij
            cp = [cp [zeros(1,p-1) hSjk(S(p,1),S(p,2)) zeros(1,number_sij-p)]' ];

            dpp=[]; dppx=[];
            for k = 1 : K        
                if k == S(p,2)
                    onesp2 = ones(M,1); onesp2(S(p,1))=0;
                    dpp = [dpp ; xi*hSjk(S(p,1),k)*onesp2];
                    onesp2x = zeros(M,1);  
                    dppx = [dppx ; onesp2x];
                else
                    dpp = [dpp ; hSjk(S(p,1),k)*ones(M,1)];
                    dppx = [dppx ; hSjk(S(p,1),k)*ones(M,1)];
                end
            end            
            dp = [dp dpp]; dpx = [dpx dppx];

            epp=[];
            for k = 1 : K        
                if k == S(p,2)
                    zerosp2 = zeros(M,1); zerosp2(S(p,1))=1;
                    epp = [epp ; zerosp2];
                else
                    epp = [epp ; hPki(S(p,2),k)/hPki(S(p,2),S(p,2))*ones(M,1)];
                end
            end            
            ep = [ep epp];
            fpp=[];
            for k = 1 : K        
                if k == S(p,2)
                    zerosp2 = zeros(M,1); zerosp2(S(p,1))=1;
                    fpp = [fpp ; zerosp2];
                else
                    fpp = [fpp ; hSjk(S(p,1),k)/hSjk(S(p,1),S(p,2))*ones(M,1)];
                end
            end            
            fp = [fp fpp];
        end

        %inititlization for both algorithms
        tempy = hSjk(S(1,1),S(1,2));
        for p =1 : number_sij
            tempy(p) = hSjk(S(p,1),S(p,2));
        end
        [x,indy]=max(tempy);
        ym = zeros(2*number_sij,1);ym(indy)=Pmax; %ybas=ym(1:number_sij);
                
        %scheduling
        % log2(1+ hPki(4,4)*rhoP/(hPki(4,4)*ybb(1)+sigma))>Rk
        ybas=ym(1:number_sij);
        ybas(indy)=min([Pmax ; -ck(S(indy,2)); -bjk(S(indy,1),S(indy,2))]);
%             (rhoP/(2^Rk-1)-sigma/hPki(S(indy,2),S(indy,2)));
%             (rhoP/(2^Rk-1)-sigma/hSjk(S(indy,1),S(indy,2)))]);
        %j=S(indy,1)  k=S(indy,2)
        
        %%%%% BB %%%%%%
        %find the bound
        b_BB = [];
        for p =1 : number_sij
            b_BB = [b_BB; cp(:,p)'*ones(number_sij,1)*Pmax/tjk(S(p,1),S(p,2))];
        end
        %initilization of BB
        kBB = 1;  Bset = [zeros(number_sij,1) b_BB]; 
        xl = zeros(number_sij,1); xu = b_BB;
        if feasibility(xl,  Pmax,S, cp,dp,ep,fp, tjk, bjk, ck,R) %if xl is in the G
            Upk = -fobj(xl ,   weight, S   );
            Lowk = -fobj(xu ,   weight, S   );
        else
            Upk = 0;Lowk = 0;
        end

        Bset = [Bset; Lowk Upk];% the last row of Bset is for low and upper bound values
        
        bb_k=1;
        Blow = Bset; low_ind = 1; %index of the set whose lower bound is highest
        while (Upk(end)-Lowk(end)>eps_BB) & (bb_k<max_itr)
            if bb_k==17
                dfd=0;
            end
            %find the set which produces the lower bound previously
            %brunching 
            vec_low = Bset(end,1:2:end);% lower bound of each set, odd positions
            [temp,low_ind] = min(vec_low);            
            B_temp = Bset(1:end-1,2*(low_ind-1)+1:2*(low_ind-1)+2);%find the corresponding set 
            Bsetx = Bset(:,2*(low_ind-1)+1:2*(low_ind-1)+2);
            Bset(:,2*(low_ind-1)+1:2*(low_ind-1)+2) = [];%remove the set
            length_temp = B_temp(1:end,2)-B_temp(1:end,1);%the last row is not used
            [tempbb, ind_edge] = max(length_temp);%find the longest edge
            
            Bset_low_half = B_temp;  
            Bset_low_half(ind_edge, 2) = sum(Bset_low_half(ind_edge, :))/2;
            if feasibility(Bset_low_half(:,1),  Pmax,S, cp,dp,ep,fp, tjk, bjk, ck,R) %if xl is in the G
%                 if number_sij>1
                    xlow1 = flow_bound(Bset_low_half,number_sij,cp, ep,R,ck,fp,bjk,Pmax,tjk,S,dp    );
%                 else
%                     xlow1 = Bset_low_half(:,2);%do nothing
%                 end           
                Bset_low_half(:,2) = xlow1; %shrink the box
                Lowk1 = -fobj(xlow1 ,   weight, S   );% additional lower bound
                
                %find the upper bound
%                 proj1 = proj(Bset_low_half,Pmax,S,cp,dp,ep,fp, tjk, bjk, ck,R);
                up1_temp=0;
                for ij = 1 : number_sij
                    up_vec1 = Bset_low_half(:,1); %let be xmin
                    up_vec1(ij) = xlow1(ij); %
                    up1_temp(ij)= -fobj(up_vec1 ,   weight, S   );
                end
                Upk1 = min(up1_temp);%-fobj(Bset_low_half(:,1) ,   weight, S   );%additional upper bound
                %Upk1 =  -fobj(Bset_low_half(:,1) ,   weight, S   );%additional upper bound
%                 Upk1 =  -fobj(proj1 ,   weight, S   );%additional upper bound
           else
                Upk1 = 0;Lowk1 = 0;
            end
            
            Bset_up_half = B_temp;  Bset_up_half(ind_edge, 1) = sum(Bset_up_half(ind_edge, :))/2;
            if feasibility(Bset_up_half(:,1),  Pmax,S, cp,dp,ep,fp, tjk, bjk, ck,R) %if xl is in the G
%                 if number_sij>1
                    xlow2 = flow_bound(Bset_up_half,number_sij,cp, ep,R,ck,fp,bjk,Pmax,tjk,S ,dp   );
%                 else
%                     xlow2 = Bset_up_half(:,2); % do nothing
%                 end
                Bset_up_half(:,2) = xlow2;%shrink boxk
                Lowk2 = -fobj(xlow2 ,   weight, S   );% additional lower bound
                
                %find the upper bound
%                 proj2 = proj(Bset_up_half,Pmax,S,cp,dp,ep,fp, tjk, bjk, ck,R);
                up2_temp=0;
                for ij = 1 : number_sij
                    up_vec2 = Bset_up_half(:,1); %let be xmin
                    up_vec2(ij) = xlow2(ij); %
                    up2_temp(ij)= -fobj(up_vec2 ,   weight, S   );
                end
                Upk2 = min(up2_temp);%-fobj(Bset_up_half(:,1) ,   weight, S   );%additional upper bound
                %Upk2 =  -fobj(Bset_up_half(:,1) ,   weight, S   );%additional upper bound
%                  Upk2 =  -fobj(proj2 ,   weight, S   );%additional upper bound
            else
                Upk2 = 0;Lowk2 = 0;
            end            
            
            if (Lowk1<=Upk1) & (Lowk2<=Upk2) % if resolution issues happen, just use the current low
                Bset = [Bset [Bset_low_half; Lowk1 Upk1]]; %add this new set
                Bset = [Bset [Bset_up_half; Lowk2 Upk2]]; %add this new set  

                %Bounding 
                Lowknew = min(Bset(end,1:2:end));%the new lower bound
                Lowk = [Lowk Lowknew]; %  need to remember this index
                Upnew = min(Bset(end,2:2:end));%the new upper bound
                Upk = [Upk Upnew]; % need to remember this index
                %[Lowk; Upk]
                %prunning
                vec_pru = Bset(end,1:2:end);% lower bound of each set, odd positions
                ind_pru = find(vec_pru>Upk(end));
                if ~isempty(ind_pru)
                    Bset(:,[2*(ind_pru-1)+1 2*(ind_pru-1)+2])=[];
                end
                bb_k = bb_k+1;
            else
                Bset = [Bset Bsetx];%restore the set
                bb_k = max_itr;
            end
        end
        
        %recover the x
        [tempbb,opt_ind] = min(Bset(end,1:2:end));%find which box contributes the lower bound
        Boptimal = Bset(1:end-1,2*(opt_ind-1)+1:2*(opt_ind-1)+2);%extract the box
%        if number_sij>1
            xlow = Boptimal(:,2);%flow_bound(Boptimal,number_sij,cp, ep,R,ck,fp,bjk,Pmax,tjk,S ,dp   );
            up_temp=0; v_all = [];
            for ij = 1 : number_sij
                up_vec = Boptimal(:,1); %let be xmin
                up_vec(ij) = xlow(ij); %
                v_all = [v_all up_vec];
                up_temp(ij)= -fobj(up_vec ,   weight, S   );
            end
            [Upk2 vind] = min(up_temp);
            x=v_all(1:end,vind);
%         else
%             x = Boptimal(:,1);
%         end
        
        %recover y from x
        xsize = length(x);
        [t,ybb]=feasibility(x,  Pmax,S, cp,dp,ep,fp, tjk, bjk, ck,R);
        if isempty(ybb)%due to the resolution of Matlab, A might be treated as singular, and our analytical method might not be working robust
            [t,ybb]=feasibility(Boptimal(:,1),  Pmax,S, cp,dp,ep,fp, tjk, bjk, ck,R);
        end
        
%         ybb
%         structure
        a = find(ybb>0.01);
        b=[];
        for ia = 1 : length(a)
            b = [b; S(a(ia),2)];
        end
        c = unique(b);
        if length(b)~=length(c)
            sum4 = sum4+1;
        end

        ratebb(ict)=frate(ybb,xsize, cp,dpx,R,tjk, S   );
        
        % % % % %%%% use fmincon for SCA
        itr=1;
        while itr<=5
            fun =@sca_obj;
            A = zeros(2*number_sij, 2*number_sij); %for now, 9c and 9d
            b = zeros(3*number_sij, 1);
            tempmat = eye(number_sij); 
            for i = 1 : number_sij %S(ind(i),1) is j, S(ind(i),2) is k, i correspons to {j,k}        
                A(3*(i-1)+1,:) = [ ep(:,i)'*R zeros(1,number_sij)];
                b(3*(i-1)+1) = -ck(S(i,2));
                A(3*(i-1)+2,:) = [fp(:,i)'*R zeros(1,number_sij)];
                b(3*(i-1)+2) = -bjk(S(i,1),S(i,2));
                cxx = (dp(:,i)'*R*ym(1:number_sij)+ tjk(S(i,1),S(i,2)));
                abtemp = dp(:,i)'*R/log(2)/ cxx;
                A(3*(i-1)+3,:) = cxx*[abtemp -tempmat(i,:)];
                b(3*(i-1)+3) = cxx*(-log2(dp(:,i)'*R*ym(1:number_sij)+ tjk(S(i,1),S(i,2)) ) +abtemp*ym(1:number_sij));
            end
            A = [A; [ones(1,number_sij) zeros(1,number_sij)]];    b = [b; Pmax];  %add the power constraint
%             for i =1 : size(A,1)
%                 xcale = max(abs(A(i,:)) );
%                 A(i,:) = A(i,:)/xcale;
%                 b(i,:) = b(i,:)/xcale;
%             end
%             bxx = max(abs(b));
%             A = A/bxx; b = b/bxx;
            Aeq = [];
            beq = [];
            lb = [zeros(number_sij,1);-inf*ones(number_sij,1)];
            ub = [];

            y0 = ym;%ones(number_sij+1,1);

            options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
            y = fmincon(@(y) -sca_obj(y, cp, dp, R, tjk,weight,S,scale),y0,A,b,Aeq,beq,lb,ub,[], options);
            ym = y;
            %ym(1:number_sij);    
            %sca_obj(y, cp, dp, R, tjk,weight,S,scale);
            itr = itr+1;
        end
        yca2=y;
        ratesca2(ict)=frate(yca2(1:number_sij),number_sij, cp,dpx,R,tjk, S   );
        ratebas(ict)=frate(ybas,number_sij, cp,dpx,R,tjk, S   );
        %[ymc yca2(1:number_sij)]
        %[ratemc  ratesca2 ratebas]
        [iM ict]
    end
    Ratebb(iM) = mean(ratebb)
     Ratesca2(iM) = mean(ratesca2)
    Ratebas(iM) = mean(ratebas)
    errorx(iM) = sum4
    
end
plot(Mvec,Ratebas, Mvec, Ratesca2,  Mvec, ratebb )
% if max(A*y-b)>0
%     display('not feasible')
% end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define the rate function frate(x)
function [ratex] = frate(y,xsize, cp,dpx,R,tjk, S   )
    ratex=0;
    for i =1 : xsize
        if y(i)>0.01
            index=0;  mat2=[];
            ratex = ratex + log2(1+cp(:,i)'*y/(dpx(:,i)'*R*y+tjk(S(i,1),S(i,2))));
        end
    end
end

%define the objective function f(x)
function [f] = fobj(xtilde,   weight, S   )
    f=0;
    xsize = length(xtilde);
    for i = 1 : length(xtilde)
        %S
        %xtilde
        f = f + log2(1+xtilde(i));
        
    end 
end
 

%define the SCA objective function
function [f] = sca_obj(yz, cp, dp, R, tjk,weight,S, scale)
    f=0;    
    xsize = length(yz)/2;
    y = yz(1:xsize); z=yz(xsize+1:end);
    for p = 1 : xsize
        f = f +   log2(exp(1))*log(scale*(cp(:,p)'+dp(:,p)'*R)*y+ scale*tjk(S(p,1),S(p,2)) ) - log(scale)-z(p);
    end  
    
end

%bisection methods to find the projection on G 
function [temp] = proj(Box,Pmax,S,cp,dp,ep,fp, tjk, bjk, ck,R) 

     zk = Box(:,2);%upper corner
     xmin = Box(:,1);%lower corner, also the lower bound
     amin = 0; 
     amax = 1;
     delta = 0.001;%1/max(zk)/1;%0.00001;
     while amax-amin>=delta
         alpha = (amin+amax)/2;
         if alpha*zk>=xmin %the projection is inside of Box
             [temp1,y] =feasibility(alpha*zk,  Pmax,S, cp,dp,ep,fp, tjk, bjk, ck,R) ;% in G                 
             if  ~temp1% not feasible  
                amax = alpha;
             else %  feasible 
                amin = alpha;
             end     
         else
             amin = alpha;
         end
    end
    temp=amin*zk;%alpha*zk;
end
 
 


%bisection methods to find the projection on G 
function [temp,y] = feasibility(x,  Pmax,S, cp,dp,ep,fp, tjk, bjk, ck,R) 
    xsize = length(x);
    A = zeros(3*xsize, xsize);
    b = zeros(3*xsize, 1);
    for i = 1 : xsize %S(ind(i),1) is j, S(ind(i),2) is k, i correspons to {j,k}        
        A(3*(i-1)+1,:) =   - (cp(:,i)'-x(i)*dp(:,i)'*R)/tjk(S(i,1),S(i,2));
        b(3*(i-1)+1) = -x(i);
        A(3*(i-1)+2,:) = ep(:,i)'*R;
        b(3*(i-1)+2) = -ck(S(i,2));
        A(3*(i-1)+3,:) = sign(x(i))*fp(:,i)'*R; % this constraint is not active if x(i)=0
        b(3*(i-1)+3) = -sign(x(i))*bjk(S(i,1),S(i,2));% this constraint is not active if x(i)=0
    end
    
    A = [A; ones(1,xsize)];
    b = [b; Pmax];    
    Aeq = [];
    beq = [];
    lb = [zeros(xsize,1)];
    ub = [];
    options = optimoptions('linprog','Display', 'off');%, 'OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
    y = linprog([],A,b,Aeq,beq,lb,ub,options);
    %options = optimoptions('fmincon','Display', 'off');%,'OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
    %y0=zeros(xsize,1);
    %y = fmincon(@(y) 1,y0,A,b,Aeq,beq,lb,ub,[], options);
           
    if isempty(y)% strcmp(cvx_status, 'Solved') & min(y)>=0
        temp = 0; %feasible
    else
        temp = 1; %infeasible
    end
end

%finding a tighter lower bound  
function [xnew] = flow_bound(xbox,number_sij,cp, ep,R,ck,fp,bjk,Pmax,tjk,S ,dp   )
    xmin = xbox(1:number_sij,1);
    xmax = xbox(1:number_sij,2);
    for p =1 :  number_sij
        if xmax(p)==0
            xnew(p) = 0;
            continue; % if xmax(p) is already zero, no need to tighten it
        end
        A = zeros(2*number_sij, number_sij);
        b = zeros(2*number_sij, 1);
        for i = 1 : number_sij %S(ind(i),1) is j, S(ind(i),2) is k, i correspons to {j,k}        
            A(2*(i-1)+1,:) = [ep(:,i)'*R ];
            b(2*(i-1)+1) = -ck(S(i,2));
            if i==p
                A(2*(i-1)+2,:) =  fp(:,i)'*R ; % because we are increasing xmin(p), so xmin(p) should not be zer 
                b(2*(i-1)+2) = - bjk(S(i,1),S(i,2)); 
            else
                A(2*(i-1)+2,:) = sign(xmin(i))*fp(:,i)'*R ; % this constraint is not active if x(i)=0
                b(2*(i-1)+2) = -sign(xmin(i))*bjk(S(i,1),S(i,2));% this constraint is not active if x(i)=0
            end
        end
        A = [A; ones(1,number_sij)  ]; %power budget constraint
        b = [b; Pmax];    
        A = [A; (cp(:,p)'-xmax(p)*dp(:,p)'*R)/tjk(S(p,1),S(p,2)) ];
        b = [b; xmax(p)];        
       
        Aeq = [] ;
        beq = []; 
        for i =1 : number_sij %S(ind(i),1) is j, S(ind(i),2) is k, i correspons to {j,k}        
            if i==p
                continue; %skip if i=p
            end
            Aeq = [Aeq; (cp(:,i)'-xmin(i)*dp(:,i)'*R)/tjk(S(i,1),S(i,2))];
            beq = [beq;xmin(i)];
        end
                
%         lb = [zeros(number_sij,1)];
%         ub = [];
%         y0 = xmin;
%         options = optimoptions('fmincon','Display', 'off');%,'OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
%         y = fmincon(@(y) -cp(:,p)'/tjk(S(p,1),S(p,2))*y/(dp(:,p)'*R/tjk(S(p,1),S(p,2))*y+1),y0,A,b,Aeq,beq,lb,ub,[], options);
%             
% 
%         if isempty(y)
%             fdfd=0;
%         end
%         xnew1(p) = cp(:,p)'/tjk(S(p,1),S(p,2))*y/(dp(:,p)'*R/tjk(S(p,1),S(p,2))*y+1);
%         

        if number_sij>1
            warning('off')
            %normalize Aeq and beq
            for ia = 1 : size(Aeq,1)
                anorm = max( abs( Aeq(ia,:)));
                if anorm ==0
                    continue;
                end
                Aeq(ia,:) = Aeq(ia,:)/ anorm;
                beq(ia,:) = beq(ia,:)/ anorm;
            end
            Aeq_tilde =Aeq;  Aeq_tilde(:,p)=[];     Aeqinv=inv(Aeq_tilde);
            svd_A = svd(Aeq_tilde);
            if min(svd_A)<0.0001%use fmincon        
                lb = [zeros(number_sij,1)];
                ub = [];
                y0 = xmin;
                options = optimoptions('fmincon','Display', 'off');%,'OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
                y = fmincon(@(y) -cp(:,p)'/tjk(S(p,1),S(p,2))*y/(dp(:,p)'*R/tjk(S(p,1),S(p,2))*y+1),y0,A,b,Aeq,beq,lb,ub,[], options);
                xnew2(p) = cp(:,p)'/tjk(S(p,1),S(p,2))*y/(dp(:,p)'*R/tjk(S(p,1),S(p,2))*y+1);
            else%use analytical way
                A_tilde =A;  A_tilde(:,p)=[];  
                app = A(:,p); aeqpp = Aeq(:,p);
                vec_temp1 = b - A_tilde*Aeqinv*beq;
                vec_temp2 = app - A_tilde*Aeqinv*aeqpp; %maybe vec_temp2 has negative elements
                index_min = find(vec_temp2>0);%there are to use min
                index_max = find(vec_temp2<0);
                upper_bound = vec_temp1(index_min)./vec_temp2(index_min);%x<= these
                lower_bound = vec_temp1(index_max)./vec_temp2(index_max);%x> these
    
                yp = min(upper_bound);
                ynew_tilde = Aeqinv*(beq-aeqpp*yp);
                ynew = [ynew_tilde(1:p-1);yp;ynew_tilde(p:end)];
                xnew2(p) =  cp(:,p)'/tjk(S(p,1),S(p,2))*ynew/(dp(:,p)'*R/tjk(S(p,1),S(p,2))*ynew+1);
            end

            xnew(p)=xnew2(p);%max(xnew1(p),xnew2(p));
        else
            ynew = min(b./A);
            xnew(p) =  cp(:,p)'/tjk(S(p,1),S(p,2))*ynew/(dp(:,p)'*R/tjk(S(p,1),S(p,2))*ynew+1);
        end
        xnew(p) = max(0,xnew(p));%we assume that xmin(p) can be improved and add the constraint fp, but it is possible that not feasible, so force it zero if negative.
    end
end 
