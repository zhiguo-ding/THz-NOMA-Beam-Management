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
xi=0;%100000000;%1.0000e+15;
taux = 0.0000001;tauy = 50; epsilon_cbv=0.01;
al=2; %path loss
c0=3*10^8; fc = 300*10^9; rl = 5*exp(-3);%molecular absorption
ct=5000;
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
        DBF = inv(Habf)*sqrt(D_power);%/sqrt(K); %normalized digital beamforming - trace(DBF*DBF')=1

        for k = 1 : K
            beam(:,k) = ABF*DBF(:,k);
        end
        %ABF=beam; %beam is the composite precoding

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
        ratebas(ict)=frate(ybas,number_sij, cp,dpx,R,tjk, S   );
                
        tempz = [];
        for p =1 : number_sij
            rhotemp = min([Pmax ; -ck(S(p,2)); -bjk(S(p,1),S(p,2))]);
            tempz(p) = log2(1+hSjk(S(p,1),S(p,2))*rhotemp/tjk(S(p,1),S(p,2)));
        end
        [x,indz]=max(tempz);
        ymz = zeros(number_sij,1);
        ymz(indz)=min([Pmax ; -ck(S(indz,2)); -bjk(S(indz,1),S(indz,2))]);
        ratebas_new(ict)=frate(ymz,number_sij, cp,dpx,R,tjk, S   ); 
        
        if ratebas_new(ict)~=ratebas(ict)
            dfd=0;
        end
        [iM ict]
    end 
    Ratebas(iM) = mean(ratebas)
    Ratebas_new(iM) = mean(ratebas_new)
    
    
end
plot(Mvec,Ratebas,Mvec,Ratebas_new  )
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
