function [ yini age_arr S0_imm initS] = make_ics_naive( pars, arrSlu, arrIlu, arrCIlu, age)
%make_ics Summary of this function goes here
% Setep the initial conditions with everyone susceptible
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Complete naive individuals
% initial values are set to be 0 <<-??
yini = zeros(1,pars.novars);
%pop = [16.42 30.29 39.4 13.2]*pars.N;
pop = pars.demographic*pars.N;
age_arr =  get_popweights_4_hk(pars,pop); %age groups ratio
%age_arr = pars.demographic; 
seed_arr = (age_arr.*pars.seed)./pars.N; 

%%%%vvvvvv 20150206
%initS = initial_prev_Xu_uniform(pars);
%initS = initial_prev_Wu_expon();

%%%Default 20150307
if isfield(pars, 'inittitres_flag')
if pars.inittitres_flag == 1 %with 2% seropositive, average 2.266 seroprevalence
    initS = initial_prev_Wu_uniform(pars);
end
if pars.inittitres_flag == 0                       %compelte susceptible population
    initS = zeros(pars.maxa, pars.maxi);
    initS(:,1) = 1;
end
if pars.inittitres_flag == 2                       %compelte susceptible population
    initS = initial_prev_Wu_uniform(pars,1);
end
if pars.inittitres_flag == 3                       %manual setting partial protected population
    initS = initial_prev_Wu_uniform(pars,1);
end
else %if there is no information of inittitres_flag
    initS = initial_prev_Wu_uniform(pars);
end


%initS = initial_naive(pars);
initS0_imm = zeros(1,pars.maxa);
S0_imm = zeros(1,pars.maxa);
for a = 1:pars.maxa
  initS0_imm(a) = initS(a,1).*pars.PUAb;
  S0_imm(a) = (age_arr(a) - seed_arr(a))*initS0_imm(a);  %S(age)*Ti(i)*PUAB
end

% initialize susceptible numbers
for a=1:pars.maxa
    for i=1:pars.maxi
        for j=1:pars.maxj
            for k=1:pars.maxk
                %yini(arrSlu(a,i,j,k)) = (age_arr(a) - seed_arr(a))*initial_prev_naive(pars,a,i,j,k); %not finished yet
                yini(arrSlu(a,i,j,k)) = (age_arr(a) - seed_arr(a))*initS(a,i); %only for 1st strain
            end
        end
    end
end




% initialize infected numbers
for a=1:pars.maxa
    for X=1:pars.maxX % this will cause total number not 1
        %X = 1;
        %yini(arrIlu(X,a,1,1,1)) = seed_arr(a)./pars.maxa; %%%should I divided by pars.maxa???
        yini(arrIlu(X,a,1,1,1)) = seed_arr(a);
    end
end


% naive seroprevalence 
function inits = initial_naive(par)
        
        inits = zeros(par.maxa,par.maxi);
        for a1=1:par.maxa
          for i1=1:par.maxi
             if i1==1
               inits(a1,i1) = 1;
             end
          end
        end
end


% naive seroprevalence 
function init_prev = initial_prev_Xu()
        init_prev = []; 
        AbBoostRate = 3;
        poi = poisspdf([0:pars.maxi-1],AbBoostRate);
        %seroconvert = [sum(poi(1:3)) sum(poi(4:end))];
        for a=1:4
            if a==4
                seroconvert_obs = 0.02;
            end
            if a~=4
                seroconvert_obs = 0.01;
            end
            ratio = seroconvert_obs/sum(poi(4:end));
            naive = 1 - sum(poi)*ratio;
            init_prev(a,:) = poi*(ratio);
            init_prev(a,1) = init_prev(a,1)+naive;        
        end
end


% naive seroprevalence 
function init_prev = initial_prev_Xu3()
        init_prev = []; 
        AbBoostRate = 3;
        poi = poisspdf([0:pars.maxi-1],AbBoostRate);
        %seroconvert = [sum(poi(1:3)) sum(poi(4:end))];
        for a=1:4
            if a==4
                seroconvert_obs = 0.02;
            end
            if a~=4
                seroconvert_obs = 0.01;
            end
            ratio = seroconvert_obs/sum(poi(4:end));
            naive = 1 - sum(poi)*ratio;
            init_prev(a,:) = poi*(ratio);
            init_prev(a,1) = init_prev(a,1)+naive;        
        end
end


function init_prev = initial_prev_johnson_2009(pa,a,i,j,k)
        init_prev = [];
        if(a==1)
            prev = 25/58;
            if i==2 
                init_prev = prev;
            elseif i==1
                init_prev = 1-prev;
            else
                init_prev = 0;
            end
        elseif(a==2)
            prev = 0.25;
            if i==2
                init_prev = prev;
            elseif i==1
                init_prev = 1-prev;
            else
                init_prev = 0;
            end
        elseif(a==3)
            prev = 0.18;
            if i==2
                init_prev = prev;
            elseif i==1
                init_prev = 1-prev;
            else
                init_prev = 0;
            end
        elseif(a==4)
            prev = 0.2;
            if i==2
                init_prev = prev;
            elseif i==1
                init_prev = 1-prev;
            else
                init_prev = 0;
            end
        elseif(a==5)
            prev = 0.30;
            if i==2
                init_prev = prev;
            elseif i==1
                init_prev = 1-prev;
            else
                init_prev = 0;
            end
        else
            error('Problem in initial_prev_johnson_2009');
        end
end

end