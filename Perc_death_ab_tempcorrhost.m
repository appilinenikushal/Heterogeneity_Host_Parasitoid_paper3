function [mean_temp_corr] = Perc_death_ab_tempcorrhost(L,l,g,eh,ep,t,p_location,p_density,h_indices,sites)
%Looking at Host parasitoid system on a lattice with side length L
%growth rate of host = l ; growth rate of parasitoid = g;
%fraction of host migrating = eh; fraction of parasitoid migrating = ep;
% t = length of time series
%p_location = location of initial parasitoid
%p_density = initial fraction of max density
%h_indices = habitable indices, set by setdifference(total,uninhabitable sites)
%absorbing boundaries, death at heterogeneities
%calculating pairwise temporal corelation for 1000 pairs 
A = zeros(L,L);
A(h_indices) = 1; %setting up the habitat landscape
%Data = cell(1,2);
%tsh = zeros(4,1000);
%tsp = zeros(4,1000);
Host = A*l; %initial host population
Parasitoid = zeros(L,L);
Parasitoid(p_location)= g*p_density; %select a habitat place and introduce parasitoid
tempcorr_ts = cell(1,100); %storing timeseries for each of the 1000
for i=1:100
    tempcorr_ts{i} = zeros(2,1000); %storing timeseries of the last 1000 for a 2 adjacent sites
end
for i=1:t
    %competition and parasitism
    H = l*min(Host,1).*exp(-Parasitoid);
    P = g*min(Host,1).*(1 - exp(-Parasitoid));

    %migration with absorbing boundary conditions
    EH = zeros(L+2,L+2);
    EP = zeros(L+2,L+2);
    EH(2:L+1,2:L+1) = H;
    EP(2:L+1,2:L+1) = P;
    for q=1:L
        for w=1:L
            Host(q,w) = (1-eh)*H(q,w) + eh*(EH(q,w+1) + EH(q+2,w+1) + EH(q+1,w) + EH(q+1,w+2))/4;
            Parasitoid(q,w) = (1-ep)*P(q,w) + ep*(EP(q,w+1) + EP(q+2,w+1) + EP(q+1,w) + EP(q+1,w+2))/4;
        end
    end

    %Reproduction and noise
    %Host = l*(Host);
    Host(Host<=0) = 0;
    Host(A==0) = 0; %no insects where there's no habitat
    %Parasitoid = g*(Parasitoid);
    Parasitoid(Parasitoid <=0) = 0;
    Parasitoid(A==0)=0; %no insects where there's no habitat


    if i>t-1000
    %     tsh(1,i+1000-t) = Host(sites(1,1),sites(1,2)); % storing host ts for 2 occupied
    %     tsh(2,i+1000-t) = Host(sites(2,1),sites(2,2)); % % storing host ts for 3 occupied
    %     tsh(3,i+1000-t) = Host(sites(3,1),sites(3,2)); % storing host ts for 4 occupied
    %     tsh(4,i+1000-t) = sum(sum(Host))/(numel(h_indices)); %overall
    % 
    %     tsp(1,i+1000-t) = Parasitoid(sites(1,1),sites(1,2)); % storing parasitoid ts for 2 occupied
    %     tsp(2,i+1000-t) = Parasitoid(sites(2,1),sites(2,2)); % % storing host ts for 3 occupied
    %     tsp(3,i+1000-t) = Parasitoid(sites(3,1),sites(3,2)); % storing host ts for 4 occupied
    %     tsp(4,i+1000-t) = sum(sum(Parasitoid))/(numel(h_indices)); %overall
        for rt=1:100
            we = tempcorr_ts{rt};
            we(1,i+1000-t) = Host(sites(rt,1),sites(rt,2)); %storing timeseries after steady state
            we(2,i+1000-t) = Host(sites(rt,3),sites(rt,4));
            tempcorr_ts{rt} = we;
        end
        
    end
end
temp_corr = zeros(1,100);
for i=1:100
    qw = tempcorr_ts{i}; %getting time series
    R = corrcoef(qw(1,:),qw(2,:));
    temp_corr(1,i) = R(1,2);
end
mean_temp_corr = mean(temp_corr);
end