function [Host_TS,Parasitoid_TS ] = Perc_ref_ref_lattice(L,l,g,eh,ep,t,p_location,p_density,h_indices)
%Looking at Host parasitoid system on a lattice with side length L
%growth rate of host = l ; growth rate of parasitoid = g;
%fraction of host migrating = eh; fraction of parasitoid migrating = ep;
% t = length of time series
%p_location = location of initial parasitoid
%p_density = initial fraction of max density
%h_indices = habitable indices, set by setdifference(total,uninhabitable sites)
%storing entire lattices for last 1000 timesteps
%reflecting boundaries, no movement/death into heterogeneities
A = zeros(L,L);
A(h_indices) = 1; %setting up the habitat landscape
%Data = cell(1,2);
% tsh = zeros(4,1000);
% tsp = zeros(4,1000);
Host = A*l; %initial host population
Parasitoid = zeros(L,L);
Parasitoid(p_location)= g*p_density; %select a habitat place and introduce parasitoid
Host_TS = cell(1,1000);
Parasitoid_TS = cell(1,1000);
for i=1:t
    %competition and parasitism
    H = l*min(Host,1).*exp(-Parasitoid);
    P = g*min(Host,1).*(1 - exp(-Parasitoid));

    %migration with reflecting boundary conditions
    EA = zeros(L+2,L+2);
    EA(2:L+1,2:L+1) = A;
    H1 = zeros(L,L);
    P1 = zeros(L,L);   %storing info about inmigration
    H2 = zeros(L,L);
    P2 = zeros(L,L);   %storing info about out inmigration
    for q=1:L
        for w=1:L
            if EA(q+1,w+1) ~=0 %habitat place
                if EA(q,w+1)*EA(q+2,w+1)*EA(q+1,w)*EA(q+1,w+2) == 1     %all neighbors occupied
                    H1(q,w) = eh*(H(q-1,w) + H(q+1,w) + H(q,w-1) + H(q,w+1))/4;
                    P1(q,w) = ep*(P(q-1,w) + P(q+1,w) + P(q,w-1) + P(q,w+1))/4; %migration into focal site
                    H2(q,w+1) = H2(q,w+1) -eh*H(q,w+1)/4;  %current - outmigration
                    H2(q,w-1) = H2(q,w-1) -eh*H(q,w-1)/4;  %current - outmigration
                    H2(q+1,w) = H2(q+1,w) -eh*H(q+1,w)/4;  %current - outmigration
                    H2(q-1,w) = H2(q-1,w) -eh*H(q-1,w)/4;  %current - outmigration
                    P2(q,w+1) = P2(q,w+1) -ep*P(q,w+1)/4;  %current - outmigration
                    P2(q,w-1) = P2(q,w-1) -ep*P(q,w-1)/4;  %current - outmigration
                    P2(q+1,w) = P2(q+1,w) -ep*P(q+1,w)/4;  %current - outmigration
                    P2(q-1,w) = P2(q-1,w) -ep*P(q-1,w)/4;  %current - outmigration
                else        %atleast one neighbor unoccupied
                   
                    if EA(q,w+1) == 1      %north neighbor               
                        H1(q,w) = H1(q,w) +  H(q-1,w)*eh/4; %current + in migration from north
                        P1(q,w) = P1(q,w) +  P(q-1,w)*ep/4; %current + in migration
                        H2(q-1,w) = H2(q-1,w) -H(q-1,w)*eh/4; %current - outmigration
                        P2(q-1,w) = P2(q-1,w) -P(q-1,w)*ep/4; %current - outmigration
                    end
                    if EA(q+2,w+1) == 1      %souths neighbor                       
                        H1(q,w) = H1(q,w) +  H(q+1,w)*eh/4; %current + in migration from south
                        P1(q,w) = P1(q,w) +  P(q+1,w)*ep/4; %current + in migration
                        H2(q+1,w) = H2(q+1,w) -H(q+1,w)*eh/4; %current - outmigration
                        P2(q+1,w) = P2(q+1,w) -P(q+1,w)*ep/4; %current - outmigration
                    end
                    if EA(q+1,w) == 1      %east neighbor              
                        H1(q,w) = H1(q,w) +  H(q,w-1)*eh/4; %current + in migration from east
                        P1(q,w) = P1(q,w) +  P(q,w-1)*ep/4; %current + in migration
                        H2(q,w-1) = H2(q,w-1) -H(q,w-1)*eh/4; %current - outmigration
                        P2(q,w-1) = P2(q,w-1) -P(q,w-1)*ep/4; %current - outmigration
                    end
                    if EA(q+1,w+2) == 1      %west neighbor
                        H1(q,w) = H1(q,w) +  H(q,w+1)*eh/4; %current + in migration from west
                        P1(q,w) = P1(q,w) +  P(q,w+1)*ep/4; %current + in migration
                        H2(q,w+1) = H2(q,w+1) -H(q,w+1)*eh/4; %current - outmigration
                        P2(q,w+1) = P2(q,w+1) -P(q,w+1)*ep/4; %current - outmigration
                    end
                end
            end
        end
    end
    Host = H + H1 +H2; %pre migration + in migration + out migration
    Parasitoid = P + P1 +P2;

    %Reproduction and noise
    %Host = l*(Host);
    Host(Host<=0) = 0;
    Host(A==0) = 0; %no insects where there's no habitat
    %Parasitoid = g*(Parasitoid);
    Parasitoid(Parasitoid <=0) = 0;
    Parasitoid(A==0)=0; %no insects where there's no habitat
    if i>t-1000
        
            Host_TS{i+1000-t} = Host;
            Parasitoid_TS{i+1000-t} = Parasitoid;
       
    end
end

end