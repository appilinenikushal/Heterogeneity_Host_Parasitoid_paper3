function [] = PH_OG_death(N,l,g,it)
op=load("all_info.mat");
initialconfiguration = op.all_info; % initial configuration
%it is the iteration number, for cluster computing
%N = lattice size
%l,g growth rates
initial = 1+ 2*96*(it-1);
final = initial + 2*96-1;
Data = cell(1,192);
    parfor i=initial:final       %no of realizations, setup the lattice with different heterogeniety 
        Subdata = cell(4,10);
        e = [0.05 0.075 0.1 0.2 0.3 0.4 0.425 0.45 0.475 0.5];  %host and parasitoid migration we look into
      %  het = [0,0.01, 0.05,0.1]; %3 heterogeneity values (0,0.01, 0.05,0.1)
        p_location = initialconfiguration{i}{1}; %location of parasitoid
        p_density = initialconfiguration{i}{2};  %initial density
        h_ind = cell(1,4); %storing h_indices for each heterogeneity fraction
        sites = cell(1,4); %storing 3 sites for each heterogeneity fraction
        op = 1:128*128;
        h_ind{1} = op';                                     %0%
        sites{1} = initialconfiguration{i}{6}; 
        h_ind{2} = setdiff(op',initialconfiguration{i}{3}); %1%
        sites{2} = initialconfiguration{i}{7}; 
        h_ind{3} = setdiff(op',initialconfiguration{i}{4}); %5%
        sites{3} = initialconfiguration{i}{8}; 
        h_ind{4} = setdiff(op',initialconfiguration{i}{5}); %10%
        sites{4} = initialconfiguration{i}{9}; 
        for j=1:4 %simulating diff heterogeneity
            for k1=1:10  %population migration
                 A = Perc_death_ab(N,l,g,e(k1),e(k1),5000,p_location,p_density,h_ind{j},sites{j}); 
                 Subdata{j,k1} = A;
            end
        end 
        Data{i} = Subdata;
    end
    j1 = sprintf('OG_death_%d_%d_%d_%d_%d.mat',N,l,g,initial,final); %save the data file for 192 iterations
    parsave(j1,Data)  %save each data file
end