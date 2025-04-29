function [] = parsave(fname,vari)
%saves a file inside parfor loop

save(fname,'vari',"-v7.3")
end