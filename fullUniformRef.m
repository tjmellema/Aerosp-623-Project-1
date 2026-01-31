% will do all uniform mesh stuff
function fullUniformRef(grifile)
    for i = 1:3
        newgrifile = [grifile(1:end-4),num2str(i),grifile(end-3:end)];
        uniformRefinement(grifile,i)
        mesh_verification(newgrifile)
    end
end