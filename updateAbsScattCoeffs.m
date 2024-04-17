function [abs_map,scatt_map] = updateAbsScattCoeffs(abs_map,scatt_map,key,abs_coeffs,scatt_coefs,Vertex1,Vertex2)
    norm = Vertex2-Vertex1;
    if any(norm<0)
        disp('Not valid vertex positions')
        return 
    else
        Idx_X1 = round(100*Vertex1(1))+1;
        Idx_Y1 = round(100*Vertex1(2))+1;
        Idx_X2 = round(100*Vertex2(1));
        Idx_Y2 = round(100*Vertex2(2));

        try
            abs_map{key}(Idx_X1:Idx_X2,Idx_Y1:Idx_Y2,:) = reshape(abs_coeffs,1,1,length(abs_coeffs)).*ones(1+Idx_X2-Idx_X1,1+Idx_Y2-Idx_Y1,length(abs_coeffs));
            scatt_map{key}(Idx_X1:Idx_X2,Idx_Y1:Idx_Y2,:) = reshape(scatt_coefs,1,1,length(scatt_coefs)).*ones(1+Idx_X2-Idx_X1,1+Idx_Y2-Idx_Y1,length(scatt_coefs));
        catch
            disp('The vertex positions exceeds the dimensions of that surface')
            return
        end
    end
end