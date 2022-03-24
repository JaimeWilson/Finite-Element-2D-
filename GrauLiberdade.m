%=========================================================================
% Rela��o entre a numera��o Elementar e Global dos graus de liberdade
%=========================================================================

function [iulDoF, ivlDoF, iplDoF] = GrauLiberdade()

    % Grau de liberdade local da velocidade (u e v)
    for i = 1 : 9
        
        iulDoF(i) = i;
        ivlDoF(i) = i + 9;

    end
    
    % Grau de liberdade da press�o
    for i = 1 : 3
        
        iplDoF(i) = i + 18;

    end
    
end

