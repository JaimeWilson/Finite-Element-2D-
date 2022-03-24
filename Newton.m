%=========================================================================
% Método de Newton para solução do sistema nao linear
%=========================================================================

function [C] = Newton(NELE, NDoF, NEX, NEY, DomNodeID, ASSMtrx, X, Y, ro, mi, p_in, p_out, C0)
    
    Erro = 0.0001;
    itermax = 10;

    [iulDoF, ivlDoF, iplDoF] = GrauLiberdade();
    C = C0;
    DC = zeros(length(C0), 1);

    % Cálculo do Vetor Resíduo Global Rv
    [Rv] = FormRv(NELE, NDoF, NEX, NEY, DomNodeID, ASSMtrx, X, Y, iulDoF, ivlDoF, iplDoF, ro, mi, p_in, p_out, C);  
    Rvnorm = norm(Rv);

    iter = 0;

    iter
    Rvnorm

    itervec(1) = 1;
    Rvnormvec(1) = Rvnorm;

    while Rvnorm > Erro & iter < itermax
        
        iter = iter + 1;
        
        % Cálculo da Matriz Jacobiana Global
        [J] = FormJ(NELE, NDoF, NEX, NEY, DomNodeID, ASSMtrx, iulDoF, ivlDoF, iplDoF, X, Y, ro, mi, C); 
        b = -Rv;
 
        [LSP, USP] = lu(J);
        YLU = LSP\b;
        DC = USP\YLU;
        C = DC + C;
   
        % Cálculo do Vetor Resíduo Global
        [Rv] = FormRv(NELE, NDoF, NEX, NEY, DomNodeID, ASSMtrx, X, Y, iulDoF, ivlDoF, iplDoF, ro, mi, p_in, p_out, C);
        Rvnorm = norm(Rv);
   
        iter
        Rvnorm
   
        itervec(iter + 1) = iter + 1;
        Rvnormvec(iter + 1) = Rvnorm;

    end

    if norm(Rv) > Erro
        
        display(' Did not converge ');

    end
    
end


   
