%=========================================================================
% Cálculo da Matriz Jacobiana Global
%=========================================================================

function [J] = FormJ(NELE, NDoF, NEX, NEY, DomNodeID, ASSMtrx, iulDoF, ivlDoF, iplDoF, X, Y, ro, mi, C)

    NDoFele = 21;
    NDIM = NDoFele*NELE;
    
    Jval = zeros(NDIM, 1);
    irow = zeros(NDIM, 1);
    jcol = zeros(NDIM, 1);
    
    icont = 1;
    
    for iele = 1 : NELE
        
        [ElemJ] = GetElemJ(iele, NEX, NEY, DomNodeID, ASSMtrx, iulDoF, ivlDoF, iplDoF, X, Y, ro, mi, C);
        
        for ilDoF = 1 : NDoFele
            
            igDoF = ASSMtrx(ilDoF, iele);
            
            for jlDoF = 1 : NDoFele
                
                jgDoF = ASSMtrx(jlDoF, iele);
                irow(icont) = igDoF;
                jcol(icont) = jgDoF;
                Jval(icont) = ElemJ(ilDoF, jlDoF);
                icont = icont + 1;
      
            end
            
        end
        
    end
    
    J = sparse(irow, jcol, Jval, NDoF, NDoF);
    
    spy(J)
    title('Padrão de Dispersão')
   
end
