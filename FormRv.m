%=========================================================================
% Cálculo do Vetor Resíduo Global
%=========================================================================

function [Rv] = FormRv(NELE, NDoF, NEX, NEY, DomNodeID, ASSMtrx, X, Y, iulDoF, ivlDoF, iplDoF, ro, mi, p_in, p_out, C)

    NDoFele = 21;
    
    Rv = zeros(NDoF, 1);
    
    for iele = 1 : NELE
        
        [ElemRv] = GetElemRv(iele, NEX, NEY, DomNodeID, ASSMtrx, X, Y, iulDoF, ivlDoF, iplDoF, ro, mi, p_in, p_out, C);
        
        for ilDoF = 1 : NDoFele
            
            igDoF = ASSMtrx(ilDoF, iele);
            
            Rv(igDoF) = Rv(igDoF) + ElemRv(ilDoF);
   
        end
        
    end
    
end