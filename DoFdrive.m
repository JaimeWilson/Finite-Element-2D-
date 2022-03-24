%=========================================================================
% Relação entre a numeração Elementar e Global dos graus de liberdade
%=========================================================================

function [NDoF, ASSMtrx] = DoFdrive(NELE, NNODES, DomNodeID)

    ASSMtrx = zeros(21, NELE);
    GlobalDoFID = zeros(5, NNODES); % Matriz de graus de liberdade por nó

    % Definition of the array GlobalDoFID(ilDoF, ignode)
    igDoFcount = 1;
    
    for iele = 1 : NELE
        
        icount = igDoFcount;
        nadd = 0;
        ilDoFcount = 0;
        
        for ilnode = 1 : 9
            
            ignode = DomNodeID(ilnode, iele);
            
            if GlobalDoFID(1, ignode) == 0
                
                nadd = nadd + 1;
        
            end
            
        end
    
        for ilnode = 1 : 9
            
            ignode = DomNodeID(ilnode, iele);
            
            if GlobalDoFID(1, ignode) == 0
                
                GlobalDoFID(1, ignode) = icount;
                GlobalDoFID(2, ignode) = icount + nadd;
                
                icount = icount + 1;
                ilDoFcount = ilDoFcount + 2;
                    
            end
            
        end
        
        ignode = DomNodeID(9, iele);
        GlobalDoFID(3, ignode) = GlobalDoFID(2, ignode) + 1;
        GlobalDoFID(4, ignode) = GlobalDoFID(3, ignode) + 1;
        GlobalDoFID(5, ignode) = GlobalDoFID(4, ignode) + 1;
        
        ilDoFcount = ilDoFcount + 3;
        igDoFcount = igDoFcount + ilDoFcount;

    end
    
    NDoF = igDoFcount - 1;
    
    % Definition of Assembly Matrix
    for iele = 1 : NELE
        
        ilDoF = 0;
        
        for ilnode = 1 : 9
            
            ignode = DomNodeID(ilnode, iele);
            ilDoF = ilDoF + 1;
            ASSMtrx(ilDoF, iele) = GlobalDoFID(1, ignode);
         
        end
        
        for ilnode = 1 : 9
            
            ignode = DomNodeID(ilnode, iele);
            ilDoF = ilDoF + 1;
            ASSMtrx(ilDoF, iele) = GlobalDoFID(2, ignode);
              
        end   
    
        ignode = DomNodeID(9, iele);
        ilDoF = ilDoF + 1;
        ASSMtrx(ilDoF, iele) = GlobalDoFID(3, ignode);
        ilDoF = ilDoF + 1;
        ASSMtrx(ilDoF, iele) = GlobalDoFID(4, ignode);
        ilDoF = ilDoF + 1;
        ASSMtrx(ilDoF, iele) = GlobalDoFID(5, ignode);

    end
      
end