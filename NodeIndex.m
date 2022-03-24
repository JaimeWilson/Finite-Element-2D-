%=========================================================================
% Relação entre a numeração Elementar e Global dos nós
%=========================================================================

function [NELE, NNODES, DomNodeID] = NodeIndex(NEX, NEY)

    NELE = NEX*NEY;
    NNODES = (2*NEX + 1)*(2*NEY + 1);
    
    DomNodeID = zeros(9, NELE);
    NodeCount = 1;
    
    for iele = 1 : NELE
        
        iWestEle = 0; % Neighbour West?
        iSouthEle = 0; % Neighbour South?
    
        % Find West element neighbour
        if iele > NEY
            
            iWestEle = iele - NEY;
    
        end
    
        % Find South element neighbour
        iSouth = mod(iele - 1, NEY);
        
        if iSouth == 0
            
            iSouthEle = 0;
    
        else
            
            iSouthEle = iele - 1;
    
        end
    
        % Set DomNodeID for nodes 1, 4, 8 if element has West Neighbour
        if iWestEle ~= 0
            
            DomNodeID(1, iele) = DomNodeID(2, iWestEle);
            DomNodeID(4, iele) = DomNodeID(3, iWestEle);
            DomNodeID(8, iele) = DomNodeID(6, iWestEle); 
    
        end
    
        % Set DomNodeID for nodes 1, 2, 5 if element has South Neighbour
        if iSouthEle ~= 0 
            
            DomNodeID(1, iele) = DomNodeID(4, iSouthEle);
            DomNodeID(2, iele) = DomNodeID(3, iSouthEle);
            DomNodeID(5, iele) = DomNodeID(7, iSouthEle);   
        
        end    
    
        for ilnode = 1 : 9
            
            if DomNodeID(ilnode, iele) == 0
                
                DomNodeID(ilnode, iele) = NodeCount;
                NodeCount = NodeCount + 1;
       
            end
            
        end
        
    end
    
end
    
    
   
        