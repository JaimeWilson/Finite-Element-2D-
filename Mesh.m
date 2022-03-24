%=========================================================================
% Cálculo de Coordenadas da Malha - problema 2D
%=========================================================================

function [X, Y] = Mesh(L, H, NEX, NEY, DomNodeID)

    DX = L/NEX;
    DY = H/NEY;
    
    for ix = 1 : NEX
        
        for jy = 1 : NEY
            
            iele = (ix - 1)*NEY + jy;
            
            xloc(1) = (ix - 1)*DX;
            xloc(2) = xloc(1) + DX;
            xloc(3) = xloc(2);
            xloc(4) = xloc(1);
            xloc(5) = xloc(1) + DX/2;
            xloc(6) = xloc(2);
            xloc(7) = xloc(5);
            xloc(8) = xloc(1);
            xloc(9) = xloc(5);
            
            yloc(1) = (3/4)*(2*(jy - 1)*(DY/2)) + (1/4)*(2*(jy - 1)*(DY/2))*cos(2*pi*xloc(1)/(L));
            yloc(2) = (3/4)*(2*(jy - 1)*(DY/2)) + (1/4)*(2*(jy - 1)*(DY/2))*cos(2*pi*xloc(2)/(L));
            yloc(3) = (3/4)*((2*jy)*(DY/2)) + (1/4)*((2*jy)*(DY/2))*cos(2*pi*xloc(3)/(L));
            yloc(4) = (3/4)*((2*jy)*(DY/2)) + (1/4)*((2*jy)*(DY/2))*cos(2*pi*xloc(4)/(L));
            yloc(5) = (3/4)*(2*(jy - 1)*(DY/2)) + (1/4)*(2*(jy - 1)*(DY/2))*cos(2*pi*xloc(5)/(L));
            yloc(6) = (3/4)*((2*jy - 1)*(DY/2)) + (1/4)*((2*jy - 1)*(DY/2))*cos(2*pi*xloc(6)/(L));
            yloc(7) = (3/4)*((2*jy)*(DY/2)) + (1/4)*((2*jy)*(DY/2))*cos(2*pi*xloc(7)/(L));
            yloc(8) = (3/4)*((2*jy - 1)*(DY/2)) + (1/4)*((2*jy - 1)*(DY/2))*cos(2*pi*xloc(8)/(L));
            yloc(9) = (3/4)*((2*jy - 1)*(DY/2)) + (1/4)*((2*jy - 1)*(DY/2))*cos(2*pi*xloc(9)/(L));
                    
            for ilnode = 1 : 9
                
                ignode = DomNodeID(ilnode, iele);
                X(ignode) = xloc(ilnode);
                Y(ignode) = yloc(ilnode);
                
            end
               
        end
           
    end

end
        