%=========================================================================
% Gráfico das Velocidades 
%=========================================================================

function Graph(L, H, NELE, DomNodeID, ASSMtrx, X, Y, C)

    for iele = 1 : NELE
        
        for ilnode = 1 : 9
            
            ignode = DomNodeID(ilnode, iele);
            VXG(ignode) = C(ASSMtrx(ilnode, iele));
            VYG(ignode) = C(ASSMtrx(ilnode + 9, iele));
   
        end
        
    end
    
    figure('Name', 'Velocidade', 'NumberTitle','off'), grid on, hold on
    axis([-0.1*L 1.1*L -0.1*H 1.1*H])
    title('VELOCIDADE')
    xlabel('X')
    ylabel('Y')
    quiver(X, Y, VXG, VYG, 'b', 'LineWidth', 1.2);

end