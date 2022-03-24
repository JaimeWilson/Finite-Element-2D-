%=========================================================================
% Cálculo da Vazão de Entrada Qe e da Vazão de Saída Qs
%=========================================================================

function [Qe, Qs] = Vazao(L, H, X, Y, NEX, NEY, NELE, DomNodeID, ASSMtrx, C)
    
    Qe = 0;
    Qs = 0;
    
    % Matriz auxiliar para localizar nós nas fronteiras
    iaux = [1 2 4 1; 5 6 7 8; 2 3 3 4];
    counts = 1;
    counte = 1;
    
    for iele = 1 : NELE
        
        % Elemento na fronteira East (4)
        if iele <= NEY
            
            for isidenode = 1 : 3
                
                ilDoF = iaux(isidenode, 4);
                Vxe(counte) = C(ASSMtrx(ilDoF, iele));
                Hxe(counte) = X(DomNodeID(ilDoF, iele));
                Hye(counte) = Y(DomNodeID(ilDoF, iele));
                Qe = Qe + (1/3)*(H/NEY)*Vxe(counte);
                counte = counte + 1;
 
            end
            
        end
        
        % Elemento na fronteira East (2)
        if iele > (NEX - 1)*NEY
            
            for isidenode = 1 : 3
                
                ilDoF = iaux(isidenode, 2);
                Vxs(counts) = C(ASSMtrx(ilDoF, iele));
                Hxs(counts) = X(DomNodeID(ilDoF, iele));
                Hys(counts) = Y(DomNodeID(ilDoF, iele));
                Qs = Qs + (1/3)*(H/NEY)*Vxs(counts);
                counts = counts + 1;
  
            end
            
        end
        
    end
    
    figure('Name', 'Vazão', 'NumberTitle','off'), grid on, hold on
    quiver(Hxe, Hye, Vxe, zeros(1, length(Vxe)),'b', 'LineWidth', 1.5)
    axis([-0.011*L 0.025*L -0.1*H 1.1*H])
    title(['VAZÃO: Q = ', num2str(Qe)])
    xlabel('u')
    ylabel('Y')
    
% %     figure(revisar)
% %     quiver(Hxs, Hys, Vxs, zeros(1, length(Vxs)))
% %     axis([revisar]) 
% %     title('VAZÃO')

end


