%=========================================================================
% Cálculo da Matriz Jacobiana Elementar
%=========================================================================

function [ElemJ] = GetElemJ(iele, NEX, NEY, DomNodeID, ASSMtrx, iulDoF, ivlDoF, iplDoF, X, Y, ro, mi, C)

    % Inicialização
    ElemJ = zeros(21, 21);
    
    % Integração Numérica - Quadratura Gaussiana
    NGP = 3; % Número de pontos de Gauss
    XGP = [-0.7745966 0.0 0.7745966]; % Posição dos pontos de Gauss
    WGP = [0.555 0.888 0.555]; % Pesos dos pontos de Gauss
    
    % Matriz auxiliar para localizar nós nas fronteiras
    iaux = [1 2 4 1; 5 6 7 8; 2 3 3 4];
    
    % Valor elementar das coordenadas dos nós, velocidade e pressão
    for ilnode = 1 : 9
        
        ignode = DomNodeID(ilnode, iele);
        XY(1, ilnode) = X(ignode);
        XY(2, ilnode) = Y(ignode);
    
    end
    
    for ilnode = 1 : 9
        
        velocity(1, ilnode) = C(ASSMtrx(ilnode, iele));
        velocity(2, ilnode) = C(ASSMtrx(ilnode + 9, iele));

    end   
    
    for ipress = 1 : 3
        
        pressure(ipress) = C(ASSMtrx(ipress + 18, iele));
    
    end
    
    % Loop nos pontos de Gauss
    for igp = 1 : NGP
        
        NCSI = XGP(igp);
        WI = WGP(igp);
        
        for jgp = 1 : NGP
            
            NETA = XGP(jgp);
            WJ = WGP(jgp);
            W = WI*WJ; 
            
            [Phi, GradPhi, Qui] = BasisFunc(NCSI, NETA);
            
            % Cálculo da Jacobiana da mudança de variáveis:
            JAC = GradPhi * XY';
            % JACinv = inversa da Jacobiana.
            JACinv = inv(JAC);
            % JACdet = determinante da Jacobiana.
            JACdet = det(JAC);
            % Cálculo do GradPhixy (derivada em relação à x e y das funções base) 
            GradPhixy = JACinv*GradPhi;
            % Cálculo de u,v (velocides) e pre(pressão)
            UV = velocity*Phi;
            Pre = pressure*Qui;
            % Derivada da velocidade em relação a xy
            dUV = velocity*GradPhixy';
            
            % Quantidade de movimento
            for ilnode = 1 : 9
                
                iu = iulDoF(ilnode);
                iv = ivlDoF(ilnode);
                
                for jlnode = 1 : 9
                    
                    ju = iulDoF(jlnode);
                    jv = ivlDoF(jlnode);
                    
                    ElemJ(iu, ju) = ElemJ(iu, ju) + ...
                                    W*JACdet*(ro*Phi(ilnode)*(Phi(jlnode)*dUV(1, 1) + UV(1)*GradPhixy(1, jlnode) + UV(2)*GradPhixy(2, jlnode)) + ...
                                              2*mi*GradPhixy(1, ilnode)*GradPhixy(1, jlnode) + ...
                                              mi*GradPhixy(2, ilnode)*GradPhixy(2, jlnode));
                    
                    ElemJ(iu, jv) = ElemJ(iu, jv) + W*JACdet*(ro*Phi(ilnode)*Phi(jlnode)*dUV(1, 2) + mi*GradPhixy(2, ilnode)*GradPhixy(1, jlnode));
               
                    ElemJ(iv, ju) = ElemJ(iv, ju) + W*JACdet*(ro*Phi(ilnode)*Phi(jlnode)*dUV(2, 1) + mi*GradPhixy(1, ilnode)*GradPhixy(2, jlnode));

                    ElemJ(iv, jv) = ElemJ(iv, jv) + ...
                                    W*JACdet*(ro*Phi(ilnode)*(Phi(jlnode)*dUV(2, 2) + UV(1)*GradPhixy(1, jlnode) + UV(2)*GradPhixy(2, jlnode)) + ...
                                              2*mi*GradPhixy(2, ilnode)*GradPhixy(2, jlnode) + ...
                                              mi*GradPhixy(1, ilnode)*GradPhixy(1, jlnode));

                end
         
                for jpress = 1 : 3
                    
                    jp = iplDoF(jpress);
                    
                    ElemJ(iu, jp) = ElemJ(iu, jp) + W*JACdet*(-GradPhixy(1, ilnode)*Qui(jpress));
            
                    ElemJ(iv, jp) = ElemJ(iv, jp) + W*JACdet*(-GradPhixy(2, ilnode)*Qui(jpress));
         
                end
                
            end
            
            % Conservação de massa
            for ipress = 1 : 3
                
                ip = iplDoF(ipress);
                
                for jlnode = 1 : 9
                    
                    ju = iulDoF(jlnode);
                    jv = ivlDoF(jlnode);
                    
                    ElemJ(ip, ju) = ElemJ(ip, ju) + W*JACdet*(GradPhixy(1, jlnode)*Qui(ipress));
            
                    ElemJ(ip, jv) = ElemJ(ip, jv) + W*JACdet*(GradPhixy(2, jlnode)*Qui(ipress));   
         
                end
                
            end
            
        end
        
    end
    
    % Condições de contorno

    % Elemento na fronteira West (4)
    if iele <= NEY
        
        % Loop nos pontos de Gauss
        NCSI = -1;
        
        for jgp = 1 : NGP
            
            NETA = XGP(jgp);
            WJ = WGP(jgp);
      
            [Phi, GradPhi, Qui] = BasisFunc(NCSI, NETA);
            
            % Cálculo da Jacobiana de mudança de variáveis:
            JAC = GradPhi*XY';
            % JACinv = inversa da Jacobiana.
            JACinv = inv(JAC);
            % JACdet = determinante da Jacobiana.
            JACdet = det(JAC);
            % Cálculo do GradPhixy 
            GradPhixy = JACinv*GradPhi;
            % Cálculo de u,v (velocides) e pre(pressão)
            UV = velocity*Phi;
            Pre = pressure*Qui;
            % Derivada da velocidade em relação a xy
            dUV = velocity*GradPhixy';
            
            for isidenode = 1 : 3
                
                ilnode = iaux(isidenode, 4);
                iu = iulDoF(ilnode);
                iv = ivlDoF(ilnode);
                 
                for jlnode = 1 : 9
                    
                    ju = iulDoF(jlnode);
                    jv = ivlDoF(jlnode);
   
                    ElemJ(iu, ju) = ElemJ(iu, ju) + WJ*(2*mi*GradPhixy(1, jlnode))*Phi(ilnode)*JAC(2, 2);
                    ElemJ(iv, ju) = ElemJ(iv, ju) + 0;
                                        
                    ElemJ(iv, ju) = ElemJ(iv, ju) + WJ*(mi*GradPhixy(2, jlnode))*Phi(ilnode)*JAC(2, 2);   
                    ElemJ(iv, jv) = ElemJ(iv, jv) + WJ*(mi*GradPhixy(1, jlnode))*Phi(ilnode)*JAC(2, 2); 
                
                end
                
            end
            
        end
        
    end
    
    % Elemento na fronteira East (2)
    if iele > (NEX - 1)*NEY
        
        % Loop nos pontos de Gauss
        NCSI = 1;
        
        for jgp = 1 : NGP
            
            NETA = XGP(jgp);
            WJ = WGP(jgp);
     
            [Phi, GradPhi, Qui] = BasisFunc(NCSI, NETA);
            
            % Cálculo da Jacobiana de mudança de variáveis:
            JAC = GradPhi*XY';
            % JACinv = inversa da Jacobiana.
            JACinv = inv(JAC);
            % JACdet = determinante da Jacobiana.
            JACdet = det(JAC);
            % Cálculo do GradPhixy 
            GradPhixy = JACinv*GradPhi;
            % Calcular de u,v (velocides) e pre(pressão)
            UV = velocity*Phi;
            Pre = pressure*Qui;
            % Derivada da velocidade em relação a xy
            dUV = velocity*GradPhixy'; 
            
            for isidenode = 1 : 3
                
                ilnode = iaux(isidenode, 2);
                iu = iulDoF(ilnode);
                iv = ivlDoF(ilnode);
                
                for jlnode = 1 : 9
                    
                    ju = iulDoF(jlnode);
                    jv = ivlDoF(jlnode);
     
                    ElemJ(iu, ju) = ElemJ(iu, ju) - WJ*(2*mi*GradPhixy(1, jlnode))*Phi(ilnode)*JAC(2, 2);
                    ElemJ(iv, ju) = ElemJ(iv, ju) - 0;
                    
                    ElemJ(iv, ju) = ElemJ(iv, ju) - WJ*(mi*GradPhixy(2, jlnode))*Phi(ilnode)*JAC(2, 2); 
                    ElemJ(iv, jv) = ElemJ(iv, jv) - WJ*(mi*GradPhixy(1, jlnode))*Phi(ilnode)*JAC(2, 2); 

                end
                
            end
            
        end
        
    end
    
    % Elemento na fronteira North (3)
    if mod(iele, NEY) == 0
        
        for isidenode = 1 : 3
            
            ilnode = iaux(isidenode, 3);
            iu = iulDoF(ilnode);
            iv = ivlDoF(ilnode);
            
            for j = 1 : 21
                
                ElemJ(iu, j) = 0;
                ElemJ(iv, j) = 0;
      
            end
            
            ElemJ(iu, iu) = 1;
            ElemJ(iv, iv) = 1;
   
        end
        
    end
    
    % Elemento na fronteira South (1)
    if mod(iele, NEY) == 1
        
        for isidenode = 1 : 3
            
            ilnode = iaux(isidenode, 1);
            iv = ivlDoF(ilnode);
            
            for j = 1 : 21
                
                ElemJ(iv, j) = 0;
      
            end
            
            ElemJ(iv, iv) = 1;
   
        end
        
    end
    
end