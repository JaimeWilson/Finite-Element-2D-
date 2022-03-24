%=========================================================================
% C�lculo do Vetor Res�duo Elementar
%=========================================================================

function [ElemRv] = GetElemRv(iele, NEX, NEY, DomNodeID, ASSMtrx, X, Y, iulDoF, ivlDoF, iplDoF, ro, mi, p_in, p_out, C)

    % Inicializa��o
    ElemRv = zeros(21, 1);
    
    % Integra��o Num�rica - Quadratura Gaussiana
    NGP = 3; % N�mero de pontos de Gauss
    XGP = [-0.7745966 0.0 0.7745966]; % Posi��o dos pontos de Gauss
    WGP = [0.555 0.888 0.555]; % Pesos dos pontos de Gauss
    
    % Matriz auxiliar para localizar n�s nas fronteiras
    iaux = [1 2 4 1; 5 6 7 8; 2 3 3 4];
    
    % Valor elementar das coordenadas dos n�s, velocidade e press�o
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
            
            % C�lculo da Jacobiana de mudan�a de vari�veis:
            JAC = GradPhi*XY';
            % JACinv = inversa da Jacobiana.
            JACinv = inv(JAC);
            % JACdet = determinante da Jacobiana.
            JACdet = det(JAC);
            % C�lculo do GradPhixy (derivada em rela��o � x e y das fun��es base)
            GradPhixy = JACinv*GradPhi;
            % C�lculo de u,v (velocides) e pre(press�o)
            UV = velocity*Phi;
            Pre = pressure*Qui;
            % Derivada da velocidade em rela��o a xy
            dUV = velocity*GradPhixy';
            
            % Res�duo elementar da velocidade - Quantidade de movimento 
            for ilnode = 1 : 9
                
                iu = iulDoF(ilnode);
                iv = ivlDoF(ilnode);
                
                ElemRv(iu) = ElemRv(iu) + W*JACdet*(ro*Phi(ilnode)*(UV(1)*dUV(1, 1) + UV(2)*dUV(1, 2)) + ...
                                                    GradPhixy(1, ilnode)*(-Pre + 2*mi*dUV(1, 1)) + ...
                                                    mi*GradPhixy(2, ilnode)*(dUV(1, 2) + dUV(2, 1)));
                            
                ElemRv(iv) = ElemRv(iv) + W*JACdet*(ro*Phi(ilnode)*(UV(1)*dUV(2, 1) + UV(2)*dUV(2, 2)) + ...
                                                    GradPhixy(2, ilnode)*(-Pre + 2*mi*dUV(2, 2)) + ...
                                                    mi*GradPhixy(1, ilnode)*(dUV(1, 2) + dUV(2, 1)));
                                                   
            end  
      
            % Res�duo elementar da press�o - Conserva��o de massa
            for ipress = 1 : 3
                
                ip = iplDoF(ipress);
                
                ElemRv(ip) = ElemRv(ip) + W*JACdet*(dUV(1, 1) + dUV(2, 2))*Qui(ipress);
      
            end
            
        end
        
    end
    
    % Condi��es de contorno 
    
    % Elemento na fronteira West (4)
    if iele <= NEY
        
        % Loop nos pontos de Gauss
        NCSI = -1;  
        
        for jgp = 1 : NGP
            
            NETA = XGP(jgp);
            WJ = WGP(jgp);
            
            [Phi, GradPhi, Qui] = BasisFunc(NCSI, NETA);
            
            % C�lculo da Jacobiana de mudan�a de vari�veis:
            JAC = GradPhi*XY';
            % JACinv = inversa da Jacobiana.
            JACinv = inv(JAC);
            % JACdet = determinante da Jacobiana.
            JACdet = det(JAC);
            % C�lculo do GradPhixy (derivada em rela��o � x e y das fun��es base)
            GradPhixy = JACinv*GradPhi;
            % C�lculo de u,v (velocides) e pre(press�o)
            UV = velocity*Phi;
            Pre = pressure*Qui;
            % Derivada da velocidade em rela��o a xy
            dUV = velocity*GradPhixy';
            
            for isidenode = 1 : 3
                
                ilnode = iaux(isidenode, 4);
                iu = iulDoF(ilnode);
                iv = ivlDoF(ilnode);
                
                ElemRv(iu) = ElemRv(iu) + WJ*(-p_in + 2*mi*dUV(1, 1))*Phi(ilnode)*JAC(2, 2);   
                
                ElemRv(iv) = ElemRv(iv) + WJ*mi*(dUV(1, 2) + dUV(2, 1))*Phi(ilnode)*JAC(2, 2); 
    
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
           
            % C�lculo da Jacobiana de mudan�a de vari�veis:
            JAC = GradPhi*XY';
            % JACinv = inversa da Jacobiana.
            JACinv = inv(JAC);
            % JACdet = determinante da Jacobiana.
            JACdet = det(JAC);
            % C�lculo do GradPhixy (derivada em rela��o � x e y das fun��es base)
            GradPhixy = JACinv*GradPhi;
            % C�lculo de u,v (velocides) e pre(press�o)
            UV = velocity*Phi;
            Pre = pressure*Qui;
            % Derivada da velocidade em rela��o a xy
            dUV = velocity*GradPhixy';
           
            for isidenode = 1 : 3
                    
                ilnode = iaux(isidenode, 2);
                iu = iulDoF(ilnode);
                iv = ivlDoF(ilnode);
               
                ElemRv(iu) = ElemRv(iu) - WJ*(-p_out + 2*mi*dUV(1, 1))*Phi(ilnode)*JAC(2, 2); 
               
                ElemRv(iv) = ElemRv(iv) - WJ*mi*(dUV(1, 2) + dUV(2, 1))*Phi(ilnode)*JAC(2, 2);
                
            end
            
        end
        
    end
    
    % Elemento na fronteira North (3)
    if mod(iele, NEY) == 0
        
        for isidenode = 1 : 3
            
            ilnode = iaux(isidenode, 3);
            iu = iulDoF(ilnode);
            iv = ivlDoF(ilnode);
            
            ElemRv(iu) = velocity(1, ilnode) - 0;
            ElemRv(iv) = velocity(2, ilnode) - 0;
   
        end
        
    end
    
    % Elemento na fronteira South (1)
    if mod(iele, NEY) == 1
        
        for isidenode = 1 : 3
            
            ilnode = iaux(isidenode, 1);
            iv = ivlDoF(ilnode);
            
            ElemRv(iv) = velocity(2, ilnode) - 0.0;
        
        end
        
    end
    
end
 