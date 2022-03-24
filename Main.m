%=========================================================================
% C�digo Baseado no MEF para Solu��o de Escoamentos Incompress�veis
% Data: 11-03-2018
% Autores: Jaime Casta�eda Barbosa
%=========================================================================
close all, clear all, clc

%=========================================================================
% PREPROCESSAMENTO

% Propriedades
ro  = 10; 
mi  = 1;

% Condi��es do escoamento
p_in = 100; % Press�o de entrada
p_out = 0; % Press�o de sa�da

% Geometria
L  = 1; 
H  = 0.1;

% Malha
NEX = 10;
NEY = 10;

% C�lculo do n�mero de elementos, n�mero de n�s e rela��o entra a numera��o local e global dos n�s
[NELE, NNODES, DomNodeID] = NodeIndex(NEX, NEY);

% C�lculo de coordenadas dos n�s da malha - problema 2D
[X, Y] = Mesh(L, H, NEX, NEY, DomNodeID);

% C�lculo do n�mero de graus de liberdade e rela��o da numeracao local e global dos graus de liberdade
[NDoF, ASSMtrx] = DoFdrive(NELE, NNODES, DomNodeID);

%=========================================================================
% RESOLU��O DO PROBLEMA

% Solu��o pelo m�todo de Newton
C0 = ones(NDoF, 1); % Chute inicial
[C] = Newton(NELE, NDoF, NEX, NEY, DomNodeID, ASSMtrx, X, Y, ro, mi, p_in, p_out, C0);

% C�lculo da vaz�o de entrada Qe e da vaz�o de saida Qs 
[Qe, Qs] = Vazao(L, H, X, Y, NEX, NEY, NELE, DomNodeID, ASSMtrx, C);

%=========================================================================
% POSPRECESSAMENTO

Graph(L, H, NELE, DomNodeID, ASSMtrx, X, Y, C);

%=========================================================================
