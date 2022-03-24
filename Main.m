%=========================================================================
% Código Baseado no MEF para Solução de Escoamentos Incompressíveis
% Data: 11-03-2018
% Autores: Jaime Castañeda Barbosa
%=========================================================================
close all, clear all, clc

%=========================================================================
% PREPROCESSAMENTO

% Propriedades
ro  = 10; 
mi  = 1;

% Condições do escoamento
p_in = 100; % Pressão de entrada
p_out = 0; % Pressão de saída

% Geometria
L  = 1; 
H  = 0.1;

% Malha
NEX = 10;
NEY = 10;

% Cálculo do número de elementos, número de nós e relação entra a numeração local e global dos nós
[NELE, NNODES, DomNodeID] = NodeIndex(NEX, NEY);

% Cálculo de coordenadas dos nós da malha - problema 2D
[X, Y] = Mesh(L, H, NEX, NEY, DomNodeID);

% Cálculo do número de graus de liberdade e relação da numeracao local e global dos graus de liberdade
[NDoF, ASSMtrx] = DoFdrive(NELE, NNODES, DomNodeID);

%=========================================================================
% RESOLUÇÃO DO PROBLEMA

% Solução pelo método de Newton
C0 = ones(NDoF, 1); % Chute inicial
[C] = Newton(NELE, NDoF, NEX, NEY, DomNodeID, ASSMtrx, X, Y, ro, mi, p_in, p_out, C0);

% Cálculo da vazão de entrada Qe e da vazão de saida Qs 
[Qe, Qs] = Vazao(L, H, X, Y, NEX, NEY, NELE, DomNodeID, ASSMtrx, C);

%=========================================================================
% POSPRECESSAMENTO

Graph(L, H, NELE, DomNodeID, ASSMtrx, X, Y, C);

%=========================================================================
