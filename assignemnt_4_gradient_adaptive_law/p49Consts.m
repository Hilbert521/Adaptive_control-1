clc;
clear all;

m = 20;
beta = 0.1;
k = 5;


lambda_1 = 4;
lambda_0 = 4;

alfa = 1;
beta_g = 0.1;

Gamma = [10 0 0;
         0 10 0;
         0 0 10];
     
     
theta_star = [m; beta; k];


a11 = 0;
a12 = 1;
a21 = -k/m;
a22 = -beta/m;    %%%%

b1 = 0;
b2 = 1/m;

A = [a11 a12;
     a21 a22];
 
B = [b1;
     b2];
 

 v = 2;
 w = [1; 2; 3];