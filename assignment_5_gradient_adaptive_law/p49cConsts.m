clc;
clear all;

m = 20;
beta = 0.1;
k = 5;


lambda_1 = 0.1;
lambda_0 = 3;

alfa = 3;
beta_g = 0.1;

     
     
theta_star = [m; beta; k];


a11 = 0;
a12 = 1;
a21 = -k/m;
a22 = -beta/m;    

b1 = 0;
b2 = 1/m;

A = [a11 a12;
     a21 a22];
 
B = [b1;
     b2];
 

P0 = [3 0 0;
      0 3 0;
      0 0 3];
  
x0 = [0;0];


R_0 = 50;            
beta = 5;
ZERO = [0 0 0;
        0 0 0;
        0 0 0];