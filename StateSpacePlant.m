classdef StateSpacePlant
properties
    A
    B
    C
    D
end
methods
    function self=StateSpacePlant(A,B,C,D)
    self.A=A;
    self.B=B;
    self.C=C;
    self.D=D;
    end
    function [y,x_new]=update(self,x)
    %Generate Noise
    w=normrnd(0,1);
    v=normrnd(0,1);
    %Update state and output
    x_new=self.A*x+self.B*w;
    y=self.C*x_new+self.D*v;

    end
end
end