function [qar,dqards,dqardt] = vtqarshape(s,t)
%VTQARSHAPE  quartic shape functions for vector s, t
%   [qar,dqards,dqardt] = vtqarshape(s,t);
%   input
%          s         first triangle coordinate   
%          t         second triangle coordinate  
%  
%   output
%          qar        shape function
%          dqards     s derivative of qar
%          dqardt     t derivative of qar
%
%    TIFISS function: QL; 17 April 2011.
% Copyright (c) 2011 D.J. Silvester and Q. Liao
      n=length(s);
      one = 1.0e0*ones(n,1); zero=0.0e0*one;
      xi(:,1) = one-s-t;
      xi(:,2) = s;
	  xi(:,3) = t;
      dxids(:,1) = -one;
      dxids(:,2) = one;
      dxids(:,3) = zero;
      dxidt(:,1) = -one;
      dxidt(:,2) = zero;
      dxidt(:,3) = one;
%
       qar(:,1) = 32/3*xi(:,1).*(xi(:,1)-1/4).*(xi(:,1)-1/2).*(xi(:,1)-3/4);
       qar(:,2) = 32/3*xi(:,2).*(xi(:,2)-1/4).*(xi(:,2)-1/2).*(xi(:,2)-3/4);
       qar(:,3) = 32/3*xi(:,3).*(xi(:,3)-1/4).*(xi(:,3)-1/2).*(xi(:,3)-3/4);
       qar(:,4) = 64*xi(:,2).*xi(:,3).*(xi(:,2)-1/4).*(xi(:,3)-1/4);
       qar(:,5) = 64*xi(:,1).*xi(:,3).*(xi(:,1)-1/4).*(xi(:,3)-1/4);
       qar(:,6) = 64*xi(:,1).*xi(:,2).*(xi(:,1)-1/4).*(xi(:,2)-1/4);
       qar(:,7) = 128/3*xi(:,2).*xi(:,3).*(xi(:,2)-1/4).*(xi(:,2)-1/2);
       qar(:,8) = 128/3*xi(:,2).*xi(:,3).*(xi(:,3)-1/4).*(xi(:,3)-1/2);
       qar(:,9) = 128/3*xi(:,1).*xi(:,3).*(xi(:,3)-1/4).*(xi(:,3)-1/2);
       qar(:,10)= 128/3*xi(:,1).*xi(:,3).*(xi(:,1)-1/4).*(xi(:,1)-1/2);
       qar(:,11)= 128/3*xi(:,1).*xi(:,2).*(xi(:,1)-1/4).*(xi(:,1)-1/2);
       qar(:,12)= 128/3*xi(:,1).*xi(:,2).*(xi(:,2)-1/4).*(xi(:,2)-1/2);
       qar(:,13)= 128*xi(:,1).*xi(:,2).*xi(:,3).*(xi(:,1)-1/4);
       qar(:,14)= 128*xi(:,1).*xi(:,2).*xi(:,3).*(xi(:,2)-1/4);
       qar(:,15)= 128*xi(:,1).*xi(:,2).*xi(:,3).*(xi(:,3)-1/4); 
       
   %dqards    
       dqards(:,1) = (128/3*xi(:,1).^3-48*xi(:,1).^2+44/3*xi(:,1)-1).*dxids(:,1);
       dqards(:,2) = (128/3*xi(:,2).^3-48*xi(:,2).^2+44/3*xi(:,2)-1).*dxids(:,2);
       dqards(:,3) = (128/3*xi(:,3).^3-48*xi(:,3).^2+44/3*xi(:,3)-1).*dxids(:,3);
       dqards(:,4) = (128*xi(:,2).*xi(:,3).^2-32*xi(:,2).*xi(:,3)-16*xi(:,3).^2+4*xi(:,3)).*dxids(:,2)...
                    +(128*xi(:,2).^2.*xi(:,3)-16*xi(:,2).^2-32*xi(:,2).*xi(:,3)+4*xi(:,2)).*dxids(:,3);
       dqards(:,5) = (128*xi(:,1).*xi(:,3).^2-32*xi(:,1).*xi(:,3)-16*xi(:,3).^2+4*xi(:,3)).*dxids(:,1)...
                    +(128*xi(:,1).^2.*xi(:,3)-16*xi(:,1).^2-32*xi(:,1).*xi(:,3)+4*xi(:,1)).*dxids(:,3);
       dqards(:,6) = (128*xi(:,1).*xi(:,2).^2-32*xi(:,1).*xi(:,2)-16*xi(:,2).^2+4*xi(:,2)).*dxids(:,1)...
                    +(128*xi(:,1).^2.*xi(:,2)-16*xi(:,1).^2-32*xi(:,1).*xi(:,2)+4*xi(:,1)).*dxids(:,2);
       dqards(:,7) = (128*xi(:,2).^2.*xi(:,3)-64*xi(:,2).*xi(:,3)+16/3.*xi(:,3)).*dxids(:,2)...
                    +(16/3*xi(:,2).*(4*xi(:,2)-1).*(2*xi(:,2)-1)).*dxids(:,3);
       dqards(:,8) = (16/3*xi(:,3).*(4*xi(:,3)-1).*(2*xi(:,3)-1)).*dxids(:,2)...
                    +(128*xi(:,2).*xi(:,3).^2-64*xi(:,2).*xi(:,3)+16/3*xi(:,2)).*dxids(:,3);
       dqards(:,9) = (16/3*xi(:,3).*(4*xi(:,3)-1).*(2*xi(:,3)-1)).*dxids(:,1)...
                    +(128*xi(:,1).*xi(:,3).^2-64*xi(:,1).*xi(:,3)+16/3*xi(:,1)).*dxids(:,3);
       dqards(:,10)= (128*xi(:,1).^2.*xi(:,3)-64*xi(:,1).*xi(:,3)+16/3*xi(:,3)).*dxids(:,1)...
                    +(16/3*xi(:,1).*(4*xi(:,1)-1).*(2*xi(:,1)-1)).*dxids(:,3);
       dqards(:,11)= (128*xi(:,1).^2.*xi(:,2)-64*xi(:,1).*xi(:,2)+16/3*xi(:,2)).*dxids(:,1)...
                    +(16/3*xi(:,1).*(4*xi(:,1)-1).*(2*xi(:,1)-1)).*dxids(:,2);
       dqards(:,12)= (16/3*xi(:,2).*(4*xi(:,2)-1).*(2*xi(:,2)-1)).*dxids(:,1)...
                    +(128*xi(:,1).*xi(:,2).^2-64*xi(:,1).*xi(:,2)+16/3*xi(:,1)).*dxids(:,2);
       dqards(:,13)= 256*xi(:,1).*xi(:,2).*xi(:,3).*dxids(:,1) + 128*xi(:,1).^2.*xi(:,3).*dxids(:,2) + 128*xi(:,1).^2.*xi(:,2).*dxids(:,3)...
                    -32*((xi(:,2).*xi(:,3).*dxids(:,1) + xi(:,3).*xi(:,1).*dxids(:,2) + xi(:,1).*xi(:,2).*dxids(:,3))); 
       dqards(:,14)= 128*xi(:,2).^2.*xi(:,3).*dxids(:,1) + 256*xi(:,1).*xi(:,2).*xi(:,3).*dxids(:,2) + 128*xi(:,1).*xi(:,2).^2.*dxids(:,3)...
                    -32*((xi(:,2).*xi(:,3).*dxids(:,1) + xi(:,3).*xi(:,1).*dxids(:,2) + xi(:,1).*xi(:,2).*dxids(:,3))); 
       dqards(:,15)= 128*xi(:,2).*xi(:,3).^2.*dxids(:,1) + 128*xi(:,1).*xi(:,3).^2.*dxids(:,2) + 256*xi(:,1).*xi(:,2).*xi(:,3).*dxids(:,3)...
                    -32*((xi(:,2).*xi(:,3).*dxids(:,1) + xi(:,3).*xi(:,1).*dxids(:,2) + xi(:,1).*xi(:,2).*dxids(:,3))); 
  %dqardt 
       dqardt(:,1) = (128/3*xi(:,1).^3-48*xi(:,1).^2+44/3*xi(:,1)-1).*dxidt(:,1);
       dqardt(:,2) = (128/3*xi(:,2).^3-48*xi(:,2).^2+44/3*xi(:,2)-1).*dxidt(:,2);
       dqardt(:,3) = (128/3*xi(:,3).^3-48*xi(:,3).^2+44/3*xi(:,3)-1).*dxidt(:,3);
       dqardt(:,4) = (128*xi(:,2).*xi(:,3).^2-32*xi(:,2).*xi(:,3)-16*xi(:,3).^2+4*xi(:,3)).*dxidt(:,2)...
                    +(128*xi(:,2).^2.*xi(:,3)-16*xi(:,2).^2-32*xi(:,2).*xi(:,3)+4*xi(:,2)).*dxidt(:,3);
       dqardt(:,5) = (128*xi(:,1).*xi(:,3).^2-32*xi(:,1).*xi(:,3)-16*xi(:,3).^2+4*xi(:,3)).*dxidt(:,1)...
                    +(128*xi(:,1).^2.*xi(:,3)-16*xi(:,1).^2-32*xi(:,1).*xi(:,3)+4*xi(:,1)).*dxidt(:,3);
       dqardt(:,6) = (128*xi(:,1).*xi(:,2).^2-32*xi(:,1).*xi(:,2)-16*xi(:,2).^2+4*xi(:,2)).*dxidt(:,1)...
                    +(128*xi(:,1).^2.*xi(:,2)-16*xi(:,1).^2-32*xi(:,1).*xi(:,2)+4*xi(:,1)).*dxidt(:,2);
       dqardt(:,7) = (128*xi(:,2).^2.*xi(:,3)-64*xi(:,2).*xi(:,3)+16/3.*xi(:,3)).*dxidt(:,2)...
                    +(16/3*xi(:,2).*(4*xi(:,2)-1).*(2*xi(:,2)-1)).*dxidt(:,3);
       dqardt(:,8) = (16/3*xi(:,3).*(4*xi(:,3)-1).*(2*xi(:,3)-1)).*dxidt(:,2)...
                    +(128*xi(:,2).*xi(:,3).^2-64*xi(:,2).*xi(:,3)+16/3*xi(:,2)).*dxidt(:,3);
       dqardt(:,9) = (16/3*xi(:,3).*(4*xi(:,3)-1).*(2*xi(:,3)-1)).*dxidt(:,1)...
                    +(128*xi(:,1).*xi(:,3).^2-64*xi(:,1).*xi(:,3)+16/3*xi(:,1)).*dxidt(:,3);
       dqardt(:,10)= (128*xi(:,1).^2.*xi(:,3)-64*xi(:,1).*xi(:,3)+16/3*xi(:,3)).*dxidt(:,1)...
                    +(16/3*xi(:,1).*(4*xi(:,1)-1).*(2*xi(:,1)-1)).*dxidt(:,3);
       dqardt(:,11)= (128*xi(:,1).^2.*xi(:,2)-64*xi(:,1).*xi(:,2)+16/3*xi(:,2)).*dxidt(:,1)...
                    +(16/3*xi(:,1).*(4*xi(:,1)-1).*(2*xi(:,1)-1)).*dxidt(:,2);
       dqardt(:,12)= (16/3*xi(:,2).*(4*xi(:,2)-1).*(2*xi(:,2)-1)).*dxidt(:,1)...
                    +(128*xi(:,1).*xi(:,2).^2-64*xi(:,1).*xi(:,2)+16/3*xi(:,1)).*dxidt(:,2);
       dqardt(:,13)= 256*xi(:,1).*xi(:,2).*xi(:,3).*dxidt(:,1) + 128*xi(:,1).^2.*xi(:,3).*dxidt(:,2) + 128*xi(:,1).^2.*xi(:,2).*dxidt(:,3)...
                    -32*((xi(:,2).*xi(:,3).*dxidt(:,1) + xi(:,3).*xi(:,1).*dxidt(:,2) + xi(:,1).*xi(:,2).*dxidt(:,3))); 
       dqardt(:,14)= 128*xi(:,2).^2.*xi(:,3).*dxidt(:,1) + 256*xi(:,1).*xi(:,2).*xi(:,3).*dxidt(:,2) + 128*xi(:,1).*xi(:,2).^2.*dxidt(:,3)...
                    -32*((xi(:,2).*xi(:,3).*dxidt(:,1) + xi(:,3).*xi(:,1).*dxidt(:,2) + xi(:,1).*xi(:,2).*dxidt(:,3))); 
       dqardt(:,15)= 128*xi(:,2).*xi(:,3).^2.*dxidt(:,1) + 128*xi(:,1).*xi(:,3).^2.*dxidt(:,2) + 256*xi(:,1).*xi(:,2).*xi(:,3).*dxidt(:,3)...
                    -32*((xi(:,2).*xi(:,3).*dxidt(:,1) + xi(:,3).*xi(:,1).*dxidt(:,2) + xi(:,1).*xi(:,2).*dxidt(:,3))); 
  
                 
return