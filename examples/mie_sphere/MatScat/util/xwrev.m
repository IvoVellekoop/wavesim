function [ rv ] = xwrev( v )
%XWREV Own implementation of wrev function.
%
%   The wrev function is only available in the Wavelet Toolbox. This
%   function provides the same functionality.
%   
%   RV = XWREV(V) returns the reversed version RV of the vector V. 
%

rv = v(end:-1:1);

end

