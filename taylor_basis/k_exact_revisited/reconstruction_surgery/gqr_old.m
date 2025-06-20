function [Q, U, R, S] = gqr_old(A, B, flag)
%gqr_old    Generalised QR factorisation.
%   [Q, U, R, S] = gqr_old(A, B) produces a generalised QR factorisation of
%   A and B such that A' = QRU and B' = QS, where Q and U are unitary,
%   and R and S are upper triangular.
%
%   [___] = gqr_old(S,B,0) produces an economy-size decomposition of R & U.

if ( nargin < 3 ), flag = {}; else, flag = {flag}; end

[Q, S] = qr(B.');
[U, R] = qr(flipud(fliplr(A*Q)), flag{:});  %#ok<FLUDLR>
R = rot90(R.',2);
U = flipud(fliplr(U.'));                    %#ok<FLUDLR>
end