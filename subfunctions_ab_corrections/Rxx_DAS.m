% function DAS_img = Rxx_DAS(Rxx)
%
% Input: 
% Rxx : 3D matrix  Rxx(emission,detection,time)
% Returns DAS_img of Rxx matrix for plotting

function DAS_img = Rxx_DAS(Rxx)

[Nx, ~, Nz] = size(Rxx);

DAS_img = zeros(Nz, Nx);

for iz = 1:Nz
    % Extract diagonal of covariance matrix slice at depth iz
    DAS_img(iz, :) = abs(diag(Rxx(:, :, iz)))';  
end
end
