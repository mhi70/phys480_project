x_emitter = 1:5;
z_emitter = 0;
xin = 1:5;
z = 1:10;


diff_x = x_emitter - xin.';
[x, y] = size(diff_x);
receive_x = ones(x, y, length(z)).*diff_x;
receive_z = meshgrid(z-z_emitter, xin,xin);
receive_d = sqrt(receive_d.^2 + (receive_z.').^2);




% Xin = meshgrid((x_emitter - xin.').^2,z);
% d_transmit = sqrt(Xin+((z - z_emitter).^2).');
% 
% Xout = meshgrid((xin - x_emitter.').^2,xin,z);
% Zout = meshgrid(z, xin);
% d_receive =sqrt(Xout + ((Zout - z_emitter).^2));
