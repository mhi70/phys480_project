% function [Rxx,DAS_img] = Rxx_DAS_focusing(R,Img)
%
%
% Inputs:
%
% Ruu : 3D matrix  Ruu(emission,detection,time). Can be single-emission, or plane-wave (PW) emission
% kgrid : kwave structure
% array : structure defining array properties (created by define_transducer_array_2D.m or something like that)
% Img : structure with image properties. It should have, for this function to work, the fields:
%       xvec & zvec [m] are vectors designating the desired image pixels.
%           zvec (depths) is WITH RESPECT TO the array depth (so should start at 0 or above, no negative numbers)
%       c & c_ref [m/s]
%
%       optional parameters
%           fnum : f number, for apodization
%           directivity : optional flag, ==1 to compensate for directivity of individual elements
%           Dx_max_factor : only required if directivity==1
%           plotflag : default is 0.
%                       Set to ==1 to plot the image and an Rxx slice for each emission and reception (slower but fun to watch)
%                       Set to ==2 to plot the image and an Rxx slice for each emission
%           zplot_ind : default is middle of the grid
%    for PW emission:
%       angles [deg], if PW emission (this is the flag for PW emission, otherwise assumes single-emission)
%       timeOrigin [s] :
%           timeOrigin must include the input pulse length, and the emission beam shape:
%               pulse_duration_s = pulse.signal_length/pulse.Fs;
%               timeOrigin = pulse_duration_s/2;  % giving Rxx_subtract
%           if PW emission,     timeOrigin = timeOrigin + abs(element_pitch/2*sind(angles)/Img.c_ref).'; %[s] The reference is set on the central piezo
%           if single emission, timeOrigin = timeOrigin*ones(array.element.num,1);
%
% Outputs:
%
% Rxx : 3D matrix Rxx(x_in,x_out,z). Complex.
% DAS_img : 2D intensity Bmode image created with delay_and_sum (DAS)

% Created by Laura Cobus, 2020
%
% Version 1.03
% Last updated: Aptil 29, 2025
%
% New for this version:
%
% - some stuff cleaned up - unnecessary inputs removed


function [Rxx, DAS_img] = Rxx_DAS_focusing(Ruu,kgrid,array,Img)

plotflag = 0; % default
if isfield(Img,'plotflag')
    plotflag=Img.plotflag;
end
if ~isfield(Img,'zplot_ind')
    Img.zplot_ind = round(length(Img.zvec)); 
end

Img.zplot_mm = Img.zplot_ind*kgrid.dx*1000;

if ~isfield(Img,'fignum')
    Img.fignum=200;
end
% set up for plotting
if plotflag>0
    nplots = 4;
    figure(101);
    subplot1 = subplot(nplots,2,[1,3]); cla(subplot1)
    subplot2 = subplot(nplots,2,[2,4]); cla(subplot2)
    subplot3 = subplot(nplots,2,[5]); cla(subplot3)
    subplot4 = subplot(nplots,2,[6]); cla(subplot4)
    subplot5 = subplot(nplots,2,[7]); cla(subplot5)
    subplot6 = subplot(nplots,2,[8]); cla(subplot6)
end

log_compression = @(x) 20*log10(x/max(x(:)));

c = Img.c;
c_ref = Img.c_ref;
xvec = Img.xvec/c_ref*c;
zvec = Img.zvec/c_ref*c;

Nz = length(zvec);
Nx = length(xvec);

% depth_plot_idx = find(zvec>=zvec(end)*4/5/c_ref*c,1,'first'); % arbitrary

[zvec,xvec] = meshgrid(zvec,xvec);
% meshgrid creates 2D matrices with [z,x]; in zvec, each row is a copy of zvec. In xvec, each column is a copy of xvec,
% matlab indexes by [row,colum]. So choosing zvec(1,:) gives zvec 1D, and choosing xvec(:,1) gives xvec 1D.

element_num = size(Ruu,2);  % # elements in array

if isfield(Img,'angles')
    timeOrigin = Img.timeOrigin;
    angles = Img.angles;
    PWflag=1;
    disp('PW focusing');
else
    PWflag=0;
    timeOrigin = zeros(element_num,1);
end


if(PWflag)
    emission_num=length(angles);
else
    emission_num=element_num; % just over emission elements
end

trans_pos_x = array.element.lateral; % [m]
trans_pos_z = array.element.depth; % [m]

%%%%%%%%%%%%%%%%%%%%%%%  DIRECTIVITY MASK/S %%%%%%%%%%%%%%%%%%%%%%%%%
% construct mask to incorporate directivity of each element, and/or plane-wave spatial extent
% also limit maximum dx = |xin - xout| with apod_mask_Dx

if isfield(Img,'directivity')
    apod_mask_receive = zeros(array.element.num,Nx,Nz);
    apod_mask_transmit = zeros(array.element.num,Nx,Nz);
    if Img.directivity==1

        % criterion:  z/2/(xin-xout)>=fnum. So receive apodization must obey:  xin-xout >= 2*fnum/z
        delta_x_fnum = 2*abs(zvec.')/2/Img.fnum;
        apod_mask_Dx = repmat(abs(Img.xvec - Img.xvec.'),1,1,Nz)<= repmat(permute(delta_x_fnum/Img.Dx_max_factor,[3,1,2]),Nx,1,1);
        % apod_mask_Dx is a binary matrix with dimensions [Dx,Dx,
        for m_d = 1:array.element.num % loop over detection elements - laziest way to construct apod_mask, to work with the main loop
            apod_mask_receive(m_d,:,:) = abs(trans_pos_x(m_d)-xvec)<= delta_x_fnum;
        end
        if PWflag==1
            for m_e = 1:emission_num
                % if angles(m_e)>0
                %     x_crit = zvec.*tand(angles(m_e))+trans_pos_x(1); % [m]
                %     x_crit = trans_pos_x(round(array.element.num/2)) - x_crit; % rel. to x = 0
                %     % ... so abs(focus point x position) > abs(x_crit)
                %     x_crit_ind = (xvec) < (x_crit); % these ones no good
                % else
                x_crit = zvec.*tand(angles(m_e))+trans_pos_x(end);
                x_crit = trans_pos_x(round(array.element.num/2)) + x_crit; % rel. to x = 0
                % ... so abs(focus point x position) > abs(x_crit)
                x_crit_ind = abs(xvec) < abs(x_crit); % these ones are good
                % end
                % x_crit = trans_pos_x(round(array.element.num/2)) - x_crit; % rel. to x = 0
                % % ... so abs(focus point x position) > abs(x_crit)
                % x_crit_ind = abs(xvec)<abs(x_crit); % these ones no good
                % dimensions x_crit_ind[emission,x,z]
                apod_mask_transmit(m_e,x_crit_ind) = 1;
            end
        else
            apod_mask_transmit = apod_mask_receive;
        end
    end
else
    Img.directivity = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FOCUSSING   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

apod_mask = ones(Nx,Nx,Nz);
Rxx = zeros(Nx,Nx,Nz); % input x ,output x ,depth z
DAS_img = zeros(Nx,Nz);
for m_e = 1:emission_num
    tic
    disp(['Emission: ' num2str(m_e) ' calculating...']);
    if(PWflag)
        transmit_distance =  abs(trans_pos_z(1)-zvec)*cosd(angles(m_e))+xvec*sind(angles(m_e)); % [m]
    else
        transmit_distance = sqrt((trans_pos_z(1)-zvec).^2+(trans_pos_x(m_e)-xvec).^2);  % [m]
    end

    if Img.directivity==1
        at = squeeze(apod_mask_transmit(m_e,:,:)); % [x,z]
        at = permute(repmat(at,1,1,Nx),[1,3,2]);   % [x,x,z]
    end

    for m_d = 1:array.element.num % loop over detection elements

        receive_distance =  sqrt((trans_pos_z(m_d)-zvec.').^2+(trans_pos_x(m_d)-xvec.').^2);  % [m]

        time_delay = (permute(repmat(transmit_distance,1,1,Nx),[1,3,2])+permute(repmat(receive_distance,1,1,Nx),[3,2,1]))/Img.c;  % [s]
        % above: [xin, xout, z]  and  [xout, xin, z]

        if Img.directivity==1
            ar = squeeze(apod_mask_receive(m_d,:,:)).';
            ar = permute(repmat(ar,1,1,Nx),[3,2,1]);
            apod_mask = apod_mask_Dx; % binary  at.*ar.*
        end



        % % showing apodization masks. They should all have dimensions [Nx,Nz,Nz].
        if plotflag>0 && m_d==1 && Img.directivity==1 % can't put it in outside loop as ar is only defined here
            imagesc(subplot4,Img.xvec*1000,Img.xvec*1000,squeeze(apod_mask_Dx(:,:,Img.zplot_ind)));
            xlabel(subplot4,'x_{in}');
            ylabel(subplot4,'x_{out}');
            axis(subplot4,'image')
            colorbar(subplot4)
            title(subplot4,['Dx apodization, z = ' num2str(Img.zplot_mm) ' mm'])

            imagesc(subplot5,Img.xvec*1000,Img.xvec*1000,squeeze(at(:,:,Img.zplot_ind)).')
            axis(subplot5,'image')
            colorbar(subplot5)
            title(subplot5,['Transmit apodization, z = ' num2str(Img.zplot_mm) ' mm'])
            xlabel(subplot5,'x_{in}');
            ylabel(subplot5,'x_{out}');

            imagesc(subplot6,Img.xvec*1000,Img.xvec*1000,squeeze(ar(:,:,Img.zplot_ind)));
            xlabel(subplot6,'x_{in}');
            ylabel(subplot6,'x_{out}');
            axis(subplot6,'image')
            colorbar(subplot6)
            title(subplot6,['Receive apodization, z = ' num2str(Img.zplot_mm) ' mm'])
            drawnow
            pause(0.2)
        end

        time_delay = time_delay + timeOrigin(m_e,:);  % original time delays applied

        signal_m = squeeze(Ruu(m_e,m_d,:));

        tmp = interp1(kgrid.t_array,signal_m,time_delay);  % [s]

        %        Vq = interp1(X,V,Xq) interpolates to find Vq, the values of the
        %                      underlying function V=F(X) at the query points Xq.
        % X and V must be the same length

        tmp(isnan(tmp))=0; % Avoid NaNs
        tmp = tmp.*apod_mask;

        Rxx = Rxx + tmp;

        tmp_1d = zeros(Nx,Nz);
        for xx=1:Nx
            tmp_1d(xx,:) = tmp(xx,xx,:);
        end
        DAS_img = DAS_img + tmp_1d;


        if plotflag==1
            if(PWflag)
                % disp('Rotating...');
                Rxx2 = rot90(Rxx,3);
                % disp('Finished rotating.');
            else
                Rxx2=Rxx;
            end
            image_xx1 = zeros(size(Rxx2,3),size(Rxx2,2));
            for z = 1:size(Rxx2,3)
                R_xx_z = squeeze(Rxx2(:,:,z));
                image_xx1(z,:) = diag(R_xx_z)';
            end
            R_xx_z = squeeze(Rxx2(:,:,Img.zplot_ind));

            % bmode, DAS
            imagesc(subplot1,xvec(:,1)*1000,zvec(1,:)*1000,log_compression(abs(DAS_img).'))
            caxis(subplot1,[-40 0]);
            axis(subplot1,'image');
            title(subplot1,'Bmode, DAS');
            yline(subplot1,Img.zplot_ind*kgrid.dx*1000,'w--','Linewidth',2)

            % bmode, Rxx
            imagesc(subplot2,xvec(:,1)*1000,zvec(1,:)*1000,log_compression(abs(image_xx1)))
            caxis(subplot2,[-40 0]);
            axis(subplot2,'image')
            title(subplot2,'Bmode, Rxx');

            imagesc(subplot3,Img.xvec*1000,Img.xvec*1000,abs(R_xx_z));
            axis(subplot3,'image')
            title(subplot3,['Rxx, z = ' num2str(Img.zplot_mm) ' mm']);
            xlabel(subplot3,'x_{in}');
            ylabel(subplot3,'x_{out}');

            sgtitle(['Emission ' num2str(m_e) ', detection ' num2str(m_d)]);

            drawnow;
            pause(0.05);

            clear Rxx2
        end

    end

    if plotflag==2
        if(PWflag)
            Rxx2 = rot90(Rxx,3);
        else
            Rxx2=Rxx;
        end
        image_xx1 = zeros(size(Rxx2,3),size(Rxx2,2));
        for z = 1:size(Rxx2,3)
            R_xx_z = squeeze(Rxx2(:,:,z));
            image_xx1(z,:) = diag(R_xx_z)';
        end
        R_xx_z = squeeze(Rxx2(:,:,Img.zplot_ind));

        % bmode, DAS
        imagesc(subplot1,Img.xvec*1000,Img.zvec*1000,log_compression(abs(DAS_img).'))
        caxis(subplot1,[-40 0]);
        axis(subplot1,'image')
        title(subplot1,'Bmode, DAS');
        yline(subplot1,Img.zplot_ind*kgrid.dx*1000,'w--','Linewidth',2)

        % bmode, Rxx
        imagesc(subplot2,Img.xvec*1000,Img.zvec*1000,log_compression(abs(image_xx1)))
        caxis(subplot2,[-40 0]);
        axis(subplot2,'image')
        title(subplot2,'Bmode, Rxx');

        imagesc(subplot3,Img.xvec*1000,Img.xvec*1000,abs(R_xx_z));
        axis(subplot3,'image')
        title(subplot3,['Rxx, z = ' num2str(Img.zplot_mm) ' mm']);
        xlabel(subplot3,'x_{in}');
        ylabel(subplot3,'x_{out}');

        sgtitle(['Emission ' num2str(m_e) ' of ' num2str(emission_num)]);

        drawnow;
        pause(0.01);
        clear Rxx2

    end

end