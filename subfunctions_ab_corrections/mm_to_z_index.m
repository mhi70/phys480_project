% Converts depth in mm as it appears on image to z_index
% Input: 
% mm: depth in mm
% Img: image parameters
% 
% Output: 
% z_ind: z index corresponding to depth 'mm'

function z_ind = mm_to_z_index(mm, Img)
dz = abs(Img.zvec(1) - Img.zvec(2));
z_length = abs(Img.zvec(1) - mm/1000);
z_ind = round(z_length/dz);
end