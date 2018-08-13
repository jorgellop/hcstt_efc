%{
DM Writing Function: Write Sinusoid
- Places a Sin Function on Mirror Surface
- Sin calculated using DE_DMMapSin.m function
- Utilizes DE_DMArrayToVect to shape map for write
*** ASSUMES MIRROR CONNECTION ALREADY PRESENT; DOES NOT CLOSE CONNECTION
*** Defaults to DE_DMMapSin output setting 7: only returns heigh
        If desired, hnm can be plotted with imagesc. Refer to DE_DMMapSin
 
******************************************************
- Arguments:
    h0          = Max poke height in nm
    q           = angle of sinusoid
    x0          = actuators per cycle
    alp         = phase delay
    drv_info    = DM info from OPEN_mutliDM
- Returns:
    hnm         = Surface map as matrix in nm
    hV          = Vector of voltage percentages for writing
******************************************************

Compiled By:    Daniel Echeverri
Last Modified:  08/04/2016
%}

function im = hcstt_TakeCamImage(setCExp,openhimg,CExp_new)

global cam img Xcr Ycr CExp s drv_inf flat itr Idat himg

if(setCExp)
    CExp = CExp_new;
    cam.Timing.Exposure.Set(CExp);
end
%Acquire image
cam.Acquisition.Freeze(true);
% %Extract image from RAM to array in Matlab
[ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
% %Reshape array into plottable matrix
img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
% %Crop data to (200x200) around user input center; Save image
Idat(:,:,1)  = img.Data(int16(Ycr-199):int16(Ycr+200),int16(Xcr-199):int16(Xcr+200));
im = Idat(:,:,1);
% %Draw cropped flat image
if(openhimg)
    himg = imshow(img.Data, 'Border', 'tight');
    title('Camera Image');
    set(himg, 'CData', Idat(:,:,1));
    axis tight equal off;
end

end
