%{
Speckling Nulling Optimization: 
- Places speckle in given location and measures inverse supression ratio
* Meant to be used with optimization code
        - Optimization code calls this function and attempts to minimize
        output

******************************************************
- Arguments:
    NONE        User Input on Image       

- Returns:
    NONE        Outputs a .FITS cube with properties:
                    WidthxHeightxImage
                    200  x 200  x (heights -> phase)
- Dependencies:
    OPEN_multiDM()  Create DM connection
    UPDATE_multiDM()Place new DM Map
    CLOSE_multiDM() Disconnect DM
    uc480DotNet.dll Library of CCD Camera Communications
    DE_DMW_Sin      Create Sinusoidal Map for DM
******************************************************


Compiled By:    Daniel Echeverri and Jerry Xuan
Last Modified:  08/18/2016
%}

function [Power] = DE_supMinimizer_2(par)

    global icount flat REtxt apd Pdat Pflat Idat ParT sb newp cam drv_inf img Xcr Ycr Samp Rawnm currentSecond startingSecond
    
    icount = icount + 1;
    
    %Store trial parameters in vector for later access
    ParT(icount,:) = par;
    
    %Rename incoming parameters for readability
    h0  = par(1);       %Amplitude
    q   = par(2);       %Angle of Sin
    x   = par(3);       %Actuators/cycle
    alp = par(4);       %Phase delay
    
    %Place Sinusoid on DM and wait to settle
    DM_Command = DE_DMMapSin(h0, q, x, alp);
    FlatCommand = NK_MultiDM_Command(DM_Command, flat);
    JR_UPDATE_MultiDM(drv_inf, FlatCommand);
    pause(.05);
    
    %______Take Measurment_______
    %Clear PM data store of any old data
    % newp.Write(4, 'PM:DS:CLear');
    
    %Acquire image
    % cam.Acquisition.Freeze(true);
    
    %Measure Time
    currentTime = clock;
    currentSecond = currentTime(6) + 60*currentTime(5) + 3600*currentTime(4);
    
    %Take Data; Enable state returns to 0 automaticaly when 50 sample done
    Pdat(icount) = APDMeasurement(apd);
    
    %_______Process and Save Image Data_______
%     %Extract image from RAM to array in Matlab
%     [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
%     %Reshape array into plottable matrix
%     img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
%     %Crop data to (200x200) around user input center
%     Idat(:,:,icount)  = img.Data(int16(Ycr-99):int16(Ycr+100),int16(Xcr-99):int16(Xcr+100));

%     REtxt = fopen(Rawnm, 'wt');
    fprintf(REtxt, 'Seconds: %09.3f     ', currentSecond - startingSecond);
    fprintf(REtxt, 'Parameters: Amp. = %09.5f, Angle = %09.5f, Act. = %09.5f, Phase = %09.5f \n', ...
        ParT(icount,1), ParT(icount, 2), ParT(icount,3), ParT(icount,4));
    fprintf(REtxt, '       Power: %05.2e      Suppression Ratio: %05.2f \n \n', Pdat(icount), Pflat/Pdat(icount));
%     fclose(REtxt);
    
    Power = Pdat(icount)
end