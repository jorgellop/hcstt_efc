function [ OUT ] = hexSegMirror_getSupport( hexMirror_struct )
%hexSegMirror_getSupport Returns the support of the pupil function defined
%by a hexagonally segmented mirror 
%   Input: hexMirror_struct - Structure with the following variables 
%   apDia - flat to flat aperture diameter (samples)
%   gapWidth - width of the gap between segments (samples)
%   numRings - number of rings in the segmented mirror (samples)
%   N - size of NxN computational grid 

apDia = hexMirror_struct.apDia; % flat to flat aperture diameter (samples)
gapWidth = hexMirror_struct.gapWidth; % samples
numRings = hexMirror_struct.numRings;% Number of rings in hexagonally segmented mirror 
N = hexMirror_struct.N;

OUT = zeros(N);

hexFlatDiam = (apDia-numRings*2*gapWidth)/(2*numRings+1);
hexSep = hexFlatDiam + gapWidth;

for ringNum = 0:numRings

    cenrow = ringNum*hexSep;
    cencol = 0;
    
    [ OUT ] = hexSegMirror_addHexagon( cenrow,cencol, hexFlatDiam, OUT );
    
    for face = 1:6
        
        step_dir = pi/6*(2*face+5);
        steprow = hexSep*sin(step_dir);
        stepcol = hexSep*cos(step_dir);
        
        stepnum = 1;

        while(stepnum<=ringNum)
            cenrow = cenrow + steprow;
            cencol = cencol + stepcol;
            if(face==6 && stepnum==ringNum)
                %disp(['Finished ring ',num2str(ringNum)]);
            else
                [ OUT ] = hexSegMirror_addHexagon( cenrow,cencol, hexFlatDiam, OUT );
            end
            stepnum = stepnum + 1;
            
        end
    end
end


end

