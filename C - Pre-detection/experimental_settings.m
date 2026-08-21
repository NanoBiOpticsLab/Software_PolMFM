function [microscope_settings, PSF_theo, Detect]=experimental_settings()

dlgTitle = 'Experimental parameters';

prompt(1) = {'Pixel-size xy [nm]'};
prompt(2) = {'Pixel-size z [nm]'};
prompt(3) = {'Refractive index'};
prompt(4) = {'Numeric aperture NA'};
prompt(5) = {'Emission wavelength'};
prompt(6) = {'Excitation wavelength'};
prompt(7) = {'Microscope'};
prompt(8) = {'Camera CCD sensitivity [e-/count]'};
prompt(9) = {'Camera EM gain'};
prompt(10) = {'Photons converstion ? (yes=1;no=0)'};

defaultValue{1} = num2str(120);
defaultValue{2} = num2str(378);
defaultValue{3} = num2str(1.518);
defaultValue{4} = num2str(1.4);
defaultValue{5} = num2str(617);
defaultValue{6} = num2str(560);
defaultValue{7} = 'widefield';
defaultValue{8} = num2str(15.4);
defaultValue{9} = num2str(200);
defaultValue{10} = num2str(0);


userValue = inputdlg(prompt,dlgTitle,1,defaultValue);

if( ~ isempty(userValue))
    microscope_settings.pixel_size.xy = str2double(userValue{1});
    microscope_settings.pixel_size.z  = str2double(userValue{2});   
    microscope_settings.RI            = str2double(userValue{3});   
    microscope_settings.NA            = str2double(userValue{4});
    microscope_settings.Em            = str2double(userValue{5});   
    microscope_settings.Ex            = str2double(userValue{6});
    microscope_settings.type          = userValue{7};   
    microscope_settings.Cam_sensitivity= str2double(userValue{8}); 
    microscope_settings.Cam_EM        = str2double(userValue{9}); 
    microscope_settings.Photons_convert=str2double(userValue{10});

   
end



%- Calculate theoretical PSF and show it 
[PSF_theo.xy_nm, PSF_theo.z_nm] = sigma_PSF_BoZhang_v1(microscope_settings);
PSF_theo.xy_pix = PSF_theo.xy_nm / microscope_settings.pixel_size.xy ;
PSF_theo.z_pix  = PSF_theo.z_nm  / microscope_settings.pixel_size.z ;

%- Calculate size of detection region and show it
Detect.region.xy = round(2*PSF_theo.xy_pix)+1;       % Size of detection zone in xy 
Detect.region.z  = round(2*PSF_theo.z_pix)+1;        % Size of detection zone in z 

end