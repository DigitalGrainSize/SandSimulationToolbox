
function simgrains(config_file)

% 
% Written by Daniel Buscombe, various times in 2011 - 2013
% while at
% School of Marine Science and Engineering, University of Plymouth, UK
% then
% Grand Canyon Monitoring and Research Center, U.G. Geological Survey, Flagstaff, AZ 
% please contact:
% dbuscombe@usgs.gov
% for lastest code version please visit:
% https://github.com/dbuscombe-usgs
% see also (project blog):
% http://dbuscombe-usgs.github.com/
% Buscombe, D. and Rubin, D.M., 2012, Advances in the Simulation and Automated Measurement of Well-Sorted Granular Material, Part 1: Simulations. 
% Journal of Geophysical Research - Earth Surface 117, F02001.
%====================================
%   This function is part of 'sand simulation toolbox' software
%   This software is in the public domain because it contains materials that originally came 
%   from the United States Geological Survey, an agency of the United States Department of Interior. 
%   For more information, see the official USGS copyright policy at 
%   http://www.usgs.gov/visual-id/credit_usgs.html#copyright
%====================================
warning off all

%==================================================================
%============ Part 1: read the config file ==============
%==================================================================
tmp.fid = fopen(config_file);
if tmp.fid<0
    fprintf(2,...
        'error: config_file can''t be read. Is it on your MATLAB path?\n');
    return
end
tmp.str=fgets(tmp.fid);
while ischar(tmp.str)
    if min([~(tmp.str(1)=='%') ~(length(tmp.str)<4)])
        eval(['Settings.' tmp.str ';']);
    end
    tmp.str=fgets(tmp.fid);
end
fclose(tmp.fid);

clear tmp ans

%==================================================================
%============ Part 2: make sure necessary tools, folders, and settings are in place ==============
%==================================================================

if ~isfield(Settings,'num_grains') && ~isfield(Settings,'input_coords_file');
    fprintf(2,...
        'error: need to supply num_grains or input coordinate file. Program exiting ... update config file\n')
    return
end

if ~Settings.print_surface && ~Settings.print_3d &&...
        ~Settings.print_3d_slice_composite && ~Settings.save_polytopes &&...
        ~Settings.save_slicestats && ~Settings.save_particle_centers;
    fprintf(2,...
        'error: no save or print options given. Is there any point running the simulation? Program exiting ... update config file\n')
    return
end

if ~isfield(Settings,'input_coords_file')
    if Settings.num_grains>=20000 && Settings.use_model && ~Settings.use_compiled;
        res=input(['Generating ',num2str(Settings.num_grains),...
            ' grains without a compiled program can take a VERY long time.',...
            'Are you sure you want to continue? (y/n) '],'s');
        if ~isempty(strmatch(res,'n')) || ~isempty(strmatch(res,'N')) ||...
                ~isempty(strmatch(res,'no'));
            fprintf(2,'You said no. Program exiting ... update config file\n')
            return
        end
    end
end

if ~isfield(Settings,'input_coords_file')
    if isfield(Settings,'model')
        % if cvt toolbox required, check matlab can find it
        if Settings.model<=3 && ~Settings.use_compiled; % i.e. requires cvt toolbox
            if ~isdir(Settings.path_to_cvt) || ~isfield(Settings,'path_to_cvt');
                fprintf(2,...
                    'error: supplied path to cvt directory does not exist. Program exiting ... update config file\n')
                return
            else
                addpath(Settings.path_to_cvt);
            end
        end
    end
end

% if mpt toolbox required, check matlab can find it
if Settings.save_polytopes || Settings.save_slicestats ||...
        Settings.print_surface ||...
        Settings.print_3d || Settings.print_3d_slice_composite;
    if ~isdir(Settings.path_to_mpt) || ~isfield(Settings,'path_to_mpt');
        fprintf(2,...
            'error: supplied path to mpt directories does not exist. Program exiting ... update config file\n')
        return
    else
        addpath(genpath(Settings.path_to_mpt));
    end
end

if isfield(Settings,'input_coords_file') && ~isempty(Settings.input_coords_file);
    if exist(Settings.input_coords_file,'file')~=2
        fprintf(2,...
            'error: supplied input_coords_file does not exist. Program exiting ... update config file\n')
        return
    end
end

if Settings.save_polytopes || Settings.save_slicestats ||...
        Settings.save_particle_centers;
    if ~isdir(Settings.save_output_direc) || ~isfield(Settings,'save_output_direc')
        mkdir(Settings.save_output_direc)
        %fprintf(2,...
        %    'error: supplied path to saved output directory does not exist. Program exiting ... update config file\n')
        return
    end
end

if Settings.print_surface || Settings.print_3d ||...
        Settings.print_3d_slice_composite;
    if ~isdir(Settings.print_output_direc) || ~isfield(Settings,'print_output_direc');
        fprintf(2,...
            'error: supplied path to output print directory does not exist. Program exiting ... update config file\n')
        return
    end
end

if Settings.print_surface || Settings.print_3d ||...
        Settings.print_3d_slice_composite;
    if ~isfield(Settings,'out_image_format')||...
            isempty(Settings.out_image_format);
        imfmt='.eps';
        imflag='-deps';
        Settings.out_image_format='eps';
    else
        imfmt=['.',Settings.out_image_format];
        imflag=['-d',Settings.out_image_format];
    end
    
    if ~isfield(Settings,'res_dpi')||...
            isempty(Settings.res_dpi);
        dpi='-r1000';
    else
        dpi=['-r',num2str(Settings.res_dpi)];
    end
end

% set defaults
if isempty(Settings.modify_throats) || ~isfield(Settings,'modify_throats');
    Settings.modify_throats=0;
    disp('modify_throats not given. Throats will not be modified')
end

if isempty(Settings.grain_conc) || ~isfield(Settings,'grain_conc');
    Settings.grain_conc=0.7;
    disp('grain_conc not given. Set to 0.7')
end

if ~isfield(Settings,'input_coords_file')
    if isfield(Settings,'model')
        % compiled program only works with models 1, 2 and 3
        if Settings.model>3 && Settings.use_compiled;
            Settings.use_compiled=0;
        end
    end
end

%==================================================================
%============ Part 3: generate or load in points ==============
%==================================================================
if Settings.gen_coords; % generate coordinates?
    disp('generating coordinates ...')
    if Settings.use_model; % use a model?
        disp('using a model to generate coordinates')
        
        if Settings.use_compiled; % use a compiled model (cvt, written in fortran -lots quicker than equivalent matlab)
            disp('using a compiled model to generate coordinates')
            if exist(Settings.compiled_file,'file')==2 &&...
                    Settings.model<=3; % does the file exist?, is it appropriate?
                % run the model:
                [status]=system([Settings.compiled_file,' ',...
                    num2str(Settings.num_grains),' ',...
                    num2str(Settings.model),' ',...
                    num2str(Settings.modeloptions_batch),' ',...
                    num2str(Settings.modeloptions_it_max),...
                    ' ''out.txt''']);
                if status==0
                    disp('Compiled program execution successful')
                    grain_centres=load('out.txt');
                    delete out.txt
                    
                    if Settings.print_surface;
                        surface_printfile=['surface_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_batch',...
                            num2str(Settings.modeloptions_batch),'_itmax',...
                            num2str(Settings.modeloptions_it_max),imfmt];
                    end
                    if Settings.print_3d;
                        threeD_printfile=['threeD_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_batch',...
                            num2str(Settings.modeloptions_batch),'_itmax',...
                            num2str(Settings.modeloptions_it_max),imfmt];
                    end
                    if Settings.print_3d_slice_composite;
                        threeDcomp_printfile=['threeDcomp_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_batch',...
                            num2str(Settings.modeloptions_batch),'_itmax',...
                            num2str(Settings.modeloptions_it_max),imfmt];
                    end
                    
                else
                    fprintf(2,...
                        'error running compiled program. Program exiting ...\n')
                    return
                end
                
            else % doesn't seem to exist
                fprintf(2,...
                    'error: use_compiled=1, but compiled_file does not exist. Program exiting ... update config file\n')
                return
            end % use compiled model? loop
        else % use matlab cvt
            disp('using cvt to generate particle centres')
            if ~isempty(Settings.model) && isnumeric(Settings.model);
                
                if Settings.model==4; % model is CP
                    disp('model is CP')
                    if ~isfield(Settings,'model4_lambda') ||...
                            isempty(Settings.model4_lambda);
                        options.lam=10;
                    else
                        options.lam=Settings.model4_lambda;
                    end
                    if ~isfield(Settings,'model4_var') ||...
                            isempty(Settings.model4_var);
                        options.var=0.1;
                    else
                        options.var=Settings.model4_var;
                    end
                    
                    if Settings.print_surface;
                        surface_printfile=['surface_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_lam',num2str(options.lam),'_var',...
                            num2str(options.var),imfmt];
                    end
                    if Settings.print_3d;
                        threeD_printfile=['threeD_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_lam',num2str(options.lam),'_var',...
                            num2str(options.var),imfmt];
                    end
                    if Settings.print_3d_slice_composite;
                        threeDcomp_printfile=['threeDcomp_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_lam',num2str(options.lam),'_var',...
                            num2str(options.var),imfmt];
                    end
                    
                elseif Settings.model==5; % model is strauss
                    disp('Model is Strauss')
                    if ~isfield(Settings,'model5_c') ||...
                            isempty(Settings.model5_c);
                        options.c=1;
                    else
                        options.c=Settings.model5_c;
                    end
                    if ~isfield(Settings,'model5_delta') ||...
                            isempty(Settings.model5_delta);
                        options.delta=Settings.num_grains/10000;
                    else
                        options.delta=Settings.model5_delta;
                    end
                    
                    if Settings.print_surface;
                        surface_printfile=['surface_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_c',num2str(options.c),'_delta',...
                            num2str(options.delta),imfmt];
                    end
                    if Settings.print_3d;
                        threeD_printfile=['threeD_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_c',num2str(options.c),'_delta',...
                            num2str(options.delta),imfmt];
                    end
                    if Settings.print_3d_slice_composite;
                        threeDcomp_printfile=['threeDcomp_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_c',num2str(options.c),'_delta',...
                            num2str(options.delta),imfmt];
                    end
                    
                else % model is 1, 2 or 3 (pvt, cvt or hal)
                    disp('model is pvt, cvt or halton')
                    if ~isfield(Settings,'modeloptions_batch') ||...
                            isempty(Settings.modeloptions_batch);
                        options.batch = 1000;
                    else
                        options.batch=Settings.modeloptions_batch;
                    end
                    
                    if ~isfield(Settings,'modeloptions_it_max') ||...
                            isempty(Settings.modeloptions_it_max);
                        options.it_max = 20;
                    else
                        options.it_max=Settings.modeloptions_it_max;
                    end
                    
                    if Settings.print_surface;
                        surface_printfile=['surface_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_batch',num2str(options.batch),'_itmax',...
                            num2str(options.it_max),imfmt];
                    end
                    if Settings.print_3d;
                        threeD_printfile=['threeD_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_batch',num2str(options.batch),'_itmax',...
                            num2str(options.it_max),imfmt];
                    end
                    if Settings.print_3d_slice_composite;
                        threeDcomp_printfile=['threeDcomp_',...
                            num2str(Settings.num_grains),'_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),'_batch',num2str(options.batch),'_itmax',...
                            num2str(options.it_max),imfmt];
                    end
                    
                end % set options loop
                % note: no options for model 6
                if Settings.model==6;
                    if Settings.print_surface;
                        surface_printfile=['surface_',num2str(Settings.num_grains),...
                            '_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),imfmt];
                    end
                    if Settings.print_3d;
                        threeD_printfile=['threeD_',num2str(Settings.num_grains),...
                            '_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),imfmt];
                    end
                    if Settings.print_3d_slice_composite;
                        threeDcomp_printfile=['threeDcomp_',num2str(Settings.num_grains),...
                            '_C',num2str(Settings.grain_conc),...
                            '_mod',num2str(Settings.model),imfmt];
                    end
                end
                
                %==========================================================
                %==========================================================
                % run function to generate particle centers
                disp('generating particle centres')
                grain_centres=gen_points(Settings.model,Settings.num_grains,options)';
                
                %==========================================================
                %==========================================================
                
            else % model is empty or not numeric
                fprintf(2,...
                    'error: use_model=1, but model is empty or non-numeric. Program exiting ... update config file\n')
                return
            end % is model empty? loop
            
        end % use compiled model? loop
        
    elseif Settings.gen_coords_f_noise; % gen according to f_noise?
        disp('generating coordinates using 1/f noise')
        if ~isempty(Settings.f_exponent) && isnumeric(Settings.f_exponent)
            disp(['generating particle distribution according to noise on f to exponent: ',...
                num2str(Settings.f_exponent)])
            CDF=noiseonf_3d(100,100,100,Settings.f_exponent);
            grain_centres=discrete_pdf_sample3d(Settings.num_grains,CDF)';
        else % no f exponent
            fprintf(2,...
                'error: gen_coords_f_noise=1, but f_exponent is empty. Program exiting ... update config file\n')
            return
        end  % f exponent? loop
        
        if Settings.print_surface;
            surface_printfile=['surface_',num2str(Settings.num_grains),...
                '_C',num2str(Settings.grain_conc),'_fnoise_exp',...
                num2str(Settings.f_exponent),imfmt];
        end
        if Settings.print_3d;
            threeD_printfile=['threeD_',num2str(Settings.num_grains),...
                '_C',num2str(Settings.grain_conc),'_fnoise_exp',...
                num2str(Settings.f_exponent),imfmt];
        end
        if Settings.print_3d_slice_composite;
            threeDcomp_printfile=['threeDcomp_',num2str(Settings.num_grains),...
                '_C',num2str(Settings.grain_conc),'_fnoise_exp',...
                num2str(Settings.f_exponent),imfmt];
        end
        
    elseif Settings.gen_coords_supplied_image; % generate points according to image?
        if ~isempty(Settings.supplied_image);
            if exist(Settings.supplied_image,'file')==2; % does the file exist?
                
                I=double(imread(Settings.supplied_image)); % read image
                if numel(size(I))==3 % if RGB
                    I=col2gray(I); % convert to greyscale
                end
                [a,b]=size(I);
                I=I(1:min(a,b),1:min(a,b)); % crop to smallest dimension
                
                disp(['generating particle distribution according to input image: ',...
                    Settings.supplied_image])
                % if supplied_image_do_filter=1
                if Settings.supplied_image_do_filter
                    
                    if ~isfield(Settings,'supplied_image_disksize') ||...
                            isempty(Settings.supplied_image_disksize);
                        disksize = ceil(size(I,1)/25);
                    else
                        disksize=Settings.supplied_image_disksize;
                    end
                    if ~isfield(Settings,'supplied_image_gain') ||...
                            isempty(Settings.supplied_image_gain);
                        gain = 10;
                    else
                        gain=Settings.supplied_image_gain;
                    end
                    if ~isfield(Settings,'supplied_image_cutoff') ||...
                            isempty(Settings.supplied_image_cutoff);
                        cutoff = 0.5;
                    else
                        cutoff =Settings.supplied_image_cutoff;
                    end
                    
                    H = fspecial('disk',disksize);
                    CDF = conv2(-I,H,'valid'); % convolve image with filter
                    scale=single(size(I,1)/size(CDF,1)); % may not be same size
                    CDF=imageresize(CDF,scale,scale); % so we rescale
                    G = adjcontrast(rescale(CDF,0,1), gain, cutoff); % adjust contrast
                    CDF=G./sum(sum(G));
                    
                    if Settings.print_surface;
                        surface_printfile=['surface_',num2str(Settings.num_grains),...
                            '_C',num2str(Settings.grain_conc),'_supp_image_gain',...
                            num2str(gain),'_cutoff',num2str(cutoff),...
                            '_disksize',num2str(disksize),imfmt];
                    end
                    if Settings.print_3d;
                        threeD_printfile=['threeD_',num2str(Settings.num_grains),...
                            '_C',num2str(Settings.grain_conc),'_supp_image_gain',...
                            num2str(gain),'_cutoff',num2str(cutoff),...
                            '_disksize',num2str(disksize),imfmt];
                    end
                    if Settings.print_3d_slice_composite;
                        threeDcomp_printfile=['threeDcomp_',num2str(Settings.num_grains),...
                            '_C',num2str(Settings.grain_conc),'_supp_image_gain',...
                            num2str(gain),'_cutoff',num2str(cutoff),...
                            '_disksize',num2str(disksize),imfmt];
                    end
                    
                else
                    CDF=I./sum(sum(I));
                    
                    if Settings.print_surface;
                        surface_printfile=['surface_',num2str(Settings.num_grains),...
                            '_C',num2str(Settings.grain_conc),...
                            '_supp_image',imfmt];
                    end
                    if Settings.print_3d;
                        threeD_printfile=['threeD_',num2str(Settings.num_grains),...
                            '_C',num2str(Settings.grain_conc),...
                            '_supp_image',imfmt];
                    end
                    if Settings.print_3d_slice_composite;
                        threeDcomp_printfile=['threeDcomp_',num2str(Settings.num_grains),...
                            '_C',num2str(Settings.grain_conc),...
                            '_supp_image',imfmt];
                    end
                    
                end
                disp('generating particle centres')
                grain_centres=discrete_pdf_sample3d(Settings.num_grains,CDF)';
                
            else % doesn't seem to exist
                fprintf(2,...
                    'error: gen_coords_supplied_image=1, but supplied_image does not exist. Program exiting ... update config file\n')
                return
            end % supplied image exist? loop
            
        else % no supplied image
            fprintf(2,...
                'error: gen_coords_supplied_image=1, but supplied_image does not exist. Program exiting ... update config file\n')
            return
        end
    end
    
    % don't generate coordinates, i.e. read them in
elseif isfield(Settings,'input_coords_file') &&...
        ~isempty(Settings.input_coords_file);
    Settings.save_particle_centers=0; % set this to zero if not already - no point saving them
    [pathstr,namestr]=fileparts(Settings.input_coords_file);
    addpath(pathstr)
    disp('reading existing particle centres from file')
    try
        grain_centres=load(Settings.input_coords_file);
    catch
        try
            grain_centres=csvread(Settings.input_coords_file);
        catch
            fprintf(2,...
                'error: cannot read input_coords_file using load or csvread. Program exiting ... update coords file\n')
            return
        end
    end
    
    if Settings.print_surface;
        surface_printfile=['surface_C',num2str(Settings.grain_conc),...
            '_supp_coords_',namestr,imfmt];
    end
    if Settings.print_3d;
        threeD_printfile=['threeD_C',num2str(Settings.grain_conc),...
            '_supp_coords_',namestr,imfmt];
    end
    if Settings.print_3d_slice_composite;
        threeDcomp_printfile=['threeDcomp_C',num2str(Settings.grain_conc),...
            '_supp_coords_',namestr,imfmt];
    end
    
else % no input coords file
    fprintf(2,...
        'error: gen_coords=0, but input_coords_file is empty. What gives? . Program exiting ... update config file\n')
    return
    
end % gen_coords? loop



%==================================================================
%============ Part 4: generate sediment structure from points ==============
%==================================================================
if Settings.save_polytopes || Settings.save_slicestats ||...
        Settings.print_surface ||...
        Settings.print_3d || Settings.print_3d_slice_composite;
    
    try
        disp('generating sediment structure from points')
        grain_centres=rescale(grain_centres,0,1);
        Options.pbound=unitbox(3,1);
        
        % Matrix: number of points times dimension
        if size(grain_centres,1)>size(grain_centres,2)
            P=mpt_voronoi(grain_centres,Options);
            [index]=find_inside_domain(P,grain_centres',0,99);
        else
            P=mpt_voronoi(grain_centres',Options);
            [index]=find_inside_domain(P,grain_centres,0,99);
        end
        clc
        disp('reducing dimensionality')
        [P,keptrows] = reduce(P);
        Pfilt=P(index);
        [Q]=P_reduceconc_nobar(Pfilt,Settings.grain_conc);
        
    catch
        fprintf(2,...
            ['error in computing polytopes from particle centers: check mpt settings \n'])
    end
    
end


%==================================================================
%============ Part 5: saving and printing ==============
%==================================================================

% ==== save particle centers ========
if Settings.save_particle_centers;
    if ~isfield(Settings,'save_output_file') ||...
            isempty(Settings.save_output_file);
        outfile='out.txt';
    else
        outfile=Settings.save_output_file;
    end
    [pathstr,namestr]=fileparts(Settings.save_output_direc);
    savedfile=[pathstr,filesep,namestr,filesep,outfile];
    disp('saving')
    savedfile=check_savedfile(savedfile,'txt');
    
    try
        csvwrite(savedfile,grain_centres);
        disp('--------------------------------------')
        disp(['particle centre coordinates written to ',savedfile])
        disp('--------------------------------------')
    catch
        fprintf(2,['error in writing particle centre coordinates file to ',...
            savedfile,'\n'])
    end
end

% ==== save particle polytopes ========
if Settings.save_polytopes;
    if ~isfield(Settings,'save_polytope_file') ||...
            isempty(Settings.save_polytope_file);
        outfile='polytopes.mat';
    else
        outfile=Settings.save_polytope_file;
    end
    [pathstr,namestr]=fileparts(Settings.save_output_direc);
    savedfile=[pathstr,filesep,namestr,filesep,outfile];
    
    savedfile=check_savedfile(savedfile,'mat');
    
    try
        disp('saving mat file')
        save(savedfile,'grain_centres','P','Pfilt','Q');
        disp('--------------------------------------')
        disp(['particle polytopes written to ',savedfile])
        disp('--------------------------------------')
    catch
        fprintf(2,['error in writing particle polytopes file to ',...
            savedfile,'\n'])
    end
    
end

% ==== save slice stats ========
if Settings.save_slicestats;
    if ~isfield(Settings,'save_slice_file') ||...
            isempty(Settings.save_slice_file);
        outfile='slicestats.mat';
    else
        outfile=Settings.save_slice_file;
    end
    [pathstr,namestr]=fileparts(Settings.save_output_direc);
    savedfile=[pathstr,filesep,namestr,filesep,outfile];
    
    savedfile=check_savedfile(savedfile,'mat');
    
    try
        disp('saving slices')
        cutlocs=[0:.01:.2,.25:.05:.8];
        [M,N]=slice_stats(Q);
        save(savedfile,'M','N','cutlocs');
        disp('--------------------------------------')
        disp(['particle slice stats written to ',savedfile])
        disp('--------------------------------------')
    catch
        fprintf(2,['error in writing particle slice stats file to ',...
            savedfile,'\n'])
    end
    
end


% ==== generate and print surface ========
if Settings.print_surface;
    
    [pathstr,namestr]=fileparts(Settings.print_output_direc);
    savedfile=[pathstr,filesep,namestr,filesep,surface_printfile];
    
    savedfile=check_savedfile(savedfile,Settings.out_image_format);
    
    try
        disp('printing 3d surface')
        print_3d_surface(Q)
        print(imflag,dpi,savedfile)
        close
        disp('--------------------------------------')
        disp(['particle surface printed to ',savedfile])
        disp('--------------------------------------')
    catch
        fprintf(2,'error in printing particle surface\n')
        return
    end
end


% ==== generate and print 3d plot========
if Settings.print_3d;
    
    [pathstr,namestr]=fileparts(Settings.print_output_direc);
    savedfile=[pathstr,filesep,namestr,filesep,threeD_printfile];
    
    savedfile=check_savedfile(savedfile,Settings.out_image_format);
    
    try
        disp('printing 3d particles')
        Options.shade=1; Options.linestyle='-'; Options.wire=0;
        plot(Q); view(3); axis tight
        print(imflag,dpi,savedfile)
        close
        disp('--------------------------------------')
        disp(['3d particles printed to ',savedfile])
        disp('--------------------------------------')
    catch
        fprintf(2,'error in printing 3d particles\n')
        return
    end
end

% ==== generate and print 3d slice composite ========
if Settings.print_3d_slice_composite;
    
    [pathstr,namestr]=fileparts(Settings.print_output_direc);
    savedfile=[pathstr,filesep,namestr,filesep,threeDcomp_printfile];
    
    savedfile=check_savedfile(savedfile,Settings.out_image_format);
    
    try
        disp('printing a 3d slice composite')
        do_3d_composite(Q)
        print('-dtiffnocompression', '-r400', savedfile)
        close
        disp('--------------------------------------')
        disp(['3d particles composite plot printed to ',savedfile])
        disp('--------------------------------------')
    catch
        fprintf(2,'error in printing 3d composite plot of particles\n')
        return
    end
end


warning on all


%==========================================================
%==========================================================
%=============   end of main script  ======================
%==========================================================
%==========================================================

function savedfile=check_savedfile(savedfile,filetype)

   % see if the requested file already exists, and if so find an
   % alternative to avoid overwriting
   counter=1;
   while exist(savedfile,'file')==2
       disp([savedfile,' already exists!'])
       savedfile=[savedfile(1:regexp(savedfile,filetype)-2),'_',...
           num2str(counter),['.',filetype]];
       counter=counter+1;
   end


%==================================================================
function print_3d_surface(Q)

    figure
    hold on
    for ii=fliplr([0:.1:1])%fliplr([(1:num).*nanmean(M)] )
        eval(['[Pcut, sliced] = slice(Q,2,',num2str(ii),');'])
        
        for jj=1:length(Pcut)
            eval(['E=extreme(Pcut(',num2str(jj),'));'])
            
            ey=E(:,2); ex=E(:,1);
            K=convhull(ex,ey);
            ex=ex(K); ey=ey(K);
            
            h=patch(ex,ey,[rand(1,3)]); set(h,'LineStyle','-','FaceAlpha',0,...
                'EraseMode', 'background')
            
        end
    end

axis tight off

%==================================================================
function img_out=imageresize(img_in, rowScale, colScale)

    [row col dep]=size(img_in);
    
    if dep<2
        img_out=zeros(floor(row*rowScale),floor(col*colScale));
        for i=1:floor(row*rowScale)
            for j=1:floor(col*colScale)
                img_out(i,j)=img_in(ceil(1/rowScale * i), ceil(1/colScale * j));
            end
        end
        
    else
        img_out=zeros(floor(row*rowScale),floor(col*colScale),dep);
        for i=1:floor(row*rowScale)
            for j=1:floor(col*colScale)
                for k=1:dep
                    img_out(i,j,k)=img_in(ceil(1/rowScale * i), ceil(1/colScale * j),k);
                end
            end
        end
    end

%==================================================================
function []=do_3d_composite(Q)

    im=zeros(501,601,length([0:.01:1]));
    counter=1;
    for i=[0:.01:1]
        Pcut=slice(Q,2,i);
        plot(Pcut)
        axis([0 1 0 1])
        axis off
        print -dtiffnocompression tmp.tiff
        I=double(rgb2gray(imread('tmp.tiff')));
        
        Ib=(rescale(I,0,1)~=1);
        Ib=Ib(200:700,400:1000);
        
        im(:,:,counter)=Ib;
        counter=counter+1;
    end
    
    close all
    delete tmp.tiff
    
    im_50percent=imageresize(im,.5,.5);
    im=im_50percent;
    
    [x,y,z]=meshgrid([1:size(im,2)],[1:size(im,1)],[1:size(im,3)]);
    xmin = min(x(:));
    ymin = min(y(:));
    zmin = min(z(:));
    xmax = max(x(:));
    ymax = max(y(:));
    zmax = max(z(:));
    
    hslice = surf(linspace(xmin,xmax,100),...
        linspace(ymin,ymax,100),...
        zeros(100));
    
    rotate(hslice,[-1,0,0],-90)
    xd = get(hslice,'XData');
    yd = get(hslice,'YData');
    zd = rescale(get(hslice,'ZData'),zmin,zmax);
    
    delete(hslice)
    
    h = slice(x,y,z,im,xd,yd,zd);
    set(h,'FaceColor','flat',...
        'EdgeColor','none',...
        'DiffuseStrength',5)
    
    colormap gray
    hold on
    
    hslice = surf(linspace(xmin,xmax,100),...
        linspace(ymin,ymax,100),...
        zeros(100));
    rotate(hslice,[0,-1,0],-90)
    xd = get(hslice,'XData');
    yd = get(hslice,'YData');
    zd = rescale(get(hslice,'ZData'),zmin,zmax);
    delete(hslice)
    
    h = slice(x,y,z,im,xd,yd,zd);
    set(h,'FaceColor','flat',...
        'EdgeColor','none',...
        'DiffuseStrength',5)
    
    hx = slice(x,y,z,im,xmax,[],[]);
    set(hx,'FaceColor','flat','EdgeColor','none')
    
    hy = slice(x,y,z,im,[],ymax,[]);
    set(hy,'FaceColor','flat','EdgeColor','none')
    
    axis tight
    box on
    camproj perspective
    lightangle(-45,45)
    set(gcf,'Renderer','zbuffer')


%==================================================================
function [M,N]=slice_stats(P)

    M=[]; N=[];
    for ii=[0:.01:.2,.25:.05:.8]
        [Pcut, sliced] = slice(P,2,ii);
        [xc,rc] = chebyball(Pcut); M=[M;mean(rc.*2)]; N=[N;length(rc)];
    end %ii


function [Q]=P_reduceconc_nobar(P,fct)

    E = pelemfun(@extreme, P);
    
    Q=[];
    for i=1:length(E)
        ex = fct*E{i}(:,1)+(1-fct)*mean(E{i}(:,1));  % expand x by factor fct
        ey = fct*E{i}(:,2)+(1-fct)*mean(E{i}(:,2));  % expand y
        ez = fct*E{i}(:,3)+(1-fct)*mean(E{i}(:,3));  % expand z
        Q = [Q,hull([ex,ey,ez])];
    end

%==================================================================
function [index]=find_inside_domain(P,r,min,max)

    index=[];
    for ip=1:length(P)
        E=extreme(P(ip));
        try
            if(any((r(1,ip)+E(:,1))<min) || any((r(2,ip)+E(:,2))<min) ||...
                    any((r(3,ip)+E(:,3))<min))
                continue
            elseif(any((r(1,ip)+E(:,1))>max) || any((r(2,ip)+E(:,2))>max) ||...
                    any((r(3,ip)+E(:,3))>max))
                continue
            else
                index=[index;ip];
            end
        catch
            continue
        end % try
    end % ip

%==================================================================
function out=rescale(in,mn,mx)
    
    m=min(in(:)); M=max(in(:));
    out=(mx-mn).*(in-m)/(M-m)+mn;

%==================================================================
function xyz = discrete_pdf_sample3d ( n , pdf)

    if numel(size(pdf))~=3
        pdf=repmat(pdf,[1 1 10]);
    end
    
    if sum(sum(sum(pdf)))~=1
        pdf=pdf./sum(sum(sum(pdf)));
    end
    
    if ( n <= 0 )
        xyz = [];
        return
    end
    %  "Integrate" the data over rows and columns of the region to get the CDF.
    cdf = set_discrete_cdf ( pdf );
    %  Choose N CDF values at random.
    u = rand ( n, 1 );
    %  Find the cell corresponding to each CDF value,
    %  and choose a random point in that cell.
    xyz = discrete_cdf_to_xyz ( cdf, n, u );


%==================================================================
function cdf = set_discrete_cdf ( pdf )

    % SET_DISCRETE_PDF sets a CDF from a discrete PDF.
    
    cdf=zeros(size(pdf));
    total = 0.0;
    for j = 1 : size(pdf,1)
        for i = 1 : size(pdf,2)
            for k=1:size(pdf,3)
                total = total + pdf(i,j,k);
                cdf(i,j,k) = total;
            end
        end
    end

%==================================================================
function   xyz = discrete_cdf_to_xyz ( cdf, n, u )

    xyz(1:3,1:n) = 0.0;
    low = 0.0;
    for j = 1 : size(cdf,1)
        for i = 1 : size(cdf,2)
            for l=1:size(cdf,3)
                high = cdf(i,j,l);
                for k = 1 : n
                    if ( low <= u(k) && u(k) <= high )
                        xyz(1,k) = ( ( i - 1 ) + rand ( 1, 1 ) ) / size(cdf,1);
                        xyz(2,k) = ( ( j - 1 ) + rand ( 1, 1 ) ) / size(cdf,2);
                        xyz(3,k) = ( ( l - 1 ) + rand ( 1, 1 ) ) / size(cdf,3);
                    end
                end
                low = high;
            end
        end
    end


%==================================================================
function newim = adjcontrast(im, gain, cutoff)

    if isa(im,'uint8');
        newim = double(im);
    else
        newim = im;
    end
    
    % rescale range 0-1
    newim = newim-min(min(newim));
    newim = newim./max(max(newim));
    
    newim =  1./(1 + exp(gain*(cutoff-newim)));  % Apply Sigmoid function
    
 %==================================================================   
 function f = fspecial(varargin)
    
    [filtertype, sze, sigma, radius, len, angle] = checkargs(varargin(:));    
    
    rows = sze(1); cols = sze(2);
    r2 = (rows-1)/2; c2 = (cols-1)/2;

    if strcmpi(filtertype, 'average')
	f = ones(sze)/(rows*cols);
    
    elseif strcmpi(filtertype, 'disk')
	[x,y] = meshgrid(-c2:c2, -r2:r2);
	rad = sqrt(x.^2 + y.^2);	
	f = rad <= radius;
	f = f/sum(f(:));	
	
    elseif strcmpi(filtertype, 'gaussian')
	[x,y] = meshgrid(-c2:c2, -r2:r2);
	radsqrd = x.^2 + y.^2;
	f = exp(-radsqrd/(2*sigma^2));
	f = f/sum(f(:));
	
    elseif strcmpi(filtertype, 'log') 
	[x,y] = meshgrid(-c2:c2, -r2:r2);
	radsqrd = x.^2 + y.^2;	
	f = -1/(pi*sigma^4)*(1-radsqrd/(2*sigma^2))...
	     .*exp(-radsqrd/(2*sigma^2));	   
	f = f-mean(f(:));    % Ensure 0 DC

    elseif strcmpi(filtertype, 'laplacian')	    	
	f = [1  1  1 
	     1 -8  1 
	     1  1  1]; 
    
    elseif strcmpi(filtertype, 'unsharp')	    		
	f = -fspecial('log') + [0 0 0 
		                0 1 0 
		                0 0 0];
	
    elseif strcmpi(filtertype, 'sobel')
	f = [1 2 1; 0 0 0; -1 -2 -1];
    
    elseif strcmpi(filtertype, 'prewitt')
	f = [1 1 1; 0 0 0; -1 -1 -1];	
	
    elseif strcmpi(filtertype, 'motion')	
	% First generate a horizontal line across the middle
	f = zeros(sze);
	f(floor(len/2)+1,1:len) = 1;   

	% Then rotate to specified angle
	f = imrotate(f,angle,'bilinear','loose');
	f = f/sum(f(:));
	
    else
	error('Unrecognized filter type');
	
    end	
	

%==================================================================
function [filtertype, sze, sigma, radius, len, angle] = checkargs(arg)
    
% Set defaults
    
    sze = [3 3];
    sigma = 0.5;
    radius = 5;
    len = 9;
    angle = 0;

    narg = length(arg);    

    filtertype = arg{1};
    if ~ischar(filtertype)
	error('filtertype must be specified as a string');
    end

    if strcmpi(filtertype, 'log')    
      if narg == 1
	  sze = [5 5];
      end
    end
    
    if strcmpi(filtertype, 'average')   || ...
       strcmpi(filtertype, 'gaussian')  || ...
       strcmpi(filtertype, 'log')
	if narg >= 2
	    sze = arg{2};
	    if isscalar(sze)
		sze = [sze sze];
	    end
	end
    end

    if strcmpi(filtertype, 'gaussian')  || ...
       strcmpi(filtertype, 'log')    
	if narg >= 3
	    sigma = arg{3};
	end
    end
    
    if strcmpi(filtertype, 'disk')  
	if narg >= 2
	    radius = arg{2};
	end
	sze = [2*radius+1 2*radius+1];
    end
    
    if strcmpi(filtertype, 'motion')      
	if narg >= 2
	    len = arg{2};
	end

	if narg >= 3
	    angle = arg{3};
	end
	
	% Ensure size is odd so that there is a middle point 
	% about which to rotate the filter by angle.
	if mod(len,2)         
	    sze = [len len];
	else
	    sze = [len+1 len+1];
	end
    end
    
    if ~isscalar(len) || len < 1
	error('length must be a scalar >= 1');
    end
    
    if ~isscalar(angle)
	error('angle must be a scalar');
    end

    if ~isscalar(radius) || radius < 1
	error('radius must be a scalar >= 1');
    end
    
    if ~isscalar(sigma) || sigma < 0.5
	error('sigma must be a scalar >= 0.5');
    end    
    
    if length(sze) > 2
	error('filter size must be a scalar or 2-vector');
    end
    
    if any(fix(sze) ~= sze)
	error('filter size must be integer');
    end	  

%==================================================================    
function a = col2gray(r,g,b)

    if nargin==0,
        error('Need input arguments.');
    end
    threeD = (ndims(r)==3); % Determine if input includes a 3-D array
    
    if nargin==1,
        if threeD,
            rgb = reshape(r(:),size(r,1)*size(r,2),3);
            a = zeros([size(r,1), size(r,2)]);
        else % Colormap
            rgb = r;
            a = zeros(size(r,1),1);
        end
        
    elseif nargin==2,
        error('Wrong number of arguments.');
        
    elseif nargin==3,
        warning(['COL2GRAY(R,G,B) is an obsolete syntax.',...
            'Use a three-dimensional array to represent RGB image.']);
        if (any(size(r)~=size(g)) || any(size(r)~=size(b))),
            error('R, G, and B must all be the same size.')
        end
        rgb = [r(:), g(:), b(:)];
        a = zeros(size(r));
    else
        error('Invalid input arguments.');
    end
    
    T = inv([1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703]);
    
    if isa(rgb, 'uint8')
        a = uint8(reshape(double(rgb)*T(1,:)', size(a)));
    elseif isa(rgb, 'uint16')
        a = uint16(reshape(double(rgb)*T(1,:)', size(a)));
    elseif isa(rgb, 'double')
        a = reshape(rgb*T(1,:)', size(a));
        a = min(max(a,0),1);
    end
    
    if ((nargin==1) && (~threeD)),    % col2gray(MAP)
        if ~isa(a, 'double')
            a = im2double(a);
        end
        a = [a,a,a];
    end

%==================================================================
function im=noiseonf_3d(rows,cols,dep,factor)

    if rem(rows,2)
        error('error: rows, cols, and dep must all be even')
    elseif rem(cols,2)
        error('error: rows, cols, and dep must all be even')
    elseif rem(dep,2)
        error('error: rows, cols, and dep must all be even')
    end
    
    imfft=fftshift(fftn(randn(cols,rows,dep)));
    mag=abs(imfft);
    phase=imfft./mag;
    
    [xi,yi,zi]=meshgrid(1:rows,1:cols,1:dep);
    radius=sqrt(xi.^2+yi.^2+zi.^2);
    
    radius(cols/2+1,rows/2+1,dep/2+1)=1;
    filter=1./(radius.^factor);
    im=real(ifftn(fftshift(filter.*phase)));
    %im=ifftn(fftshift(filter.*phase),'symmetric');
    im=rescale(im,0,1);
    im=im./sum(im(:));
    im=imresize(im(1:2:cols,1:2:rows,1:dep),2);

%==================================================================
function r=gen_points(model,num_grains,options)

    dim_num=3;
    it_fixed = 1;
    seed = 123456789;
    r = [];
    
    if model == 1
        init = -1; %pvt
        disp('Using the PVT model to generate particle centers ...')
        disp(['grains: ',num2str(num_grains),...
            ', batch: ',num2str(options.batch),', it_max: ',...
            num2str(options.it_max)])
    elseif model == 2 % cvt
        init = 0; %cvt
        disp('Using the CVT-Uniform model to generate particle centers ...')
        disp(['grains: ',num2str(num_grains),...
            ', batch: ',num2str(options.batch),', it_max: ',...
            num2str(options.it_max)])
    elseif model ==3
        init = 1; %hal
        disp('Using the CVT-Halton model to generate particle centers ...')
        disp(['grains: ',num2str(num_grains),...
            ', batch: ',num2str(options.batch),', it_max: ',...
            num2str(options.it_max)])
    end
    
    if model<=3
        sample = init;
        sample_num = num_grains*10;
        
        [ r, seed, it_num, it_diff, energy ] = cvt ( dim_num,...
            num_grains, options.batch, init, ...
            sample, sample_num, options.it_max, it_fixed, seed, r );
        
    elseif model==4 % CP
        
        Vert = [ ...
            -.5 -.5 -.5 ;
            .5 -.5 -.5 ;
            .5 .5 -.5 ;
            -.5 .5 -.5 ;
            -.5 -.5 .5 ;
            .5 -.5 .5 ;
            .5 .5 .5 ;
            -.5 .5 .5];
        
        xp=Vert(:,1); yp=Vert(:,2); zp=Vert(:,3);
        [xi,yi,zi] = csbinproc3d(xp,yp,zp,num_grains);
        
        disp('Using the CP model to generate particle centers ...')
        disp(['grains: ',num2str(num_grains),', intensity: ',...
            num2str(options.lam),', variance: ',num2str(options.var)])
        r = gen_varsig_3d(options.lam,num_grains,options.var,xi,yi,zi,xp,yp,zp);
        r=r';
        
    elseif model==5  % strauss
        
        Vert = [ ...
            -.5 -.5 -.5 ;
            .5 -.5 -.5 ;
            .5 .5 -.5 ;
            -.5 .5 -.5 ;
            -.5 -.5 .5 ;
            .5 -.5 .5 ;
            .5 .5 .5 ;
            -.5 .5 .5];
        
        xp=Vert(:,1); yp=Vert(:,2); zp=Vert(:,3);
        
        disp('Using the CP model to generate particle centers ...')
        disp(['grains: ',num2str(num_grains),', inhibition: ',...
            num2str(options.c),', delta: ',num2str(options.delta)])
        
        r = zeros(num_grains,3);
        % Generate the first point.
        r(1,:) = csbinproc3d(xp,yp,zp,1);
        
        % The following code is similar to the SSI process, except that we now have a
        % mechanism for accepting points that are closer than the inhibition distance.
        i = 1;
        while i<num_grains
            
            [sx, sy, sz] = csbinproc3d(xp,yp,zp,1);
            xt = [sx sy sz; r(1:i,:)];
            % Find the distance between the events.
            dist = ipdm(xt);
            % Find the distance between the candidate event
            % and the others that have been generated already.
            ind = find(dist(1:i) <= options.delta);
            m = length(ind);
            if m == 0
                % Then ok to keep the point - nothing is close.
                i = i+1;
                r(i,:) = [sx, sy, sz];
            elseif rand(1) <= options.c^m
                % Then ok to keep the point.
                i = i+1;
                r(i,:) = [sx, sy, sz];
            end
        end
        
        r=r';
        
    elseif model==6
        
        Vert = [ ...
            -.5 -.5 -.5 ;
            .5 -.5 -.5 ;
            .5 .5 -.5 ;
            -.5 .5 -.5 ;
            -.5 -.5 .5 ;
            .5 -.5 .5 ;
            .5 .5 .5 ;
            -.5 .5 .5];
        
        xp=Vert(:,1); yp=Vert(:,2); zp=Vert(:,3);
        
        disp('Using the homog model to generate particle centers ...')
        disp(['grains: ',num2str(num_grains)])
        
        r = csbinproc3d(xp, yp, zp, num_grains);
        r=r';
        
    end
    
%==================================================================
function [x,y,z] = csbinproc3d(xp, yp, zp, n)
    % CSBINPROC Generate homogeneous 3-D Poisson process.
    x = zeros(n,1);
    y = zeros(n,1);
    z = zeros(n,1);
    i = 1;
    % find the maximum and the minimum for a 'box' around
    % the region. Will generate uniform on this, and throw
    % out those points that are not inside the region.
    minx = min(xp);
    maxx = max(xp);
    miny = min(yp);
    maxy = max(yp);
    minz = min(zp);
    maxz = max(zp);
    cx = maxx-minx;
    cy = maxy - miny;
    cz = maxz - minz;
    
    while i <= n
        xt = rand(1)*cx + minx;
        yt = rand(1)*cy + miny;
        zt = rand(1)*cz + minz;
        k = inhull([xt,yt,zt],[xp,yp,zp]);
        if k == 1
            % it is in the region
            x(i) = xt;
            y(i) = yt;
            z(i) = zt;
            i = i+1;
        end
    end

%==================================================================
function [X] = gen_varsig_3d(lam,num_grains,r,xi,yi,zi,xp,yp,zp)
    
    % Get the number of children per parent.
    nchild = poissrnd(lam,1,num_grains);
    
    X = [];
    sig = r*eye(3);
    
    for i = 1:num_grains
        xc = randn(nchild(i),3)*sig + repmat([xi(i) yi(i) zi(i)],nchild(i),1);
        X = [X; xc];
    end
    
    % Find the ones that are in the region of interest.
    ind = find(inhull(X,[xp,yp,zp]));
    
    % Those are the children for the sample.
    x = X(ind,1);
    y = X(ind,2);
    z = X(ind,3);
    X=[x(:),y(:),z(:)];

%==================================================================%=======

function in = inhull(testpts,xyz,tess,tol)
% inhull: tests if a set of points are inside a convex hull
% usage: in = inhull(testpts,xyz)
% usage: in = inhull(testpts,xyz,tess)
% usage: in = inhull(testpts,xyz,tess,tol)
%
% arguments: (input)
%  testpts - nxp array to test, n data points, in p dimensions
%       If you have many points to test, it is most efficient to
%       call this function once with the entire set.
%
%  xyz - mxp array of vertices of the convex hull, as used by
%       convhulln.
%
%  tess - tessellation (or triangulation) generated by convhulln
%       If tess is left empty or not supplied, then it will be
%       generated.
%
%  tol - (OPTIONAL) tolerance on the tests for inclusion in the
%       convex hull. You can think of tol as the distance a point
%       may possibly lie outside the hull, and still be perceived
%       as on the surface of the hull. Because of numerical slop
%       nothing can ever be done exactly here. I might guess a
%       semi-intelligent value of tol to be
%
%         tol = 1.e-13*mean(abs(xyz(:)))
%
%       In higher dimensions, the numerical issues of floating
%       point arithmetic will probably suggest a larger value
%       of tol.
%
%       DEFAULT: tol = 0
%
% arguments: (output)
%  in  - nx1 logical vector
%        in(i) == 1 --> the i'th point was inside the convex hull.
%  
% Example usage: The first point should be inside, the second out
%
%  xy = randn(20,2);
%  tess = convhulln(xy);
%  testpoints = [ 0 0; 10 10];
%  in = inhull(testpoints,xy,tess)
%
% in = 
%      1
%      0
%
% A non-zero count of the number of degenerate simplexes in the hull
% will generate a warning (in 4 or more dimensions.) This warning
% may be disabled off with the command:
%
%   warning('off','inhull:degeneracy')
%
% See also: convhull, convhulln, delaunay, delaunayn, tsearch, tsearchn
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 3.0
% Release date: 10/26/06

% get array sizes
% m points, p dimensions
p = size(xyz,2);
[n,c] = size(testpts);
if p ~= c
    error 'testpts and xyz must have the same number of columns'
end
if p < 2
    error 'Points must lie in at least a 2-d space.'
end

% was the convex hull supplied?
if (nargin<3) || isempty(tess)
    tess = convhulln(xyz);
end
[nt,c] = size(tess);
if c ~= p
    error 'tess array is incompatible with a dimension p space'
end

% was tol supplied?
if (nargin<4) || isempty(tol)
    tol = 0;
end

% build normal vectors
switch p
    case 2
        % really simple for 2-d
        nrmls = (xyz(tess(:,1),:) - xyz(tess(:,2),:)) * [0 1;-1 0];
        
        % Any degenerate edges?
        del = sqrt(sum(nrmls.^2,2));
        degenflag = (del<(max(del)*10*eps));
        if sum(degenflag)>0
            warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
                ' degenerate edges identified in the convex hull'])
            
            % we need to delete those degenerate normal vectors
            nrmls(degenflag,:) = [];
            nt = size(nrmls,1);
        end
    case 3
        % use vectorized cross product for 3-d
        ab = xyz(tess(:,1),:) - xyz(tess(:,2),:);
        ac = xyz(tess(:,1),:) - xyz(tess(:,3),:);
        nrmls = cross(ab,ac,2);
        degenflag = repmat(false,nt,1);
    otherwise
        % slightly more work in higher dimensions,
        nrmls = zeros(nt,p);
        degenflag = repmat(false,nt,1);
        for i = 1:nt
            % just in case of a degeneracy
            nullsp = null(xyz(tess(i,2:end),:) - repmat(xyz(tess(i,1),:),p-1,1))';
            if size(nullsp,1)>1
                degenflag(i) = true;
                nrmls(i,:) = NaN;
            else
                nrmls(i,:) = nullsp;
            end
        end
        if sum(degenflag)>0
            warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
                ' degenerate simplexes identified in the convex hull'])
            
            % we need to delete those degenerate normal vectors
            nrmls(degenflag,:) = [];
            nt = size(nrmls,1);
        end
end

% scale normal vectors to unit length
nrmllen = sqrt(sum(nrmls.^2,2));
nrmls = nrmls.*repmat(1./nrmllen,1,p);

% center point in the hull
center = mean(xyz,1);

% any point in the plane of each simplex in the convex hull
a = xyz(tess(~degenflag,1),:);

% ensure the normals are pointing inwards
dp = sum((repmat(center,nt,1) - a).*nrmls,2);
k = dp<0;
nrmls(k,:) = -nrmls(k,:);

% We want to test if:  dot((x - a),N) >= 0
% If so for all faces of the hull, then x is inside
% the hull. Change this to dot(x,N) >= dot(a,N)
aN = sum(nrmls.*a,2);

% test, be careful in case there are many points
in = repmat(false,n,1);

% if n is too large, we need to worry about the
% dot product grabbing huge chunks of memory.
memblock = 1e6;
blocks = max(1,floor(n/(memblock/nt)));
aNr = repmat(aN,1,length(1:blocks:n));
for i = 1:blocks
    j = i:blocks:n;
    if size(aNr,2) ~= length(j),
        aNr = repmat(aN,1,length(j));
    end
    in(j) = all((nrmls*testpts(j,:)' - aNr) >= -tol,1)';
end

%==================================================================%=======


function d = ipdm(data1,varargin)
% ipdm: Inter-Point Distance Matrix
% usage: d = ipdm(data1)
% usage: d = ipdm(data1,data2)
% usage: d = ipdm(data1,prop,value)
% usage: d = ipdm(data1,data2,prop,value)
%
% Arguments: (input)
%  data1 - array of data points, each point is one row. p dimensional
%          data will be represented by matrix with p columns.
%          If only data1 is provided, then the distance matrix
%          is computed between all pairs of rows of data1.
%
%          If your data is one dimensional, it MUST form a column
%          vector. A row vector of length n will be interpreted as
%          an n-dimensional data set.
%
%  data2 - second array, supplied only if distances are to be computed
%          between two sets of points.
%
%
% Class support: data1 and data2 are assumed to be either
% single or double precision. I have not tested this code to
% verify its success on integer data of any class.
%
%
% Additional parameters are expected to be property/value pairs.
% Property/value pairs are pairs of arguments, the first of which
% (properties) must always be a character string. These strings
% may be shortened as long as the shortening is unambiguous.
% Capitalization is ignored. Valid properties for ipdm are:
%
%  'Metric', 'Subset', 'Limit', 'Result'
%
%  'Metric' - numeric flag - defines the distance metric used
%          metric = 2 --> (DEFAULT) Euclidean distance = 2-norm
%                         The standard distance metric.
%
%          metric = 1 --> 1-norm = sum of absolute differences
%                         Also sometimes known as the "city block
%                         metric", since this is the sum of the
%                         differences in each dimension.
%
%          metric = inf --> infinity-norm = maximum difference
%                         over all dimensions. The name refers
%                         to the limit of the p-norm, as p
%                         approaches infinity.
%
%          metric = 0 --> minimum difference over all dimensions.
%                         This is not really a useful norm in
%                         practice.
%
%          Note: while other distance metrics exist, IMHO, these
%          seemed to be the common ones.
%
%
%  'Result' - A string variable that denotes the style of returned
%          result. Valid result types are 'Array', 'Structure'.
%          Capitalization is ignored, and the string may be
%          shortened if you wish.
%
%          result = 'Array' --> (DEFAULT) A matrix of all
%                         interpoint distances will be generated.
%                         This array may be large. If this option
%                         is specified along with a minimum or
%                         maximum value, then those elements above
%                         or below the limiting values will be
%                         set as -inf or +inf, as appropriate.
%
%                         When any of 'LargestFew', 'SmallestFew',
%                         or 'NearestNeighbor' are set, then the
%                         resulting array will be a sparse matrix
%                         if 'array' is specified as the result.
%
%          result = 'Structure' --> A list of all computed distances,
%                         defined as a structure. This structure
%                         will have fields named 'rowindex',
%                         'columnindex', and 'distance'.
%
%                         This option will be useful when a subset
%                         criterion for the distances has been
%                         specified, since then the distance matrix
%                         may be very sparsely populated. Distances
%                         for pairs outside of the criterion will
%                         not be returned.
%
%
%  'Subset' - Character string, any of:
%
%          'All', 'Maximum', 'Minimum', 'LargestFew', 'SmallestFew',
%          'NearestNeighbor', 'FarthestNeighbor', or empty
%
%          Like properties, capitalization is ignored here, and
%          any unambiguous shortening of the word is acceptable.
%
%          DEFAULT = 'All'
%
%          Some interpoint distance matrices can be huge. Often
%          these matrices are too large to be fully retained in
%          memory, yet only the pair of points with the largest
%          or smallest distance may be needed. When only some
%          subset of the complete set of distances is of interest,
%          these options allow you to specify which distances will
%          be returned.
%
%          If 'result' is defined to be an array, then a sparse
%          matrix will be returned for the 'LargestFew', 'SmallestFew',
%          'NearestNeighbor', and 'FarthestNeighbor' subset classes.
%          'Minimum' and 'Maximum' will yield full matrices by
%          default. If a structure is specified, then only those
%          elements which have been identified will be returned.
%
%          Where a subset is specified, its limiting value is
%          specified by the 'Limit' property. Call that value k.
%
%
%          'All' -->     (DEFAULT) Return all interpoint distances
%
%          'Minimum' --> Only look for those distances above
%                        the cutoff k. All other distances will
%                        be returned as -inf.
%
%          'Maximum' --> Only look for those distances below
%                        the cutoff k. All other distances will
%                        be returned as +inf.
%
%          'SmallestFew' --> Only return the subset of the k
%                        smallest distances. Where only one data
%                        set is provided, only the upper triangle
%                        of the inter-point distance matrix will
%                        be generated since that matrix is symmetric.
%
%          'LargestFew' --> Only return the subset of the k
%                        largest distances. Where only one data
%                        set is provided, only the upper triangle
%                        of the inter-point distance matrix will
%                        be generated since that matrix is symmetric.
%
%          'NearestNeighbor' --> Only return the single nearest
%                        neighbor in data2 to each point in data1.
%                        No limiting value is required for this
%                        option. If multiple points have the same
%                        nearest distance, then return the first
%                        such point found. With only one input set,
%                        a point will not be its own nearest
%                        neighbor.
%
%                        Note that exact replicates in a single set
%                        will cause problems, since a sparse matrix
%                        is returned by default. Since they will have
%                        a zero distance, they will not show up in
%                        the sparse matrix. A structure return will
%                        show those points as having a zero distance
%                        though.
%
%          'FarthestNeighbor' --> Only return the single farthest
%                        neighbor to each point. No limiting value
%                        is required for this option. If multiple
%                        points have the same farthest distance,
%                        then return the first such point found.
%
%
%  'Limit' - scalar numeric value or []. Used only when some
%           Subset is specified.
%
%          DEFAULT = []
%
%
%  'ChunkSize' - allows a user with lower RAM limits
%          to force the code to only grab smaller chunks of RAM
%          at a time (where possible). This parameter is specified
%          in bytes of RAM. The default is 32 megabytes, or 2^22
%          elements in any piece of the distance matrix. Only some
%          options will break the problem into chunks, thus as long
%          as a full matrix is expected to be returned, there seems
%          no reason to break the problem up into pieces.
%
%          DEFAULT = 2^25
%
%
% Arguments: (output)
%  d     - array of interpoint distances, or a struct wth the
%          fields {'rowindex', 'columnindex', 'distance'}.
%
%          d(i,j) represents the distance between point i
%          (from data1) and point j (from data2).
%
%          If only one (n1 x p) array is supplied, then d will
%          be an array of size == [n1,n1].
%
%          If two arrays (of sizes n1 x p and n2 x p) then d
%          will be an array of size == [n1,n2].
%
%
% Efficiency considerations:
%  Where possible, this code will use bsxfun to compute its
%  distances.
%
%
% Example:
%  Compute the interpoint distances between all pairs of points
%  in a list of 5 points, in 2 dimensions and using Euclidean
%  distance as the distance metric.
%
%  A = randn(5,2);
%  d = ipdm(A,'metric',2)
%  d =
%            0       2.3295       3.2263       2.0263       2.8244
%       2.3295            0       1.1485      0.31798       1.0086
%       3.2263       1.1485            0       1.4318       1.8479
%       2.0263      0.31798       1.4318            0       1.0716
%       2.8244       1.0086       1.8479       1.0716            0
%
% (see the demo file for many other examples)
%
% See also: pdist
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/26/08

% Default property values
params.Metric = 2;
params.Result = 'array';
params.Subset = 'all';
params.Limit = [];
params.ChunkSize = 2^25;

% untangle the arguments
if nargin<1
    % if called with no arguments, then the user probably
    % needs help. Give it to them.
    help ipdm
    return
end

% were two sets of data provided?
pvpairs = {};
if nargin==1
    % only 1 set of data provided
    dataflag = 1;
    data2 = [];
else
    if ischar(varargin{1})
        dataflag = 1;
        data2 = [];
        pvpairs = varargin;
    else
        dataflag = 2;
        data2 = varargin{1};
        if nargin>2
            pvpairs = varargin(2:end);
        end
    end
end

% get data sizes for later
[n1,dim] = size(data1);
if dataflag == 2
    n2 = size(data2,1);
end

% Test the class of the input variables
if ~(isa(data1,'double') || isa(data1,'single')) || ...
        ((dataflag == 2) && ~(isa(data2,'double') || isa(data2,'single')))
    error('data points must be either single or double precision variables.')
end

% do we need to process any property/value pairs?
if nargin>2
    params = parse_pv_pairs(params,pvpairs);
    
    % check for problems in the properties
    
    % was a legal Subset provided?
    if ~isempty(params.Subset) && ~ischar(params.Subset)
        error('If provided, ''Subset'' must be character')
    elseif isempty(params.Subset)
        params.Subset = 'all';
    end
    valid = {'all','maximum','minimum','largestfew','smallestfew', ...
        'nearestneighbor','farthestneighbor'};
    ind = find(strncmpi(params.Subset,valid,length(params.Subset)));
    if (length(ind)==1)
        params.Subset = valid{ind};
    else
        error(['Invalid Subset: ',params.Subset])
    end
    
    % was a limit provided?
    if ~ismember(params.Subset,{'all','nearestneighbor','farthestneighbor'}) && ...
            isempty(params.Limit)
        error('No limit provided, but a Subset that requires a limit value was specified')
    end
    % check the limit values for validity
    if length(params.Limit)>1
        error('Limit must be scalar or empty')
    end
    
    switch params.Subset
        case {'largestfew', 'smallestfew'}
            % must be at least 1, and an integer
            if (params.Limit<1) || (round(params.Limit)~=params.Limit)
                error('Limit must be a positive integer for LargestFew or NearestFew')
            end
    end
    
    % was a legal Result provided?
    if isempty(params.Result)
        params.result = 'Array';
    elseif ~ischar(params.Result)
        error('If provided, ''Result'' must be character or empty')
    end
    valid = {'array','structure'};
    ind = find(strncmpi(params.Result,valid,length(params.Result)));
    if (length(ind)==1)
        params.Result = valid{ind};
    else
        error(['Invalid Result: ',params.Subset])
    end
    
    % check for the metric
    if isempty(params.Metric)
        params.Metric = 2;
    elseif (length(params.Metric)~=1) || ~ismember(params.Metric,[0 1 2 inf])
        error('If supplied, ''Metric'' must be a scalar, and one of [0 1 2 inf]')
    end
end % if nargin>2

% If Metric was given as 2, but the dimension is only 1, then it will
% be slightly faster (and equivalent) to use the 1-norm Metric.
if (dim == 1) && (params.Metric == 2)
    params.Metric = 1;
end

% Can we use bsxfun to compute the interpoint distances?
% Older Matlab releases will not have bsxfun, but if it is
% around, it will ne both faster and less of a memory hog.
params.usebsxfun = (5==exist('bsxfun','builtin'));

% check for dimension mismatch if 2 sets
if (dataflag==2) && (size(data2,2)~=dim)
    error('If 2 point sets provided, then both must have the same number of columns')
end

% Total number of distances to compute, in case I must do it in batches
if dataflag==1
    n2 = n1;
end
ntotal = n1*n2;

% FINALLY!!! Compute inter-point distances
switch params.Subset
    case 'all'
        % The complete set of interpoint distances. There is no need
        % to break this into chunks, since we must return all distances.
        % If that is too much to compute in memory, then it will fail
        % anyway when we try to store the result. bsxfun will at least
        % do the computation efficiently.
        
        % One set or two?
        if dataflag == 1
            d = distcomp(data1,data1,params);
        else
            d = distcomp(data1,data2,params);
        end
        
        % Must we return it as a struct?
        if params.Result(1) == 's'
            [rind,cind] = ndgrid(1:size(d,1),1:size(d,2));
            ds.rowindex = rind(:);
            ds.columnindex = cind(:);
            ds.distance = d(:);
            d = ds;
        end
        
    case {'minimum' 'maximum'}
        % There is no reason to break this into pieces if the result
        % sill be filled in the end with +/- inf. Only break it up
        % if the final result is a struct.
        if ((ntotal*8)<=params.ChunkSize) || (params.Result(1) == 'a')
            % its small enough to do it all at once
            
            % One set or two?
            if dataflag == 1
                d = distcomp(data1,data1,params);
            else
                d = distcomp(data1,data2,params);
            end
            
            % Must we return it as a struct?
            if params.Result(1) == 'a'
                % its an array, fill the unwanted distances with +/- inf
                if params.Subset(2) == 'i'
                    % minimum
                    d(d<=params.Limit) = -inf;
                else
                    % maximum
                    d(d>=params.Limit) = +inf;
                end
            else
                % a struct will be returned
                if params.Subset(2) == 'i'
                    % minimum
                    [dist.rowindex,dist.columnindex] = find(d>=params.Limit);
                else
                    % maximum
                    [dist.rowindex,dist.columnindex] = find(d<=params.Limit);
                end
                dist.distance = d(dist.rowindex + n1*(dist.columnindex-1));
                d = dist;
            end
            
        else
            % we need to break this into chunks. This branch
            % will always return a struct.
            
            % this is the number of rows of data1 that we will
            % process at a time.
            bs = floor(params.ChunkSize/(8*n2));
            bs = min(n1,max(1,bs));
            
            % Accumulate the result into a cell array. Do it this
            % way because we don't know in advance how many elements
            % that we will find satisfying the minimum or maximum
            % limit specified.
            accum = cell(0,1);
            
            % now loop over the chunks
            batch = 1:bs;
            while ~isempty(batch)
                
                % One set or two?
                if dataflag == 1
                    dist = distcomp(data1(batch,:),data1,params);
                else
                    dist = distcomp(data1(batch,:),data2,params);
                end
                
                % big or small as requested
                if ('i'==params.Subset(2))
                    % minimum value specified
                    [I,J,V] = find(dist>=params.Limit);
                else
                    % maximum limit
                    [I,J] = find(dist<=params.Limit);
                    I = I(:);
                    J = J(:);
                    V = dist(I + (J-1)*length(batch));
                    I = I + (batch(1)-1);
                end
                
                % and stuff them into the cell structure
                if ~isempty(V)
                    accum{end+1,1} = [I,J,V(:)]; %#ok
                end
                
                % increment the batch
                batch = batch + bs;
                if batch(end)>n1
                    batch(batch>n1) = [];
                end
                
            end
            
            % convert the cells into one flat array
            accum = cell2mat(accum);
            
            if isempty(accum)
                d.rowindex = [];
                d.columnindex = [];
                d.distance = [];
            else
                % we found something
                
                % sort on the second column, to put them in a reasonable order
                accum = sortrows(accum,[2 1]);
                
                d.rowindex = accum(:,1);
                d.columnindex = accum(:,2);
                d.distance = accum(:,3);
            end
            
        end
        
    case {'smallestfew' 'largestfew'}
        % find the k smallest/largest distances. k is
        % given by params.Limit
        
        % if only 1 set, params.Limit must be less than n*(n-1)/2
        if dataflag == 1
            params.Limit = min(params.Limit,n1*(n1-1)/2);
        end
        
        % is this a large problem?
        if ((ntotal*8) <= params.ChunkSize)
            % small potatoes
            
            % One set or two?
            if dataflag == 1
                dist = distcomp(data1,data1,params);
                % if only one data set, set the diagonal and
                % below that to +/- inf so we don't find it.
                temp = find(tril(ones(n1,n1),0));
                if params.Subset(1) == 's'
                    dist(temp) = inf;
                else
                    dist(temp) = -inf;
                end
            else
                dist = distcomp(data1,data2,params);
            end
            
            % sort the distances to find those we need
            if ('s'==params.Subset(1))
                % smallestfew
                [val,tags] = sort(dist(:),'ascend');
            else
                % largestfew
                [val,tags] = sort(dist(:),'descend');
            end
            val = val(1:params.Limit);
            tags = tags(1:params.Limit);
            
            % recover the row and column index from the linear
            % index returned by sort in tags.
            [d.rowindex,d.columnindex] = ind2sub([n1,size(dist,2)],tags);
            
            % create the matrix as a sparse one or a struct?
            if params.Result(1)=='a'
                % its an array, so make the array sparse.
                d = sparse(d.rowindex,d.columnindex,val,n1,size(dist,2));
            else
                % a structure
                d.distance = val;
            end
            
        else
            % chunks
            
            % this is the number of rows of data1 that we will
            % process at a time.
            bs = floor(params.ChunkSize/(8*n2));
            bs = min(n1,max(1,bs));
            
            % We need to find the extreme cases. There are two possible
            % algorithms, depending on how many total elements we will
            % search for.
            % 1. Only a very few total elements.
            % 2. A relatively large number of total elements, forming
            %    a significant fraction of the total set.
            %
            % Case #1 would suggest to retain params.Limit numberr of
            % elements from each batch, then at the end, sort them all
            % to find the best few. Case #2 will result in too many
            % elements to retain, so we must distinguish between these
            % alternatives.
            if (8*params.Limit*n1/bs) <= params.ChunkSize
                % params.Limit is small enough to fall into case #1.
                
                % Accumulate the result into a cell array. Do it this
                % way because we don't know in advance how many elements
                % that we will find satisfying the minimum or maximum
                % limit specified.
                accum = cell(0,1);
                
                % now loop over the chunks
                batch = (1:bs)';
                while ~isempty(batch)
                    % One set or two?
                    if dataflag == 1
                        dist = distcomp(data1(batch,:),data1,params);
                        k = find(tril(ones(length(batch),n2),batch(1)-1));
                        if ('s'==params.Subset(1))
                            dist(k) = inf;
                        else
                            dist(k) = -inf;
                        end
                    else
                        dist = distcomp(data1(batch,:),data2,params);
                    end
                    
                    % big or small as requested, keeping only the best
                    % params.Limit number of elements
                    if ('s'==params.Subset(1))
                        % minimum value specified
                        [tags,tags] = sort(dist(:),1,'ascend');
                        tags = tags(1:bs);
                        [I,J] = ndgrid(batch,1:n2);
                        ijv = [I(tags),J(tags),dist(tags)];
                    else
                        % maximum limit
                        [tags,tags] = sort(dist(:),1,'descend');
                        tags = tags(1:bs);
                        [I,J] = ndgrid(batch,1:n2);
                        ijv = [I(tags),J(tags),dist(tags)];
                    end
                    % and stuff them into the cell structure
                    accum{end+1,1} = ijv; %#ok
                    
                    % increment the batch
                    batch = batch + bs;
                    if batch(end)>n1
                        batch(batch>n1) = [];
                    end
                end
                
                % convert the cells into one flat array
                accum = cell2mat(accum);
                
                % keep only the params.Limit best of those singled out
                accum = sortrows(accum,3);
                if ('s'==params.Subset(1))
                    % minimum value specified
                    accum = accum(1:params.Limit,:);
                else
                    % minimum value specified
                    accum = accum(end + 1 - (1:params.Limit),:);
                end
                d.rowindex = accum(:,1);
                d.columnindex = accum(:,2);
                d.distance = accum(:,3);
                
                % create the matrix as a sparse one or a struct?
                if params.Result(1)=='a'
                    % its an array, so make the array sparse.
                    d = sparse(d.rowindex,d.columnindex,d.distance,n1,size(dist,2));
                end
                
            else
                % params.Limit forces us into the domain of case #2.
                % Here we cannot retain params.Limit elements from each chunk.
                % so we will grab each chunk and append it to the best elements
                % found so far, then filter out the best after each chunk is
                % done. This may be slower than we want, but its the only way.
                ijv = zeros(0,3);
                
                % loop over the chunks
                batch = (1:bs)';
                while ~isempty(batch)
                    % One set or two?
                    if dataflag == 1
                        dist = distcomp(data1(batch,:),data1,params);
                        k = find(tril(ones(length(batch),n2),batch(1)-1));
                        if ('s'==params.Subset(1))
                            dist(k) = inf;
                        else
                            dist(k) = -inf;
                        end
                    else
                        dist = distcomp(data1(batch,:),data2,params);
                    end
                    
                    [I,J] = ndgrid(batch,1:n2);
                    ijv = [ijv;[I(:),J(:),dist(:)]]; %#ok
                    
                    % big or small as requested, keeping only the best
                    % params.Limit number of elements
                    if size(ijv,1) > params.Limit
                        if ('s'==params.Subset(1))
                            % minimum value specified
                            [tags,tags] = sort(ijv(:,3),1,'ascend');
                        else
                            [tags,tags] = sort(ijv(:,3),1,'ascend');
                        end
                        ijv = ijv(tags(1:params.Limit),:);
                    end
                    
                    % increment the batch
                    batch = batch + bs;
                    if batch(end)>n1
                        batch(batch>n1) = [];
                    end
                end
                
                % They are fully trimmed down. stuff a structure
                d.rowindex = ijv(:,1);
                d.columnindex = ijv(:,2);
                d.distance = ijv(:,3);
                
                % create the matrix as a sparse one or a struct?
                if params.Result(1)=='a'
                    % its an array, so make the array sparse.
                    d = sparse(d.rowindex,d.columnindex,d.distance,n1,size(dist,2));
                end
                
            end
            
        end
        
    case {'nearestneighbor' 'farthestneighbor'}
        % find the closest/farthest neighbor for every point
        
        % is this a large problem? Or a 1-d problem?
        if dim == 1
            % its a 1-d nearest/farthest neighbor problem. we can
            % special case these easily enough, and all the distance
            % metric options are the same in 1-d.
            
            % first split it into the farthest versus nearest cases.
            if params.Subset(1) == 'f'
                % farthest away
                
                % One set or two?
                if dataflag == 1
                    [d2min,minind] = min(data1);
                    [d2max,maxind] = max(data1);
                else
                    [d2min,minind] = min(data2);
                    [d2max,maxind] = max(data2);
                end
                
                d.rowindex = (1:n1)';
                d.columnindex = repmat(maxind,n1,1);
                d.distance = repmat(d2max,n1,1);
                
                % which endpoint was further away?
                k = abs((data1 - d2min)) >= abs((data1 - d2max));
                if any(k)
                    d.columnindex(k) = minind;
                    d.distance(k) = d2min;
                end
                
            else
                % nearest. this is mainly a sort and some fussing around.
                d.rowindex = (1:n1)';
                d.columnindex = ones(n1,1);
                d.distance = zeros(n1,1);
                
                % One set or two?
                if dataflag == 1
                    % if only one data point, then we are done
                    if n1 == 2
                        % if exactly two data points, its trivial
                        d.columnindex = [2 1];
                        d.distance = repmat(abs(diff(data1)),2,1);
                    elseif n1>2
                        % at least three points. do a sort.
                        [sorted_data,tags] = sort(data1);
                        
                        % handle the first and last points separately
                        d.columnindex(tags(1)) = tags(2);
                        d.distance(tags(1)) = sorted_data(2) - sorted_data(1);
                        d.columnindex(tags(end)) = tags(end-1);
                        d.distance(tags(end)) = sorted_data(end) - sorted_data(end-1);
                        
                        ind = (2:(n1-1))';
                        
                        d1 = sorted_data(ind) - sorted_data(ind-1);
                        d2 = sorted_data(ind+1) - sorted_data(ind);
                        
                        k = d1 < d2;
                        d.distance(tags(ind(k))) = d1(k);
                        d.columnindex(tags(ind(k))) = tags(ind(k)-1);
                        k = ~k;
                        d.distance(tags(ind(k))) = d2(k);
                        d.columnindex(tags(ind(k))) = tags(ind(k)+1);
                    end % if n1 == 2
                else
                    % Two sets of data. still really a sort and some fuss.
                    if n2 == 1
                        % there is only one point in data2
                        d.distance = abs(data1 - data2);
                        % d.columnindex is already set correctly
                    else
                        % At least two points in data2
                        % We need to sort all the data points together, but also
                        % know which points from each set went where. ind12 and
                        % bool12 will help keep track.
                        ind12 = [1:n1,1:n2]';
                        bool12 = [zeros(n1,1);ones(n2,1)];
                        [sorted_data,tags] = sort([data1;data2]);
                        
                        ind12 = ind12(tags);
                        bool12 = bool12(tags);
                        
                        % where did each point end up after the sort?
                        loc1 = find(~bool12);
                        loc2 = find(bool12);
                        
                        % for each point in data1, what is the (sorted) data2
                        % element which appears most nearly to the left of it?
                        cs = cumsum(bool12);
                        leftelement = cs(loc1);
                        
                        % any points which fell below the minimum element in data2
                        % will have a zero for the index of the element on their
                        % left. fix this.
                        leftelement = max(1,leftelement);
                        
                        % likewise, any point greater than the max in data2 will
                        % have an n2 in left element. this too will be a problem
                        % later, so fix it.
                        leftelement = min(n2-1,leftelement);
                        
                        % distance to the left hand element
                        dleft = abs(sorted_data(loc1) - sorted_data(loc2(leftelement)));
                        dright = abs(sorted_data(loc1) - sorted_data(loc2(leftelement+1)));
                        
                        % find the points which are closer to the left element in data2
                        k = (dleft < dright);
                        d.distance(ind12(loc1(k))) = dleft(k);
                        d.columnindex(ind12(loc1(k))) = ind12(loc2(leftelement(k)));
                        k = ~k;
                        d.distance(ind12(loc1(k))) = dright(k);
                        d.columnindex(ind12(loc1(k))) = ind12(loc2(leftelement(k)+1));
                        
                    end % if n2 == 1
                end % if dataflag == 1
            end % if params.Subset(1) == 'f'
            
            % create the matrix as a sparse one or a struct?
            if params.Result(1)=='a'
                % its an array, so make the array sparse.
                d = sparse(d.rowindex,d.columnindex,d.distance,n1,n2);
            end
            
        elseif (ntotal>1000) && (((params.Metric == 0) && (params.Subset(1) == 'n')) || ...
                ((params.Metric == inf) && (params.Subset(1) == 'f')))
            % nearest/farthest neighbour in n>1 dimensions, but for an
            % infinity norm metric. Reduce this to a sequence of
            % 1-d problems, each of which will be faster in general.
            % do this only if the problem is moderately large, since
            % we must overcome the extra overhead of the recursive
            % calls to ipdm.
            
            % do the first dimension
            if dataflag == 1
                d = ipdm(data1(:,1),data1(:,1),'subset',params.Subset,'metric',params.Metric,'result','struct');
            else
                d = ipdm(data1(:,1),data2(:,1),'subset',params.Subset,'metric',params.Metric,'result','struct');
            end
            
            % its slightly different for nearest versus farthest here
            % now, loop over dimensions
            for i = 2:dim
                if dataflag == 1
                    di = ipdm(data1(:,i),data1(:,i),'subset',params.Subset,'metric',params.Metric,'result','struct');
                else
                    di = ipdm(data1(:,i),data2(:,i),'subset',params.Subset,'metric',params.Metric,'result','struct');
                end
                
                % did any of the distances change?
                if params.Metric == 0
                    % the 0 norm, with nearest neighbour, so take the
                    % smallest distance in any dimension.
                    k = d.distance > di.distance;
                else
                    % inf norm. so take the largest distance across dimensions
                    k = d.distance < di.distance;
                end
                
                if any(k)
                    d.distance(k) = di.distance(k);
                    d.columnindex(k) = di.columnindex(k);
                end
            end
            
            % create the matrix as a sparse one or a struct?
            if params.Result(1)=='a'
                % its an array, so make the array sparse.
                d = sparse(d.rowindex,d.columnindex,d.distance,n1,n2);
            end
            
        elseif ((ntotal*8) <= params.ChunkSize)
            % None of the other special cases apply, so do it using brute
            % force for the small potatoes problem.
            
            % One set or two?
            if dataflag == 1
                dist = distcomp(data1,data1,params);
            else
                dist = distcomp(data1,data2,params);
            end
            
            % if only one data set and if a nearest neighbor
            % problem, set the diagonal to +inf so we don't find it.
            if (dataflag==1) && (n1>1) && ('n'==params.Subset(1))
                diagind = (1:n1) + (0:n1:(n1^2-1));
                dist(diagind) = +inf;
            end
            
            if ('n'==params.Subset(1))
                % nearest
                [val,j] = min(dist,[],2);
            else
                % farthest
                [val,j] = max(dist,[],2);
            end
            
            % create the matrix as a sparse one or a struct?
            if params.Result(1)=='a'
                % its an array, so make the array sparse.
                d = sparse((1:n1)',j,val,n1,size(dist,2));
            else
                % a structure
                d.rowindex = (1:n1)';
                d.columnindex = j;
                d.distance = val;
            end
            
        else
            
            % break it into chunks
            bs = floor(params.ChunkSize/(8*n2));
            bs = min(n1,max(1,bs));
            
            % pre-allocate the result
            d.rowindex = (1:n1)';
            d.columnindex = zeros(n1,1);
            d.distance = zeros(n1,1);
            
            % now loop over the chunks
            batch = 1:bs;
            while ~isempty(batch)
                
                % One set or two?
                if dataflag == 1
                    dist = distcomp(data1(batch,:),data1,params);
                else
                    dist = distcomp(data1(batch,:),data2,params);
                end
                
                % if only one data set and if a nearest neighbor
                % problem, set the diagonal to +inf so we don't find it.
                if (dataflag==1) && (n1>1) && ('n'==params.Subset(1))
                    diagind = 1:length(batch);
                    diagind = diagind + (diagind-2+batch(1))*length(batch);
                    dist(diagind) = +inf;
                end
                
                % big or small as requested
                if ('n'==params.Subset(1))
                    % nearest
                    [val,j] = min(dist,[],2);
                else
                    % farthest
                    [val,j] = max(dist,[],2);
                end
                
                % and stuff them into the result structure
                d.columnindex(batch) = j;
                d.distance(batch) = val;
                
                % increment the batch
                batch = batch + bs;
                if batch(end)>n1
                    batch(batch>n1) = [];
                end
                
            end
            
            % did we need to return a struct or an array?
            if params.Result(1) == 'a'
                % an array. make it a sparse one
                d = sparse(d.rowindex,d.columnindex,d.distance,n1,n2);
            end
            
        end % if dim == 1
        
end  % switch params.Subset

% End of mainline

% ======================================================
% begin subfunctions
% ======================================================
    function d = distcomp(set1,set2,params)
        % Subfunction to compute all distances between two sets of points
        dim = size(set1,2);
        % can we take advantage of bsxfun?
        % Note: in theory, there is no need to loop over the dimensions. We
        % could Just let bsxfun do ALL the work, then wrap a sum around the
        % outside. In practice, this tends to create large intermediate
        % arrays, especially in higher numbers of dimensions. Its also when
        % we might gain here by use of a vectorized code. This will only be
        % a serious gain when the number of points is relatively small and
        % the dimension is large.
        if params.usebsxfun
            % its a recent enough version of matlab that we can
            % use bsxfun at all.
            n1 = size(set1,1);
            n2 = size(set2,1);
            if (dim>1) && ((n1*n2*dim)<=params.ChunkSize)
                % its a small enough problem that we might gain by full
                % use of bsxfun
                switch params.Metric
                    case 2
                        d = sum(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim])).^2,3);
                    case 1
                        d = sum(abs(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim]))),3);
                    case inf
                        d = max(abs(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim]))),[],3);
                    case 0
                        d = min(abs(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim]))),[],3);
                end
            else
                % too big, so that the ChunkSize will have been exceeded, or just 1-d
                if params.Metric == 2
                    d = bsxfun(@minus,set1(:,1),set2(:,1)').^2;
                else
                    d = abs(bsxfun(@minus,set1(:,1),set2(:,1)'));
                end
                for i=2:dim
                    switch params.Metric
                        case 2
                            d = d + bsxfun(@minus,set1(:,i),set2(:,i)').^2;
                        case 1
                            d = d + abs(bsxfun(@minus,set1(:,i),set2(:,i)'));
                        case inf
                            d = max(d,abs(bsxfun(@minus,set1(:,i),set2(:,i)')));
                        case 0
                            d = min(d,abs(bsxfun(@minus,set1(:,i),set2(:,i)')));
                    end
                end
            end
        else
            % Cannot use bsxfun. Sigh. Do things the hard (and slower) way.
            n1 = size(set1,1);
            n2 = size(set2,1);
            if params.Metric == 2
                % Note: While some people might use a different Euclidean
                % norm computation based on expanding the square of the
                % difference of two numbers, that computation is inherantly
                % inaccurate when implemented in floating point arithmetic.
                % While it might be faster, I won't use it here. Sorry.
                d = (repmat(set1(:,1),1,n2) - repmat(set2(:,1)',n1,1)).^2;
            else
                d = abs(repmat(set1(:,1),1,n2) - repmat(set2(:,1)',n1,1));
            end
            for i=2:dim
                switch params.Metric
                    case 2
                        d = d + (repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1)).^2;
                    case 1
                        d = d + abs(repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1));
                    case inf
                        d = max(d,abs(repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1)));
                    case 0
                        d = min(d,abs(repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1)));
                end
            end
        end
        % if 2 norm, then we must sqrt at the end
        if params.Metric==2
            d = sqrt(d);
        end

% ==============================================================
%    end main ipdm
%    begin included function - parse_pv_pairs
% ==============================================================
function params=parse_pv_pairs(params,pv_pairs)
% parse_pv_pairs: parses sets of property value pairs, allows defaults
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.
%
% Example usage:
% First, set default values for the parameters. Assume we
% have four parameters that we wish to use optionally in
% the function examplefun.
%
%  - 'viscosity', which will have a default value of 1
%  - 'volume', which will default to 1
%  - 'pie' - which will have default value 3.141592653589793
%  - 'description' - a text field, left empty by default
%
% The first argument to examplefun is one which will always be
% supplied.
%
%   function examplefun(dummyarg1,varargin)
%   params.Viscosity = 1;
%   params.Volume = 1;
%   params.Pie = 3.141592653589793
%
%   params.Description = '';
%   params=parse_pv_pairs(params,varargin);
%   params
%
% Use examplefun, overriding the defaults for 'pie', 'viscosity'
% and 'description'. The 'volume' parameter is left at its default.
%
%   examplefun(rand(10),'vis',10,'pie',3,'Description','Hello world')
%
% params =
%     Viscosity: 10
%        Volume: 1
%           Pie: 3
%   Description: 'Hello world'
%
% Note that capitalization was ignored, and the property 'viscosity'
% was truncated as supplied. Also note that the order the pairs were
% supplied was arbitrary.

npv = length(pv_pairs);
n = npv/2;

if n~=floor(n)
    error 'Property/value pairs must come in PAIRS.'
end
if n<=0
    % just return the defaults
    return
end

if ~isstruct(params)
    error 'No structure for defaults was supplied'
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for i=1:n
    p_i = lower(pv_pairs{2*i-1});
    v_i = pv_pairs{2*i};
    
    ind = strmatch(p_i,lpropnames,'exact');
    if isempty(ind)
        ind = find(strncmp(p_i,lpropnames,length(p_i)));
        if isempty(ind)
            error(['No matching property found for: ',pv_pairs{2*i-1}])
        elseif length(ind)>1
            error(['Ambiguous property name: ',pv_pairs{2*i-1}])
        end
    end
    p_i = propnames{ind};
    
    % override the corresponding default in params.
    % Use setfield for comptability issues with older releases.
    params = setfield(params,p_i,v_i); %#ok
    
end



