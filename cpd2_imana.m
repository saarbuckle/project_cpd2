function varargout = cpd2_imana(what,varargin)
% saarbuckle 08/2019
% for analysis of cpd1 dataset

%% ------------------------- Directories ----------------------------------
cwd        = cd; % get current directory when called, and return at end of script
% paths to project directories
codeDir         = '/Users/sarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/projects/project_passivePatterns';
baseDir         = '/Users/sarbuckle/DATA/cpd2';   % base directory for analysis
behavDir        = [baseDir '/data'];             % behavioural data directory           
anatomicalDir   = [baseDir '/anatomicals'];  
freesurferDir   = [baseDir '/surfaceFreesurfer'];  
wbDir           = [baseDir '/surfaceWB'];
atlasDir        = '/Users/sarbuckle/DATA/Atlas_templates/FS_LR_164';
caretDir        = [baseDir '/surfaceCaret'];     
gpCaretDir      = [caretDir '/fsaverage_sym'];
imagingDir      = [baseDir '/imaging_data'];
regDir          = [baseDir '/RegionOfInterest'];   
pcmDir          = [baseDir '/PCM_models'];
glmDir          = [baseDir '/GLM_firstlevel_fast2'];  
%% ------------------------- Subj info ------------------------------------
subj_name       = {'p02','p03','p05','p06','p07','p08','p09','p10'};
%% ------------------------- Exp info -------------------------------------
%glmToUse = 6; % FAST without high-pass filtering
% Analysis parameters
cpd_chords=[eye(5); 1,1,0,0,0; 1,0,1,0,0; 1,0,0,1,0; 1,0,0,0,1; 0,1,1,0,0;...
    0,1,0,1,0; 0,1,0,0,1; 0,0,1,1,0; 0,0,1,0,1; 0,0,0,1,1;...
    0,0,1,1,1; 0,1,0,1,1; 0,1,1,0,1; 0,1,1,1,0; 1,0,0,1,1;...
    1,0,1,0,1; 1,0,1,1,0; 1,1,0,0,1; 1,1,0,1,0; 1,1,1,0,0;...
    0,1,1,1,1; 1,0,1,1,1; 1,1,0,1,1; 1,1,1,0,1; 1,1,1,1,0;...
    1,1,1,1,1;];

chords = pp1_imana('chords');

numDigits = sum(chords,2);
%% ------------------------- ROI info -------------------------------------
% Freesurfer & ROI parameters
hemName    = {'LeftHem','RightHem'};   % freesurfer hemisphere folder names 
atlasA     = 'x';
atlasname  = 'fsaverage_sym';
hem        = {'lh','rh'};                                                   % left & right hemi folder names/prefixes
hemLetter   = {'L','R'};
%hemNum     = [1 2];
regname    = {'ba3A','ba3B','ba1','ba2','rM1','cM1','SII','S1','M1','SPLa','SPLp','VPL','MGN','LGN'};        % Cortical ROIs, 5 = S1, 6 = M1; Thalamus = 12                                             % roi names, independent of hemisphere    
cortical   = repmat([1 1 1 1 1 1 1 1 1 1 1 0 0 0],1,2);
thalamic   = repmat([0 0 0 0 0 0 0 0 0 0 0 1 1 1],1,2);
regSide    = [ones(size(regname)),...                                       % Hemisphere of the roi
                ones(size(regname)).*2];                                    % [1 = left hemi (contra), 2 = right hemi (ipsi)]
regType    = [1:length(regname),...                                         % roi # (within hemisphere)
                1:length(regname)];
numregions = max(regType);                                   


%% ------------------------- ANALYSES -------------------------------------
switch what
    case 'SPM_remapDir'                                                     % remap the image files to SPM structures for each subject
        D = cpd2_imana('LIST'); % get subject info
        % remaps to imagingDir
        for s = unique(D.sn)'
            subjName = D.origSN{s};
            fprintf('%s : ',subjName);
            SPM_fname  = fullfile(glmDir,subjName,'SPM.mat');      % where SPM folder currently is  
            rawDataDir = fullfile(imagingDir,subjName);            % where associated volume files are saved
            spmj_move_rawdata(SPM_fname,rawDataDir);
            fprintf('done\n');
        end
    case 'renumber_chords'
        % The order of the chords in my experiment (pp1) are different from
        % that in Naveed's (cpd1). Thus, we need to renumber the chords
        % from Naveed's to match mine. We edit the SPM_info.mat file.
        
        % The single finger and two-finger chords are in the correct order.
        % The third & fourth finger chords are different condition numbers.
        
        error('this case has been run on this dataset already. DO NOT RE-RUN, else you will mess up condition labels!')
        %***************************************************************%
        % % % DANGER: THIS CASE HAS BEEN ALREADY RUN. DO NOT RE-RUN % % % 
        %***************************************************************%
        
        
%         I = cpd2_imana('LIST_subjs'); % get subject info
%         pp_chords = pp1_imana('chords');
%         [~,~,ppNum]=intersect(cpd_chords,pp_chords,'rows','stable'); % find correspondence between chord condition numbers from the two experiments


%         % re-label the glm info file
%         for s = unique(I.sn)'
%             D = load(fullfile(glmDir,subj_name{s},'SPM_info.mat'));
%             D.chord_OLD = D.chord;
%             X = indicatorMatrix('identity',D.chord);
%             D.chord = X*ppNum;
%             save(fullfile(glmDir,subj_name{s},'SPM_info.mat'),'-struct','D');
%             fprintf('s%02d done\n',s);
%         end

%         % re-label the betas file
%         pp_chords = pp1_imana('chords');
%         [~,~,ppNum]=intersect(cpd_chords,pp_chords,'rows','stable');
%         T = load(fullfile(regDir,'glm1_reg_betas.mat'));
%         for i=1:numel(T.sn)
%             X = indicatorMatrix('identity',T.tt{i});
%             T.tt_OLD{i,1} = T.tt{i};
%             T.tt{i} = X*ppNum;
%         end
%         save(fullfile(regDir,'glm1_reg_betas.mat'),'-struct','T');
    
    case 'LIST_subjs'
        % reads txt file from basedir that contains subject info
        D = dload(fullfile(baseDir,'subj_info.txt'));
        
        if nargout==0
            fprintf('\nSN\torigSN\tAge\tGender\tHandedness\tID');
            fprintf('\n--\t------\t---\t------\t----------\t--');
            for s = unique(D.sn)'
                S = getrow(D,ismember(D.sn,s));
                fprintf('\n%02d\t%s\t%d\t%d\t%d\t\t%s',S.sn,S.origSN{1},S.age,S.gender,S.handedness,S.ID{1});
            end
            fprintf('\n');
        else
            D = rmfield(D,{'ID'});
            varargout = {D};
        end
    case 'getSubjs'
        % returns vector of subject numbers for which we have fully
        % analyzed fmri data available.
        I = cpd2_imana('LIST_subjs');
        sn = unique(I.sn)';
        varargout = {sn};      
    case 'LIST_rois'
        fprintf('\nROI#\tName\tHemipshere\tCortical\tThalamic');
        fprintf('\n----\t----\t----------\t--------\t-------\n');
        for r = 1:length(regSide)
            fprintf('%d\t%s\t%s\t\t%d\t\t%d\n',r,regname{r-(regSide(r)-1)*length(regname)},hem{regSide(r)},cortical(r-(regSide(r)-1)*length(regname)),thalamic(r-(regSide(r)-1)*length(regname)));
        end  
    
    case 'SURF_WBresample'   % Reslice indiv surfaces into fs_lr standard mesh
        % This reslices from the individual surfaces into the the fs_lr
        % standard mesh - This replaces calls to freesurfer_registerXhem,
        % freesurfer_mapicosahedron_xhem, & caret_importfreesurfer. It
        % requires connectome wb to be installed, added to the bash_profile
        % (on terminal), and updated on the startup.m file
        vararginoptions(varargin, {'sn'});
        fprintf('reslicing %s...',subj_name{sn});
        surf_resliceFS2WB(subj_name{sn}, freesurferDir, wbDir,'resolution','164k'); 
        fprintf('done\n');    
        
    case 'makeGrayMask'
        % makes gray-matter FOV mask 
        I = cpd2_imana('LIST_subjs'); % get subject info
        sn = cpd2_imana('getSubjs');
        for s=sn
            nam     = {};
            nam{1}  = fullfile(anatomicalDir,I.origSN{s},'rmask_noskull.nii');
            nam{2}  = fullfile(anatomicalDir,I.origSN{s},['c1' I.origSN{s} '_anatomical.nii']);
            out     = fullfile(anatomicalDir,I.origSN{s},'rmask_gray.nii');
            spm_imcalc(nam,out,'i1>0 & i2>0.4');            
        end
    case 'SURF_segThalamus'
        % segment the thalamus nuclei
        I = cpd2_imana('LIST_subjs');
        sn = cpd2_imana('getSubjs');
        vararginoptions(varargin,{'sn'});
        % do segmentation
        % % ***SEGMENTATION MUST BE DONE VIA TERMINAL
%         cmd = ['segmentThalamicNuclei.sh ' subj_name{sn} ' ' freesurferDir];
%         fprintf('%s\n',cmd);
%         [status,result] = system(cmd);
%         if status; error(result); end
        for s=sn
            % convert segmentation to nii and move segmentation file to
            % anatomicals directory
            inName  = fullfile(freesurferDir,I.origSN{s},'mri','ThalamicNuclei.v10.T1.mgz');
            outName = fullfile(anatomicalDir,I.origSN{s},[I.origSN{s} '_thalamicNuclei_raw.nii']);
            cmd = sprintf('mri_convert %s %s', inName, outName);
            fprintf('%s\n',cmd);
            [status,result] = system(cmd);
            if status; error(result); end
            % reslice the thalamic masks into same space as epis (this will
            % also mask them to be within bounds of epi FOV)
            mask = fullfile(anatomicalDir,I.origSN{s},['rmask_gray.nii']);
            T = spm_vol(outName);
            outName = fullfile(anatomicalDir,I.origSN{s},[I.origSN{s} '_thalamicNuclei.nii']);
            M = spm_vol(mask);
            spmj_reslice_vol(T,M.dim,M.mat,outName);
            fprintf('%s...done.\n',I.origSN{s});     
        end
    case 'ROI_define_OLD'                                                       % Define rois: BA rois, M1/S1 cut to hand area, BA rois cut to hand area
        % Define the ROIs of the group fsaverage atlas for each subject's
        % surface reconstruction. 
        % Output saved for each subject ('s#_regions.mat').
        % The save variable for each subject is a cell array of size
        % {1,#rois}. 
        I = cpd2_imana('LIST_subjs'); % get subject info
        sn = cpd2_imana('getSubjs');
        vararginoptions(varargin,{'sn'});

        linedef = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
        
        thalamicMarker = [8133,8115,8109; 8233,8215,8209]; % values of voxels assigned to each of the 3 thalamic regions (vpl,mgn,lgn) split by hemisphere
        corticalMarker = [1,2,3,4,5,6,1,2,7,8]; % values of surface nodes assigned to the cortical rois in their respective paint files
        corticalFile   = {'D1','D1','D1','D1','D1','D1','D2','D2','D2','D2'}; % different files for coritcal rois
        
        for s = sn % for each subject
            R = {};
            j = 1; % overall region ticker
            subjName = I.origSN{s};
            fprintf('%s...',subjName);
            caretSubjDir   = fullfile(caretDir,[atlasA subjName]);
            for h = 1:2 % per hemisphere
                t = 1; % thalamic ticker
                D1 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI_pp1.paint'])); % ba3A, ba3B, ba1, ba2, rM1, cM1
                D2 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI.paint']));     % premade ROIs from fsaverage_sym: M1 and S1
                for r = 1:numregions % make regions
                    % make R region structure for participant
                    R{j}.name  = [subjName '_' regname{r} '_' hem{h}];
                    % Mapping of rois is different if cortical or thalamic:
                    if cortical(r)==1 && thalamic(r)==0 % cortical rois
                        eval(['D=' corticalFile{r} ';']);
                        R{j}.type     = 'surf_nodes';
                        R{j}.location = find(D.data(:,1)==corticalMarker(r));
                        R{j}.white    = fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                        R{j}.pial     = fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                        R{j}.topo     = fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                        R{j}.linedef  = linedef;
                        R{j}.flat     = fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                        R{j}.image    = fullfile(glmDir,subjName,'mask.nii,1'); % load glm mask from first glm
                        if strcmp(corticalFile{r},'D1')
                            R{j}.origFile = 'ROI_pp1.paint';
                        elseif strcmp(corticalFile{r},'D2')
                            R{j}.origFile = 'ROI.paint';
                        end
                    elseif thalamic(r)==1 && cortical(r)==0 % thalamic rois
                        % make temporary functional mask of thalamic nuclei
                        R{j}.type     = 'roi_image';
                        R{j}.value    = thalamicMarker(h,t);
                        R{j}.file     = fullfile(anatomicalDir,subjName,sprintf('%s_thalamicNuclei.nii',subjName));
                        t = t+1;
                    end
                    % add indexing to structure R
                    R{j}.regNum   = j;
                    R{j}.regType  = r;
                    R{j}.cortical = cortical(r);
                    R{j}.thalamic = thalamic(r);
                    R{j}.hemi     = h;
                    j = j+1;
                end    
            end
            % threshold across CS, both hemispheres
            exculdePairs = [1,5; 1,6; 2,5; 2,6; 3,5; 3,6; 7,8];
            excludePairs = [exculdePairs; exculdePairs+numregions]; 
            R = region_calcregions(R,'exclude',excludePairs,'exclude_thres',0.75);
            cd(regDir);
            save(['regions_' subjName '.mat'],'R');
            fprintf('..done\n');
            clear R
        end
    case 'ROI_define'                                                       % Define rois: BA rois, M1/S1 cut to hand area, BA rois cut to hand area
        % Define the ROIs of the group fsaverage atlas for each subject's
        % surface reconstruction. 
        % Output saved for each subject ('s#_regions.mat').
        % The save variable for each subject is a cell array of size
        % {1,#rois}. 
        I = cpd2_imana('LIST_subjs');
        sn = cpd2_imana('getSubjs');
        vararginoptions(varargin,{'sn'});

        linedef = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
        
        thalamicMarker = [8133,8115,8109,8233,8215,8209]; % values of voxels assigned to each of the 3 thalamic regions (vpl,mgn,lgn) split by hemisphere
        corticalMarker = [1,2,3,4,5,6,7,1,2,7,8]; % values of surface nodes assigned to the cortical rois
        corticalFile   = {'D1','D1','D1','D1','D1','D1','D1','D2','D2','D2','D2'}; % different files for coritcal rois
        
        for s = sn % for each subject
            R = {};
            j = 1; % overall region ticker
            t = 1; % thalamic ticker
            subjName = I.origSN{s};
            fprintf('%s...',subjName);
            mask      = fullfile(glmDir,subjName,'mask.nii,1');  % load mask file now 
            for h = 1:2 % per hemisphere
                wbFilePrefix = fullfile(wbDir,subjName,[subjName '.' hemLetter{h}]);
                regFile1 = fullfile(atlasDir,['ROI_pp1.164k.' hemLetter{h} '.label.gii']);
                regFile2 = fullfile(atlasDir,['ROI.164k.' hemLetter{h} '.label.gii']);
                D1       = gifti(regFile1);
                D2       = gifti(regFile2);
                for r = 1:numregions % make regions
                    % make R region structure for participant
                    R{j}.name     = [subjName '_' regname{r} '_' hem{h}];
                    R{j}.regNum   = j;
                    R{j}.regType  = r;
                    R{j}.cortical = cortical(r);
                    R{j}.thalamic = thalamic(r);
                    R{j}.hemi     = h;
                    % Mapping of rois is different if cortical or thalamic:
                    if cortical(r)==1 && thalamic(r)==0 % cortical surface rois
                        eval(['D=' corticalFile{r} ';']);
                        R{j}.type     = 'surf_nodes_wb';
                        R{j}.location = find(D.cdata(:,1)==corticalMarker(r));
                        R{j}.white    = [wbFilePrefix '.white.164k.surf.gii'];
                        R{j}.pial     = [wbFilePrefix '.pial.164k.surf.gii'];
                        R{j}.linedef  = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
                        R{j}.image    = mask;    % functional mask
                        if strcmp(corticalFile{r},'D1')
                            R{j}.origFile = regFile1;
                        elseif strcmp(corticalFile{r},'D2')
                            R{j}.origFile = regFile2;
                        end
                    elseif thalamic(r)==1 && cortical(r)==0 % thalamic rois
                        % make temporary functional mask of thalamic nuclei
                        R{j}.type   = 'roi_image';
                        R{j}.value  = thalamicMarker(t);
                        R{j}.file   = fullfile(anatomicalDir,subjName,sprintf('%s_thalamicNuclei.nii',subjName));
                        t = t+1;
                    end
                    % update roi ticker
                    j = j+1;
                end    
            end
            % threshold across CS, both hemispheres
            exculdePairs = [1,5; 1,6; 2,5; 2,6; 3,5; 3,6; 8,9];
            excludePairs = [exculdePairs; exculdePairs+numregions]; % do same exclusion for rois in both hemispsheres
            R = region_calcregions(R,'exclude',excludePairs,'exclude_thres',0.75);
            cd(regDir);
            save(['regions_' subjName '.mat'],'R');
            fprintf('..done\n');
            clear R
        end
    case 'ROI_make_nii'                                                     % OPTIONAL   :  Convert ROI def (.mat) into multiple .nii files (to check!)
        I = pp1_imana('LIST_subjs');
        sn = pp1_imana('getSubjs');
        glm = 4;
        vararginoptions(varargin,{'sn'});        
        
        for s=sn
            subjName = I.origSN{s};
            glmSubjDir = fullfile(glmDir{glm},subjName);

            cd(glmSubjDir);
            % load ROI definition
            load(fullfile(regDir,['regions_' subjName '.mat']));

            % loop over rois
            for rr = [1:28]
                % mask volume
                mask = fullfile(glmSubjDir,'mask.nii');           
                % Save region file as nifti
                cd(regDir);
                region_saveasimg(R{rr},mask);      
            end
            
        end        

    case 'ROI_defineOLD'                                                      
        % Define the ROIs of the group fsaverage atlas for each subject's
        % surface reconstruction. 
        % Output saved for each subject ('s#_regions.mat').
        % The save variable for each subject is a cell array of size
        % {1,#rois}. 

        I = cpd2_imana('LIST'); % get subject info

        linedef = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
        tmpMaskFile = {};  % for volumetric sub-cortical regions (these tmp masks are deleted after mapping to subject)
        for s = unique(I.sn)'
            
            caretSubjDir = fullfile(caretDir,[atlasA I.origSN{s}]);
            mask         = fullfile(glmDir,I.origSN{s},'mask.nii,1');  
            
            for h = 1:2 
                
                D1 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI_pp1.paint'])); % ba3A, ba3B, ba1, ba2, rM1, cM1
                D2 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI.paint']));     % premade ROIs from fsaverage_sym: M1 and S1
                
                for r = 1:numregions
                    if r<7; D = D1; rr = r;         % pp1 rois
                    elseif (r==7 || r==7+numregions); D = D2; rr = 1;    % S1
                    elseif (r==8 || r==8+numregions); D = D2; rr = 2;    % M1
                    elseif (r==9 || r==9+numregions); D = D2; rr = 7;    % SPLa
                    elseif (r==10 || r==10+numregions); D = D2; rr = 8;  % SPLp
                    end
                    idx = r+(h-1)*numregions;
                    % make R region structure for participant
                    R{idx}.name  = [I.origSN{s} '_' regname{r} '_' hem{h}];
                    if r<11 % cortical surface regions
                        R{idx}.type     = 'surf_nodes';
                        R{idx}.location = find(D.data(:,1)==rr);
                        R{idx}.white    = fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                        R{idx}.pial     = fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                        R{idx}.topo     = fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                        R{idx}.linedef  = linedef;
                        R{idx}.flat     = fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                        R{idx}.image    = mask; % functional mask
%                     elseif r==11 % thalamus volumetric region
%                         if h==1; rr = 10; elseif h==2; rr = 49; end % 10 is left thalamus, 49 is right thalamus
%                         R{idx}.type     = 'roi_image';
%                         R{idx}.file     = fullfile(scAnatDir,subj_name{s},sprintf('%s_%d_%d.nii',subj_name{s},h,rr));
%                         R{idx}.value    = 1;
%                         R{idx}.image    = fullfile(glmDir{1},subj_name{s},sprintf('mask_tmp%d%d.nii',h,rr)); % functional mask
%                         % make temporary functional mask of thalamus
%                         % This is done b/c region_calcregions assumes mask
%                         % is same dimension as roi image. Thus,
%                         % region_calcregions will basically redo the mask,
%                         % but that's okay.
%                         F      = spm_vol(R{idx}.file);
%                         F.data = spm_read_vols(F);
%                         M      = spm_vol(mask);
%                         M.data = spm_read_vols(M);
%                         V(1)   = F;
%                         V(2)   = M;
%                         spm_imcalc(V,R{idx}.image,'i1.*i2');
%                         tmpMaskFile{end+1} = R{idx}.image;
%                         clear V
                    end
                end    
            end
            %R = region_calcregions(R,'exclude',[2,6; 5,6],'exclude_thres',0.75);
            exculdeRoi = [1,5; 1,6; 2,5; 2,6; 3,5; 3,6; 7,8];
            R = region_calcregions(R,'exclude',[exculdeRoi; exculdeRoi+numregions],'exclude_thres',0.75);
            cd(regDir);
            save(['regions_' I.origSN{s} '.mat'],'R');
            fprintf('\n %s done\n',I.origSN{s})
            clear R
        end
        for i = 1:numel(tmpMaskFile) 
            delete(tmpMaskFile{i});
        end
        
    case 'ROI_getBetas'                                                     % STEP 5.3   :  Harvest activity patterns from specified rois
        glm = 1;
        I   = cpd2_imana('LIST_subjs'); % get subject info
        sn  = cpd2_imana('getSubjs');
        roi = 1:28;
        append = 0; % just add betas to currently existing datastructure
        vararginoptions(varargin,{'sn','glm','roi','append'});

        T=[];
        fprintf('extracting betas\n');
        if append
            fprintf('...and appending them to already existing file\n');
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        end
        % harvest
        for s=sn % for each subj
            subjName = I.origSN{s};
            fprintf('%s...',subjName);
            % load files
            load(fullfile(glmDir, subjName,'SPM.mat'));  % load subject's SPM data structure (SPM struct)
            D = load(fullfile(glmDir, subjName,'SPM_info.mat'));
            load(fullfile(regDir,sprintf('regions_%s.mat',subjName)));          % load subject's region parcellation & depth structure (R)
            
            % TR img info
            V = SPM.xY.VY; 

            % remove run means from patterns
            D.part = D.run + 8*(D.sess-1);
            C0   = indicatorMatrix('identity',D.part); 
            ofInterest = 1:size(C0,1); % indicies for regressors of interest
            
            for r = roi % for each region
                % get raw data/psc for voxels in region
                Y = region_getdata(V,R{r});  % Data Y is N x P (P is in order of transpose of R{r}.depth)
                %P = region_getdata(Q,R{r});
                % estimate region betas
                [betaW,resMS,~,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','runwise');
                % toss stuff into output structure
                S.sn                 = s;
                S.roi                = r;
                S.tt                 = {D.chord};
                S.run                = {D.run};
                S.sess               = {D.sess};
                S.part               = {D.part};
                S.numDigits          = {numDigits(D.chord)};
                % remove nuisance regressor betas
                betaUW               = bsxfun(@rdivide,beta,sqrt(resMS));
                betaUW               = betaUW(ofInterest,:);
                betaW                = betaW(ofInterest,:);
                raw_beta             = beta(ofInterest,:);
                % add data to output structure
                S.betaW_noRunMean    = {betaW-C0*pinv(C0)*betaW};
                S.betaUW_noRunMean   = {betaUW-C0*pinv(C0)*betaUW};
                S.betaW              = {betaW};        
                S.betaUW             = {betaUW};  
                S.raw_beta           = {raw_beta};
                %S.psc                = {P};
                S.resMS              = {resMS};
                S.xyzcoord           = {R{r}.data'}; % excl already applied
                S.depth = {[NaN]};     try S.depth = {R{r}.depth(R{r}.excl==0,:)'}; catch end
                S.flatcoord = {[NaN]}; try S.flatcoord = {R{r}.flatcoord(R{r}.excl==0,1:2)'}; catch end
                T = addstruct(T,S);
                fprintf('%d.',r)
            end
            fprintf('..done\n');
        end
        % save T
        save(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)),'-struct','T'); 
        fprintf('\n')
    case 'ROI_pattConsist'                                                
        % Crossvalidated Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 
        % This stat is useful for determining which GLM model yeilds least
        % variable beta estimates. That GLM should be the one you use for
        % future analysis cases.

        glm = 1;
        I = cpd2_imana('LIST'); % get subject info
        sn = unique(I.sn)';
        roi = 1:8;
        removeMean = 1;
        conds = 1:31;
        vararginoptions(varargin,{'sn','glm','roi','removeMean'});
        
        % % Calculate pattern consistency for each roi, each subj.
        % Do so separately per session per subject.
        R = []; % output structure
        for g = glm
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',g))); % loads in struct 'T'
            for r = roi
                for s = sn
                    % get subj data
                    S      = getrow(T,(T.sn==s & T.roi==r));
                    b      = [];
                    b.beta = S.raw_beta{1};
                    b.run  = S.run{1};
                    b.part = S.part{1};
                    b.tt   = S.tt{1};
                    b.sess = S.sess{1};
                    b      = getrow(b,ismember(b.tt,conds)); % restrict to specific conditions
                    %for ii = unique(b.sess)' % per session
                    bs = b;    
                    %bs = getrow(b,b.sess==ii);
                        % calculate the pattern consistency
                        rs.r2              = rsa_patternConsistency(bs.beta,bs.part,bs.tt,'removeMean',removeMean);
                        [rs.r2_cv,rs.r_cv] = rsa_patternConsistency_crossval(bs.beta,bs.part,bs.tt,'removeMean',removeMean);
                        rs.sn              = s;
                        rs.roi             = r;
                        rs.glm             = g;
                        rs.numConds        = numel(conds);
                        rs.passive         = 0;
                        %rs.sess            = ii;
                        rs.removeMean      = removeMean;
                        R = addstruct(R,rs);
                    %end
                end
            end
        end
        %pivottable(R.glm,R.sn,R.r2,'mean','numformat','%0.4f');
        %pivottable(R.glm,R.sn,R.r2_cv,'mean','numformat','%0.4f');
        varargout = {R};
        % output arranged such that each row is an roi, each col is subj
    case 'ROI_stats'
        % housekeeping
        glm       = 1;
        numDigits = sum(chords,2);
        I         = cpd2_imana('LIST_subjs'); % get subject info
        To        = []; % output structures
        Td        = [];
        % get data
        T   = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        roi = unique(T.roi)';
        % do stats
        for s = unique(I.sn)' % for each subject
            fprintf('\nSubject: %d\n',s)
            D           = load(fullfile(glmDir, I.origSN{s}, 'SPM_info.mat'));   % load subject's trial structure
            D.numDigits = arrayfun(@(x) numDigits(x),D.chord);
            D.part      = (D.sess-1)*8+ D.run;
            C0          = indicatorMatrix('identity',D.part);
            numRun      = size(C0,2);
            numCond     = numel(unique(D.chord));
            ofInterest  = 1:(numCond*numRun); % indicies for regressors of interest
            
            for r = roi % for each region
                S = getrow(T,(T.sn==s & T.roi==r)); % subject's region data
                betaW       = S.betaW{1}; 
                betaW_nmean = betaW(ofInterest,:)-C0*pinv(C0)*betaW(ofInterest,:); % run mean subtraction  
                % % Toverall structure stats
                % crossval second moment matrix
                [G,Sig]      = pcm_estGCrossval(betaW_nmean(ofInterest,:),D.part,D.chord);
                So.sig       = rsa_vectorizeIPM(Sig);
                So.G         = rsa_vectorizeIPM(G);
                So.G_wmean   = rsa_vectorizeIPM(pcm_estGCrossval(betaW(ofInterest,:),D.part,D.chord));
                % squared dissimilarities
                So.ldc_wmean = rsa.distanceLDC(betaW,D.part,D.chord);        % rdm crossvalidated, on patterns without run mean patterns removed
                So.ldc       = rsa.distanceLDC(betaW_nmean,D.part,D.chord);  % rdm crossvalidated, patterns with run means removed
                % PSC
%                 So.psc       = nanmean(S.psc{1},2)';
%                 if glm==2 % only regressors for chords
%                     So.psc_chord = [1:31]; 
%                     So.psc_numD  = numDigits;
%                 elseif glm==3 % regressors for chords and thumb response
%                     So.psc_chord = [1:32]; 
%                     So.psc_numD  = [numDigits,99];
%                 end
                % Calculate avg. betas for each condition
                Q            = [];
                Q.raw_beta   = S.raw_beta{1};
                Q.tt         = D.chord;
                Q            = tapply(Q,{'tt'},{'raw_beta','mean'});
                So.avg_betas = mean(Q.raw_beta,2)';
                So.avg_tt    = Q.tt';
                % indexing fields
                So.sn       = s;
                So.roi      = r;
                So.numVox   = size(betaW,2);
                So.regSide  = regSide(r);
                So.regType  = regType(r);
                To          = addstruct(To,So);
                % calc avg. chord patterns for each number of digits 
                d           = [];
                d.betaW     = betaW_nmean;
                d.numDigits = D.numDigits;
                d.run       = D.run;
                d.part      = D.part;
                d.roi       = ones(size(d.run)).*r;
                d.chord     = D.chord;
                d           = getrow(d,d.chord<32);
                d5          = getrow(d,d.numDigits==5);
                d           = getrow(d,d.numDigits<5 & d.numDigits>0);
                d           = tapply(d,{'numDigits','part','roi'},{'betaW','mean'});
                d           = addstruct(d,d5);
                d           = rmfield(d,{'chord'});
                % calc distance between avg patterns for 1 finger up to
                % 5 finger chords:
                td.ldc = rsa.distanceLDC(d.betaW,d.part,d.numDigits)';
                td.distPair  = [1:10]';
                td.digitDiff = [1;2;3;4;1;2;3;1;2;1];
                td.roi       = ones(10,1).*r;
                td.sn        = ones(10,1).*s;
                Td           = addstruct(Td,td);
                fprintf('%d.',r)
            end % each region
        end % each subject

        % % save
        save(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)),'-struct','To');
        save(fullfile(regDir,sprintf('glm%d_reg_TnumDigits.mat',glm)),'-struct','Td');
        fprintf('done.\n')
    case 'ROI_getSingleFingerTuning'
        glm = 1;
        roi = 1:8;
        I = cpd2_imana('LIST'); % get subject info
        sn = unique(I.sn);
        vararginoptions(varargin,{'glm','roi','sn'});
        conds = 1:5;
        numCond = numel(conds);
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        D = []; % output structure
        v = ones(2,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b           = [];
            b.beta      = T.raw_beta{ii};
            b.tt        = T.tt{ii};
            b           = getrow(b,ismember(b.tt,conds));
            b           = tapply(b,{'tt'},{'beta','mean'});
            Gtt         = cov(b.beta');
            % rescale betaNoMax by setting the lowest value to 0,
            % calculating distance from maxB to each value.
            b.betaScaled= b.beta + abs(min(b.beta,[],1));
            maxBs       = max(b.betaScaled,[],1);
            avgDists    = sum(repmat(maxBs,numCond,1) - b.betaScaled,1)./(numCond-1);
            sftBeta     = mean(avgDists./maxBs); % normalize distances by max distance, such that values = 1 indicate exact single finger preference
            sftEV       = saSims('SFT:calcExpectedValue',Gtt,numel(unique(T.run{ii})),size(b.beta,2)); % expected value of the null
            % add to output structure
            d.passive   = v;
            d.sn        = v.*T.sn(ii);
            d.roi       = v.*T.roi(ii);
            d.glm       = v.*glm;
            d.sft       = [sftBeta;sftEV];
            d.isEV      = [0;1];
            D = addstruct(D,d);
        end
        varargout = {D};
    case 'ROI_rdmStability'
        glm = 1;
        roi = 1:8;
        I = cpd2_imana('LIST'); % get subject info
        sn = unique(I.sn);
        vararginoptions(varargin,{'glm','roi','sn'});
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        % housekeeping
        T = getrow(T,ismember(T.sn,sn));
        D = []; % output structure
        take = logical(tril(ones(numel(sn)),-1));
        for r = roi
            t = getrow(T,T.roi==r);
            R = corr(t.ldc');
            d = [];
            d.numSN = numel(sn); 
            d.roi   = r;
            d.corrs = R(take)';
            d.corr  = mean(R(take));
            % calc confidence bounds
            rz        = fisherz(d.corrs)'; 
            d.is_mean = fisherinv(mean(rz));
            d.is_LB   = fisherinv(d.is_mean - 1.96*stderr(rz));
            d.is_UB   = fisherinv(d.is_mean + 1.96*stderr(rz));
            D = addstruct(D,d);
        end
        varargout = {D};
    
    case 'ROI_calcSFT'
        % wrapper to perform SFT analysis
        I = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        glm = [];
        roi = 1:28;
        vararginoptions(varargin,{'glm','roi','sn'});
        
        numSim = 100; % # simulated datasets per model per participant (mvnrnd and sparse tuning models)
        conds  = 1:5; % conditions to analyze (single finger conditions)
        
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        
        fprintf('\nsubj\troi\tsvar\t\tevar\t\ts/e ratio\n');
        fprintf('----\t---\t----\t\t----\t\t---------\n');
        D = []; % output structure
        v = ones(3,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b          = [];
            b.beta     = T.raw_beta{ii}; 
            b.tt       = T.tt{ii};
            b.run      = T.run{ii};
            b          = getrow(b,ismember(b.tt,conds));

            % some simulation params needed:
            numVox  = ceil(size(b.beta,2)/5)*5; % round up so equal # of voxels per condition (for sparse patterns)
            numRun  = numel(unique(T.run{ii}));
            numCond = numel(conds);
            % explicitly estimate run mean covariances
            %[evar,svar] = pp1_imana('estErrVarMeans',b.beta,b.tt,b.run);
            
            % remove run means
            C0 = indicatorMatrix('identity',b.run);
            b.beta = b.beta -C0*pinv(C0)*b.beta;
            [evar,svar] = pp1_imana('estErrVar',b.beta,b.tt,b.run);
            
            % accounting for mean pattern separately, for noisy data the
            % signal variance can become negative. It's a bummer. To
            % counteract this, we set signal sigma to be realmin, so as to
            % avoid weird selectivity values.
            if svar<0
                svar = realmin;
            end
            
            b       = tapply(b,{'tt'},{'beta','mean'});
            Gtt     = cov(b.beta'); % single-finger G
            
            % 1. calc tuning of actual voxel data
            sftBeta = pp1_imana('estSingleFingerTuning',b.beta);
            sftBeta = mean(sftBeta);
            
            % 2. calc expected tuning of voxels with ~N(0,G)
            [sftEV,sftDistEV] = pp1_imana('SFT:calcExpectedValueG',evar,svar,Gtt,numVox,numRun,numSim); % expected value of the null
            
            % 3. calc expected tuning of voxels with ~sparse tuning, but where signal-to-noise equals that voxel data
            [sftSp,sftDistSp] = pp1_imana('SFT:calcExpectedValueSparse',evar,svar,1,numVox,numRun,numSim); % 1= tuned to one condition (perfectly sparse)
            
            % 4. do prob test for each subject on their gauss sft distribution
            pEV = sum(sftBeta<=sftDistEV)/length(sftDistEV);
            if isempty(pEV)
                pEV = realmin;
            end
            
            % 5. do prob test for each subject on their sparse sft distribution
            pSp = sum(sftBeta<=sftDistSp)/length(sftDistSp);
            if isempty(pSp)
                pSp = realmin;
            end
            
            % add to output structure
            d.passive   = v;
            d.sn        = v.*T.sn(ii);
            d.roi       = v.*T.roi(ii);
            d.glm       = v.*glm;
            d.evar      = v.*evar;
            d.svar      = v.*svar;
            d.sft       = [sftBeta;sftEV;sftSp];
            d.sftProb   = [0;pEV;pSp];
            d.isEV      = [0;1;2];
            D = addstruct(D,d);
            fprintf('%s\t%02d\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),svar,evar,svar/evar);
        end
        save(fullfile(regDir,sprintf('sft_glm%d',glm)),'-struct','D');
        varargout = {D};
    case 'ROI_calcSFT_fthres'
        % wrapper to perform SFT analysis, specifically focused on voxels
        % that show significant finger resposnes (voxel F-stat > F-crit)
        I   = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        glm = [];
        roi = 1:12;
        vararginoptions(varargin,{'glm','roi','sn'});
        
        numSim = 100; % # simulated datasets per model per participant (mvnrnd and sparse tuning models)
        conds  = 1:5; % conditions to analyze (single finger conditions)
        
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        
        fprintf('\nsubj\troi\tsig voxels (%%)\t\tsvar\t\tevar\t\ts/e ratio');
        fprintf('\n----\t---\t--------------\t\t----\t\t----\t\t---------\n');
        D = []; % output structure
        v = ones(3,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b          = [];
            b.beta     = T.raw_beta{ii};
            b.tt       = T.tt{ii};
            b.run      = T.run{ii};
            b          = getrow(b,ismember(b.tt,conds));
            [F,Fcrit]  = pp1_imana('calcFstat',b.beta,b.tt,b.run);
            numVoxOrig = size(b.beta,2);
            b.beta     = b.beta(:,F>=Fcrit); % using non-cv f stat as thresholder
            numVoxSig  = sum(F>=Fcrit);
            if numVoxSig==0 % if no voxels meet significance threshold, boot
                d.passive   = v;
                d.sn        = v.*T.sn(ii);
                d.roi       = v.*T.roi(ii);
                d.numVoxF   = v.*numVoxSig;
                d.numVoxTot = v.*size(T.raw_beta{ii},2);
                d.avgF      = v.*nan;
                d.glm       = v.*glm;
                d.sft       = [nan;nan;nan];
                d.sftProb   = [nan;nan;nan];
                d.isEV      = [0;1;2];
                D = addstruct(D,d);
                fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,nan,nan,nan);
                continue
            end
            % some simulation params needed:
            %numVox  = ceil(size(b.beta,2)/5)*5; % round up so equal # of voxels per condition (for sparse patterns)
            numVox  = ceil(numVoxSig/5)*5; % round up so equal # of voxels per condition
            numRun  = numel(unique(T.run{ii}));
            numCond = numel(conds);
            % explicitly estimate run mean covariances
            %[evar,svar] = pp1_imana('estErrVarMeans',b.beta,b.tt,b.run);
            
            % remove run means from these specific voxels
            C0 = indicatorMatrix('identity',b.run);
            b.beta = b.beta -C0*pinv(C0)*b.beta;
            [evar,svar] = pp1_imana('estErrVar',b.beta,b.tt,b.run);
            
            % accounting for mean pattern separately, for noisy data the
            % signal variance can become negative. It's a bummer. To
            % counteract this, we set signal sigma to be realmin, so as to
            % avoid weird selectivity values.
            if svar<0
                svar = realmin;
            end
            
            b       = tapply(b,{'tt'},{'beta','mean'});
            Gtt     = cov(b.beta'); % single-finger G
            
            % 1. calc tuning of actual voxel data
            sftBeta = pp1_imana('estSingleFingerTuning',b.beta);
            sftBeta = mean(sftBeta);
            
            % 2. calc expected tuning of voxels with ~N(0,G)
            [sftEV,sftDistEV] = pp1_imana('SFT:calcExpectedValueG',evar,svar,Gtt,numVox,numRun,numSim); % expected value of the null
            
            % 3. calc expected tuning of voxels with ~sparse tuning, but where signal-to-noise equals that voxel data
            [sftSp,sftDistSp]   = pp1_imana('SFT:calcExpectedValueSparse',evar,svar,1,numVox,numRun,numSim); % 1= tuned to one condition (perfectly sparse)
            %[sftSp2,sftDistSp2] = pp1_imana('SFT:calcExpectedValueSparse',evar,svar,2,numVox,numRun,numSim); % 2= tuned to two conditions
            
            % 4. do prob test for each subject on their gauss sft distribution
            pEV = sum(sftBeta<=sftDistEV)/length(sftDistEV);
            if isempty(pEV)
                pEV = realmin;
            end
            
            % 5. do prob test for each subject on their sparse sft
            % distributions
            pSp = sum(sftBeta<=sftDistSp)/length(sftDistSp);
            if isempty(pSp)
                pSp = realmin;
            end
%             pSp2 = sum(sftBeta<=sftDistSp2)/length(sftDistSp2);
%             if isempty(pSp2)
%                 pSp2 = realmin;
%             end
            
            % add to output structure
            d.passive   = v;
            d.sn        = v.*T.sn(ii);
            d.roi       = v.*T.roi(ii);
            d.numVoxF   = v.*numVoxSig;
            d.numVoxTot = v.*size(T.raw_beta{ii},2);
            d.avgF      = v.*mean(F(F>=Fcrit));
            d.glm       = v.*glm;
%             d.sft       = [sftBeta;sftEV;sftSp;sftSp2];
            d.sft       = [sftBeta;sftEV;sftSp];
            d.sftProb   = [0;pEV;pSp];
            d.isEV      = [0;1;2];
            D = addstruct(D,d);
            fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,svar,evar,svar/evar);
        end
        save(fullfile(regDir,sprintf('sft_glm%d_fthres',glm)),'-struct','D');
        varargout = {D};    
        
    case 'ROI_getChordActivity'
        % makes a nearly-universal plotting structure from ROI_stats output
        glm = 1;
        vararginoptions(varargin,{'glm'});
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        D = [];
        d = [];
        v = ones(size(T.avg_betas,2),1);
        numDigits = sum(chords,2);
        for i = 1:size(T.sn,1)
            d = [];
            d.avgBeta   = T.avg_betas(i,:)';
            d.chord     = T.avg_tt(i,:)';
            d.numDigits = arrayfun(@(x) numDigits(x),d.chord);
            d.roi       = v.*T.roi(i);
            d.sn        = v.*T.sn(i);
            D = addstruct(D,d);
        end
        varargout = {D};    
    case 'PLOT_activityPerNumDigits'
        % plot distance between chords (split by number of fingers
        % stimulated in chord) across M1 and S1 subregions.
        glm = 1;
        roi = 1:4;
        I = cpd2_imana('LIST'); % get subject info
        sn = unique(I.sn);
        vararginoptions(varargin,{'glm','roi','sn'});
        
        D  = cpd2_imana('ROI_getChordActivity','glm',glm);
        D  = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi));
        D  = getrow(D,D.numDigits<6);
        Dr = tapply(D,{'sn','roi','numDigits'},{'avgBeta','mean'});
        % plot
        warning off
        style.use('numDigits');
        plt.bar([Dr.roi Dr.numDigits],Dr.avgBeta,'split',Dr.numDigits);
        xtick = {};
        for r = roi
            xtick = {xtick{:} '','',regname{r},'',''};
        end
        plt.set('xticklabel',xtick,'xticklabelrotation',45);
        plt.legend('east',{'1 digit','2 digits','3 digits','4 digits','5 digits'});
        plt.labels('region','avg. raw beta');
        drawline(0,'dir','horz');
        warning on
        varargout = {Dr,D};    
    
    case 'GLM_contrast'                                                 % STEP 3.3   :  Make t-stat contrasts for specified GLM estimates.
        % Make t-stat contrasts for specified GLM estimates.
        % enter sn, glm #
        % models each chord and also error trials
        vararginoptions(varargin,{'sn'});
        % Go to subject's directory
        cd(fullfile(glmDir, subj_name{sn}));
        load SPM;
        SPM = rmfield(SPM,'xCon');
        T   = load('SPM_info.mat');
        %_____t contrast for each chord condition
        for d = 1:31
            con               = zeros(1,size(SPM.xX.X,2));
            con(:,T.chord==d) = 1;
            con               = con/sum(con);
            SPM.xCon(d) = spm_FcUtil('Set',sprintf('chord_%d',d), 'T', 'c',con',SPM.xX.xKXs);
        end

        %____do the constrasts
        SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save('SPM.mat','SPM');
    
    case 'fingerpics'                                                       % Makes jpegs of finger activity patterns on cortical surface M1/S1
        sn  = 1;
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        
        for s=sn
            for g=glm
                cpd2_imana('surf_mapFingers','sn',s,'glm',g)
                %cpd2_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[1,5]);
                %cpd2_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[2]);
                %cpd2_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[3]);
                %cpd2_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[4]);
            end
        end
    case 'surf_mapFingers'                                                % Map locations of finger patterns- run after glm estimation
        % map volume images to metric file and save them in individual surface folder
        sn  = 1;
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        hemisphere = 1;   % left hemi

        for c = 1:31 % contrast #
            fileList{c} = sprintf('spmT_00%02d.nii',c); % see case 'contrast' for contrast number index
        end;
        for s = sn
            for h = hemisphere
                caretSDir = fullfile(caretDir,[atlasA,subj_name{s}],hemName{h});
                %specname  = fullfile(caretSDir,[atlasA,subj_name{s} '.' hem{h}   '.spec']);
                white     = fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                pial      = fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                topo      = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                
                C1 = caret_load(white);
                C2 = caret_load(pial);
                
                images = {};
                for f=1:length(fileList)
                    images{f} = fullfile(glmDir,subj_name{s},fileList{f});
                end
                
                M  = caret_vol2surf_own(C1.data,C2.data,images,'topo',topo,'exclude_thres',0.75,'ignore_zeros',1);
                caret_save(fullfile(caretSDir,sprintf('%s_glm%d_hemi%d_finger.metric',subj_name{s},glm,h)),M);
            end
        end;   
    case 'surf_fingerpatterns'             % Make finger pattern jpegs
        %close all;
        sn  = 1;
        glm = 3;
        numDigits = [1,5];
        vararginoptions(varargin,{'sn','glm','numDigits'});
        
        h = 1; % left hemi
        groupDir = [gpCaretDir filesep hemName{h} ];
        cd(groupDir);
        switch(h)
            case 1 % left hemi
                coord  = 'lh.FLAT.coord';
                topo   = 'lh.CUT.topo';
                %data  = 'lh.surface_shape';  
%                 xlims=[-4 15]; % may need to adjust locations for pics
%                 ylims=[-10 9];
                xlims=[-4 18]; % may need to adjust locations for pics
                ylims=[-9 20];
            case 2 % right hemi
                coord  = 'rh.FLAT.coord';
                topo   = 'rh.CUT.topo';
                %data  = 'rh.surface_shape';
                xlims  = [-10 20];
                ylims  = [-15 30];
        end
        
        % load Central Sulcus border line (to plot as dashed line in pics)
        border = fullfile(caretDir,'fsaverage_sym',hemName{h},['BA_borders.border']);
        B      = caret_load(border);
        % set path to caret surface patterns
        data = fullfile(caretDir,['x' subj_name{sn}],hemName{h},sprintf('%s_glm%d_hemi%d_finger.metric',subj_name{sn},glm,h));

        % plot topographic image of surface reconstruction (w/out patterns)
        figure('Color',[1 1 1]); % make figure
        sshape = fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.surface_shape']);
        M      = caret_plotflatmap('col',2,'data',sshape,'border',B.Border,...
                 'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims,'bordersize',10);
        colormap('bone');
        %close gcf;
        
        % plot patterns for single finger stimulation
        figure('Color',[1 1 1]); % make figure 
        digits    = sum(chords,2);
        digitCols = find(ismember(digits,numDigits));
        j = 1;
        for i = digitCols'
            subplot(1,length(digitCols),j);
            [M,d]   = caret_plotflatmap('M',M,'col',i,'data',data,...
                        'border',B.Border,'bordersize',10,'topo',topo,'coord',coord,'bordercolor',{'k.'});
            colormap('parula');
            j = j+1;
        end;
        
        mm = 4; % force colour scaling on patterns
        % loop through both figures and all subplots to: force colour
        % scaling and label conditions
        for i = 1:length(digitCols)
            subplot(1,length(digitCols),i);
            caxis([-mm/2 mm]);   % scale color across plots
            set(gca,'XTick',[]); % remove X and Y axis ticks
            set(gca,'YTick',[]);
            %axis equal;
            box on
            ax = get(gca);
            ax.XAxis.LineWidth = 4;
            ax.YAxis.LineWidth = 4;
        end % for each subplot

    case '0' % ------------ PCM: pcm analyses. ----------------------------    
    case 'PCM_getData'
        % Get betas for roi from subjects in PCM-friendly format.
        % Betas do not have run means removed.
        sn  = 1;
        glm = 1;
        roi = 2; % only one roi supported
        %depth = [-Inf Inf];
        vararginoptions(varargin,{'sn','glm','roi','depth'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        % load betas
        B = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));
        B = getrow(B,B.roi==roi);
        % outputs
        Y = {};
        partVec = {};
        condVec = {};
        for i = 1:length(sn)
            % get subject data
            s = sn(i);
            b = getrow(B,B.sn==s);
            bb = [];
%             depthIdx = b.depth{1}>=depth(1) & b.depth{1}<=depth(2);
            bb.run   = cell2mat(b.part); % since runs defined per session, use partitions (sess*run split)
            bb.chord = cell2mat(b.tt);
            bb.betas = cell2mat(b.betaW);
%             bb.betas = bb.betas(:,depthIdx);
            bb.betas = bb.betas;
            bb = getrow(bb,ismember(bb.chord,1:31)); % restrict to passive patterns only
            % put subj data into pcm variables
            Y{i}         = bb.betas;
            partVec{i}   = bb.run;
            condVec{i}   = bb.chord;
            G_hat(:,:,i) = pcm_estGCrossval(Y{i},partVec{i},condVec{i});
        end
        varargout = {Y,partVec,condVec,G_hat};
    case 'PCM_fitModels_oneROI'
        % fits pcm models to data from one region
        sn     = [];
        glm    = [];
        roi    = [];
        depth  = [-Inf Inf]; % depth of voxels to include in model analysis
        plotit = 0; % plot fits
        saveit = 1; % save fits
        runEffect = 'random';
        vararginoptions(varargin,{'sn','glm','roi','plotit','saveit','depth'});
        
        % get data
        [Y,partVec,condVec,G_hat_all] = cpd2_imana('PCM_getData','sn',sn,'roi',roi,'glm',glm,'depth',depth);
        G_hat = mean(G_hat_all,3);
        
        % get model structure
        M = pp1_imana('PCM_defineModels',G_hat); 
        
        % choose proper way to deal with run effects
        runEffect = 'random';

        % fit all models
        [T,theta,G_pred]        = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'isCheckDeriv',0,'verbose',1); % perform derivative checks (to ensure I've done correct deriv implementation for nonlinear models)
        [Tcv,theta_cv,Gcv_pred] = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta,'fitScale',1,'verbose',1,'isCheckDeriv',0);
        
        % do group CV gainExp model fittings:
        [M,T,theta,G_pred,Tcv,Gcv_pred] = pp1_imana('PCM_cvGainExpFits',Y,partVec,condVec,G_hat(1:5,1:5),M,T,theta,G_pred,Tcv,Gcv_pred);
        
        % save fits?
        if saveit
           outfile = fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d_00',glm,roi));
           save(outfile,'M','T','theta','Tcv','theta_cv','G_hat_all','G_pred','Gcv_pred'); 
        end
        
        %keyboard
        varargout = {Tcv,T,M,theta_cv,G_pred,Y,partVec,condVec};
    case 'PCM_fitGroup'
        % Does group-level pcm fitting for multiple rois
        sn  = cpd2_imana('getSubjs');
        roi = [1:6]; % no thalamic coverage in 5/8 subjects
        glm = 1;
        vararginoptions(varargin,{'sn','roi','glm'});
        fprintf('PCM GROUP fitting.\nsubjs included in model fitting : ');
        for s=sn
            fprintf('%d.',s);
        end
        for r = roi
            fprintf('\nhemi: %s  |  roi: %s\n',hemLetter{regSide(r)},regname{r-(regSide(r)-1)*length(regname)});
            cpd2_imana('PCM_fitModels_oneROI','sn',sn,'roi',r,'glm',glm,'saveit',1);
        end
        
    case 'PCM_getFits'
        % Gets models fits across regions & arranges into plotting
        % structure.
        % Assumes null model is model 1 and noiseceiling is last model.
        glm = 1;
        roi = [];
        nNull = 1; % which model is null?
        nCeil = 5; % which model is noise ceiling ?
        vararginoptions(varargin,{'glm','roi','nNull','nCeil'});
        D   = []; % output structure
        for r = roi
            Tcv = [];
            % if exists, load pcm fits for region (otherwise, skip region)
            fprintf('\nroi %d...',r);
            try
                load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d',glm,r)));
            catch
                fprintf('no file.');
                continue
            end
            % scale likelihoods
            Tcv.likelihood_norm = bsxfun(@minus,Tcv.likelihood,Tcv.likelihood(:,nNull));
            % arrange into plotting structure
            numSubjs   = size(Tcv.SN,1);
            numModels  = numel(M);
            nameModels = {};
            Q = [];
            for m = 1:numModels
                % get model names
                nameModels{end+1,1} = M{m}.name;
                % get thetas (if model has any)
                if M{m}.numGparams==0
                    q.thetaCV = num2cell(nan(numSubjs,1),2);
                else
                    q.thetaCV = num2cell(theta_cv{m}',2);
                end
                q.model = ones(numSubjs,1).*m;
                q.sn    = [1:numSubjs]';
                Q = addstruct(Q,q);
            end
            v = ones(numModels,1);
            for j = 1:numSubjs
                d.sn  = v.*Tcv.SN(j);
                d.roi = v.*r;
                d.model = [1:numModels]';
                d.modelName  = nameModels;
                d.likeNormCV = Tcv.likelihood_norm(j,:)';
                d.likeCV     = Tcv.likelihood(j,:)';
                d.likeNorm   = v.*nan;
                d.likeNorm(nCeil) = T.likelihood(j,nCeil) - Tcv.likelihood(j,nNull); % upper noise ceiling
                d.thetaCV    = Q.thetaCV(Q.sn==j);
                D = addstruct(D,d);
            end
            fprintf('done.');
        end
        fprintf('\n');
        varargout = {D};
    case 'PCM_plotFits'
        type = varargin{1}; % bar, box, or r2
        glm = 1;
        roi = [1:6];
        sn  = 1:8;
        nNull = 1;
        nCeil = 8;
        modelsToPlot = [nNull, 4:nCeil];%[nNull:nCeil];
        
        
        inputs = {'glm',glm,'roi',roi,'sn',sn,'nNull',nNull,'nCeil',nCeil,'modelsToPlot',modelsToPlot};
        switch type
            case 'bar'
                cpd2_imana('PCM_plotFitsBar',inputs{:});
            case 'box'
                cpd2_imana('PCM_plotFitsBox',inputs{:});
            case 'r2'
                cpd2_imana('PCM_plotPseudoR2',inputs{:});
        end
        
    case 'PCM_plotPseudoR2'
        % plots pseudo R2 value (0 = null, 1 = upper noise ceiling)
        glm = [];
        roi = [];
        sn  = [];
        nNull = [];
        nCeil = [];
        modelsToPlot = [];
        vararginoptions(varargin,{'glm','roi','sn','nNull','nCeil','modelsToPlot'});
        % get pcm fits
        D = cpd2_imana('PCM_getFits','glm',glm,'roi',roi,'nNull',nNull,'nCeil',nCeil);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.sn,sn) & ismember(D.model,modelsToPlot));
        % calc pseudo R2
        upperNoiseCeil = kron(D.likeNorm(D.model==nCeil),ones(numel(unique(D.model)),1));
        D.pseudoR2 = D.likeNormCV./upperNoiseCeil;
        % check: [D.roi D.sn upperNoiseCeil]
        % % each subject should be scaled by a different upper noise ceiling model per roi
        
        % plot pcm fits
        numModels  = length(unique(D.model));
        nameModels = D.modelName(1:numModels);
        numPlots   = numel(roi);
        for i = 1:numPlots
            r = roi(i);
            subplot(1,numPlots,i);
            sty = style.custom(plt.helper.get_shades(numModels-2,'jet','descend'));
            plt.box(D.model,D.pseudoR2,'subset',D.model~=nNull & D.model~=nCeil & D.roi==r,'split',D.model,'style',sty);
            plt.set('xticklabel',{nameModels{2:numModels-1}},'xticklabelrotation',45);
            plt.labels('','pseudo model R2',regname{r});
            % plot noise ceilings
            drawline(median(D.pseudoR2(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-.'); % lower noise ceiling
            drawline(1,'dir','horz','linestyle','-');
            drawline(0,'dir','horz','linestyle','-','color',[0.7 0.7 0.7]);
            legend off
            ylim([0 1]);
        end
        plt.match('y');
        varargout = {D};
    case 'PCM_plotFitsBox'  
        % loads fit results per roi and plots them.
        glm = [];
        roi = [];
        sn  = [];
        nNull = [];
        nCeil = [];
        modelsToPlot = [];
        vararginoptions(varargin,{'glm','roi','sn','nNull','nCeil','modelsToPlot'});
        % get pcm fits
        D = cpd2_imana('PCM_getFits','glm',glm,'roi',roi,'nNull',nNull,'nCeil',nCeil);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.sn,sn) & ismember(D.model,modelsToPlot));
        % plot pcm fits
        numModels  = length(unique(D.model));
        nameModels = D.modelName(1:numModels);
        numPlots   = numel(roi);
        for i = 1:numPlots
            r = roi(i);
            subplot(1,numPlots,i);
            sty = style.custom(plt.helper.get_shades(numModels-2,'jet','descend'));
            plt.box(D.model,D.likeNormCV,'subset',D.model~=nNull & D.model~=nCeil & D.roi==r,'split',D.model,'style',sty);
            plt.set('xticklabel',{nameModels{2:numModels-1}},'xticklabelrotation',45);
            plt.labels('','relative log likelihood',regname{r});
            % plot noise ceilings
            drawline(mean(D.likeNormCV(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-.');
            drawline(mean(D.likeNorm(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-');
            legend off
            ylims = ylim;
            %ylim([0 ylims(2)]);
        end
        plt.match('y');
        varargout = {D};
    case 'PCM_plotFitsBar'  
        % loads fit results per roi and plots them.
        glm = [];
        roi = [];
        sn  = [];
        nNull = [];
        nCeil = [];
        modelsToPlot = [];
        vararginoptions(varargin,{'glm','roi','sn','nNull','nCeil','modelsToPlot'});
        numPlots = numel(roi);
        figure('Color',[1 1 1]);
        for i = 1:numPlots
            r = roi(i);
            load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d',glm,r)));
            T   = rmfield(T,{'reg'}); 
            T   = getrow(T,ismember(T.SN,sn));
            Tcv = rmfield(Tcv,{'reg'});
            Tcv = getrow(Tcv,ismember(Tcv.SN,sn));
            
            % % for plotting specified models
            Tp.SN           = T.SN;
            Tp.noise        = T.noise(:,modelsToPlot);
            Tp.scale        = T.scale(:,modelsToPlot);
            Tp.run          = T.run(:,modelsToPlot);
            Tp.likelihood   = T.likelihood(:,modelsToPlot);
                
            Tpcv.SN         = Tcv.SN;
            Tpcv.noise      = Tcv.noise(:,modelsToPlot);
            Tpcv.scale      = Tcv.scale(:,modelsToPlot);
            Tpcv.run        = Tcv.run(:,modelsToPlot);
            Tpcv.likelihood = Tcv.likelihood(:,modelsToPlot);
            
            Mp = {};
            for m = modelsToPlot
                Mp{end+1} = M{m};
            end
            % %
            Tpcv.likelihood_norm = bsxfun(@minus,Tpcv.likelihood,Tpcv.likelihood(:,1));
            % plot fits (errorbars are stderr)
            subplot(1,numPlots,i);
            pcm_plotModelLikelihood(Tpcv,Mp,'upperceil',Tp.likelihood(:,end),'style','bar','Nceil',length(modelsToPlot));
            set(gca,'xticklabelrotation',45);
            ylabel('relative log-likelihood')
            title(regname{r});
            box off;
        end
        plt.match('y');
        varargout = {Tpcv};
    
    case 'PCM_plotGpred'  
        % loads fit results per roi and plots them.
        glm   = 1;
        roi   = 2;
        nNull = 1;
        nCeil = 5;
        modelsToPlot = [nNull:nCeil];
        vararginoptions(varargin,{'glm','roi','sn'});
        if length(roi)>1
           error('can only call case with one roi'); 
        end
        % load fits
        load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d.mat',glm,roi)));
        % plot
        numPlots = numel(modelsToPlot);
        j = 1;
        for i = 1:numPlots
            % plot fits
            %subplot(2,ceil(numPlots/2),j);
            subplot(1,numPlots,j);
            imagesc(G_pred{i});
            %title(sprintf('%s : %d params',M{i}.name,M{i}.numGparams));
            title(sprintf('%s : %s',regname{roi},M{i}.name));
            axis square,
            colorbar
            % drawlines for chord distinctions
            drawline(5.5,'dir','horz');
            drawline(5.5,'dir','vert');
            drawline(15.5,'dir','horz');
            drawline(15.5,'dir','vert');
            drawline(25.5,'dir','horz');
            drawline(25.5,'dir','vert');
            drawline(30.5,'dir','horz');
            drawline(30.5,'dir','vert');
            j = j+1;
        end    
        varargout = {G_pred};
    case 'PCM_plotThetas'
        % loads fit results per roi and plots them.
        glm   = 3;
        roi   = [1:4,6];
        sn    = 1:5;
        model = 7;
        vararginoptions(varargin,{'glm','roi','sn'});
        % load plotting-friendly data structure
        D = cpd2_imana('PCM_getFits','glm',glm);
        D = getrow(D,D.model==model & ismember(D.roi,roi) & ismember(D.sn,sn));
        D.thetaCV = cell2mat(D.thetaCV); 
        % plot
        sty = style.custom(plt.helper.get_shades(length(roi),'parula','descend'));
        plt.trace([1:14],D.thetaCV(:,1:14),'split',D.roi,'style',sty);
        plt.labels('parameter no.','theta_cv value','single finger theta parameters');

    otherwise
        fprintf('%s: no such case.\n',what)
   
end



end