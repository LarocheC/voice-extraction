function ltfatmex(varargin)
%LTFATMEX   Compile Mex/Oct interfaces
%   Usage:  ltfatmex;
%           ltfatmex(...);
%
%   LTFATMEX compiles the C backend in order to speed up the execution of
%   the toolbox. The C backend is linked to Matlab and Octave through Mex
%   and Octave C++ interfaces.
%
%   The action of LTFATMEX is determined by one of the following flags:
%
%     'compile'  Compile stuff. This is the default.
%
%     'clean'    Removes the compiled functions.
%
%     'test'     Run some small tests that verify that the compiled
%                functions work.
%
%   The target to work on is determined by on of the following flags:
%
%     'lib'      Perform action on the LTFAT C library.
%
%     'mex'      Perform action on the mex / oct interfaces.
%
%     'gpc'      Perform action on the GPC code for use with MULACLAB
%
%     'playrec'  Perform action on the playrec code for use with real-time
%                block streaming framework.
%
%     'auto'     Choose automatically which targets to work on based on
%                the operation system etc. This is the default.
%
%
%   Url: http://ltfat.sourceforge.net/doc/ltfatmex.php

% Copyright (C) 2005-2013 Peter L. Søndergaard <soender@users.sourceforge.net>.
% This file is part of LTFAT version 1.4.2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
%   AUTHOR : Peter L. Søndergaard.
%   TESTING: NA
%   REFERENCE: NA

% Verify that comp_pgauss is in path
if ~exist('comp_pgauss','file')
  disp(' ');
  disp('--- LTFAT - The Linear Time Frequency Analysis toolbox. ---');
  disp(' ')
  disp('To start the toolbox, call LTFATSTART as the first command.');
  disp(' ');
  return;
end;

bp=mfilename('fullpath');
bp=bp(1:end-length(mfilename));

definput.flags.target={'auto','lib','mex','gpc','playrec','java'};
definput.flags.command={'compile','clean','test'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Remember the current directory.
curdir=pwd;

do_lib  = flags.do_lib || flags.do_auto;
do_mex  = flags.do_mex || flags.do_auto;
do_gpc  = flags.do_gpc || (flags.do_auto && ~isoctave);
do_playrec  = flags.do_playrec && (~isoctave);
do_java  = flags.do_java;

if isoctave
  extname='oct';
  ext='oct';
else
  extname='mex';
  ext=mexext;
end;

if ispc
   makefilename='Makefile_mingw';
   make_exe = 'mingw32-make';
end;

if isunix
   makefilename='Makefile_unix';
   make_exe = 'make';
end;

if ismac
   makefilename='Makefile_mac';
   make_exe = 'make';
end;


portaudioLib = 'libportaudio';
if ispc
   sharedExt = 'dll';
else
   sharedExt = 'so';
end

fftw_lib_names = {'libfftw3-3', 'libfftw3f-3' };
clear mex;
% -------------- Handle cleaning --------------------------------
if flags.do_clean
    
  if ~isoctave
    % Delete files permanently (bypass trashcan on Windows/Mac
    % but remember the old state
    oldstate = recycle('off');
  end;
  
  if do_lib
    disp('========= Cleaning libltfat ===============');
    cd([bp,'src']);
    [status,result]=system([make_exe, ' -f ',makefilename,' clean']);
    %disp('Done.');    
  end;
  
  if do_mex
    fprintf('========= Cleaning %s interfaces =========\n', upper(extname));
    if isoctave
      deletefiles([bp,'oct'],'*.oct');
      deletefiles([bp,'oct'],'*.o');
    else
     cd([bp,'mex']);
     [status,result]=system([make_exe, ' -f ',makefilename,' clean',...
                     ' EXT=',mexext]); 
    end;
  end;
  
  if do_gpc
    disp('========= Cleaning GPC ====================');
    %deletefiles([bp,'thirdparty',filesep,'PolygonClip'],['PolygonClip.',mexext]);
    cd([bp,'thirdparty',filesep,'PolygonClip']);
    clear mex; 
    [status,result]=system([make_exe, ' -f ',makefilename,' clean',' EXT=',mexext]);
  end;
  
  if do_playrec
    disp('========= Cleaning PLAYREC ================');
    %deletefiles([bp,'thirdparty',filesep,'PolygonClip'],['PolygonClip.',mexext]);
    cd([bp,'thirdparty',filesep,'Playrec']);
    clear mex; 
    [status,result]=system([make_exe, ' -f ',makefilename,' clean',' EXT=',mexext]);
  end;
  
  if do_java
    disp('========= Cleaning JAVA ================');
    %deletefiles([bp,'thirdparty',filesep,'PolygonClip'],['PolygonClip.',mexext]);
    cd([bp,'blockproc',filesep,'java']);
    %clear java; 
    [status,result]=system([make_exe, ' -f ',makefilename,' clean']);
  end;
  
  if ~isoctave
    recycle(oldstate);
  end;
  
  cd(curdir);
end;

% -------------- Handle compiling  --------------------------------

if flags.do_compile
  if do_lib
    disp('========= Compiling libltfat ==============');
    cd([bp,'src']);
    if isoctave
        %[status,output]=system('mkoctfile -p CC');
      system([make_exe,' -f ',makefilename])
    else
        if ispc
            filesExist(cellfun(@(fEl) [bp,'mex',filesep,fEl,'.dll'], ...
                               fftw_lib_names,'UniformOutput',false));
        end;
      cd([bp,'src']);
      clear mex; 
      [status,result]=system([make_exe, ' -f ',makefilename,...
                  ' MATLABROOT=','"',matlabroot,'"',...
                  ' ARCH=',computer('arch')]);
      if(~status)
        disp('Done.');
      else
        error('Failed to build LTFAT libs:\n %s',result);
      end
    end;
  end;
  
  if do_mex
    fprintf('========= Compiling %s interfaces ========\n', upper(extname));
    if ispc
        filesExist(cellfun(@(fEl) ['mex',filesep,fEl,'.dll'], ...
                           fftw_lib_names,'UniformOutput',false));
    end;
    clear mex; 

    if isoctave
        cd([bp,'oct']);
        [status,result]=system([make_exe, ' -f ',makefilename])
    else
        cd([bp,'mex']);
        [status,result]=system([make_exe, ' -f ',makefilename,...
                            ' MATLABROOT=','"',matlabroot,'"',...
                            ' EXT=',mexext,...
                            ' ARCH=',computer('arch')]);
    end;
    if(~status)
      disp('Done.');
    else
      error('Failed to build %s interfaces: %s \n',upper(extname),result);
    end
  end;
  
  if do_gpc
    disp('========= Compiling GPC ===================');
    % Compile the PolygonClip interface to GPC for use with mulaclab
    cd([bp,'thirdparty',filesep,'PolygonClip']);
    clear mex; 
    [status,result]=system([make_exe, ' -f ',makefilename,...
                     ' MATLABROOT=','"',matlabroot,'"',...
                     ' EXT=',mexext,...
                     ' ARCH=',computer('arch')]);
    if(~status)
      disp('Done.');
    else
      error('Failed to build GPC:\n %s',result);
    end
  end;
  if do_playrec
    disp('========= Compiling PLAYREC ===============');
    % Compile the Playrec (interface to portaudio) for the real-time block-
    % stream processing
    % Check if portaudio library is present in the Matlab instalation
    binArchPath = [matlabroot,filesep,'bin',filesep,computer('arch')];
    playrecPath = [bp,'thirdparty',filesep,'Playrec'];
    foundPAuser = dir([playrecPath,filesep,'*portaudio*',sharedExt]);
    foundPAmatlab = dir([binArchPath,filesep,'*portaudio*',sharedExt]);

    if ~isempty(foundPAuser)
       if numel(foundPAuser)>2
          error('Ambiguous portaudio libraries in %s. Please leave just one.',playrecPath);
       end
       portaudioLib = foundPAuser(1).name;
       fprintf('Using %s.\n',portaudioLib);
    elseif ~isempty(foundPAmatlab)
       if numel(foundPAmatlab)>2
          error('Ambiguous portaudio libraries in %s.',binArchPath);
       end
       portaudioLib = foundPAmatlab(1).name;
       fprintf('Using %s.\n',portaudioLib);
    else
       error(['Portaudio not found. Please download Portaudio http://www.portaudio.com\n',...
              'and build it as a shared library and copy it to the\n',...
              '%s directory. \n'],playrecPath);
    end

    cd([bp,'thirdparty',filesep,'Playrec']);
    clear mex; 
    
    % Crop library filename
    if numel(portaudioLib)>3
      if strcmp(portaudioLib(1:3),'lib')
        portaudioLib = portaudioLib(4:end);
      end
      dotsInPath = strfind(portaudioLib,'.');
      if ~isempty(dotsInPath) 
        portaudioLib = portaudioLib(1:dotsInPath(1)-1);
      end
    end
    
    [status,result]=system([make_exe, ' -f ',makefilename,...
                     ' MATLABROOT=','"',matlabroot,'"',...
                     ' EXT=',mexext,...
                     ' PORTAUDIO=',portaudioLib,...
                     ' ARCH=',computer('arch')]);
    if(~status)
      disp('Done.');
    else
      error('Failed to build PLAYREC:\n %s',result);
    end
  end;
  
  if do_java
    disp('========= Compiling JAVA classes ===================');
    % Compile the JAVA classes
    cd([bp,'blockproc',filesep,'java']);
    clear mex; 
    [status,result]=system([make_exe, ' -f ',makefilename]);
    if(~status)
      disp('Done.');
    else
      error('Failed to build JAVA classes:\n %s',result);
    end
  end;
end;

% -------------- Handle testing ---------------------------------------

if flags.do_test
  
  if do_mex
    
    fprintf('========= Testing %s interfaces ==========\n', extname);
    fprintf('1.: Test if comp_pgauss.%s was compiled: ',ext);
    fname=['comp_pgauss.',ext];
    if exist(fname,'file')
      disp('SUCCESS.');
    else
      disp('FAILED.');
    end;
    
    fprintf('2.: Test if pgauss executes:              ');
    pgauss(100);
    % If the execution of the script makes it here, we know that pgauss
    % did not crash the system, so we can just print success. Same story
    % with the following entries.
    disp('SUCCESS.');

    fprintf('3.: Test if fftreal executes:             ');
    fftreal(randn(10,1),10);
    disp('SUCCESS.');

    fprintf('4.: Test if dgt executes:                 ');
    dgt(randn(12,1),randn(12,1),3,4);
    disp('SUCCESS.');

    
  end;
  
end;

% Jump back to the original directory.
cd(curdir);


function status = filesExist(filenames)
   if(~iscell(filenames))
      filenames={filenames};
   end
   for ii=1:length(filenames)
      filename = filenames{ii};
      if(~exist(filename,'file'))
         error('%s: File %s not found.',mfilename,filename);
      end
   end

function deletefiles(base,files)

L=dir([base,filesep,files]);
for ii=1:numel(L)
    s=[base,filesep,L(ii).name];
    delete(s);
end;


function status=compile_ltfat(bp)

  uses_lapack = {'comp_gabdual_long','comp_gabtight_long'};
  
% If we exit early, it is because of an error, so set status=1
status=1;

% Determine the name of the ltfat library
s=[bp,'lib',filesep,'libltfat-nomem.a'];

if ispc && ~isoctave
    if strcmp(mexext,'mexw64')
        s=[bp,'mex',filesep,'ltfat.dll'];
    else
        s=[bp,'lib',filesep,'libltfat-nomem.lib'];
    end;
end;

if ~exist(s,'file')
    disp(['The LTFAT C library cannot be found. Please compile ' ...
          'as described in the INSTALL file, or download ' ...
          'a binary release.']);
             
    return;
end;
    
if ispc && ~isoctave
    
    if ~exist([bp,'mex\libfftw3-3.dll'],'file')
        disp(['Please download the FFTW binary release for Windows as ' ...
               'described in the INSTALL file.']);
        return;
    end;
    
    % Create the .lib file for the Matlab FFTW3 lib.
    if ~exist([bp,'mex\libfftw3-3.lib'],'file')
        system(['"',matlabroot,'\sys\lcc\bin\lcc_implib.exe" -u "',...
                bp,'mex\LIBFFTW3-3.DLL"']);
    end;
        
    if ~exist([bp,'mex\libfftw3f-3.dll'],'file')
        disp(['Please download the FFTW binary release for Windows as ' ...
               'described in the INSTALL file.']);
        return;
    end;
        
    if ~exist([bp,'mex\libfftw3f-3.lib'],'file')
        system(['"',matlabroot,'\sys\lcc\bin\lcc_implib.exe" -u "',...
                bp,'mex\libfftw3f-3.dll"']);
    end;
        
end;
    
if isoctave
    cd([bp,'oct']);
    
    ext='oct';
    
    % Get the list of files.
    L=dir('*.cc');
    
    endchar=2;
else
    
    cd([bp,'mex']);
    
    ext=mexext;
    
    % Get the list of files.
    L=dir('comp_*.c');
    
    endchar=1;
end;


for ii=1:numel(L)
    filename = L(ii).name;
    objname  = [filename(1:end-endchar),ext];
    objdirinfo = dir(objname);
    
    % Make-like behaviour: build only the files where the src file is
    % newer than the object file, or the object file is missing.
    
    if ~isoctave && strcmp(mexext,'mexa64')
      
      % We don't know how to call LAPACK properly for this platform, so
      % if a mex-file uses LAPACK, skip its compilation
      if any(strcmp(uses_lapack,filename(1:end-endchar-1)))
        disp(['Skipping ',filename]);
        continue
      end;
    end;
    
    if isempty(objdirinfo) || (objdirinfo.datenum<L(ii).datenum)
        
        fprintf('Compiling %s\n',filename);
        
        if isoctave
            % Octave dynamically links to FFTW, Blas and Lapack, so they
            % are not included on the compilation line.
            mkoctfile('-c','oct-memalloc.c');
            mkoctfile('-I../thirdparty',...
                      '-I.','-I../src','-L../src','-L../lib',...
                      filename,'oct-memalloc.o',...
                      '-lltfat-nomem');            
        else
            mex('-c','mex-memalloc.c');
            if ispc
              if strcmp(mexext,'mexw64')
                mex('-I../thirdparty',...
                    '-I.','-I../src','-L../lib','-L../mex',...
                    ['-L',matlabroot,'\extern\lib\win64\microsoft'],...
                    filename,...
                    '-lltfat',...
                    '-lfftw3-3',...
                    '-lfftw3f-3');
                else  
                  mex('-I../thirdparty',...
                      '-I.','-I../src','-L../lib','-L../mex',...
                      ['-L',matlabroot,'\extern\lib\win32\lcc'],...
                      filename,'mex-memalloc.obj',...
                      '-lltfat-nomem',...
                      '-lfftw3-3',...
                      '-lfftw3f-3',...
                      '-lmwlapack',...
                      '-lmwblas');
                end;
                
            else
              
                mex('-I../thirdparty',...
                    '-I.','-I../src','-L../src','-L../lib',...
                    filename,'mex-memalloc.o',...
                    '-lltfat-nomem',...
                    '-lfftw3',...
                    '-lfftw3f',...
                    '-llapack',...
                    '-lblas');
                
            end;
            
        end;
        
        
    end;        
    
end;


status=0;


