% Create vtk-files (beam in different resolutions for REFERENCE solution in Abaqus

function createFEMInput(nxMax)
    for nx = 2:1:nxMax
        % Declaration of the Test Object: Beam:
        % x=0.0...0.03; y=0.0...0.06; z=0.0...0.2;

        xMin = 0.0;  xMax = 0.03;     % minimal and maximal x-value
        yMin = 0.0;  yMax = 0.06;     % minimal and maximal y-value
        zMin = 0.0;  zMax = 0.2;      % minimal and maximal z-value

        % Definition of the Cutting Plane (Ex.1):
        xTip = 0.015;                    % x-Position of Crack Tip
        yTip = 0.06;                     % not actually needed!  % y-Position of Crack 'Tip'
        zCrack = 0.1;                   % y-Position of Crack

        h = (xMax-xMin)/nx
        ny = ceil((yMax - yMin) / h)
        nz = ceil((zMax - zMin) / h)
        % if (mod(nz,2))
        %     nz = nz - 1
        % end

        outfileDim = [num2str(nx) 'x' num2str(ny) 'x' num2str(nz)];
        outfileTet4VTK = ['FEM_Tet4_6ts' outfileDim '.vtk'];

        % print file
        printVTKandTXTfile(outfileTet4VTK,xMax,xMin,yMax,yMin,zMin,zMax,nx,ny,nz,xTip,yTip,zCrack)
    end
end

function printVTKandTXTfile(outfileTet4VTK,xmax,xmin,ymax,ymin,zmin,zmax,nx,ny,nz,xTip,yTip,zCrack) % outfileCube8,

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate Mesh, vertices and connectivity - without cut

    % generate indices
    coords = zeros((nx+1)*(ny+1)*(nz+1),3);
    for k=1:(nz+1)
        for j=1:(ny+1)
            for i=1:(nx+1)
                coords( (k-1)*((nx+1)*(ny+1)) + (j-1)*(nx+1) + (i-1) + 1 , 1 ) = (xmax-xmin)/nx*(i-1);
                coords( (k-1)*((nx+1)*(ny+1)) + (j-1)*(nx+1) + (i-1) + 1 , 2 ) = (ymax-ymin)/ny*(j-1);
                coords( (k-1)*((nx+1)*(ny+1)) + (j-1)*(nx+1) + (i-1) + 1 , 3 ) = (zmax-zmin)/nz*(k-1);
            end
        end
    end
    % in case the object is not lying in the 1st quadrant (i.e. positive coords only), shift it accordingly:
    %coords(:,1) = coords(:,1) + xmin;
    %coords(:,2) = coords(:,2) + ymin;
    %coords(:,3) = coords(:,3) + zmin;

    % generate connections
    origconnect = zeros(nx*ny*nz,8); % for Cube8 elements
    for k=1:nz 
        for j=1:ny 
            for i=1:nx
                origconnect( i+(j-1)*nx+(k-1)*nx*ny , 1 ) = i + (j-1)*(nx+1) + (k-1)*(nx+1)*(ny+1);
                origconnect( i+(j-1)*nx+(k-1)*nx*ny , 2 ) = i+1 + (j-1)*(nx+1) + (k-1)*(nx+1)*(ny+1);
                origconnect( i+(j-1)*nx+(k-1)*nx*ny , 3 ) = i+1 + j*(nx+1) + (k-1)*(nx+1)*(ny+1);
                origconnect( i+(j-1)*nx+(k-1)*nx*ny , 4 ) = i + j*(nx+1) + (k-1)*(nx+1)*(ny+1);
                origconnect( i+(j-1)*nx+(k-1)*nx*ny , 5 ) = i + (j-1)*(nx+1) + k*(nx+1)*(ny+1);
                origconnect( i+(j-1)*nx+(k-1)*nx*ny , 6 ) = i+1 + (j-1)*(nx+1) + k*(nx+1)*(ny+1);
                origconnect( i+(j-1)*nx+(k-1)*nx*ny , 7 ) = i+1 + j*(nx+1) + k*(nx+1)*(ny+1);
                origconnect( i+(j-1)*nx+(k-1)*nx*ny , 8 ) = i + j*(nx+1) + k*(nx+1)*(ny+1);
            end
        end
    end


    % find intersected nodes:
    a=0;
    Pos = []
    for i=1:length(coords(:,1)) % for all nodes
        if (coords(i,1) < xTip && coords(i,3) == zCrack)
            a = a+1;
            Pos(a) = i;
        end
    end

    % duplicate these nodes, and append them to the "coords"-matrix:
    for j=1:length(Pos)
        coords( (nz+1)*(nx+1)*(ny+1) + j , : ) = coords( Pos(j) , : );
    end

    % update connections in "origconnect"-matrix:
    for i=1:length(Pos)
        Pos1 = origconnect(:,1)==Pos(i);
        origconnect(Pos1,1) = (nz+1)*(nx+1)*(ny+1)+i;
        
        Pos2 = find(origconnect(:,2)==Pos(i));
        if(isempty(Pos2)==0)
            origconnect(Pos2,2) = (nz+1)*(nx+1)*(ny+1)+i;
        end
        
        Pos3 = find(origconnect(:,3)==Pos(i));              % ?!?!?!
        if(isempty(Pos3)==0)
            origconnect(Pos3,3) = (nz+1)*(nx+1)*(ny+1)+i;
        end
        
        Pos4 = find(origconnect(:,4)==Pos(i));              % ?!?!?!
        if(isempty(Pos4)==0)
            origconnect(Pos4,4) = (nz+1)*(nx+1)*(ny+1)+i;
        end
    end


    nPoints = size(coords,1)
    nCells = size(origconnect,1);


    % now, we need to care about setting 6 tet4-elements into each cube:

    % coords: can remain the same!
    % origconnect: must be changed!

    % dim(origconnect) = (nCells,8);
    tetconnect = zeros(6*nCells,4);

    for oldindex = 1:nCells % number of (old) cells
        %for interimindex = 1:6 % number ot Tet4-elements in each Cube8-element
            tetconnect((1+6*(oldindex-1)),1) = origconnect(oldindex,1);
            tetconnect((1+6*(oldindex-1)),2) = origconnect(oldindex,2);
            tetconnect((1+6*(oldindex-1)),3) = origconnect(oldindex,4);
            tetconnect((1+6*(oldindex-1)),4) = origconnect(oldindex,5);
            
            tetconnect((2+6*(oldindex-1)),1) = origconnect(oldindex,2);
            tetconnect((2+6*(oldindex-1)),2) = origconnect(oldindex,3);
            tetconnect((2+6*(oldindex-1)),3) = origconnect(oldindex,4);
            tetconnect((2+6*(oldindex-1)),4) = origconnect(oldindex,5);
            
            tetconnect((3+6*(oldindex-1)),1) = origconnect(oldindex,4);
            tetconnect((3+6*(oldindex-1)),2) = origconnect(oldindex,3);
            tetconnect((3+6*(oldindex-1)),3) = origconnect(oldindex,8);
            tetconnect((3+6*(oldindex-1)),4) = origconnect(oldindex,5);
            
            tetconnect((4+6*(oldindex-1)),1) = origconnect(oldindex,2);
            tetconnect((4+6*(oldindex-1)),2) = origconnect(oldindex,5);
            tetconnect((4+6*(oldindex-1)),3) = origconnect(oldindex,6);
            tetconnect((4+6*(oldindex-1)),4) = origconnect(oldindex,3);
            
            tetconnect((5+6*(oldindex-1)),1) = origconnect(oldindex,3);
            tetconnect((5+6*(oldindex-1)),2) = origconnect(oldindex,7);
            tetconnect((5+6*(oldindex-1)),3) = origconnect(oldindex,8);
            tetconnect((5+6*(oldindex-1)),4) = origconnect(oldindex,6);
            
            tetconnect((6+6*(oldindex-1)),1) = origconnect(oldindex,6);
            tetconnect((6+6*(oldindex-1)),2) = origconnect(oldindex,5);
            tetconnect((6+6*(oldindex-1)),3) = origconnect(oldindex,8);
            tetconnect((6+6*(oldindex-1)),4) = origconnect(oldindex,3);
        %end
    end

    nTetCells = size(tetconnect,1); % = 6*nCells;

    % at this point, we are done with generating the Tet4 mesh (connections and points).


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate Tet4-VTK-File:

    P = coords';
    T = tetconnect';
    writeVTKMesh(outfileTet4VTK, P, T);

return;
end

function writeVTKMesh(filename, P, T)
%
% works for 3D (linear or quadratic) tetrahedra only, no acceptance of cubical elements yet!!
%
% The arguments 'P' and 'T' represent the coords-matrix and the connect-matrix.
%

% the points that are put into the output are going to have the following
% precision - in digits
precision = 10;

% open file
fid = fopen(filename,'w','b');
if(fid == -1)
   disp('Could not open file in writeVTKMesh')
end


% nCoords = size(P,1);
nPoints = size(P,2);
nCells = size(T,2);
nelnodes = size(T,1);

% Create HEADER for the vtk-file: -------------------------------------------------
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'vtk output\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

% Print POINTS Header: ------------------------------------------------------------
fprintf(fid, 'POINTS ');
fprintf(fid, num2str(nPoints));
fprintf(fid, ' float\n');

% Print POINTS as given in array P each into a newline of 'fid': ------------------
for pId = 1:nPoints
   for dof = 1:3
       fprintf(fid,num2str(P(dof,pId)));
       fprintf(fid,' ')
   end
   if(mod(pId,3) == 0)
    fprintf(fid,'\n');
   end
end

% Print CELLS Header: -------------------------------------------------------------
fprintf(fid, ['\nCELLS ' num2str(nCells) ' ' num2str((nelnodes+1)*nCells) '\n']);

% Print CELLS as given in arry T each into a newline of 'fid',
% including an identifying number which indicates the number of element nodes: ----
TData = [nelnodes*ones(1,size(T,2));T-1];

if(size(TData,1)==5)
fprintf(fid, '%d %d %d %d %d\n', TData);
elseif (size(TData,1)==9)
fprintf(fid, '%d %d %d %d %d %d %d %d %d\n', TData);
end

% Prepare for last Part of vtk-file: Cell-type identification: -------------------
if(nelnodes ==4)
   CELLTYPE = 10;
elseif (nelnodes ==10)
   CELLTYPE = 24;
elseif (nelnodes == 3)
   CELLTYPE = 5;
elseif (nelnodes == 8)
   CELLTYPE = 12;
   %     elseif (nelnodes == 20)
   %         CELLTYPE = ;
else
   return;
end

% Print CELL_TYPES as figured above: ---------------------------------------------
fprintf(fid, ['CELL_TYPES ' num2str(nCells) '\n']);

% for i=1:nCells
%     fprintf(fid, [num2str(CELLTYPE) '\n']);
% end
CELLTYPEData = CELLTYPE*ones(1,nCells);
fprintf(fid, '%d\n', CELLTYPEData);

fclose(fid);

return;
end
