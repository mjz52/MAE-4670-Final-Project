% Michael Zakoworotny
% Produces a mesh for a composite specimen. Currently only supports Q4
% elements due to J integral computation concerns

function [nodalcoords,connectivity,BoundaryStruct,matList] = CompositeGrid_2D(nsd,xl,xr,yb,yt,wf,Vf,nnpe,nxF,nxM,ny)

% COMPOSITE NODAL COORDINATE GENERATION
% Define composite parameters
W = xr - xl;
w_fib = wf; % width of fiber
N_fib = floor(Vf*W/wf); % number of fibers
w_gap = W/(N_fib+1); % spacing between centerlines of each fiber

nx_per_fib = max(ceil(nxF/N_fib),1); % rounds up, with minimum of 1 (should never activate tho)
nx_per_gap = max(floor(nxM/(N_fib+1)),1); % MAY NEED TO UPDATE NUMBER OF ELEMS ON SIDE REGIONS

nelemx = nx_per_fib*N_fib + nx_per_gap*(N_fib+1);
nelemy = ny;

switch nnpe
    case 4
        nel = nelemx*nelemy;               % number of total elements
    otherwise
        error('This element type is not implemented\n');
end

nno = (nelemx+1)*(nelemy+1);            % number of total nodes
                              
nodalcoords = zeros(nno,nsd);           % initialize arrays
connectivity = zeros(nnpe,nel);

% collect the coordinates of the nodes
x_fib = linspace(-w_fib/2,w_fib/2,nx_per_fib+1); % coords relative to fiber
                                                 % centerlines
x_side = linspace(0,w_gap-w_fib/2,nx_per_gap+1); % x's from side
x_cent = linspace(0,w_gap-w_fib,nx_per_gap+1);

x = [x_side(1:end-1)]; % add left side
for i = 1:(N_fib-1)
    x = [x, x_fib(1:end-1)+w_gap*i]; % add fiber
    x = [x, x_cent(1:end-1)+w_gap*i+w_fib/2]; % add center spacing
end
x = [x, x_fib(1:end-1)+w_gap*N_fib]; % add last fiber
x = [x, x_side+w_gap*N_fib+w_fib/2]; % add right side
x = unique(x);

y = linspace(yb,yt,nelemy+1);          % equal distribution of the y nodes 

%%% vv BELOW HERE ALL THE SAME FROM BEFORE vv %%%
% Removed

end
