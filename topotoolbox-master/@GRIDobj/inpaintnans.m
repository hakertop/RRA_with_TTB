function DEM = inpaintnans(DEM,varargin)

%INPAINTNANS Interpolate or fill missing values in a grid (GRIDobj)
%
% Syntax
%
%     DEMf = inpaintnans(DEM,type)
%     DEMf = inpaintnans(DEM,type,k)
%     DEMf = inpaintnans(DEM,type,k,conn)
%     DEMf = inpaintnans(DEM,DEM2)
%     DEMf = inpaintnans(DEM,DEM2,method)
%     DEMf = inpaintnans(DEM,DEM2,'tt','fit',1,'eps',20)
%
% Description
% 
%     inpaintnans fills gaps in a grid (GRIDobj) generated by measurement
%     errors or missing values. The user may choose between different
%     techniques to fill the gaps. Note that the algorithm fills only
%     pixels not connected to the DEM grid boundaries.
%
%     inpaintnans(DEM,type) or inpaintnans(DEM,type,k) fills missing values
%     in the DEM by interpolating inward from the pixels surrounding the
%     void. There are different ways to interpolate with 'laplace'
%     interpolation being the default (see also regionfill). 
%
%     inpaintnans(DEM,DEM2) or inpaintnans(DEM,DEM2,method) fills missing
%     values using interpolation from another grid. An example is that 
%     missing values in a SRTM DEM could be filled with values derived from
%     an ASTER GDEM. This technique is also referred to as the Fill and
%     Feather method (Grohman et al., 2006).
%
%     inpaintnans(DEM,DEM2,'tt') uses a hybrid method that uses both
%     laplacian interpolation and filling using a second DEM. The methods
%     accounts for potential vertical offsets between both DEMs by fitting
%     a regression surface to the boundary pixels which is used to adjust 
%     the second DEM. Missing values are then calculated by the weighted 
%     average of both techniques whereas higher weights are assigned to 
%     the laplacian technique if pixels are close to the boundaries of the
%     voids. Note that this technique will not fill voids connected to the
%     DEM boundaries. The technique is similar to the Delta Surface Fill
%     Method (DSF) described by Grohman et al. (2006).
%
%     inpaintnans(DEM,'interactive') starts an interactive tool to map a
%     region to be filled by laplacian interpolation.
%
% Input
%
%     DEM      digital elevation model with missing values
%               indicated by nans (GRIDobj)
%     type      fill algorithm 
%               'laplace' (default): laplace interpolation 
%                     as implemented in regionfill
%               'fill': elevate all values in each connected
%                     region of missing values to the minimum
%                     value of the surrounding pixels (same as 
%                     the function nibble in ArcGIS Spatial Analyst)
%               'nearest': nearest neighbor interpolation 
%                     using bwdist
%               'neighbors': this option does not close all nan-regions. It
%                     adds a one-pixel wide boundary to the valid values in
%                     the DEM and derives values for these pixels by a
%                     distance-weighted average from the valid neighbor
%                     pixels. This approach does not support the third input 
%                     argument k.
%               'interactive' (no further arguments): opens new figure and
%                     enables drawing a polygon which is going to be
%                     filled. The resulting DEM will be again displayed as
%                     hillshade.
%     k         if supplied, only connected components with 
%               less or equal number of k pixels are filled. Others
%               remain nan. Set to inf if all pixels enclosed by non-nan
%               pixels should be filled.
%     conn      Connectivity, specified as scalar 4 or 8
%     DEM2      if the second input argument is a GRIDobj, inpaintnans will
%               interpolate from DEM2 to locations of missing values in
%               DEM. 
%     method    interpolation method if second input argument is a GRIDobj.
%               {'linear'},'nearest','spline','pchip', 'cubic', or 'tt'.
%   
%     if method is 'tt' then several parameters can be applied
%
%     'eps'     parameter of the gaussian radial basis function (pixels). 
%               Larger values will give higher weight to the Laplacian
%               interpolation near the boundaries.
%     'fit'     0, 1, or 2. If 1, the function will resolve a potential 
%               elevation bias between DEM and DEM2 by raising/lowering
%               DEM2 to the elevations of the rim of nan-regions. If 2, the
%               function will fit a linear least squares model between the
%               two DEMs which will balance potential vertical offsets and
%               mismatches in inclination.
%
%     Note that unlike the other interpolation methods (as defined by the 
%     option 'method', 'tt' will not fill nan-regions connected to the DEM 
%     boundaries.
%
% Output
%
%     DEM      processed digital elevation model (GRIDobj)
%
% Example 1
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEM.Z(300:400,300:400) = nan;
%     subplot(1,2,1)
%     imageschs(DEM,[],'colorbar',false)
%     DEMn = inpaintnans(DEM);
%     subplot(1,2,2);
%     imageschs(DEMn,[],'colorbar',false)
%
% 
% See also: ROIFILL, FILLSINKS, BWDIST, STREAMobj/inpaintnans
%
% References: Grohman, Greg, George Kroenung, and John Strebeck. "Filling 
%             SRTM voids: The delta surface fill method." Photogrammetric 
%             Engineering and Remote Sensing 72.3 (2006): 213-216.
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 27. January, 2023

if nargin == 1
    DEM.Z = deminpaint(DEM.Z,varargin{:});
elseif ischar(varargin{1})
    if strcmpi(varargin{1},'neighbors')
        DEM = interpneighborpixels(DEM);
        return
    elseif strcmpi(varargin{1},'interactive')
        MASK = createmask(DEM,true);
        DEM.Z(MASK.Z) = nan;
        DEM = inpaintnans(DEM);
        DEM.Z(isnan(DEM.Z) & (~MASK.Z)) = nan;
        imageschs(DEM)
        return
    end
    DEM.Z = deminpaint(DEM.Z,varargin{:});
elseif isa(varargin{1},'GRIDobj')
    if nargin == 2
        method = 'linear';
    else
        method = varargin{2};
        method = validatestring(method,...
        {'linear','nearest','spline','pchip','cubic','tt'},'GRIDobj/inpaintnans','method',3);
    end
    
    switch method
        case 'tt'
            p   = inputParser;
            addRequired(p,'DEM2',@(x) isa(x,'GRIDobj'))
            addRequired(p,'method',@(x) strcmpi(x,'tt'))
            addParameter(p,'fit',1)
            addParameter(p,'eps',20);
            parse(p,varargin{:})
            DEM = ttinpaint(DEM,varargin{1},p.Results.eps,p.Results.fit);
        otherwise
            INAN = isnan(DEM);
            IX   = find(INAN.Z);
            [x,y] = ind2coord(DEM,IX);
            znew  = interp(varargin{1},x,y,method);
            DEM.Z(IX) = znew;
    end
   
end

end

function dem = deminpaint(dem,type,k,conn)
if nargin == 1
    type = 'laplace';
    k    = inf;
    conn = 8;
elseif nargin == 2
    k    = inf;
    conn = 8;
elseif nargin == 3
    conn = 8;
end
   
% error checking    
% clean boundary
I = isnan(dem);
% imclearborder
% I = imclearborder(I,conn);
marker = I;
marker(2:end-1,2:end-1) = false;
I = ~imreconstruct(marker,I) & I; 
clear marker


if ~isinf(k)
    I = xor(bwareaopen(I,k+1,conn),I);
end

% 
if numel(dem) < 10000^2 || ~strcmpi(type,'laplace')

% interpolation
switch lower(type)
    case 'nearest'
        % nearest neighbor interpolation
        [~,L] = bwdist(~I);
        dem = dem(L);
    case 'laplace'
        % -- use roifill (Code before 2015a)    
        % dem = roifill(dem,imdilate(I,ones(3)));
        % -- use regionfill (Code after and including 2015a)
        dem = regionfill(dem,I);
    case 'fill'
        % fill to lowest surrounding neighbor
        marker = inf(size(dem),class(dem));
        markerpixels = imdilate(I,ones(3)) & ~I;
        marker(markerpixels) = dem(markerpixels);
        mask = dem;
        mask(I | isnan(dem)) = -inf;
        marker = -marker;
        mask   = -mask;
        demrec = imreconstruct(marker,mask,conn);
        dem(I) = -demrec(I);
    otherwise
        error('type unknown')
end

else
    CC = bwconncomp(I,conn);
    STATS = regionprops(CC,'SubarrayIdx','Image');
   
    for r = 1:numel(STATS)
        rows = STATS(r).SubarrayIdx{1};
        rows = [min(rows)-1 max(rows)+1];
        cols = STATS(r).SubarrayIdx{2};
        cols = [min(cols)-1 max(cols)+1];
        
        demtemp = dem(rows(1):rows(2),cols(1):cols(2));
        inatemp = padarray(STATS(r).Image,[1 1],false);
        % -- Code before 2015a        
        % demtemp = roifill(demtemp,imdilate(inatemp,ones(3)));   
        % -- Code after and including 2015a
        demtemp = regionfill(demtemp,inatemp);  
        dem(rows(1):rows(2),cols(1):cols(2)) = demtemp;
    end
end
end



function DEM = ttinpaint(DEM,DEM2,shapeparam,fit)

INAN   = isnan(DEM);
INAN.Z = imclearborder(INAN.Z);

D      = bwdist(~INAN.Z,'euclidean');
D      = exp(-(D/shapeparam).^2);

dem    = DEM.Z;
[X,Y]  = getcoordinates(DEM,'matrix');

DEM2res = resample(DEM2,DEM);
demres  = DEM2res.Z;

CC = bwconncomp(INAN.Z,8);
STATS = regionprops(CC,'SubarrayIdx','Image');

for r = 1:numel(STATS)
    % Extract subimage
    rows = STATS(r).SubarrayIdx{1};
    rows = [min(rows)-1 max(rows)+1];
    cols = STATS(r).SubarrayIdx{2};
    cols = [min(cols)-1 max(cols)+1];

    z  = dem(rows(1):rows(2),cols(1):cols(2));
    z2 = demres(rows(1):rows(2),cols(1):cols(2));
    x  = X(rows(1):rows(2),cols(1):cols(2));
    y  = Y(rows(1):rows(2),cols(1):cols(2));
    d  = D(rows(1):rows(2),cols(1):cols(2));


    % Get subimage of the nan-area
    I       = padarray(STATS(r).Image,[1 1],false);

    % Predict in nan area using laplacian interpolation
    zlaplace = regionfill(z,I);

    % fit DEM2 to remove potential offsets
    if fit > 0
        % Get subimage with valid boundary
        B       = bwperim(imdilate(I,ones(3)));
        % number of boundary pixels
        n       = nnz(B);

        if fit == 1
            % Adjust the mean value
            b = ones(n,1)\(double(z(B)-z2(B)));
            z(I) = b + z2(I);
        elseif fit == 2
            % Fit a LS-surface to the boundary pixels
            b     = [ones(n,1) x(B) y(B)]\(double(z(B)-z2(B)));
            z(I)  = b(1) + b(2)*x(I) + b(3)*y(I) + z2(I);
        end

    else
        z(I) = z2(I);
    end
    % Calculate weighted average of the both
    z(I)    = d(I).*zlaplace(I) + (1-d(I)).*z(I);
    
    % Write back to DEM
    dem(rows(1):rows(2),cols(1):cols(2)) = cast(z,class(dem));
end
DEM.Z = dem;
end


function DEM = interpneighborpixels(DEM)

%INTERPNEIGHBORPIXELS Interpolate pixels from their neighbor pixels
%
% Syntax
%
%     DEMi = interpneighborpixels(DEM)
%
% Description
%
%     INTERPNEIGHBORPIXELS uses distance weighted averages to calculate
%     missing values (nans) for pixels that border pixels with valid
%     values. 
%
% Input arguments
%
%     DEM      GRIDobj
%
% Output arguments
%
%     DEMi     GRIDobj with
%


I = isnan(DEM.Z);
sq2 = sqrt(2);
w = 1./[sq2 1 sq2; ...
     1   0   1; ...
     sq2 1 sq2];
w = w(:);

Z = nlfilter(DEM.Z,[3 3],@fun);
DEM.Z = Z;

function b = fun(a)

a = a(:);    
if ~isnan(a(5))
    b = a(5);
    return
end

I = isnan(a);
if all(I)
    b = nan;
    return
end

I = ~I;
b = sum(a(I).*(w(I)./sum(w(I))));
end
end