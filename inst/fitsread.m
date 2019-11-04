## Copyright (C) 2019 Dirk Schmidt <fs@dirk-schmidt.net>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License version 3 as published by the 
## Free Software Foundation.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {[@var{image}, @var{header}] =} fitsread (@var{filename},[@var{hdu} [, @var{convert}]])
##
## Read image data from FITS file @var{filename}.
##
## The first two dimensions of the image array get transposed similar to Matlab's fitsread().
## This function is not fully compatible with the parameters in Matlab's fitsread().
##
## @seealso{read_fits_image} for the parameters of this function.
## @end deftypefn

function [img,header] = fitsread( filename, varargin )

  if( nargin < 1 || nargin > 3 )
    print_usage ();
  elseif( ! ischar(filename) )
    error( "fitsread: FILENAME must be a string" )
  endif

  bitpix = [];
  hdu = [];
  for idx = 1:numel(varargin)
    opt = varargin{idx};
    if ( isscalar (opt) )
      hdu = opt;
    elseif( ischar(opt) )
      bitpix = opt;
    else
      error ("fitsread: unrecognized option of class %s", class (opt))
    end
  end

  if( ! isempty(hdu) && isempty(bitpix) )
    [img,header] = read_fits_image( filename, hdu );
  elseif( ! isempty(hdu) && ! isempty(bitpix) )
    [img,header] = read_fits_image( filename, hdu, bitpix );
  elseif( isempty(hdu) && ! isempty(bitpix) )
    [img,header] = read_fits_image( filename, bitpix );
  else
    [img,header] = read_fits_image( filename );
  end

  dims = 1:ndims(img);
  dims(1:2) = dims(2:-1:1);
  img = permute( img, dims );
end
