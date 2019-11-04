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
## @deftypefn  {Function File} {} fitswrite (@var{image}, @var{filename}, [ @var{bits_per_pixel}, @var{header}] )
##
## Write @var{image} to FITS file @var{filename}.
##
## The first two dimensions of the image array get transposed, similar to Matlab's fitswrite().
## This function is not fully compatible with the parameters in Matlab's fitswrite().
##
## @seealso{save_fits_image}
## @end deftypefn

function fitswrite( img, filename, varargin )

  if (nargin < 2 || nargin > 4)
    print_usage ();
  elseif( ! ( isreal(img) && ! issparse(img) && ! isempty(img) ) )
    error( "fitswrite: IMAGE must be a real n-dimensional array" )
  elseif( ! ischar(filename) )
    error( "fitswrite: FILENAME must be a string" )
  end

  bitpix = [];
  header = [];
  for idx = 1:numel (varargin)
    opt = varargin{idx};
    if (ischar (opt) || isscalar (opt))
      bitpix = opt;
    elseif (iscell (opt) )
      header = opt;
    else
      error ("fitswrite: unrecognized option of class %s", class (opt))
    end
  end


  dims = 1:ndims(img);
  dims(1:2) = dims(2:-1:1);
  writebuf = permute( img, dims );
  save_fits_image( filename, writebuf, bitpix, header );
end

