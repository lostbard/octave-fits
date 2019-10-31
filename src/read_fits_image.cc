// Copyright (C) 2009-2015 Dirk Schmidt <fs@dirk-schmidt.net>
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <sstream>
#include <octave/oct.h>
#include <octave/version.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

extern "C"
{
#include "fitsio.h"
}

static bool any_bad_argument( const octave_value_list& args );

DEFUN_DLD( read_fits_image, args, nargout,
"-*- texinfo -*-\n\
@deftypefn {Function File} {[@var{image},@var{header}]} = read_fits_image(@var{filename},@var{hdu})\n\
Read FITS file @var{filename} and return image data in @var{image}, and the image header in @var{header}.\n\
\n\
size(@var{image}) will return NAXIS1 NAXIS2 ... NAXISN.\n\
\n\
@var{filename} can be concatenated with filters provided by libcfitsio. See:\
<http://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node81.html>\
\n\n\
Examples:\n\
\n\
1. If the file contains a single image, read_fits_image( \"filename\" ) will store the data into a 2d double matrix.\n\
\n\
2. If the file contains a data cube (continuous data, no extensions!), read_fits_image( \"filename\" ) will store the whole data cube into a 3d array.\n\
\n\
3. If the file contains a data cube, then read_fits_image( \"filename[*,*,2:5]\" ) will read the 2nd, 3rd, 4th, and 5th image, and store them into a 3d array.\n\
\n\
4. If the file contains multiple image extensions, then read_fits_image( \"filename[5]\" ) will read the 5th image. This is equivalent to read_fits_image( \"filename\", 5 ).\n\
\n\
NOTE: It's only possible to read one extension (HDU) at a time, i.e. multi-extension files need to be read in a loop.\n\
\n\
@seealso{save_fits_image, save_fits_image_multi_ext}\n \
@end deftypefn")
{
  if ( any_bad_argument(args) )
    return octave_value_list();

  octave_value fitsimage; // the octave container for the image data to be read by this function
  std::string infile = args(0).string_value ();

  if ( args.length()==2 )
  {
    std::ostringstream stream;
    stream << infile << "[" << int(args(1).scalar_value()) << "]";
    infile = stream.str();
  }

  int status=0; // must be initialized with zero (I consider this to be a bug in libcfitsio).
                // status seems not to be set to zero after successful API calls

  // Open FITS file and position to first HDU containing an image
  fitsfile *fp;
  if ( fits_open_image( &fp, infile.c_str(), READONLY, &status) > 0 )
  {
      fprintf( stderr, "Could not open file %s.\n", infile.c_str() );
      fits_report_error( stderr, status );
      return fitsimage = -1;
  }

  // Gather information about the image
  int bits_per_pixel, num_axis;
  int const MAXDIM=999; // max number supported by FITS standard
  std::vector<long> sz_axes(MAXDIM,0);
  if( fits_get_img_param( fp, sz_axes.size(), &bits_per_pixel, &num_axis, sz_axes.data(), &status) > 0 )
  {
      fprintf( stderr, "Could not get image information.\n" );
      fits_report_error( stderr, status );
      return fitsimage = -1 ;
  }
  if( 2 == num_axis )
    sz_axes[2] = 1;
  if( 1 == num_axis )
    sz_axes[1] = 1;

  #ifdef DEBUG
  std::cerr << bits_per_pixel << " " << num_axis << " ";
  for( int i=0; i<num_axis; i++ )
    std::cerr << sz_axes[i] << " ";
  std::cerr << std::endl;  
  #endif

  // Read image header
  int num_keys, key_pos;
  char card[FLEN_CARD];   // standard string lengths defined in fitsioc.h
  string_vector header;
  if( fits_get_hdrpos( fp, &num_keys, &key_pos, &status) > 0 ) // get number of keywords
  {
    fprintf( stderr, "Could not get number of header keywords\n" ) ;
    fits_report_error( stderr, status );
  }
  for( int i = 1; i <= num_keys; i++ ) // get the keywords
  {
    if ( fits_read_record( fp, i, card, &status ) )
    {
      fprintf( stderr, "Could not read header keyword\n" );
      fits_report_error( stderr, status );
    }
    else
    {
      header.append( std::string(card) );
    }
  }
  header.append( std::string("END\n") );  /* terminate listing with END */

  // Read image data and write it to an octave MArrayN type
  dim_vector dims(1,1);
  dims.resize( num_axis );
  int read_sz=sz_axes[0];
  for( int i=0; i<num_axis; i++ )
  {
    dims(i) = sz_axes[i];
    #ifdef DEBUG
      std::cerr << i << " " << sz_axes[i]  << std::endl;
    #endif
    if(i>0)
      read_sz *= sz_axes[i];
  }
  //std::cerr << "read_sz: " << read_sz << std::endl;

  #ifdef OCTAVE_API_VERSION_NUMBER
    #if OCTAVE_API_VERSION_NUMBER < 45
      MArrayN<double> image_data( dims ); // a octave double-type array
    #else
      MArray<double> image_data( dims ); // a octave double-type array
    #endif
  #else
    MArray<double> image_data( dims ); // a octave double-type array
  #endif

  int type = TDOUBLE; // convert read data to double (done by libcfitsio)
  std::vector<long> fpixel(num_axis,1); // start at first pixel in all axes

  int  anynul;
  if( fits_read_pix( fp, type, fpixel.data(), read_sz, NULL, image_data.fortran_vec(), 
                      &anynul, &status ) > 0 )
  {
       fprintf( stderr, "Could not read image.\n" );
       fits_report_error( stderr, status );
       return fitsimage = -1;
  }

  // Close FITS file
  if( fits_close_file(fp, &status) > 0 )
  {
      fprintf( stderr, "Could not close file %s.\n", infile.c_str() );
      fits_report_error( stderr, status );
  }

  octave_value_list retlist;
  retlist(0) =  image_data;
  retlist(1) =  header;

  return retlist;
}

static bool any_bad_argument( const octave_value_list& args )
{
  if ( args.length() < 1 || args.length() > 2 )
  {
    error( "read_fits_image: number of arguments - expecting read_fits_image( filename ) or read_fits_image( filename, extension )" );
    return true;
  }

  if( !args(0).is_string() )
  {
    error( "read_fits_image: filename (string) expected for first argument" );
    return true;
  }

  if( 2 == args.length() )
  {
    if( !args(1).is_scalar_type() )
    {
      error( "read_fits_image: second argument must be a non-negative scalar integer value" );
      return true;
    }
    double val = args(1).double_value();
    if( (OCTAVE__D_NINT( val ) !=  val) || (val < 0) )
    {
      error( "read_fits_image: second argument must be a non-negative scalar integer value" );
      return true;
    }

  }

  return false;
}

#if 0
%!shared testfile
%! # WFPC II 800 x 800 x 4 primary array data cube containing the 4 CCD images,
%! # plus a table extension containing world coordinate parameters.
%! # The sample file has been trimmed to 200 x 200 x 4 pixels to save disk space
%! testfile = urlwrite ( ...
%!   'https://fits.gsfc.nasa.gov/samples/WFPC2u5780205r_c0fx.fits', ...
%!   tempname() );

%!error <read_fits_image: number of arguments> read_fits_image()

%!error <read_fits_image: filename> read_fits_image(1)

%!test
%! rd=read_fits_image(testfile);
%! assert(!isempty(rd));
%! assert(numel(size(rd)), 3);
%! assert(size(rd, 1), 200);
%! assert(size(rd, 2), 200);
%! assert(size(rd, 3), 4);

%!test
%! rd=read_fits_image(sprintf("%s[0]", testfile));
%! assert(!isempty(rd));
%! assert(size(rd, 1), 200);
%! assert(size(rd, 2), 200);
%! assert(size(rd, 3), 4);

%! if exist (testfile, 'file')
%!   delete (testfile);
%! endif
#endif

