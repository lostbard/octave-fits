// Copyright (C) 2009-2015 Dirk Schmidt <fs@dirk-schmidt.net>
// Copyright (C) 2019 John Donoghue <john.donoghue@ieee.org>
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

extern "C"
{
#include "fitsio.h"
}

static bool any_bad_argument( const octave_value_list& args );

DEFUN_DLD( save_fits_image_multi_ext, args, nargout,
"-*- texinfo -*-\n\
     @deftypefn {Function File}  save_fits_image_multi_ext(@var{filename}, @var{image}, @var{bit_per_pixel})\n\
     Write @var{IMAGE} to FITS file @var{filename}.\n\n\
     Datacubes will be saved as multi-image extensions.\n\n\
     The optional parameter @var{bit_per_pixel} specifies the data type of the pixel values. Accepted string values are BYTE_IMG, SHORT_IMG, LONG_IMG, LONGLONG_IMG, FLOAT_IMG, and DOUBLE_IMG (default). Alternatively, corresponding numbers may be passed, i.e. 8, 16, 32, 64, -32, and -64.\n\n\
     Use a preceding exclamation mark (!) in the filename to overwrite an existing file.\n\n\
     Lossless file compression can be used by adding the suffix '.gz' to the filename.\n\n\
     @seealso{save_fits_image, read_fits_image}\n\
     @end deftypefn")
{
  bool verbose = false;

  if ( any_bad_argument(args) )
    return octave_value_list();

  octave_value fitsimage;
  std::string outfile = args(0).string_value ();

  int num_axis = 2;
  const NDArray image = args(1).array_value();
  long int sz_axes[2]; 
  sz_axes[0] = image.dim1();
  sz_axes[1] = image.dim2();

  int num_images;
  if (image.ndims() < 3)
    num_images = 1;
  else
    num_images = image.dim3();

  if(verbose)
    std::cerr << "num_images " <<  num_images << std::endl;

  int bitperpixel = DOUBLE_IMG;
  if( 3 == args.length() )
  {
    if( args(2).is_string() )
    {
      if( args(2).string_value() == "BYTE_IMG" )
        bitperpixel = BYTE_IMG;
      else if( args(2).string_value() == "SHORT_IMG" )
        bitperpixel = SHORT_IMG;
      else if( args(2).string_value() == "LONG_IMG" )
        bitperpixel = LONG_IMG;
      else if( args(2).string_value() == "LONGLONG_IMG" )
        bitperpixel = LONGLONG_IMG;
      else if( args(2).string_value() == "FLOAT_IMG" )
        bitperpixel = FLOAT_IMG;
      else if( args(2).string_value() == "DOUBLE_IMG" )
        bitperpixel = DOUBLE_IMG;
      else
      {
        error ("Invalid string value for 'bit_per_pixel': %s", args(2).string_value().c_str() );
        return octave_value_list();
      }
    }
    else if( args(2).is_scalar_type() )
    {
      double val = args(2).double_value();
      if( (OCTAVE__D_NINT( val ) ==  val) )
      {
        if( BYTE_IMG == val )
          bitperpixel = BYTE_IMG;
        else if( SHORT_IMG == val )
          bitperpixel = SHORT_IMG;
        else if( LONG_IMG == val )
          bitperpixel = LONG_IMG;
        else if( LONGLONG_IMG == val )
          bitperpixel = LONGLONG_IMG;
        else if( FLOAT_IMG == val )
          bitperpixel = FLOAT_IMG;
        else if( DOUBLE_IMG == val )
          bitperpixel = DOUBLE_IMG;
        else
        {
          error ("Invalid numeric value for 'bit_per_pixel': %f", val );
          return octave_value_list();
        }
      }
    }
    else
    {
      error ("Third parameter must be a valid string or a valid scalar value.\nSee 'help save_fits_image' for valid values." );
      return octave_value_list();
    }
  }

  int status=0; // must be initialized with zero (I consider this to be a bug in libcfitsio).
                // status seems not to be set to zero after successful API calls

  // Open FITS file
  fitsfile *fp;
  if ( fits_create_file( &fp, outfile.c_str(), &status) > 0 )
  {
      fits_report_error( stderr, status );
      error ("Could not open file %s.", outfile.c_str() );
      return fitsimage = -1;
  }

  long fpixel = 1;
  for( int i=0; i<num_images; i++ )
  {
    if(verbose)
      std::cerr << "image: " << i << std::endl;
    if( fits_create_img( fp, bitperpixel, num_axis, sz_axes, &status ) > 0 )
    {
      fits_report_error( stderr, status );
      error ("Could not create HDU." );
      return octave_value_list();
    }
    char keyname[9] = "XTENSION";
    char value[6] = "IMAGE";
    char comment[16] = "IMAGE extension";
    if( fits_write_key_str( fp, keyname, value, comment, &status ) > 0 )
    {
      fits_report_error( stderr, status );
      error ("Could not write XTENSION to HDU." );
      return octave_value_list();
    }
    double * datap = const_cast<double*>( image.fortran_vec() );
    if( fits_write_img( fp, TDOUBLE, fpixel, sz_axes[0]*sz_axes[1], datap + i*sz_axes[0]*sz_axes[1], &status ) > 0 )
    {
      fits_report_error( stderr, status );
      error ("Could not write image data." );
      return octave_value_list();
    }
  }

  // Close FITS file
  status = 0;
  if( fits_close_file(fp, &status) > 0 )
  {
      fits_report_error( stderr, status );
      error("Could not close file %s.", outfile.c_str() );
  }

  return octave_value_list();
}
static bool any_bad_argument( const octave_value_list& args )
{
  if ( args.length() < 2 || args.length() > 3 )
  {
    error( "save_fits_image_multi_ext: number of arguments - expecting save_fits_image_multi_ext( filename, image ) or save_fits_image_multi_ext( filename, image, bitsperpixel )" );
    return true;
  }

  if( !args(0).is_string() )
  {
    error( "save_fits_image_multi_ext: filename (string) expected for first argument" );
    return true;
  }


  return false;
}


#if 0
%!error <save_fits_image_multi_ext: number of arguments> save_fits_image_multi_ext()

%!error <save_fits_image_multi_ext: number of arguments> save_fits_image_multi_ext(1)

%!error <save_fits_image_multi_ext: filename> save_fits_image_multi_ext(1, 1)

%!test
%! testfile = tempname();
%! data = [ 1, 2, 3; 4, 5, 6; 1, 7, 11 ];
%! save_fits_image_multi_ext(testfile, data);
%! rd=read_fits_image(testfile);
%! assert(size(rd, 1), 3);
%! assert(size(rd, 2), 3);
%! assert(size(rd, 3), 1);
%! assert(data, rd)
%! if exist (testfile, 'file')
%!   delete (testfile);
%! endif

%!test
%! testfile = tempname();
%! data = [ 1, 2, 3; 4, 5, 6; 1, 7, 11 ];
%! save_fits_image_multi_ext(testfile, data, 8);
%! rd=read_fits_image(testfile);
%! assert(size(rd, 1), 3);
%! assert(size(rd, 2), 3);
%! assert(data, rd)
%! if exist (testfile, 'file')
%!   delete (testfile);
%! endif

!test
%! testfile = tempname();
%! data = [ 1, 2, 3; 4, 5, 6; 1, 7, 11 ];
%! save_fits_image_multi_ext(testfile, data, 16);
%! rd=read_fits_image(testfile);
%! assert(size(rd, 1), 3);
%! assert(size(rd, 2), 3);
%! assert(data, rd)
%! if exist (testfile, 'file')
%!   delete (testfile);
%! endif

%!test
%! testfile = tempname();
%! im0 = [ 1, 2, 3; 4, 5, 6; 1, 7, 11 ];
%! data = resize(im0, [size(im0), 4]);
%! save_fits_image_multi_ext(testfile, data);
%! rd=read_fits_image(testfile);
%! assert(size(rd, 1), 3);
%! assert(size(rd, 2), 3);
%! assert(im0, rd)
%! rd=read_fits_image(testfile, 3);
%! assert(zeros(3,3), rd)
%! if exist (testfile, 'file')
%!   delete (testfile);
%! endif
#endif
