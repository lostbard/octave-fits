#include <iostream>
#include <sstream>
#include <octave/oct.h>

extern "C"
{
#include "fitsio.h"
}

static bool any_bad_argument( const octave_value_list& args );

DEFUN_DLD( save_fits_image_multi_ext, args, nargout,
"-*- texinfo -*-\n\
     @deftypefn {Loadable Function}  save_fits_image_multi_ext(@var{filename}, @var{image}, @var{bit_per_pixel})\n\
     Write @var{IMAGE} to FITS file @var{filename}.\n\n\
     Datacubes will be saved as multi-image extensions.\n\n\
     The optional parameter @var{bit_per_pixel} specifies the data type of the pixel values. Accepted string values are BYTE_IMG, SHORT_IMG, LONG_IMG, LONGLONG_IMG, FLOAT_IMG, and DOUBLE_IMG (default). Alternatively, corresponding numbers may be passed, i.e. 8, 16, 32, 64, -32, and -64.\n\n\
     Use a preceding exclamation mark (!) in the filename to overwirte an existing file.\n\n\
     Lossless file compression can be used by adding the suffix '.gz' to the filename.\n\n\
     @seealso{save_fits_image, read_fits_image}\
     Copyright (c) 2009-2010, Dirk Schmidt <fs@@dirk-schmidt.net>\
     @end deftypefn")
{
  if ( any_bad_argument(args) )
    return octave_value_list();

  octave_value fitsimage;
  std::string outfile = args(0).string_value ();

  int num_axis = 2;
  NDArray image = args(1).array_value();
  long int sz_axes[2]; 
  sz_axes[0] = image.dim1();
  sz_axes[1] = image.dim2();
  int num_images = image.dim3();

  //std::cerr << "num_images " <<  num_images << std::endl;

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
        fprintf( stderr, "Invalid string value for 'bit_per_pixel': %s\n", args(2).string_value().c_str() );
        return octave_value_list();
      }
    }
    else if( args(2).is_scalar_type() )
    {
      double val = args(2).double_value();
      if( (D_NINT( val ) ==  val) )
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
          fprintf( stderr, "Invalid numeric value for 'bit_per_pixel': %f\n", val );
          return octave_value_list();
        }
      }
    }
    else
    {
      fprintf( stderr, "Third parameter must be a valid string or a valid scalar value.\nSee 'help save_fits_image' for valid values.\n" );
      return octave_value_list();
    }
  }

  int status=0; // must be initialized with zero (I consider this to be a bug in libcfitsio).
                // status seems not to be set to zero after successful API calls

  // Open FITS file
  fitsfile *fp;
  if ( fits_create_file( &fp, outfile.c_str(), &status) > 0 )
  {
      fprintf( stderr, "Could not open file %s.\n", outfile.c_str() );
      fits_report_error( stderr, status );
      return fitsimage = -1;
  }

  long fpixel = 1;
  for( int i=0; i<num_images; i++ )
  {
    std::cerr << "image: " << i << std::endl;
    if( fits_create_img( fp, bitperpixel, num_axis, sz_axes, &status ) > 0 )
    {
      fprintf( stderr, "Could not create HDU.\n" );
      fits_report_error( stderr, status );
      return octave_value_list();
    }
    char keyname[9] = "XTENSION";
    char value[6] = "IMAGE";
    char comment[16] = "IMAGE extension";
    if( fits_write_key_str( fp, keyname, value, comment, &status ) > 0 )
    {
      fprintf( stderr, "Could not write XTENSION to HDU.\n" );
      fits_report_error( stderr, status );
      return octave_value_list();
    }
    if( fits_write_img( fp, TDOUBLE, fpixel, sz_axes[0]*sz_axes[1], image.fortran_vec() + i*sz_axes[0]*sz_axes[1], &status ) > 0 )
    {
      fprintf( stderr, "Could not write image data.\n" );
      fits_report_error( stderr, status );
      return octave_value_list();
    }
  }

  // Close FITS file
  if( fits_close_file(fp, &status) > 0 )
  {
      fprintf( stderr, "Could not close file %s.\n", outfile.c_str() );
      fits_report_error( stderr, status );
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
