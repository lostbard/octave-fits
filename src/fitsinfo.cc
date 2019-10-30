#include <iostream>
#include <sstream>
#include <ctype.h>
#include <octave/oct.h>
#include <octave/version.h>
#include <octave/file-info.h>

extern "C"
{
#include "fitsio.h"
}

static int
get_bin_format (const std::string &coltype, std::string &type, int &len)
{
  // numbers at start
  if(::isdigit (coltype[0]))
    len = atoi (coltype.c_str());
  else if (coltype [0] == 'P')
    len = 2;
  else
    len = 1;

  // skip the numbers to get to data type
  int idx = 0;
  for (idx=0; idx<coltype.length (); idx++)
    {
      if (! ::isdigit(coltype[idx]) && coltype[idx] != 'P') break;
    }

  std::string ftype = coltype.substr (idx);
  if (ftype == "") ftype = "A";

  switch (ftype[0])
    {
      case 'B':
        type = "uint8";
        break;
      case 'D':
        type = "double";
        break;
      case 'J':
        type = "int32";
        break;
      case 'I':
        type = "int16";
        if (coltype[0] == 'P')
          type = "int32";
        break;
      case 'E':
        type = "single";
        break;
      case 'X':
        if (len <= 8)
          type = "bit8";
        else if (len <= 16)
          type = "bit16";
        else if (len <= 32)
          type = "bit32";
        else
          type = "bit64";
        break;
      case 'C':
        type = "single complex";
        break;
      case 'M':
        type = "double complex";
        break;
      default:
        type = "char";
    }

  return ftype[0];
}


static int
get_ascii_format (const std::string &coltype, std::string &type, int &len)
{

  switch(coltype[0])
    {
      case 'I':
      case 'O':
      case 'Z':
        type = "integer";
        break;
      case 'F':
      case 'E':
      case 'G':
        type = "single";
        break;
      case 'D':
        type = "double";
        break;
      default:
        type = "char";
    }

  len = atoi (coltype.substr(1).c_str());

  return coltype[0];
}

DEFUN_DLD( fitsinfo, args, nargout,
"-*- texinfo -*-\n \
@deftypefn {Function File} {[@var{info}]} = fitsinfo(@var{filename})\n \
Read information about fits format file\n \
@end deftypefn")
{
  octave_value_list retval;  // create object to store return values

  if ( args.length() == 0)
    {
      print_usage ();
      return octave_value();
    }
  if ( args.length() != 1 || !args(0).is_string() )
    {
      error( "fitsinfo: filename (string) expected as only argument" );
      return octave_value();
    }

  std::string infile = args(0).string_value ();


  // Open FITS file and position to first HDU containing an image
  int status=0;
  fitsfile *fp;
  if ( fits_open_image( &fp, infile.c_str(), READONLY, &status) > 0 )
    {
      fits_report_error( stderr, status );
      error("Could not open file %s.", infile.c_str());
      return octave_value();
    }

  // map value that will be the info result
  octave_map om;

  // get file info
  octave::file_info stat(infile);

  om.assign("filename", octave_value(infile));
  om.assign("filesize", octave_value(stat.size()));
  om.assign("filemoddate", octave_value(stat.timestamp().ctime()));

  string_vector contents;
  int hdupos, hdutype;
  double offset;

  fits_get_hdu_num(fp, &hdupos);


  for(; !status; hdupos++)
    {
      octave_quit();
  
      octave_map hdu;

      fits_get_hdu_type(fp, &hdutype, &status);

      OFF_T headstart;
      OFF_T datastart;
      OFF_T dataend;
      fits_get_hduoff (fp, &headstart, &datastart, &dataend, &status);

      double datasize = (dataend - datastart);

      // read keywords (common for all types)
      int num_keys, key_pos;
      char card[FLEN_CARD]; 
      char keyname[FLEN_KEYWORD], value[FLEN_VALUE], comment[FLEN_VALUE];

      if( fits_get_hdrpos( fp, &num_keys, &key_pos, &status) > 0 )
        {
          fits_report_error( stderr, status );
          error( "Could not get number of header keywords" ) ;
          return octave_value();
        }

      Cell key_matrix(num_keys, 3);
      int key_count =  0;
      double slope = 1.0;
      double intercept = 0.0;
      double axis1 = 1;
  
      for( int i = 1; i <= num_keys; i++ ) // get the keywords
        { 
          if ( fits_read_record( fp, i, card, &status ) )
            {
              fits_report_error( stderr, status );
              error( "Could not read header keyword" );
              return octave_value ();
            }
          else
            {
              int len;
              char ktype;
              fits_get_keytype(card, &ktype, &status);
              fits_get_keyname(card, keyname, &len, &status);
              fits_parse_value(card, value, comment, &status);
       
              if (status == 0)
                {
                  if(strcmp(keyname, "BZERO") == 0) intercept = atof(value);
                  if(strcmp(keyname, "BSCALE") == 0) slope = atof(value);
                  if(strcmp(keyname, "NAXIS1") == 0) axis1 = atof(value);
                } 
              else if(status == VALUE_UNDEFINED)
                { 
                  status = 0;
                  keyname[0] = '\0';
                  value[0] = '\0';
                  comment[0] = '\0';
                }

              key_matrix(key_count, 0) = octave_value(std::string(keyname));
              key_matrix(key_count, 1) = octave_value(std::string(value));
              key_matrix(key_count, 2) = octave_value(std::string(comment));
              key_count ++;
            }

        }

      hdu.assign("intercept", octave_value(intercept));
      hdu.assign("slope", octave_value(slope));
      hdu.assign("keywords", octave_value(key_matrix));

      if(hdutype == IMAGE_HDU)
        {
          int bits_per_pixel, num_axis;
          int const MAXDIM=999; // max number supported by FITS standard

          std::vector<long> sz_axes(MAXDIM,0);
          if( fits_get_img_param( fp, sz_axes.size(), &bits_per_pixel,
                &num_axis, sz_axes.data(), &status) > 0 )
            {
              fits_report_error( stderr, status );
              error("Could not get image information");
            }

          std::string datatype = "";
          if(bits_per_pixel == BYTE_IMG) datatype = "uint8";
          if(bits_per_pixel == SHORT_IMG) datatype = "uint16";
          if(bits_per_pixel == LONG_IMG) datatype = "uint32";
          if(bits_per_pixel == LONGLONG_IMG) datatype = "uint64";
          if(bits_per_pixel == FLOAT_IMG) datatype = "single";
          if(bits_per_pixel == DOUBLE_IMG) datatype = "double";
          hdu.assign("datatype", octave_value(datatype));

          if(bits_per_pixel < 0)
            datasize = -bits_per_pixel;
          else
            datasize = bits_per_pixel;

          Matrix dv(1, num_axis);
          for ( int i=0; i<num_axis; i++ )
            {
              dv(0, i) = sz_axes[i];
              datasize *= sz_axes[i];
            }
          if (num_axis >= 2)
            {
              // matlab shows axis as Y, X, .... so need swap the first 2
              dv(0, 0) = sz_axes[1];
              dv(0, 1) = sz_axes[0];
            }

          datasize /= 8;

          hdu.assign("size", octave_value(dv));
          Matrix missing(0,0);
          hdu.assign("missingdatavalue", octave_value(missing));

        }
      else
      {
        // ascii or bin table data
        long nrows;
        int ncols;
        fits_get_num_rows(fp, &nrows, &status);
        fits_get_num_cols(fp, &ncols, &status);

        hdu.assign("rows", octave_value(nrows));
        hdu.assign("nfields", octave_value(ncols));


        if(hdutype == ASCII_TBL)
          {
            // ascii table method of getting size
            int linesize = 1;
            Matrix fieldwidth(1, ncols);
            Matrix fieldpos(1, ncols);
            Cell fieldformat(1, ncols);
            Cell fieldprecision(1, ncols);
            Cell fieldintercept(1, ncols);
            Cell fieldslope(1, ncols);
            Cell fieldmissing(1, ncols);
 
            char keyname[FLEN_KEYWORD], colname[FLEN_VALUE];
            char coltype[FLEN_VALUE];

            for (int i = 1; i <= ncols; i++)
              {
                double fc;
                fits_make_keyn ("TBCOL", i, keyname, &status);
                fits_read_key (fp, TDOUBLE, keyname, &fc, NULL, &status);
                if(status == KEY_NO_EXIST)
                  {
                    status = 0;
                    fc = 0;
                  }

                double tmp;
                fits_make_keyn ("TSCAL", i, keyname, &status);
                fits_read_key (fp, TDOUBLE, keyname, &tmp, NULL, &status);
                if (status == KEY_NO_EXIST)
                  {
                    status = 0;
                    tmp = 1;
                  }
                fieldslope(0,i-1) = octave_value(tmp);

                fits_make_keyn ("TZERO", i, keyname, &status);
                fits_read_key (fp, TDOUBLE, keyname, &tmp, NULL, &status);
                if (status == KEY_NO_EXIST)
                  {
                    status = 0;
                    tmp = 0;
                  }
                fieldintercept(0,i-1) = octave_value(tmp);

                fits_make_keyn ("TNULL", i, keyname, &status);
                fits_read_key (fp, TSTRING, keyname, &coltype, NULL, &status);
                if(status == KEY_NO_EXIST)
                  {
                    status = 0;
                    coltype[0] = '\0';
                    fieldmissing(0,i-1) = octave_value(Matrix());
                  }
                else
                  {
                    int j = 0;
                    // we have any chars
                    for(j=0;coltype[j]!='\0' && coltype[j] == ' ';j++) {}

                    if(j == strlen(coltype))
                      coltype[0] = '\0';

                    fieldmissing(0,i-1) = octave_value(coltype);
                  }

                fits_make_keyn("TFORM", i, keyname, &status);
                fits_read_key(fp, TSTRING, keyname, coltype, NULL, &status);

                std::string prec;
                int fw; 
                get_ascii_format(coltype, prec, fw);

                fieldwidth(0,i-1) = double(fw);
                fieldformat(0,i-1) = octave_value(coltype);
                fieldprecision(0,i-1) = octave_value(prec);
                fieldpos(0,i-1) = fc;
                linesize = fc + fw;

              }

            datasize = linesize * nrows;
            hdu.assign("rowsize", octave_value(linesize));
            hdu.assign("fieldwidth", octave_value(fieldwidth));
            hdu.assign("fieldprecision", octave_value(fieldprecision));
            hdu.assign("fieldformat", octave_value(fieldformat));
            hdu.assign("fieldpos", octave_value(fieldpos));
            hdu.assign("intercept", octave_value(fieldintercept));
            hdu.assign("slope", octave_value(fieldslope));
            hdu.assign("missingdatavalue", octave_value(fieldmissing));
          }
        else if(hdutype == BINARY_TBL)
          {
            Matrix fieldsize(1, ncols);
            Cell fieldformat(1, ncols);
            Cell fieldprecision(1, ncols);
            Cell fieldintercept(1, ncols);
            Cell fieldslope(1, ncols);
            Cell fieldmissing(1, ncols);

            char keyname[FLEN_KEYWORD], colname[FLEN_VALUE];
            char coltype[FLEN_VALUE];

            for (int i = 1; i <= ncols; i++)
              {

                double tmp;
                fits_make_keyn("TSCAL", i, keyname, &status);
                fits_read_key(fp, TDOUBLE, keyname, &tmp, NULL, &status);
                if (status == KEY_NO_EXIST)
                  {
                    status = 0;
                    tmp = 1.0;
                  }
                fieldslope(0,i-1) = octave_value(tmp);

                fits_make_keyn("TZERO", i, keyname, &status);
                fits_read_key(fp, TDOUBLE, keyname, &tmp, NULL, &status);
                if (status == KEY_NO_EXIST)
                  {
                    status = 0;
                    tmp = 0;
                  }
                fieldintercept(0,i-1) = octave_value(tmp);

                fits_make_keyn("TNULL", i, keyname, &status);
                fits_read_key(fp, TSTRING, keyname, &coltype, NULL, &status);
                if(status == KEY_NO_EXIST)
                  {
                    status = 0;
                    coltype[0] = '\0';
                    fieldmissing(0,i-1) = octave_value(Matrix());
                  }
                else
                  {
                    int j = 0;
                    // we have any chars
                    for (j=0;coltype[j]!='\0' && coltype[j] == ' ';j++) {}

                    if (j == strlen(coltype))
                      coltype[0] = '\0';

                    fieldmissing(0,i-1) = octave_value(coltype);
                  }

                fits_make_keyn("TFORM", i, keyname, &status);
                fits_read_key(fp, TSTRING, keyname, coltype, NULL, &status);

                std::string prec;
                int fw; 
                get_bin_format(coltype, prec, fw);

                fieldsize(0,i-1) = double(fw);
                fieldformat(0,i-1) = octave_value(coltype);

                fieldprecision(0,i-1) = octave_value(prec);
              }

            datasize = axis1 * nrows;
            hdu.assign("rowsize", octave_value(axis1));
            hdu.assign("fieldsize", octave_value(fieldsize));
            hdu.assign("fieldprecision", octave_value(fieldprecision));
            hdu.assign("fieldformat", octave_value(fieldformat));
            hdu.assign("intercept", octave_value(fieldintercept));
            hdu.assign("slope", octave_value(fieldslope));
            hdu.assign("missingdatavalue", octave_value(fieldmissing));
          }
        }  

      // missingdatavalue
      hdu.assign("datasize", octave_value(datasize));
      hdu.assign("offset", octave_value(double(datastart)));

      if (hdupos == 1)
        {
          om.assign("primarydata", octave_value(hdu));
          contents.append(std::string("primary"));
        }
      else
        {
          if(hdupos > 1 && hdutype == IMAGE_HDU)
            {
              char  value[FLEN_VALUE];
              fits_read_key(fp, TSTRING, "XTENSION", value, NULL, &status);
              if(status == KEY_NO_EXIST)
                status = 0;
              else if (strcmp(value, "IMAGE") != 0)
                hdutype = -1;
            }

          if(hdutype == IMAGE_HDU)
            {
              om.assign("image",  octave_value(hdu));
              contents.append(std::string("image"));
            }
          else if(hdutype == BINARY_TBL)
            {
              om.assign("binarytable",  octave_value(hdu));
              contents.append(std::string("binary table"));
            }
          else if(hdutype == ASCII_TBL)
            {
              om.assign("asciitable",  octave_value(hdu));
              contents.append(std::string("ascii table"));
            }
          else
            {
              om.assign("unknown",  octave_value(hdu));
              contents.append(std::string("unknown"));
            }
        }

      fits_movrel_hdu(fp, 1, NULL, &status);
    }

  if (status == END_OF_FILE) status = 0;

  // Close FITS file
  status = 0;
  if( fits_close_file(fp, &status) > 0 )
    {
      fits_report_error( stderr, status );
      error( "Could not close file %s.", infile.c_str() );
      return octave_value ();
    }

  Cell contents_val (1, contents.numel());
  for(int i=0; i<contents.numel(); i++)
    {
      contents_val(0,i) = octave_value(contents[i]);
    }
  om.assign("contents", octave_value(contents_val));

  retval(0) = om;

  return retval;
}

#if 0
%!shared testfile
%! testfile = urlwrite ( ...
%!   'https://fits.gsfc.nasa.gov/nrao_data/tests/pg93/tst0012.fits', ...
%!   tempname() );

%!fail("fitsinfo")

%!fail("fitsinfo(1)")

%!test
%! s=fitsinfo(testfile);
%! assert(s.filename, testfile);
%! assert(s.contents, {"primary", "binary table", "unknown", ...
%!   "image", "ascii table"});
%!
%! assert(s.primarydata.size, [109 102]);
%! assert(s.primarydata.offset, 2880);
%! assert(s.primarydata.datatype, 'single');
%! assert(s.primarydata.intercept, 0);
%! assert(s.primarydata.slope, 1);
%! assert(s.primarydata.datasize, 44472);
%!
%! assert(s.asciitable.offset, 103680);
%! assert(s.asciitable.rows, 53);
%! assert(s.asciitable.rowsize, 59);
%! assert(s.asciitable.nfields, 8);
%! assert(s.asciitable.datasize, 3127);

%!test
%! if exist (testfile, 'file')
%!   delete (testfile);
%! endif

#endif
