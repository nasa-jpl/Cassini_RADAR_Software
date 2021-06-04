static const char rcs_id_sinctable_c[] =
  "@(#) $Id: SincTable.cpp,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <cmath>
#include <complex>
#include "Config.h"
#include "Constants.h"
#include "DebugInfo.h"
#include "Error.h"
#include "SimpleArray.h"
#include "SincTable.h"
#include "Utils.h"

using std::cerr;
using std::complex;
using std::cout;
using std::endl;

SincTable::SincTable() 
  : beta_(0), config_set_(false), dec_factor_(0), enable_weighting_(0),
    pedestal_(0), rel_length_(0), intp_length_(0), table_len_(0),
    table_(NULL), table_scratch_(NULL)
{
}

//----------------------------------------------------------------------
// computeTable()
//
// Compute the weighted sinc interpolator values for a single filter.
// The algorithm (and some of the variable names) are taken directly
// from S. Hensley's GeoSAR processor documentation and code.
//
// Note that the GeoSAR JurassicProk processor documentation uses the
// convention that sinc(x) = sin(x)/x, whereas the Utils::sinc method
// called here uses the convention that sinc(x) = sin(pi*x)/(pi*x),
// where x != 0 for both cases.  The documentation also contains an
// error in the equation for the rearrangement of the coefficients;
// the index to table_ is supposed to be i+j*intp_length_, not
// j*dec_factor_+i.  The GeoSAR code contains the correct equation.
//----------------------------------------------------------------------
void SincTable::computeTable()
{
  int i, j;

  double wgthgt = (1.0 - pedestal_) / 2.0;  // eta_w in GeoSAR documentation
  double offset = (table_len_ - 1.0) / 2.0; // k_0 in GeoSAR documentation
  double wa;                                // (k-k_0) in GeoSAR documentation
  double wgt;

  if (! config_set_)
  {
    ErrorMessage e("SincTable::computeTable:  config file has not been read");
    e.throwMe();
  }
  for (i = 0; i < table_len_; i++)
  {
    wa = i - offset;
    wgt = (1.0 - wgthgt) + wgthgt*cos(pi * wa / offset);
    table_scratch_[i] = sinc(wa * beta_ / dec_factor_);
    if (enable_weighting_)
    {
      table_scratch_[i] *= wgt;
    }
  }
  // Rearrange the coefficients in memory so that the values needed for
  // a particular interpolation are contiguous.
  double wbig=0;
  for (j = 0; j < dec_factor_; j++)
    {
      double wsum=0;
    for (i = 0; i < intp_length_; i++)
      {
      table_[i+j*intp_length_] = table_scratch_[j+i*dec_factor_];
      wsum=wsum+(table_scratch_[j+i*dec_factor_])*(table_scratch_[j+i*dec_factor_]);
      }
    for (i = 0; i < intp_length_; i++)
      {
	// table_[i+j*intp_length_] /= wsum*wsum;
      }
    cout << j << " " << wsum << endl;
    wbig+=wsum;
    }
  cout << wbig << endl;
}

//----------------------------------------------------------------------
// config()
//
// Read the configuration file parameters that characterize this sinc
// interpolation filter and compute the filter coefficient values.
//
// filter_id is a string used to uniquely identify the set of parameters
// (keywords) to be used for this filter.  The current set of recognized
// keywords is
//
// SINC_INTERP_FILTER_XXX_BETA                [required]
// SINC_INTERP_FILTER_XXX_DEC_FACTOR          [required]
// SINC_INTERP_FILTER_XXX_ENABLE_WEIGHTING    [required]
// SINC_INTERP_FILTER_XXX_PEDESTAL            [required]
// SINC_INTERP_FILTER_XXX_REL_LENGTH          [required]
//
// where XXX is the value of filter_id.
//----------------------------------------------------------------------

void SincTable::config(Config& cfg, const string& filter_id)
{
  string keyword_prefix = "SINC_INTERP_FILTER_" + filter_id + "_";

  beta_ = cfg.getDouble(keyword_prefix + "BETA");
  dec_factor_ = cfg.getInt(keyword_prefix + "DEC_FACTOR");
  enable_weighting_ = cfg.getInt(keyword_prefix + "ENABLE_WEIGHTING");
  pedestal_ = cfg.getDouble(keyword_prefix + "PEDESTAL");
  rel_length_ = cfg.getInt(keyword_prefix + "REL_LENGTH");

  // Perform rudimentary error-checking on the values read from the
  // config file.
  if (beta_ <= 0.0 || beta_ > 1.0)
  {
    ErrorMessage e("Illegal value " + toStr(beta_) + " for keyword " +
      keyword_prefix + "BETA:  value must be in range (0.0, 1.0]");
    e.throwMe();
  }
  if (dec_factor_ < 1)
  {
    ErrorMessage e("Illegal value " + toStr(dec_factor_) + " for keyword " +
      keyword_prefix + "DEC_FACTOR:  value must be >= 1");
    e.throwMe();
  }
  if (pedestal_ < 0.0 || pedestal_ > 1.0)
  {
    ErrorMessage e("Illegal value " + toStr(pedestal_) + " for keyword " +
      keyword_prefix + "PEDESTAL:  value must be in range [0.0, 1.0]");
    e.throwMe();
  }
  if (rel_length_ < 1)
  {
    ErrorMessage e("Illegal value " + toStr(rel_length_) + " for keyword " +
      keyword_prefix + "REL_LENGTH:  value must be >= 1");
    e.throwMe();
  }

  // Compute the tables of weighted sinc interpolator values.
  intp_length_ = round_double(rel_length_/beta_);
  delay_ = (double) intp_length_ / 2.0;
  table_len_ = intp_length_ * dec_factor_;
  table_scratch_ = (double *) make_array(sizeof(double), 1, table_len_);
  if (table_scratch_ == (double *) NULL)
  {
    ErrorMessage e("SincTable::config:  cannot allocate space for table_scratch_");
    e.throwMe();
  }
  table_ = (double *) make_array(sizeof(double), 1, table_len_);
  if (table_ == (double *) NULL)
  {
    ErrorMessage e("SincTable::config:  cannot allocate space for table_");
    e.throwMe();
  }
  config_set_ = true;

  computeTable(); // this need config_set_ = true !
}

SincTable::~SincTable()
{
  free_array((void *) table_, 1, table_len_);
  free_array((void *) table_scratch_, 1, table_len_);
}

//----------------------------------------------------------------------
// complex<float> interpolate()
//
// Given a function F sampled at values x_min + (k-1)*x_delta, calculate
// an approximation at an arbitrary value x by interpolating the samples
// weighted by the sinc interpolator values calculated for this filter.
//
// This method does not return an error value; it should only be called
// if a sufficient number of samples (i.e., intp_length_/2) exist around
// the value to be interpolated on both sides.
//
// Note that the documentation for the GeoSAR JurassicProk processor
// has some errors in the equations given in Section 2.5.  The code in
// this method has been adapted from the GeoSAR code itself, which is
// correct.  See comments below.
//----------------------------------------------------------------------

complex<float> SincTable::interpolate(double x, double x_min, double x_delta,
  int index_min, int index_max, complex<float> samples[])
{
  int i, j;

  double rho = (x - x_min)/x_delta + delay_;
  int k_0 = (int) rho;
  double r_frac = rho - k_0;
  // Do not use the equation for i_lag from the documentation.
  int i_lag = (int) (dec_factor_ * r_frac);
  int offset = i_lag * intp_length_;

  complex<float> sum = complex<float> (0.0, 0.0);

  if (k_0 >= index_min+intp_length_-1 && k_0 <= index_max)
  {
    for (i = 0; i < intp_length_; i++)
    {
      // The documentation uses (k_0 + i) as the sample increment, but
      // (k_0 - i) is correct when k_0 is defined as the sum of
      // (x - x_min)/x_delta and intp_length_/2 (as opposed to the
      // difference).
      j = k_0 - i;
      sum = sum + (float)table_[i + offset]*samples[j];
    }
  }
  return sum;
}

//----------------------------------------------------------------------
// double interpolate()
//
// Given a function F sampled at values x_min + (k-1)*x_delta, calculate
// an approximation at an arbitrary value x by interpolating the samples
// weighted by the sinc interpolator values calculated for this filter.
//
// This method does not return an error value; it should only be called
// if a sufficient number of samples (i.e., intp_length_/2) exist around
// the value to be interpolated on both sides.
//
// Note that the documentation for the GeoSAR JurassicProk processor
// has some errors in the equations given in Section 2.5.  The code in
// this method has been adapted from the GeoSAR code itself, which is
// correct.  See comments below.
//----------------------------------------------------------------------
double SincTable::interpolate(double x, double x_min, double x_delta,
  int index_min, int index_max, double samples[])
{
  int i, j;

  double rho = (x - x_min)/x_delta + delay_;
  int k_0 = (int) rho;
  double r_frac = rho - k_0;
  // Do not use the equation for i_lag from the documentation.
  int i_lag = (int) (dec_factor_ * r_frac);
  int offset = i_lag * intp_length_;

  double sum = 0;
  if (k_0 >= index_min+intp_length_-1 && k_0 <= index_max)
  {
    for (i = 0; i < intp_length_; i++)
    {
      // The documentation uses (k_0 + i) as the sample increment, but
      // (k_0 - i) is correct when k_0 is defined as the sum of
      // (x - x_min)/x_delta and intp_length_/2 (as opposed to the
      // difference).
      j = k_0 - i;
      sum += table_[i + offset]*samples[j];
    }
  }
  return sum;
}

//----------------------------------------------------------------------
// maxIndex()
//
// Return the highest index of the 1-D array whose values will be
// interpolated.
//----------------------------------------------------------------------
int SincTable::maxIndex(double x, double x_min, double x_delta, int index_max)
{
  double rho = (x - x_min)/x_delta + delay_;
  int k_0 = (int) rho;

  return MIN(k_0, index_max);
}

//----------------------------------------------------------------------
// minIndex()
//
// Return the lowest index of the 1-D array whose values will be
// interpolated.
//----------------------------------------------------------------------
int SincTable::minIndex(double x, double x_min, double x_delta, int index_min)
{
  int j;

  double rho = (x - x_min)/x_delta + delay_;
  int k_0 = (int) rho;

  j = k_0 - intp_length_ + 1;
  return MAX(j, index_min);
}

//----------------------------------------------------------------------
// write()
//
// Write the original sinc filter coefficients followed by the rearranged
// set of filter coefficients to the specified file.  (The latter set is
// the one actually used by the processor for interpolation.)  If the
// file exists, it will be overwritten.  The final output file size will
// be 2 x sizeof(double) x table_len_.
//----------------------------------------------------------------------
void SincTable::write(const string& filename)
{
  FileMgr f(filename, "w");
  int i;

  if (! config_set_)
  {
    ErrorMessage e("SincTable::write:  config file has not been read");
    e.throwMe();
  }
  for (i = 0; i < table_len_; i++)
  {
    f.write(table_scratch_[i]);
  }
  for (i = 0; i < table_len_; i++)
  {
    f.write(table_[i]);
  }
  f.close();
}

//----------------------------------------------------------------------
// configSet()
//
// Return the value of the flag that indicates whether the config file
// has been read.
//----------------------------------------------------------------------
bool SincTable::configSet()
{
  return(config_set_);
}

//----------------------------------------------------------------------
// decimationFactor()
//
// Return the decimation factor for this sinc interpolation filter.
//----------------------------------------------------------------------
int SincTable::decimationFactor()
{
  if (! config_set_)
  {
    ErrorMessage e("SincTable::decimationFactor:  config file has not been read");
    e.throwMe();
  }
  return(dec_factor_);
}

//----------------------------------------------------------------------
// interpolatorLength()
//
// Return the length of the sinc interpolation filter.
//----------------------------------------------------------------------
int SincTable::interpolatorLength()
{
  if (! config_set_)
  {
    ErrorMessage e("SincTable::interpolatorLength:  config file has not been read");
    e.throwMe();
  }
  return(intp_length_);
}

//----------------------------------------------------------------------
// bool inBounds()
//
// Check whether a sufficient number of samples around x exist to perform
// the interpolation.
//----------------------------------------------------------------------
bool SincTable::inBounds(double x, double x_min, double x_delta,
  int index_min, int index_max)
{
  double rho = (x - x_min)/x_delta + delay_;
  int k_0 = (int) rho;

  if (k_0 < index_min+intp_length_-1 || k_0 > index_max)
  {
    return false;
  }
  return true;
}

