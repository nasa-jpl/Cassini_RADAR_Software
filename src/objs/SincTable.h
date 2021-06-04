#ifndef SINCTABLE_H
#define SINCTABLE_H

static const char rcs_id_sinctable_h[] =
  "@(#) $Id: SincTable.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include "Config.h"

class SincTable{
  public:
    SincTable();
    void computeTable();
    void config(Config& cfg);
    void config(Config& cfg, const string& filter_id);
    complex<float> interpolate(double x, double x_min, double x_delta,
      int index_min, int index_max, complex<float> samples[]);
    double interpolate(double x, double x_min, double x_delta,
      int index_min, int index_max, double samples[]);
    int maxIndex(double x, double x_min, double x_delta, int index_max);
    int minIndex(double x, double x_min, double x_delta, int index_min);
    ~SincTable();
    void write(const string& filename);
    bool configSet();
    int decimationFactor();
    int interpolatorLength();
    bool inBounds(double x, double x_min, double x_delta,
      int index_min, int index_max);

  private:
                              // Definitions below are from S. Hensley's
                              //   GeoSAR processor documentation
    double beta_;             // Beta parameter
    bool config_set_;
    int dec_factor_;          // Decimation factor d_f
    int enable_weighting_;    // 0 => no weighting of sinc function values
                              //   Otherwise, use weighting
    double pedestal_;         // w_f, pedestal height for weighting function
    int rel_length_;          // Relative filter length L_rel
    int intp_length_;         // Interpolator length L
    int table_len_;           // Total number of weighted sinc interpolator
                              //   values N_I
    double delay_;            // L/2
    double *table_;           // Table of weighted sinc interpolator values
                              //   used by processor
    double *table_scratch_;   // Table of weighted sinc interpolator values
                              //   prior to rearrangement; for diagnostics only
};
#endif
