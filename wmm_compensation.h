#ifndef COMPONENTS_GEOMETRY_WMM_COMPENSATION_H_
#define COMPONENTS_GEOMETRY_WMM_COMPENSATION_H_

#include "vector.h"

void wmm_compensate(Vector vecA, Vector gravity, const double dec, const double dip);

void wmm_compensate_dec_only(Vector vecA, Vector gravity, const double dec);

void wmm_compensate_dip_only(Vector vecA, Vector gravity, const double dip);

#endif // COMPONENTS_GEOMETRY_WMM_COMPENSATION_H_