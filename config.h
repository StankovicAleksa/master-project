#ifndef CONFIG_H
#define CONFIG_H

//#define BOOST_FLOAT

#ifdef BOOST_FLOAT
#include <boost/multiprecision/cpp_dec_float.hpp>
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<12>> Real;
#else
typedef double Real;
#endif
inline Real rmax(Real a, Real b) { return a>b? a:b;}
inline Real rmin(Real a, Real b) { return a<b? a:b;}
inline Real rabs(Real a) { return a>0? a:-a;}
//typedef double Real;
//number<cpp_dec_float<50> >
#endif
