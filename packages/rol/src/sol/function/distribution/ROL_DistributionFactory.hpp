// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_DISTRIBUTIONFACTORY_HPP
#define ROL_DISTRIBUTIONFACTORY_HPP

#include "Teuchos_ParameterList.hpp"

#include "ROL_Dirac.hpp"
#include "ROL_Gaussian.hpp"
#include "ROL_TruncatedGaussian.hpp"
#include "ROL_Uniform.hpp"
#include "ROL_Logistic.hpp"
#include "ROL_Triangle.hpp"
#include "ROL_Parabolic.hpp"
#include "ROL_RaisedCosine.hpp"
#include "ROL_Laplace.hpp"
#include "ROL_Cauchy.hpp"
#include "ROL_Smale.hpp"
#include "ROL_Arcsine.hpp"
#include "ROL_Kumaraswamy.hpp"
#include "ROL_Exponential.hpp"

namespace ROL {

  enum EDistribution {
    DISTRIBUTION_DIRAC = 0, 
    DISTRIBUTION_GAUSSIAN,
    DISTRIBUTION_TRUNCATEDGAUSSIAN,
    DISTRIBUTION_UNIFORM, 
    DISTRIBUTION_LOGISTIC,
    DISTRIBUTION_TRIANGLE,
    DISTRIBUTION_PARABOLIC,
    DISTRIBUTION_RAISEDCOSINE,
    DISTRIBUTION_LAPLACE,
    DISTRIBUTION_CAUCHY,
    DISTRIBUTION_SMALE,
    DISTRIBUTION_ARCSINE,
    DISTRIBUTION_KUMARASWAMY,
    DISTRIBUTION_EXPONENTIAL,
    DISTRIBUTION_LAST
  };

  inline std::string EDistributionToString(EDistribution ed) {
    std::string retString;
    switch(ed) {
      case DISTRIBUTION_DIRAC:             retString = "Dirac";                 break;
      case DISTRIBUTION_GAUSSIAN:          retString = "Gaussian";              break;
      case DISTRIBUTION_TRUNCATEDGAUSSIAN: retString = "Truncated Gaussian";    break;
      case DISTRIBUTION_UNIFORM:           retString = "Uniform";               break;
      case DISTRIBUTION_LOGISTIC:          retString = "Logistic";              break;
      case DISTRIBUTION_TRIANGLE:          retString = "Triangle";              break;
      case DISTRIBUTION_PARABOLIC:         retString = "Parabolic";             break;
      case DISTRIBUTION_RAISEDCOSINE:      retString = "Raised Cosine";         break;
      case DISTRIBUTION_LAPLACE:           retString = "Laplace";               break;
      case DISTRIBUTION_CAUCHY:            retString = "Cauchy";                break;
      case DISTRIBUTION_SMALE:             retString = "Smale";                 break;
      case DISTRIBUTION_ARCSINE:           retString = "Arcsine";               break;
      case DISTRIBUTION_KUMARASWAMY:       retString = "Kumaraswamy";           break;
      case DISTRIBUTION_EXPONENTIAL:       retString = "Exponential";           break;
      case DISTRIBUTION_LAST:              retString = "Last Type (Dummy)";     break;
      default:                             retString = "INVALID EDistribution"; break;
    }
    return retString;
  }

  inline int isValidDistribution(EDistribution ed) {
    return( (ed == DISTRIBUTION_DIRAC) ||
            (ed == DISTRIBUTION_GAUSSIAN) ||
            (ed == DISTRIBUTION_TRUNCATEDGAUSSIAN) ||
            (ed == DISTRIBUTION_UNIFORM) ||
            (ed == DISTRIBUTION_LOGISTIC) ||
            (ed == DISTRIBUTION_TRIANGLE) ||
            (ed == DISTRIBUTION_PARABOLIC) ||
            (ed == DISTRIBUTION_RAISEDCOSINE) ||
            (ed == DISTRIBUTION_LAPLACE) ||
            (ed == DISTRIBUTION_CAUCHY) ||
            (ed == DISTRIBUTION_SMALE) ||
            (ed == DISTRIBUTION_ARCSINE) ||
            (ed == DISTRIBUTION_KUMARASWAMY) || 
            (ed == DISTRIBUTION_EXPONENTIAL) );
  }

  inline EDistribution & operator++(EDistribution &type) {
    return type = static_cast<EDistribution>(type+1);
  }

  inline EDistribution operator++(EDistribution &type, int) {
    EDistribution oldval = type;
    ++type;
    return oldval;
  }

  inline EDistribution & operator--(EDistribution &type) {
    return type = static_cast<EDistribution>(type-1);
  }

  inline EDistribution operator--(EDistribution &type, int) {
    EDistribution oldval = type;
    --type;
    return oldval;
  }

  inline EDistribution StringToEDistribution(std::string s) {
    s = removeStringFormat(s);
    for ( EDistribution tr = DISTRIBUTION_DIRAC; tr < DISTRIBUTION_LAST; tr++ ) {
      if ( !s.compare(removeStringFormat(EDistributionToString(tr))) ) {
        return tr;
      }
    }
    return DISTRIBUTION_UNIFORM;
  }

  template<class Real>
  inline Teuchos::RCP<Distribution<Real> > DistributionFactory(Teuchos::ParameterList &parlist) {
    std::string dist = parlist.sublist("SOL").sublist("Distribution").get("Name","Dirac");
    EDistribution ed = StringToEDistribution(dist);
    switch(ed) {
      case DISTRIBUTION_DIRAC:             return Teuchos::rcp(new Dirac<Real>(parlist));
      case DISTRIBUTION_GAUSSIAN:          return Teuchos::rcp(new Gaussian<Real>(parlist));
      case DISTRIBUTION_TRUNCATEDGAUSSIAN: return Teuchos::rcp(new TruncatedGaussian<Real>(parlist));
      case DISTRIBUTION_UNIFORM:           return Teuchos::rcp(new Uniform<Real>(parlist));
      case DISTRIBUTION_LOGISTIC:          return Teuchos::rcp(new Logistic<Real>(parlist));
      case DISTRIBUTION_TRIANGLE:          return Teuchos::rcp(new Triangle<Real>(parlist));
      case DISTRIBUTION_PARABOLIC:         return Teuchos::rcp(new Parabolic<Real>(parlist));
      case DISTRIBUTION_RAISEDCOSINE:      return Teuchos::rcp(new RaisedCosine<Real>(parlist));
      case DISTRIBUTION_LAPLACE:           return Teuchos::rcp(new Laplace<Real>(parlist));
      case DISTRIBUTION_CAUCHY:            return Teuchos::rcp(new Cauchy<Real>(parlist));
      case DISTRIBUTION_SMALE:             return Teuchos::rcp(new Smale<Real>(parlist));
      case DISTRIBUTION_ARCSINE:           return Teuchos::rcp(new Arcsine<Real>(parlist));
      case DISTRIBUTION_KUMARASWAMY:       return Teuchos::rcp(new Kumaraswamy<Real>(parlist));
      case DISTRIBUTION_EXPONENTIAL:       return Teuchos::rcp(new Exponential<Real>(parlist));
      default:                             return Teuchos::null;
    }
  }
}
#endif
