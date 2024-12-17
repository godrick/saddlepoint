#ifndef SUBUNITARYMULTINOMIALCGF_H_INCLUDED
#define SUBUNITARYMULTINOMIALCGF_H_INCLUDED

#include "baseCGF.h"
#include "extractor.h"
#include "saddlepoint_types.h"
#include "CGF_Defaults.h"
#include "parametric_submodelCGF.h"
#include "multinomialCGF.h"

namespace saddlepoint {
namespace CGFs_via_templates {
class SubunitaryMultinomialCGF : public MultinomialCGF{
public:
    template <class t_vector_type, class N_type, class prob_vector_type>
    static auto K(const t_vector_type& tvec, const N_type& N, const prob_vector_type& prob_vec) {
      // The CGF for a modified Multinomial distribution.
      // It is useful in scenarios where a multinomial category are known to be 0.
      
      // For a random vector Y = (Y_1,...,Y_{d+1}),
      // prob_vec represents the probability vector for the first 'd' categories in a Multinomial distribution.
      // These probabilities can have a sum less than 1, indicating the presence of an additional zero-outcome category.
      
      // tvec is the vector of CGF arguments corresponding to each category, with the last category (d+1) assumed to be zero.
      // N is the parameter of the Multinomial distribution, representing the number of trials.
      
      // The function computes the CGF by adjusting the probability vector to account for the zero-outcome category.
      // The adjusted probabilities are calculated as p_i = prob_vec[i] / sum(prob_vec), for i = 1 to d.
      // K = N*log(sum_i prob_vec_i) + N*log(sum_i p_i e^{t_i})

      // zm1, the vector whose entries are zm1[i] = exp(tvec[i])-1
      auto zm1 = (tvec.array().exp() - 1);

      auto prob_sum = prob_vec.sum();
      auto adjusted_prob_vec = prob_vec.array() / (prob_sum);
      return N*log(prob_sum) + N*log1p((adjusted_prob_vec.array() * zm1.array()).sum() );
    }
};
}

namespace CGFs_with_AD {


class SubunitaryMultinomialExtractor {
public:
  template <class vector_type>
  auto operator()(const vector_type& parameter_vector) const {
    auto pvec = parameter_vector.tail(parameter_vector.size() - 1);
    if((pvec.array() < 0).any() || pvec.sum() > 1)
      throw std::invalid_argument("Each element of probability vector must be non-negative and the sum must not exceed 1. ");
    // if(pvec.sum() >= 1 || pvec.sum() <= 0) throw std::invalid_argument("The sum of prob_vec must be non-negative and strictly less than 1.");
    return std::make_pair(parameter_vector[0], pvec);
  }

};


using SubunitaryMultinomialCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::CGFwithExtractor<saddlepoint::CGFs_via_templates::SubunitaryMultinomialCGF, SubunitaryMultinomialExtractor>>;
// This CGF expects its arguments as a single vector, in the order specified by SubunitaryMultinomialExtractor

class SubunitaryMultinomialModelCGF : public SeparateParametersCGF<saddlepoint::CGFs_via_templates::SubunitaryMultinomialCGF, const ScalarAdaptor*, const VectorAdaptor*> {
private:
  using TemplateBaseCGF = saddlepoint::CGFs_via_templates::SubunitaryMultinomialCGF;
  using Base = SeparateParametersCGF<saddlepoint::CGFs_via_templates::SubunitaryMultinomialCGF, const ScalarAdaptor*, const VectorAdaptor*>;
public:
  SubunitaryMultinomialModelCGF(const ScalarAdaptor* n_adaptor, const VectorAdaptor* prob_vec_adaptor)
    : Base(TemplateBaseCGF(), n_adaptor, prob_vec_adaptor) {}

  SubunitaryMultinomialModelCGF() = delete; // otherwise adaptor pointers are uninitialised
};
} // namespace CGFs_with_AD

} // namespace saddlepoint

#endif // SUBUNITARYMULTINOMIALCGF_H_INCLUDED
