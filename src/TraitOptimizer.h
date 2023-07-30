//
// File: JointTraitChromosomeLikelihood.h
// Authors:
//   Anat Shafir
// Created: 17/04/2023
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef CHROMEVOL_TRAITOPTIMIZER_H
#define CHROMEVOL_JTRAITOPTIMIZER_H

#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/AbstractNumericalDerivative.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>

#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include "LikelihoodUtils.h"
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/JointPhyloLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
// From Seqlib:
#include <vector>
#include <map>
#include <utility>
#include <string>
#include <omp.h>
#include <time.h>
namespace bpp
{
/**
 * @brief This class is used to optimize trait model based likelihood functions
 *
 *
 */ 

class TraitOptimizer

{
    public:

    void printLikelihoodParameters(JointPhyloLikelihood* lik, unsigned int optimized, vector<string> paramsNames) const;
    void printLikelihoodParameters(SingleProcessPhyloLikelihood* lik, unsigned int optimized, vector<string> paramsNames) const;

    template <typename T>
    typename std::enable_if<(std::is_same<T, SingleProcessPhyloLikelihood>::value) || (std::is_same<T, JointPhyloLikelihood>::value), void>::type
    optimizeTraitModel(T* lik, double tol, uint numOfIterations, std::vector<string> &traitParamNames, std::vector<string> &fixedParamsTrait){
        DerivableSecondOrder* f = lik;
        BrentOneDimension* optimizer = new BrentOneDimension(f);
        optimizer->setVerbose(1);
        optimizer->setProfiler(0);
        optimizer->setMessageHandler(0);
        optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        optimizer->setMaximumNumberOfEvaluations(100);
        // setting bracketing for Brent optimization
        optimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);
        ParameterList params;
        double currentLikelihood = lik->getValue();
        double prevLikelihood;
        std::cout << "Optimizing trait model ...." << std::endl;
        for (size_t i = 0; i < numOfIterations; i++){
            std::cout << "Iteration #" << i << std::endl;
            for (size_t j = 0; j < traitParamNames.size(); j++){
                prevLikelihood = currentLikelihood;
                params = lik->getParameters();
                const string nameOfParam = traitParamNames[j];
                std::cout << "Previous value of " << nameOfParam  << " is: "+ std::to_string(params.getParameter(nameOfParam).getValue()) << std::endl;
                if (std::count(fixedParamsTrait.begin(), fixedParamsTrait.end(), nameOfParam)){
                    continue;
                }
                Parameter param = params.getParameter(nameOfParam);
                std::cout << "Parameter name is: " << nameOfParam << std::endl;
                std::shared_ptr<IntervalConstraint> bounds = std::dynamic_pointer_cast<IntervalConstraint>(lik->getParameter(nameOfParam).getConstraint());
                 
                if ((i == 1) & (numOfIterations > 2)){
                    optimizer->getStopCondition()->setTolerance(tol* 2);
                }else{
                    optimizer->getStopCondition()->setTolerance(tol);
                }
                optimizer->setInitialInterval(bounds->getLowerBound()+0.00001, bounds->getUpperBound()-0.00001); 
           
                optimizer->init(params.createSubList(param.getName()));
                currentLikelihood = optimizer->optimize();
                std::cout << "parameter value after optimization "+ std::to_string(lik->getParameter(param.getName()).getValue()) << std::endl;
                std::cout << "***" << std::endl;

            }
            printLikelihoodParameters(lik, 1, traitParamNames);
            if (std::abs(prevLikelihood-currentLikelihood) < tol){
                break;
            }
        }
        delete optimizer;

    }




};
} // end of namespace bpp.
#endif // CHROMEVOL_TRAITOPTIMIZER_H
