// 
//
// File: ChromosomeSubstitutionModel.cpp
// Created by: Anat Shafir
// Created on: 2020
//


#include "ChromosomeNumberDependencyFunction.h"

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;
/*****************************************************************************/
// Functions *** *** *** **** *** *** *** *** **** *** *** *** *** **** ***
/*****************************************************************************/
double ConstantDependencyFunction::getRate(std::vector<Parameter*> params, double state) const{
  return params[0]->getValue();
}
/**************************************************************************************/
double ConstantDependencyFunction::getRate(std::vector<double> params, double state) const{
  return params[0];

}

/**************************************************************************************/
double LinearDependencyFunction::getRate(std::vector<Parameter*> params, double state) const{
  double func_res = params[0]->getValue() + ((state-1)*params[1]->getValue());
  if ((func_res < 0) && (std::abs(func_res) > epsilon_neg_rate)){
    std::shared_ptr<IntervalConstraint> interval0 = dynamic_pointer_cast<IntervalConstraint>(params[0]->getConstraint());
    std::shared_ptr<IntervalConstraint> interval1 = dynamic_pointer_cast<IntervalConstraint>(params[1]->getConstraint());
    std::cout << "Intervals of param0 are " << interval0->getLowerBound() << " " << interval0->getUpperBound() << std::endl;
    std::cout << "Intervals of param1 are " << interval1->getLowerBound() << " " << interval1->getUpperBound() << std::endl;
    std::cout << "Current state is " << state << std::endl;
    std::cout << "function is in the form of p0 + p1*(s-1)" << std::endl;
    std::cout << "p0 is " << params[0]->getValue() << std::endl;
    std::cout << "p1 is " << params[1]->getValue() << std::endl;
    std::cout << "min domain is " << domainMin_<< std::endl;
    std::cout << "max domain is " << domainMax_ << std::endl;
    string paramName  = params[0]->getName();
    //return 0;

    throw Exception("LinearDependencyFunction::getRate(): Negative rate!!! Parameters are " + paramName + "\n");
  }else if (std::abs(func_res) < epsilon_neg_rate){
    return 0;
  }
  return func_res;

}
/**************************************************************************************/
double LinearDependencyFunction::getRate(std::vector<double> params, double state) const{
  double func_res = params[0] + ((state-1)*params[1]);
  if ((func_res < 0) && (std::abs(func_res) > epsilon_neg_rate)){
    throw Exception("LinearDependencyFunction::getRate(): Negative rate!!!");
  }else if (std::abs(func_res) < epsilon_neg_rate){
    return 0;
  }
  return func_res;

}
/**************************************************************************************/
double LinearDependencyFunction::getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, double minChrNum, double maxChrNum){
  if (index == 0){
    return parsimonyBound;
  }else if (index == 1){
    if (continuous_){
      return params[0]-(params[0]*(domainMax_+domainMax_)/2) + parsimonyBound;

    }else{
      return params[0]-(params[0]*(maxChrNum+minChrNum)/2) + parsimonyBound;

    }
    
  }
  throw Exception("LinearDependencyFunction::getParsimonyBound(): No such index!");
}
/**************************************************************************************/
void LinearDependencyFunction::updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, double maxChrNum){
  if (index == 0){
    if (continuous_){
      *lowerBound = std::max(-params.getParameter(paramsNames[1]).getValue()*(domainMax_), -params.getParameter(paramsNames[1]).getValue()*(domainMin_));

    }else{
      *lowerBound = std::max(lowerBoundOfRateParam, -params.getParameter(paramsNames[1]).getValue()*(maxChrNum-1));
      
    }
    *upperBound = upperBoundOfRateParam;

    
  }else if (index == 1){
    if (continuous_){
      if (domainMax_ < 0){
        *lowerBound = -upperBoundLinearRateParam;
        *upperBound = std::min(-params.getParameter(paramsNames[0]).getValue() / domainMax_, -params.getParameter(paramsNames[0]).getValue() / domainMin_);

      }else{
        if (domainMin_ < 0){ //signs of x_min and x_max are not the same
          if (params.getParameter(paramsNames[0]).getValue() < 0){
            throw Exception("LinearDependencyFunction::updateBounds(): if p0 is negative and min and max have not the same sign, it means that the rate is negative in some cases!");
          }else{
            *lowerBound = std::min(-params.getParameter(paramsNames[0]).getValue() / domainMax_, -params.getParameter(paramsNames[0]).getValue() / domainMin_);
          }
          *upperBound = std::max(-params.getParameter(paramsNames[0]).getValue() / domainMax_, -params.getParameter(paramsNames[0]).getValue() / domainMin_);

        }else{
          *upperBound = upperBoundLinearRateParam;
          *lowerBound = std::max(-params.getParameter(paramsNames[0]).getValue() / domainMax_, -params.getParameter(paramsNames[0]).getValue() / domainMin_);
        }
      }

    }else{
      *lowerBound = -params.getParameter(paramsNames[0]).getValue()/(maxChrNum-1);
      *upperBound = upperBoundLinearRateParam; 

    }
      
  }else{
    throw Exception("LinearDependencyFunction::updateBounds(): index out of bounds!!");
    
  }
  std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(paramsNames[index]).getConstraint());
  interval->setLowerBound(*lowerBound, interval->strictLowerBound());
}
/**************************************************************************************/
void LinearDependencyFunction::updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound){
  std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>((&(f->getParameter(paramName)))->getConstraint());
  interval->setLowerBound(lowerBound, interval->strictLowerBound());

}
/**************************************************************************************/
void LinearDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber){
  if (index == 0){
    //if (continuous_){
      //*lowerBound = std::max(std::max(-upperBoundLinearRateParam*(domainMax_), -upperBoundLinearRateParam*(domainMin_)), lowerBoundOfRateParam);
      // if (domainMin_ > 0){
      //   *lowerBound = std::min(std::min(-upperBoundLinearRateParam*(domainMax_), -upperBoundLinearRateParam*(domainMin_)), lowerBoundOfRateParam);

      // }else{
        //*lowerBound = lowerBoundOfRateParam;
      //}
      
      
    //}else{
    *lowerBound = lowerBoundOfRateParam;
    //}
    *upperBound = upperBoundOfRateParam;

    
  }else if (index == 1){
    if (continuous_){
      if (domainMax_ < 0){
        *lowerBound = -upperBoundLinearRateParam;
        *upperBound = std::min(-paramValues[0] / domainMax_, -paramValues[0] / domainMin_);

      }else{
        if (domainMin_ < 0){ //signs of x_min and x_max are not the same
          if (paramValues[0] < 0){
            throw Exception("LinearDependencyFunction::updateBounds(): if p0 is negative and min and max have not the same sign, it means that the rate is negative in some cases!");
          }else{
            *lowerBound = std::min(-paramValues[0] / domainMax_, -paramValues[0] / domainMin_);
          }
          *upperBound = std::max(-paramValues[0] / domainMax_, -paramValues[0] / domainMin_);

        }else{
          *upperBound = upperBoundLinearRateParam;
          *lowerBound = std::max(-paramValues[0] / domainMax_, -paramValues[0] / domainMin_);
        }
      }

    }else{
      *lowerBound = -paramValues[0]/(maxChrNumber-1);
      *upperBound = upperBoundLinearRateParam;   

    }
    
  }else{
    throw Exception("LinearDependencyFunction::updateBounds(): index out of bounds!!");
    
  }
}
/**************************************************************************************/
void LinearDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber){
  if (index == 0){
    if (continuous_){
      //*lowerBound = std::max(-upperBoundLinearRateParam*(domainMax_), -upperBoundLinearRateParam*(domainMin_));
      if (domainMin_ > 0){
        *lowerBound = std::min(-upperBoundLinearRateParam*(domainMax_), -upperBoundLinearRateParam*(domainMin_));
      }else{
        *lowerBound = lowerBoundOfRateParam;

      }
    }else{
      *lowerBound = -upperBoundLinearRateParam*(maxChrNumber - 1);
    }
    *upperBound = upperBoundOfRateParam;

  }else if (index == 1){
    if (continuous_){
      if (domainMax_ < 0){
        *lowerBound = -upperBoundLinearRateParam;
        *upperBound = std::min(-upperBoundOfRateParam / domainMax_, -upperBoundOfRateParam / domainMin_);

      }else{
        if (domainMin_ < 0){ //signs of x_min and x_max are not the same
          *lowerBound = std::min(-upperBoundOfRateParam / domainMax_, -upperBoundOfRateParam / domainMin_);
          *upperBound = upperBoundLinearRateParam;

        }else{
          *upperBound = upperBoundLinearRateParam;
          *lowerBound = std::max(-upperBoundOfRateParam / domainMax_, -upperBoundOfRateParam / domainMin_);
        }
      }

    }else{
      *lowerBound = -upperBoundOfRateParam/(maxChrNumber -1);
      *upperBound = upperBoundLinearRateParam;

    }
    
    
  }else{
    throw Exception("LinearDependencyFunction::getAbsoluteBounds(): index out of bounds!!");
  }
}
/**************************************************************************************/
double LinearBDDependencyFunction::getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, double minChrNum, double maxChrNum){
  if (index != 0){
    throw Exception("LinearDependencyFunction::getParsimonyBound(): index out of bounds!!");
  }
  return (parsimonyBound * 2)/(minChrNum + maxChrNum);
}

/**************************************************************************************/
double LinearBDDependencyFunction::getRate(std::vector<Parameter*> params, double state) const{
  return params[0]->getValue() * state;
}
/**************************************************************************************/
double LinearBDDependencyFunction::getRate(std::vector<double> params, double state) const{
  auto res = params[0] * state;
  if (res < 0){
    throw Exception("LinearBDDependencyFunction::getRate(): Negative rate!!!");
  }
  return res;
}
/**************************************************************************************/
double ExponentailDependencyFunction::getRate(std::vector<Parameter*> params, double state) const{
  return params[0]->getValue() * std::exp((state-1)*params[1]->getValue());
}
/**************************************************************************************/

double ExponentailDependencyFunction::getRate(std::vector<double> params, double state) const{
  auto res = params[0] * std::exp((state-1) * params[1]);
  if (res < 0){
    throw Exception("ExponentailDependencyFunction::getRate(): Negative rate!!!");
  }
  return res;
}
/**************************************************************************************/
void ExponentailDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);

}
/**************************************************************************************/
void ExponentailDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber){
  if (index > 1){
    throw Exception("ExponentailDependencyFunction::getAbsoluteBounds(): Too many parameters!!!");
  }
  if (index == 0){
    *lowerBound = lowerBoundOfRateParam;
    *upperBound = upperBoundOfRateParam;

  }else if (index == 1){
    *lowerBound = lowerBoundOfExpParam;
    *upperBound = upperBoundExpParam/(maxChrNumber-1);

  }
 
}
/**************************************************************************************/
double ExponentailDependencyFunction::getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, double minChrNum, double maxChrNum){
  if (index == 0){
    return parsimonyBound;
  }else if (index == 1){
    return (std::log(parsimonyBound) + std::log(params[0]));
  }else{
    throw Exception("ExponentailDependencyFunction::getParsimonyBound(): ERROR! such parameter does not exist!");
  }
}
/**************************************************************************************/
double PolynomialDependencyFunction::getRate(std::vector<Parameter*> params, double state) const{
  return (params[0]->getValue()) * pow(state + params[1]->getValue(), params[2]->getValue());

}
/**************************************************************************************/
double PolynomialDependencyFunction::getRate(std::vector<double> params, double state) const{
  auto res = (params[0]) * pow(state + params[1], params[2]);
  if (res < 0){
    throw Exception("PolynomialDependencyFunction::getRate(): Negative rate!!!");
  }
  return res;

}
/**************************************************************************************/

void PolynomialDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);

}
void PolynomialDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber){
  if (index == 0){
      *lowerBound = 0;
      *upperBound = upperBoundOfRateParam;

  }else if (index == 1){
    *lowerBound = -domainMin_;
    *upperBound = upperBoundOfRateParam;
  }else if(index == 2){
    *lowerBound = lowerBoundOfExpParam;
    *upperBound = upperBoundExpParam;

  }else{
    throw Exception("PolynomialDependencyFunction::getAbsoluteBounds: index out of bounds!!");
    
  }

}
double LognormalDependencyFunction::getRate(std::vector<Parameter*> params, double state) const{
  auto rangeFactor = params[0]->getValue();
  double scalingFactor = domainMax_/logNormalDomainFactor;
  double transformedState = state/scalingFactor;
  auto mu = params[1]->getValue();
  auto sigma = params[2]->getValue();
  double pi = 2 * acos(0.0);
  auto eq_part_1 = 1/(transformedState*sigma*sqrt(2 * pi));
  auto eq_part_2 = std::exp(-(pow(log(transformedState)-mu, 2)/(2*pow(sigma, 2))));
  return rangeFactor *eq_part_1 * eq_part_2;

}
/*************************************************************************************/
double LognormalDependencyFunction::getRate(std::vector<double> params, double state) const{
  auto rangeFactor = params[0];
  double scalingFactor = domainMax_/logNormalDomainFactor;
  double transformedState = state/scalingFactor;
  auto mu = params[1];
  auto sigma = params[2];
  double pi = 2 * acos(0.0);
  auto eq_part_1 = 1/(transformedState*sigma*sqrt(2 * pi));
  auto eq_part_2 = std::exp(-(pow(log(transformedState)-mu, 2)/(2*pow(sigma, 2))));
  auto res =  rangeFactor *eq_part_1 * eq_part_2;
  if (res < 0){
    throw Exception("LognormalDependencyFunction::getRate(): Negative rate!!!");
  }
  return res;

}
/**************************************************************************************/
void LognormalDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);

}
/*************************************************************************************/
void LognormalDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber){
  *lowerBound = lowerBoundOfRateParam;
  if (index == 0){  // for the range parameter   
    *upperBound = upperBoundOfRateParam;
  }else if (index == 1){ // mu
    *upperBound = upperBoundLinearRateParam;
  }else if (index == 2){  // sigma
    *upperBound = upperBoundLinearRateParam*2;

  }else{
    throw Exception("LognormalDependencyFunction::getAbsoluteBounds(): index out of bounds!!");
    
  }
}

/**************************************************************************************/
double RevSigmoidDependencyFunction::getRate(std::vector<Parameter*> params, double state) const{
  // p1 is the range parameter
  auto p1 = params[0]->getValue();
  // p2 is the exponent multiplier parameter
  auto p2 = params[1]->getValue();
  // p3 is the shift parameter (should manipulate the cut of the reverse sigmoid tail)
  auto p3 = params[2]->getValue();
  // f(x) = p1* (e^-p2(x-p3)/(1+(e^-p2(x-p3))))
  return p1/(1+(std::exp(p2-(p3*state))));// (std::exp(-p2*(x-p3))/(1+(std::exp(-p2*(x-p3)))));

}
/**************************************************************************************/
double RevSigmoidDependencyFunction::getRate(std::vector<double> params, double state) const{
  // p1 is the range parameter
  auto p1 = params[0];
  // p2 is the exponent multiplier parameter
  auto p2 = params[1];
  // p3 is the shift parameter (should manipulate the cut of the reverse sigmoid tail)
  auto p3 = params[2];
  // f(x) = p1* (e^-p2(x-p3)/(1+(e^-p2(x-p3))))
  auto res =  p1/(1+(std::exp(p2-(p3*state))));// (std::exp(-p2*(x-p3))/(1+(std::exp(-p2*(x-p3)))));
  if (res < 0){
    throw Exception("RevSigmoidDependencyFunction::getRate(): Negative rate!!!");
  }
  return res;

}
/**************************************************************************************/
void RevSigmoidDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber){
  
  if (index == 0){  // for the range parameter   
    *lowerBound = lowerBoundOfRateParam;
    *upperBound = upperBoundOfRateParam;
  }else if (index == 1){ // for the exponent parameter
    *lowerBound = lowerBoundOfRateParam;
    *upperBound = (double)(domainMax_-domainMin_+1);
  }else if (index == 2){  // the shift parameter
    *lowerBound = -revSigmoidExpRateParam;
    *upperBound = lowerBoundOfRateParam;

  }else{
    throw Exception("RevSigmoidDependencyFunction::getAbsoluteBounds(): index out of bounds!!");
    
  }
}
/*****************************************************************************/
void RevSigmoidDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);

}
/*****************************************************************************/
double LogitnormalDependencyFunction::getRate(std::vector<Parameter*> params, double state) const{
  auto rangeFactor = params[0]->getValue();
  double transformedState = (state-domainMin_+1)/(domainMax_-domainMin_+2);
  auto mu = params[1]->getValue();
  auto sigma = params[2]->getValue();
  double pi = 2 * acos(0.0);
  double expr_1 = 1/(sigma*sqrt(2 * pi));
  double expr_2 = 1/(transformedState*(1-transformedState));
  double logit_expr = log(transformedState/(1-transformedState));
  double expr_3 = std::exp(-(pow(logit_expr-mu, 2)/(2*pow(sigma, 2))));
  return rangeFactor * expr_1 * expr_2 * expr_3;

}
/*****************************************************************************/
double LogitnormalDependencyFunction::getRate(std::vector<double> params, double state) const{
  auto rangeFactor = params[0];
  double transformedState = (state-domainMin_+1)/(domainMax_-domainMin_+2);
  auto mu = params[1];
  auto sigma = params[2];
  double pi = 2 * acos(0.0);
  double expr_1 = 1/(sigma*sqrt(2 * pi));
  double expr_2 = 1/(transformedState*(1-transformedState));
  double logit_expr = log(transformedState/(1-transformedState));
  double expr_3 = std::exp(-(pow(logit_expr-mu, 2)/(2*pow(sigma, 2))));
  auto res = rangeFactor * expr_1 * expr_2 * expr_3;
  if (res < 0){
    throw Exception("LogitnormalDependencyFunction::getRate(): Negative rate!!!");

  }
  return res;

}
/**************************************************************************************/
void LogitnormalDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);

}
/*************************************************************************************/
void LogitnormalDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber){
  *lowerBound = lowerBoundOfRateParam;
  if (index == 0){  // for the range parameter   
    *upperBound = upperBoundOfRateParam;
  }else if (index == 1){ // mu
    *upperBound = upperBoundLinearRateParam;
  }else if (index == 2){  // sigma
    *upperBound = upperBoundLinearRateParam*2;

  }else{
    throw Exception("LogitnormalDependencyFunction::getAbsoluteBounds(): index out of bounds!!");
    
  }
}
