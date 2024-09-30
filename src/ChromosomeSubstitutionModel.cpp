//
// File: ChromosomeSubstitutionModel.cpp
// Created by: Anat Shafir
// Created on: 2020
//


#include "ChromosomeSubstitutionModel.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

using namespace std;
int lowerLimitBaseNumber = 0;
/*****************************************************************************/
// Functions *** *** *** **** *** *** *** *** **** *** *** *** *** **** ***
/*****************************************************************************/
double ConstantDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  return params[0]->getValue();
}

/**************************************************************************************/
double LinearDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  double func_res = params[0]->getValue() + ((double)(state-1)*params[1]->getValue());
  if (func_res < 0){
    return 0;
  }
  return func_res;

}
/**************************************************************************************/
double LinearDependencyFunction::getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, int minChrNum, int maxChrNum){
  if (index == 0){
    return parsimonyBound;
  }else if (index == 1){
    return params[0]-(params[0]*(maxChrNum+minChrNum)/2) + parsimonyBound;
  }
  throw Exception("LinearDependencyFunction::getParsimonyBound(): No such index!");
}
/**************************************************************************************/
void LinearDependencyFunction::updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, int maxChrNum){
  if (index == 0){
    *lowerBound = std::max(lowerBoundOfRateParam, -params.getParameter(paramsNames[1]).getValue()*(maxChrNum-1));
    *upperBound = upperBoundOfRateParam;
    
  }else if (index == 1){
    *lowerBound = -params.getParameter(paramsNames[0]).getValue()/(maxChrNum-1);
    *upperBound = upperBoundLinearRateParam;       
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
void LinearDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  if (index == 0){
    *lowerBound = lowerBoundOfRateParam;
    *upperBound = upperBoundOfRateParam;
    
  }else if (index == 1){
    *lowerBound = -paramValues[0]/(maxChrNumber-1);
    *upperBound = upperBoundLinearRateParam;       
  }else{
    throw Exception("LinearDependencyFunction::updateBounds(): index out of bounds!!");
    
  }
}
/**************************************************************************************/
void LinearDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
  if (index == 0){
    *lowerBound = -upperBoundLinearRateParam*(maxChrNumber - 1);
    *upperBound = upperBoundOfRateParam;
  }else if (index == 1){
    *lowerBound = -upperBoundOfRateParam/(maxChrNumber -1);
    *upperBound = upperBoundLinearRateParam;
  }else{
    throw Exception("LinearDependencyFunction::getAbsoluteBounds(): index out of bounds!!");
  }
}
/**************************************************************************************/
double LinearBDDependencyFunction::getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, int minChrNum, int maxChrNum){
  if (index != 0){
    throw Exception("LinearDependencyFunction::getParsimonyBound(): index out of bounds!!");
  }
  return (parsimonyBound * 2)/(minChrNum + maxChrNum);
}

/**************************************************************************************/
double LinearBDDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  return params[0]->getValue() * (double)state;
}
/**************************************************************************************/
double ExponentailDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  return params[0]->getValue() * std::exp((double)(state-1)*params[1]->getValue());
}
/**************************************************************************************/

void ExponentailDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);

}
/**************************************************************************************/
void ExponentailDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
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
double ExponentailDependencyFunction::getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, int minChrNum, int maxChrNum){
  if (index == 0){
    return parsimonyBound;
  }else if (index == 1){
    return (std::log(parsimonyBound) + std::log(params[0]));
  }else{
    throw Exception("ExponentailDependencyFunction::getParsimonyBound(): ERROR! such parameter does not exist!");
  }
}
/**************************************************************************************/
double PolynomialDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  return (params[0]->getValue()) * pow((double)(state) + params[1]->getValue(), params[2]->getValue());

}
/**************************************************************************************/

void PolynomialDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);

}
void PolynomialDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
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
double LognormalDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  auto rangeFactor = params[0]->getValue();
  double scalingFactor = (double)domainMax_/logNormalDomainFactor;
  double transformedState = (double)state/scalingFactor;
  auto mu = params[1]->getValue();
  auto sigma = params[2]->getValue();
  double pi = 2 * acos(0.0);
  auto eq_part_1 = 1/((double)(transformedState)*sigma*sqrt(2 * pi));
  auto eq_part_2 = std::exp(-(pow(log(transformedState)-mu, 2)/(2*pow(sigma, 2))));
  return rangeFactor *eq_part_1 * eq_part_2;

}
/**************************************************************************************/
void LognormalDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);

}
/*************************************************************************************/
void LognormalDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
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
double RevSigmoidDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  // p1 is the range parameter
  auto p1 = params[0]->getValue();
  // p2 is the exponent multiplier parameter
  auto p2 = params[1]->getValue();
  // p3 is the shift parameter (should manipulate the cut of the reverse sigmoid tail)
  auto p3 = params[2]->getValue();
  // f(x) = p1* (e^-p2(x-p3)/(1+(e^-p2(x-p3))))
  auto x = static_cast<double>(state);
  return p1/(1+(std::exp(p2-(p3*x))));// (std::exp(-p2*(x-p3))/(1+(std::exp(-p2*(x-p3)))));

}
/**************************************************************************************/
void RevSigmoidDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
  
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
void RevSigmoidDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);

}
/*****************************************************************************/
double LogitnormalDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  auto rangeFactor = params[0]->getValue();
  double transformedState = ((double)state-(double)domainMin_+1)/((double)domainMax_-(double)domainMin_+2);
  auto mu = params[1]->getValue();
  auto sigma = params[2]->getValue();
  double pi = 2 * acos(0.0);
  double expr_1 = 1/(sigma*sqrt(2 * pi));
  double expr_2 = 1/(transformedState*(1-transformedState));
  double logit_expr = log(transformedState/(1-transformedState));
  double expr_3 = std::exp(-(pow(logit_expr-mu, 2)/(2*pow(sigma, 2))));
  return rangeFactor * expr_1 * expr_2 * expr_3;

}
/**************************************************************************************/
void LogitnormalDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);

}
/*************************************************************************************/
void LogitnormalDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
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


/*****************************************************************************/
// CompositeParameter *** *** *** **** *** *** *** *** **** *** *** *** *** 
/******************************************************************************/
double compositeParameter::getRate(size_t state) const{
  return func_->getRate(params_, state);
} 
/******************************************************************************/
ChromosomeNumberDependencyFunction* compositeParameter::setDependencyFunction(ChromosomeNumberDependencyFunction::FunctionType funcType){
  switch (funcType)
  {
  case ChromosomeNumberDependencyFunction::CONSTANT:
    return new ConstantDependencyFunction();
  case ChromosomeNumberDependencyFunction::LINEAR:
    return new LinearDependencyFunction();
  case ChromosomeNumberDependencyFunction::EXP:
    return new ExponentailDependencyFunction();
  case ChromosomeNumberDependencyFunction::LINEAR_BD:
    return new LinearBDDependencyFunction();
  case ChromosomeNumberDependencyFunction::LOGNORMAL:
    return new LognormalDependencyFunction();
  case ChromosomeNumberDependencyFunction::POLYNOMIAL:
    return new PolynomialDependencyFunction();
  case ChromosomeNumberDependencyFunction::REVERSE_SIGMOID:
    return new RevSigmoidDependencyFunction();
  case ChromosomeNumberDependencyFunction::LOGITNORMAL:
    return new LogitnormalDependencyFunction();
  default:
    throw Exception("compositeParameter::getDependencyFunction(): No such function!!");
  }
}

/****************************************************************************/
std::vector<double> compositeParameter::getParameterValues() const{
  std::vector<double> values;
  for (size_t i = 0; i < getSize(); i++){
    values.push_back(params_[i]->getValue());
  }
  return values;
 
}
/*****************************************************************************/
std::vector<std::string> compositeParameter::getRelatedParameterNames(ParameterList &params, std::string pattern){
  std::vector<std::string> paramNames = params.getParameterNames();
  std::vector<std::string> matchingNames;
  for (size_t i = 0; i < paramNames.size(); i++){
    std::string fullParamName = paramNames[i];
    if (fullParamName.find(pattern) != string::npos){
      matchingNames.push_back(fullParamName);
    }
  }
  return matchingNames;
}

/*****************************************************************************/
// Chromosome model ///////////////////////////////////////////////////////////
/******************************************************************************/
ChromosomeSubstitutionModel :: ChromosomeSubstitutionModel(
  const ChromosomeAlphabet* alpha, 
  std::vector<double> gain, 
  std::vector<double> loss, 
  std::vector<double> dupl, 
  std::vector<double> demi,
  int baseNum,
  std::vector<double> baseNumR,
  unsigned int chrRange, 
  rootFreqType freqType,
  std::vector<int> rateChangeType,
  bool demiOnlyForEven,
  bool simulated):
    AbstractParameterAliasable("Chromosome."),
    AbstractSubstitutionModel(alpha, std::shared_ptr<const StateMap>(new CanonicalStateMap(alpha, false)), "Chromosome."),
    gain_(0),
    loss_(0),
    dupl_(0),
    demiploidy_(0),
    baseNum_(baseNum),
    baseNumR_(0),
    maxChrRange_(chrRange),
    freqType_(freqType),
    ChrMinNum_(alpha->getMin()),
    ChrMaxNum_(alpha->getMax()),
    firstNormQ_(0),
    pijtCalledFromDeriv_(false),
    gainFunc_(),
    lossFunc_(),
    duplFunc_(),
    demiFunc_(),
    baseNumRFunc_(),
    simulated_(simulated),
    demiOnlyForEven_(demiOnlyForEven),
    vPowExp_()
    
{
    size_t startNonComposite = getNumberOfNonCompositeParams();
    for (size_t i = startNonComposite; i < ChromosomeSubstitutionModel::paramType::NUM_OF_CHR_PARAMS; i++){
      switch (i)
      {
      case ChromosomeSubstitutionModel::GAIN:
        gainFunc_ = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
        break;
      case ChromosomeSubstitutionModel::LOSS:
        lossFunc_ = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
        break;
      case ChromosomeSubstitutionModel::DUPL:
        duplFunc_ = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
        break;
      case ChromosomeSubstitutionModel::DEMIDUPL:
        demiFunc_ = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
        break;
      case ChromosomeSubstitutionModel::BASENUMR:
        baseNumRFunc_ =  static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
        break;
      default:
        throw Exception("ChromosomeSubstitutionModel :: ChromosomeSubstitutionModel(): No such parameter!");
        break;
      }

    }

    updateParameters(gain, loss, dupl, demi, baseNumR);
    computeFrequencies(false);
    isScalable_ = false;    //in ChromEvol the matrix should be not normalized
    updateMatrices();

}
/******************************************************************************/
ChromosomeSubstitutionModel::ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
  std::map<int, vector<double>> mapOfParamValues,
  int baseNum,
  unsigned int chrRange, 
  rootFreqType freqType,
  vector<int> rateChangeType,
  bool demiOnlyForEven,
  bool simulated):
    AbstractParameterAliasable("Chromosome."),
    AbstractSubstitutionModel(alpha, std::shared_ptr<const StateMap>(new CanonicalStateMap(alpha, false)), "Chromosome."),
    gain_(0),
    loss_(0),
    dupl_(0),
    demiploidy_(0),
    baseNum_(baseNum),
    baseNumR_(0),
    maxChrRange_(chrRange),
    freqType_(freqType),
    ChrMinNum_(alpha->getMin()),
    ChrMaxNum_(alpha->getMax()),
    firstNormQ_(0),
    pijtCalledFromDeriv_(false),
    gainFunc_(),
    lossFunc_(),
    duplFunc_(),
    demiFunc_(),
    baseNumRFunc_(),
    simulated_(simulated),
    demiOnlyForEven_(demiOnlyForEven),
    vPowExp_()
    
{
  size_t startNonComposite = getNumberOfNonCompositeParams();
  for (size_t i = startNonComposite; i < ChromosomeSubstitutionModel::paramType::NUM_OF_CHR_PARAMS; i++){
    switch (i)
    {
    case ChromosomeSubstitutionModel::GAIN:
      gainFunc_ = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
      break;
    case ChromosomeSubstitutionModel::LOSS:
      lossFunc_ = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
      break;
    case ChromosomeSubstitutionModel::DUPL:
      duplFunc_ = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
      break;
    case ChromosomeSubstitutionModel::DEMIDUPL:
      demiFunc_ = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
      break;
    case ChromosomeSubstitutionModel::BASENUMR:
      baseNumRFunc_ =  static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
      break;
    default:
      throw Exception("ChromosomeSubstitutionModel :: ChromosomeSubstitutionModel(): No such parameter!");
      break;
    }

  }

  updateParameters(mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::GAIN)], 
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::LOSS)], 
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DUPL)],
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)],
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)]);
  computeFrequencies(false);
  isScalable_ = false;    //in ChromEvol the matrix should be not normalized
  updateMatrices();


}


/******************************************************************************/
ChromosomeSubstitutionModel* ChromosomeSubstitutionModel::initRandomModel(
  const ChromosomeAlphabet* alpha,
  int &baseNumber,
  map<int, vector<double>> initParams,
  unsigned int chrRange,
  rootFreqType rootFrequenciesType,
  vector<int> rateChangeType,
  vector<int>& fixedParams,
  bool demiOnlyForEven,
  double parsimonyBound)
{
  std::map<int, vector<double>> mapRandomParams;
  int newBaseNumber = baseNumber;
  size_t startCompositeParams = getNumberOfNonCompositeParams();

  for (int i = 0; i < ChromosomeSubstitutionModel::paramType::NUM_OF_CHR_PARAMS; i++){
    if (std::find(fixedParams.begin(), fixedParams.end(), i) != fixedParams.end()){
      if (static_cast<ChromosomeSubstitutionModel::paramType>(i) != ChromosomeSubstitutionModel::BASENUM){
        mapRandomParams[i] = initParams[i];

      }
    }else{
      vector<double> paramValues;
      double lowerBound;
      double upperBound;
      if (static_cast<ChromosomeSubstitutionModel::paramType>(i) == ChromosomeSubstitutionModel::BASENUM){
        if (baseNumber == IgnoreParam){
          continue;
        }
        lowerBound = lowerLimitBaseNumber;
        upperBound = std::max((int)chrRange, lowerLimitBaseNumber+1);
        newBaseNumber = static_cast<int>(lowerBound + RandomTools::giveIntRandomNumberBetweenZeroAndEntry((int)(upperBound-lowerBound)));
        continue;

      }else{
        //int compositeParamType = compositeParameter::getCompositeRateType(i);
        if (initParams[i].size() != 0){
          ChromosomeNumberDependencyFunction::FunctionType funcType = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startCompositeParams]);
          ChromosomeNumberDependencyFunction* functionOp = compositeParameter::setDependencyFunction(funcType);
          functionOp->setDomainsIfNeeded(alpha->getMin(), alpha->getMax());

          auto numOfParameters = functionOp->getNumOfParameters();
          for (size_t j = 0; j < numOfParameters; j++){
            functionOp->getBoundsForInitialParams(j, paramValues, &lowerBound, &upperBound, alpha->getMax());
            double upperBoundCandidate = functionOp->getParsimonyBound(paramValues, parsimonyBound, j, alpha->getMin(), alpha->getMax());
            //compositeParameter::getBoundsForInitialParams(func, j, paramValues, &lowerBound, &upperBound, alpha->getMax(), true);
            if (parsimonyBound > 0){
              if (upperBoundCandidate >= lowerBound){
                upperBound = std::min(upperBound, upperBoundCandidate);
              }         
            }
            double randomValue = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBound, upperBound);
            paramValues.push_back(randomValue);

          }
          delete functionOp;
          
        }
      }
      mapRandomParams[i] = paramValues;

    }
  }

  ChromosomeSubstitutionModel* model = new ChromosomeSubstitutionModel(alpha, mapRandomParams, newBaseNumber, chrRange, rootFrequenciesType, rateChangeType, demiOnlyForEven);//, useExtendedFloat);
  return model;

}
/******************************************************************************/

std::vector<Parameter*> ChromosomeSubstitutionModel::createCompositeParameter(ChromosomeNumberDependencyFunction::FunctionType &func, std::string paramName, vector<double> &vectorOfValues){
  std::vector<Parameter*> params;
  if (func == ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    return params;
  }
  for (size_t i = 0; i < vectorOfValues.size(); i++){
    auto paramValue = vectorOfValues[i];
    if (paramValue == IgnoreParam){
      throw Exception("ChromosomeSubstitutionModel::createCompositeParameter(): Function is not supposed to be defined as a legal function name if the parameter should be ignored!");
    }
    double lowerBound;
    double upperBound;
    ChromosomeNumberDependencyFunction* functionOp =  compositeParameter::setDependencyFunction(func);
    functionOp->setDomainsIfNeeded(ChrMinNum_, ChrMaxNum_);
    functionOp->getAbsoluteBounds(i, &lowerBound, &upperBound, ChrMaxNum_);
    if ((simulated_) && (lowerBound > paramValue)){
      lowerBound = paramValue - EPSILON;

    }
    delete functionOp;

    // in simulations it sometimes happens when the upper bound is lower than the value itself,
    // because the max number in the simulating function is much larger. In these cases it is important to
    // change the upper bound, such that it will be  >= parameter value
    if (simulated_ && upperBound <= paramValue){
      upperBound = paramValue + EPSILON;
    }
    std::shared_ptr<IntervalConstraint> interval = make_shared<IntervalConstraint>(lowerBound, upperBound, false, true);
    Parameter* param = new Parameter("Chromosome."+ paramName + std::to_string(i), paramValue, interval);
    params.push_back(param);
    //addParameter_(param);
  }
  return params;


}
/******************************************************************************/
void ChromosomeSubstitutionModel::addCompositeParameter(std::vector<Parameter*> parameters){
  for (size_t i = 0; i < parameters.size(); i++){
    addParameter_(parameters[i]);
  }
}

/******************************************************************************/
void ChromosomeSubstitutionModel::updateParameters(vector<double> &gain, vector<double> &loss, vector<double> &dupl, vector<double> &demi, vector<double> &baseNumR){
  if (baseNum_ != IgnoreParam){
    std::shared_ptr<IntervalConstraint> interval_baseNum = make_shared<IntervalConstraint>(lowerLimitBaseNumber, (int)maxChrRange_, true, true);
    addParameter_(new Parameter("Chromosome.baseNum", baseNum_, interval_baseNum));
    
  }
  std::vector<size_t> numOfParamsVector;
  auto baseNumParams = createCompositeParameter(baseNumRFunc_, "baseNumR", baseNumR);
  numOfParamsVector.push_back(baseNumParams.size());

  auto duplParams = createCompositeParameter(duplFunc_, "dupl", dupl);
  numOfParamsVector.push_back(duplParams.size());

  auto lossParams = createCompositeParameter(lossFunc_, "loss", loss);
  numOfParamsVector.push_back(lossParams.size());

  auto gainParams = createCompositeParameter(gainFunc_, "gain", gain);
  numOfParamsVector.push_back(gainParams.size());
  
  auto demiParams = createCompositeParameter(demiFunc_, "demi", demi);
  //numOfParamsVector.push_back(demiParams.size());
  size_t maxSize = *max_element(numOfParamsVector.begin(), numOfParamsVector.end());
  for (size_t i = 0; i < maxSize; i++){
    if (i < baseNumParams.size()){
      addParameter_(baseNumParams[i]);
    }
    if (i < duplParams.size()){
      addParameter_(duplParams[i]);
    }

    if (i < lossParams.size()){
      addParameter_(lossParams[i]);
    }
    if (i < gainParams.size()){
      addParameter_(gainParams[i]);
    }

  }
  for (size_t i = 0; i < demiParams.size(); i++){
    addParameter_(demiParams[i]);
  }
  if (gainFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    gain_ = new compositeParameter(gainFunc_, "gain", gainParams);
    gain_->func_->setDomainsIfNeeded(ChrMinNum_, ChrMaxNum_);

  }
  if (lossFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    loss_ = new compositeParameter(lossFunc_, "loss", lossParams);
    loss_->func_->setDomainsIfNeeded(ChrMinNum_, ChrMaxNum_);
  }
  if (duplFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    dupl_ = new compositeParameter(duplFunc_, "dupl", duplParams);
    dupl_->func_->setDomainsIfNeeded(ChrMinNum_, ChrMaxNum_);
  }
  if (demiFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    demiploidy_ = new compositeParameter(demiFunc_, "demi", demiParams);
    demiploidy_->func_->setDomainsIfNeeded(ChrMinNum_, ChrMaxNum_);
  }
  if (baseNumRFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    baseNumR_ = new compositeParameter(baseNumRFunc_, "baseNumR", baseNumParams);
    baseNumR_->func_->setDomainsIfNeeded(ChrMinNum_, ChrMaxNum_);
  }


}

/******************************************************************************/
void ChromosomeSubstitutionModel::getCompositeParametersValues(std::string paramName, compositeParameter* param){
  for (size_t i = 0; i < param->getSize(); i++){   
    getParameterValue(paramName + std::to_string(i)); //do I really need it?

  }
}
/******************************************************************************/
void ChromosomeSubstitutionModel::getParametersValues(){
    if (gain_ != 0){
      getCompositeParametersValues("gain", gain_);
    }
    if (loss_ != 0){
      getCompositeParametersValues("loss", loss_);
    }
    if (dupl_ != 0){
      getCompositeParametersValues("dupl", dupl_);
    }
    if (demiploidy_ != 0){// && (demiploidy_ != DemiEqualDupl)){
      getCompositeParametersValues("demi", demiploidy_);
    }


    if(baseNum_ != IgnoreParam){
      baseNum_ = (int)getParameterValue("baseNum");
    }
    if (baseNumR_ != 0){
      getCompositeParametersValues("baseNumR", baseNumR_);
    }  
    //checkParametersBounds();
}

/*******************************************************************************/
void ChromosomeSubstitutionModel::checkParametersBounds() const{
  std::cout << "All bounds" <<endl;
  const ParameterList params = getParameters();
  for (size_t i = 0; i < params.size(); i++){
    std::cout << params[i].getName() << " bound: "<< dynamic_pointer_cast<IntervalConstraint>(params[i].getConstraint())->getLowerBound() <<  ", value: " << params[i].getValue() << std::endl;

  }
}
/*******************************************************************************/
void ChromosomeSubstitutionModel::updateMatrices(){
    //update model parameters
    getParametersValues();

    //update generator matrix
    size_t maxChrNum = (size_t)(getMax());
    size_t minChrNum = (size_t)(getMin()); 
    MatrixTools::fill(generator_, 0);
 
    // updating Q matrix
    for (size_t i = minChrNum; i < maxChrNum+1; i++){
        // gain
        if (i + 1 < maxChrNum+1){
            updateQWithGain(i, minChrNum);
        //loss
        }if (i-1 >= minChrNum){
            updateQWithLoss(i, minChrNum);
        //duplication         
        }if (2*i <= maxChrNum){
            updateQWithDupl(i, minChrNum);
        }else if (i != maxChrNum){
            updateQWithDupl(i, minChrNum, maxChrNum);
        }
        //demi-ploidy
        updateQWithDemiDupl(i, minChrNum, maxChrNum);   

        if (i < maxChrNum){
          if (baseNum_ != IgnoreParam){
            updateQWithBaseNumParameters(i, minChrNum, maxChrNum);

          }
        }
        
    }
    setDiagonal();  //sets Qii to -sigma(Qij)
    updateEigenMatrices();
    firstNormQ_ = getFirstNorm();

}
/*******************************************************************************/
void ChromosomeSubstitutionModel::correctBaseNumForSimulation(int maxChrNumInferred){
    //update generator matrix
    size_t maxChrNum = (size_t)(getMax());
    size_t minChrNum = (size_t)(getMin()); 
 
    // updating Q matrix
    for (size_t i = minChrNum; i < maxChrNum; i++){
      if (baseNumR_->getRate(i) < 0){
        throw Exception("ChromosomeSubstitutionModel::correctBaseNumForSimulation():Negative base number rate!");
      }
      for (size_t j = i + 1; j < maxChrNum + 1; j ++){
        if (j == maxChrNum){
          if (((j-i) <= maxChrRange_) && ((int)(j-i) > baseNum_)){
            generator_(i-minChrNum, maxChrNum-minChrNum) -= baseNumR_->getRate(i);
          }


        }else{
          if ((j-i) % baseNum_ == 0){
            if (i > (size_t)maxChrNumInferred){
              if (((j-i) <= maxChrRange_) && ((int)(j-i) > baseNum_)){
                generator_(i - minChrNum, j - minChrNum) -= baseNumR_->getRate(i);
              }
            }else{
              if (((int)j-maxChrNumInferred) >= baseNum_ ){
                if ((j-i) <= maxChrRange_){
                  if (((int)i != maxChrNumInferred) || ((int)j-maxChrNumInferred != baseNum_)){
                    generator_(i - minChrNum, j - minChrNum) -= baseNumR_->getRate(i);

                  }                 
                }
              }
            }
          } 
        }
      }
    }
    setDiagonal();  //sets Qii to -sigma(Qij)
    updateEigenMatrices();
    firstNormQ_ = getFirstNorm();

}
/*******************************************************************************/
void ChromosomeSubstitutionModel::updateQWithDemiDupl(size_t i, size_t minChrNum, size_t maxChrNum){
  //double demiploidy;
  
  if (demiploidy_ != 0){
    if (demiploidy_->getRate(i) < 0){
      throw Exception("ChromosomeSubstitutionModel::updateQWithDemiDupl(): Negative demiploidy rate!");
    }
    if (i % 2 == 0 && (double)i * 1.5 <= (double)maxChrNum){

      generator_(i-minChrNum, (size_t)((double)i * 1.5)-minChrNum) += demiploidy_->getRate(i);

                        
    }else if (i % 2 != 0 && (size_t)ceil((double)i*1.5) <= maxChrNum){
      if (!demiOnlyForEven_){
        if (i == 1){
          generator_(i-minChrNum, (size_t)ceil((double)i * 1.5)-minChrNum) += demiploidy_->getRate(i);
        }else{
          generator_(i-minChrNum, (size_t)ceil((double)i * 1.5)-minChrNum) += demiploidy_->getRate(i)/2;
          generator_(i-minChrNum, (size_t)floor((double)i * 1.5)-minChrNum) += demiploidy_->getRate(i)/2;

        }

      }


    }else{
      if (i != maxChrNum){
        generator_(i-minChrNum, maxChrNum-minChrNum) += demiploidy_->getRate(i);
      }

    }

  }
}
/*******************************************************************************/
// double ChromosomeSubstitutionModel::getRate (size_t state, double constRate, double changeRate) const{
//   if ((constRate == IgnoreParam) && (changeRate == IgnoreParam)){
//     return IgnoreParam;
//   }
//   double totalRate;
//   if (constRate == IgnoreParam){
//     // a birth-death-like model
//     totalRate = changeRate;
//   }else{
//     //const rate is not to be ignored
//     totalRate = constRate;
//   }
//   if (changeRate == IgnoreParam){
//     return totalRate; //only const rate
//   }else{
//     if (rateChangeFuncType_ == rateChangeFunc::LINEAR){
//       totalRate += (changeRate* (double)(state-1));
//     }else if (rateChangeFuncType_ == rateChangeFunc::EXP){
//       totalRate *= (exp(changeRate* (double)(state-1)));
//     }
//   }
//   return totalRate;
// }
/*******************************************************************************/
void ChromosomeSubstitutionModel::updateQWithGain(size_t i, size_t minChrNum){
  if (gain_ == 0){
    return;
  }
  double gainRate = gain_->getRate(i);
  if (gainRate < 0){
    throw Exception ("ChromosomeSubstitutionModel::updateQWithGain(): negative gain rate!");
  }
  generator_(i-minChrNum, i+1-minChrNum) += gainRate;


}
/*******************************************************************************/
void ChromosomeSubstitutionModel::updateQWithLoss(size_t i, size_t minChrNum){
  //generator_(i-minChrNum, i-1-minChrNum) = loss_ + (lossR_* i);
  if (loss_ == 0){
    return;
  }
  double lossRate = loss_->getRate(i);
  if (lossRate < 0){
    throw Exception ("ChromosomeSubstitutionModel::updateQWithLoss(): negative loss rate!");
  }
  generator_(i-minChrNum, i-1-minChrNum) += lossRate;


}
/*******************************************************************************/
void ChromosomeSubstitutionModel::updateQWithDupl(size_t i, size_t minChrNum, size_t maxChrNum){
  if (dupl_ == 0){
    return;
  }
  // if the transition is not to maxChr
  double duplRate = dupl_->getRate(i);
  if (duplRate < 0){
    throw Exception("ChromosomeSubstitutionModel::updateQWithDupl(): Negative dupl rate!");
  }

  if (maxChrNum == 0){
    generator_(i-minChrNum, (2 * i)-minChrNum) += duplRate;

  }else{
     generator_(i-minChrNum, maxChrNum-minChrNum) += duplRate;

  }
}


/********************************************************************************/
void ChromosomeSubstitutionModel::updateQWithBaseNumParameters(size_t currChrNum, size_t minChrNum, size_t maxChrNum){
  if ( baseNumR_->getRate(currChrNum) < 0){
    throw Exception("ChromosomeSubstitutionModel::updateQWithBaseNumParameters():Negative base number rate!");
  }
  for (size_t j = currChrNum + 1; j < maxChrNum + 1; j ++){
    if (j == maxChrNum){
      if ((j-currChrNum) <= maxChrRange_){
        generator_(currChrNum-minChrNum, maxChrNum-minChrNum) += baseNumR_->getRate(currChrNum);
      }
    }else{
      if ((j-currChrNum) % baseNum_ == 0){
        if ((j-currChrNum) <= maxChrRange_){
          generator_(currChrNum - minChrNum, j - minChrNum) += baseNumR_->getRate(currChrNum);
        }
      }
    }
  }
}

/******************************************************************************/
void ChromosomeSubstitutionModel::setFreq(map<int, double>& freqs)
{
  for (size_t i = 0; i < size_; ++i)
  {
    freq_[i] = freqs[static_cast<int>(i)];
  }

}
/******************************************************************************/
void ChromosomeSubstitutionModel::updateEigenMatrices()
{
  // Compute eigen values and vectors:
  if (enableEigenDecomposition())
  {
    // Look for null lines (such as stop lines)
    // ie null diagonal elements

    size_t nbStop=0;
    size_t salph = getNumberOfStates();
    vector<bool> vnull(salph); // vector of the indices of lines with
                               // only zeros

    for (size_t i = 0; i < salph; i++)
    {
      if (abs(generator_(i, i)) < NumConstants::TINY())
      {
        nbStop++;
        vnull[i]=true;
      }
      else
        vnull[i]=false;
    }
        
    if (nbStop != 0)
    {
      size_t salphok=salph - nbStop;
      
      RowMatrix<double> gk(salphok, salphok);
      size_t gi = 0, gj = 0;

      for (size_t i = 0; i < salph; i++)
      {
        if (!vnull[i])
        {
          gj = 0;
          for (size_t j = 0; j < salph; j++)
          {
            if (!vnull[j])
            {
              gk(i - gi, j - gj) = generator_(i, j);
            }
            else
              gj++;
          }
        }
        else
          gi++;
      }

      EigenValue<double> ev(gk);
      eigenValues_ = ev.getRealEigenValues();
      iEigenValues_ = ev.getImagEigenValues();

      for (size_t i = 0; i < nbStop; i++)
      {
        eigenValues_.push_back(0);
        iEigenValues_.push_back(0);
      }

      RowMatrix<double> rev = ev.getV();
      rightEigenVectors_.resize(salph, salph);
      gi = 0;
      for (size_t i = 0; i < salph; i++)
      {
        if (vnull[i])
        {
          gi++;
          for (size_t j = 0; j < salph; j++)
          {
            rightEigenVectors_(i, j) = 0;
          }

          rightEigenVectors_(i, salphok + gi - 1) = 1;
        }
        else
        {
          for (size_t j = 0; j < salphok; j++)
          {
            rightEigenVectors_(i, j) = rev(i - gi, j);
          }

          for (size_t j = salphok; j < salph; j++)
          {
            rightEigenVectors_(i, j) = 0;
          }
        }
      }
    }
    else
    {
      EigenValue<double> ev(generator_);
      rightEigenVectors_ = ev.getV();
      eigenValues_ = ev.getRealEigenValues();
      iEigenValues_ = ev.getImagEigenValues();
      nbStop = 0;
    }

    /// Now check inversion and diagonalization
    try
    {
      MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);

      // is it diagonalizable ?
      isDiagonalizable_ = true;

      if (!dynamic_cast<ReversibleSubstitutionModel*>(this))
      {
        for (auto& vi : iEigenValues_)
        {
          if (abs(vi) > NumConstants::TINY())
          {
            isDiagonalizable_ = false;
            break;
          }
        }
      }
      
      // looking for the vector of 0 eigenvalues

      vector<size_t> vNullEv;
      for (size_t i = 0; i< salph - nbStop; i++)
        if ((abs(eigenValues_[i]) < NumConstants::SMALL()) && (abs(iEigenValues_[i]) < NumConstants::SMALL()))
          vNullEv.push_back(i);
      

      // pb to find unique null eigenvalue      
      isNonSingular_=(vNullEv.size()==1);

      size_t nulleigen;
      
      double val;
      if (!isNonSingular_)
      {
        //look or check which non-stop right eigen vector elements are
        //equal.
        for (auto cnull : vNullEv)
        {
          size_t i = 0;
          while (vnull[i])
            i++;
          
          val = rightEigenVectors_(i, cnull);
          i++;
          
          while (i < salph)
          {
            if (!vnull[i])
            {
              if (abs(rightEigenVectors_(i, cnull) - val) > NumConstants::SMALL())
                break;
            }
            i++;
          }
          
          if (i >= salph)
          {
            isNonSingular_ = true;
            nulleigen=cnull;
            break;
          }
        }
      }
      else
        nulleigen=vNullEv[0];
      
      if (isNonSingular_)
      {
        eigenValues_[nulleigen] = 0; // to avoid approximation errors on long long branches
        iEigenValues_[nulleigen] = 0; // to avoid approximation errors on long long branches


      }
      else
      {
        //ApplicationTools::displayMessage("AbstractSubstitutionModel::updateMatrices : Unable to find eigenvector for eigenvalue 0. Taylor series used instead.");
        isDiagonalizable_ = false;
      }
    }
    // if rightEigenVectors_ is singular
    catch (ZeroDivisionException& e)
    {
      //ApplicationTools::displayMessage("AbstractSubstitutionModel::updateMatrices : Singularity during diagonalization. Taylor series used instead.");
      isNonSingular_ = false;
      isDiagonalizable_ = false;
    }

    if (vPowExp_.size() == 0){
      vPowExp_.resize(30);
    }
    MatrixTools::getId(salph, vPowExp_[0]);
    MatrixTools::Taylor(generator_, 30, vPowExp_);

  }

}
/******************************************************************************/
void ChromosomeSubstitutionModel::calculatePijtUsingEigenValues(double t) const{
  if (isDiagonalizable_){
    MatrixTools::mult<double>(rightEigenVectors_, VectorTools::exp(eigenValues_ * (rate_ * t)), leftEigenVectors_, pijt_);
  }else{
    std::vector<double> vdia(size_);
    std::vector<double> vup(size_ - 1);
    std::vector<double> vlo(size_ - 1);
    double c = 0, s = 0;
    double l = rate_ * t;
    for (size_t i = 0; i < size_; i++){
      vdia[i] = std::exp(eigenValues_[i] * l);
      if (iEigenValues_[i] != 0){
        s = std::sin(iEigenValues_[i] * l);
        c = std::cos(iEigenValues_[i] * l);
        vup[i] = vdia[i] * s;
        vlo[i] = -vup[i];
        vdia[i] *= c;
        vdia[i + 1] = vdia[i]; // trick to avoid computation
        i++;
      }else{
        if (i < size_ - 1){
          vup[i] = 0;
          vlo[i] = 0;
        }
      }
    }
    MatrixTools::mult<double>(rightEigenVectors_, vdia, vup, vlo, leftEigenVectors_, pijt_);
  }
}


/******************************************************************************/

const Matrix<double>& ChromosomeSubstitutionModel::getPij_t(double t) const
{
  if (t == 0)
  {
    MatrixTools::getId(size_, pijt_);
  }
  else if (isNonSingular_)
  {
    calculatePijtUsingEigenValues(t);
    
  }
  else
  {
    RowMatrix<double> pijt_temp;
    MatrixTools::getId(size_, pijt_temp);
    double s = 1.0;
    double v = rate_ * t;
    double norm = v * firstNormQ_;
    size_t m = 0;
    bool converged = false;
    //while (v > 0.5)    // exp(r*t*A)=(exp(r*t/(2^m) A))^(2^m)
    while (norm > 0.5)
    {
      m += 1;
      v /= 2;
      norm /= 2;
    }
    for (size_t iternum = 2; iternum <  vPowExp_.size(); iternum++){
      calculateExp_Qt(iternum, &s, v);

      if (iternum > 2){
        converged = checkIfReachedConvergence(pijt_, pijt_temp);
        if (converged){
          break;
        }
      }
      MatrixTools::copy(pijt_, pijt_temp);
      if (iternum > 250){
        //std :: cout << "ERROR: Pijt did not reach convergence for t = "<< t <<"!"<<endl;
        throw Exception("ChromosomeSubstitutionModel: Taylor series did not reach convergence!");
        break;
      }
      if (iternum == vPowExp_.size()-1 && !converged){  //need to add more powers to the matrix
        RowMatrix<double> new_pow;
      //new_pow.resize(size_, size_);
        MatrixTools :: mult(vPowExp_[vPowExp_.size()-1], generator_, new_pow);
        vPowExp_.push_back(new_pow);

      }

    }
    while (m > 0){  // recover the 2^m
      
      MatrixTools::mult(pijt_, pijt_, tmpMat_);
      MatrixTools::copy(tmpMat_, pijt_);

      m--;
    }
  }
  
  //just for test/////////////////////
  // bool correct = true;
  if (!pijtCalledFromDeriv_){
    for (size_t i = 0; i < size_; i++){
      for (size_t j = 0; j < size_; j++){
        if (pijt_(i,j) < 0){
          pijt_(i,j) = NumConstants::VERY_TINY(); // trying to do it exactly as in ChromEvol. Maybe the "nan" problem will be solved
          //pijt_(i,j) = 0;
        }
        else if (pijt_(i, j) > 1){
          pijt_(i,j) = 1.0;
        }

      }
    }

  }

  /////////////////////////////
  return pijt_;
}
double ChromosomeSubstitutionModel::getFirstNorm() const{
  double norm = 0;
  for (size_t i = 0; i < size_; i++){
    for (size_t j = 0; j < size_; j++){
      norm += fabs(generator_(i,j));

    }
  }
  return norm;
}

/******************************************************************************/

bool ChromosomeSubstitutionModel::checkIfReachedConvergence(const Matrix<double>& pijt, const Matrix<double>& mt_prev) const{
    for (size_t i = 0; i < pijt.getNumberOfRows(); i++){
        for (size_t j = 0; j < pijt.getNumberOfColumns(); j++){
            double diff = fabs(pijt(i,j) - mt_prev(i,j));
            if (diff > get_epsilon()){
                return false;
            }else if ((pijt(i,j) + get_epsilon() < 0) || (pijt(i,j) > 1 + get_epsilon())){
              return false;
            }
        }
    }
    return true;
}

/******************************************************************************/
void ChromosomeSubstitutionModel::calculateExp_Qt(size_t pow, double s, size_t m, double v) const{
  MatrixTools::getId(size_, pijt_);
  for (size_t i = 1; i <= pow; i++){
    s *= v / static_cast<double>(i);// the initial value of v is rt/(2^m)
    MatrixTools::add(pijt_, s, vPowExp_[i]);

  }
  while (m > 0)  // recover the 2^m
  {
    MatrixTools::mult(pijt_, pijt_, tmpMat_);
    MatrixTools::copy(tmpMat_, pijt_);

    m--;
  }

}

/******************************************************************************/
const Matrix<double>& ChromosomeSubstitutionModel::getdPij_dt  (double d) const{
  pijtCalledFromDeriv_ = true;
  RowMatrix<double> pijt;
  //if (!pijt_calculated_){
  pijt = getPij_t(d);
  MatrixTools::mult(pijt, generator_, dpijt_);
  MatrixTools::scale(dpijt_, rate_);
  pijtCalledFromDeriv_ = false;
  return dpijt_;

}

/******************************************************************************/
const Matrix<double>& ChromosomeSubstitutionModel::getd2Pij_dt2(double d) const{
  pijtCalledFromDeriv_ = true;
  RowMatrix<double> pijt;
  //if (!pijt_calculated_){
  pijt = getPij_t(d);
  
  MatrixTools::mult(vPowExp_[2], pijt, d2pijt_);
  MatrixTools::scale(d2pijt_, rate_ * rate_);
  pijtCalledFromDeriv_ = false;
  return d2pijt_;
}

/******************************************************************************/
const Matrix<double>& ChromosomeSubstitutionModel::getPij_t_func2(double d) const{
  RowMatrix<double> pijt_temp;
  MatrixTools::getId(size_, pijt_temp);
  double s = 1.0;
  double v = rate_ * d;
  size_t m = 0;
  bool converged = false;
  while (v > 0.5)    // exp(r*t*A)=(exp(r*t/(2^m) A))^(2^m)
  {
    m += 1;
    v /= 2;
  }
  for (size_t iternum = 2; iternum <  vPowExp_.size(); iternum++){
    calculateExp_Qt(iternum, s, m, v);

    if (iternum > 1){
      converged = checkIfReachedConvergence(pijt_, pijt_temp);
      if (converged){
        break;
      }
    }
    MatrixTools::copy(pijt_, pijt_temp);
    if (iternum == vPowExp_.size()-1 && !converged){  //need to add more powers to the matrix
      RowMatrix<double> new_pow;
      //new_pow.resize(size_, size_);
      MatrixTools :: mult(vPowExp_[vPowExp_.size()-1], generator_, new_pow);
      vPowExp_.push_back(new_pow);

    }

  }
  return pijt_;

}
/******************************************************************************/
const Matrix<double>& ChromosomeSubstitutionModel::getPij_t_func3(double d) const{
  MatrixTools::getId(size_, pijt_);
  double s = 1.0;
  double v = rate_ * d;
  size_t m = 0;
  while (v > 0.5)    // exp(r*t*A)=(exp(r*t/(2^m) A))^(2^m)
  {
    m += 1;
    v /= 2;
  }
  for (size_t i = 1; i < vPowExp_.size(); i++)
  {
    s *= v / static_cast<double>(i);
    MatrixTools::add(pijt_, s, vPowExp_[i]);
  }
  while (m > 0)  // recover the 2^m
  {
    MatrixTools::mult(pijt_, pijt_, tmpMat_);
    MatrixTools::copy(tmpMat_, pijt_);
    m--;
  }

//  MatrixTools::print(pijt_);
  return pijt_;

}

/******************************************************************************/
void ChromosomeSubstitutionModel::calculateExp_Qt(size_t pow, double* s, double v) const{
  if (pow == 2){
    MatrixTools::getId(size_, pijt_);
    for (size_t i = 1; i <= pow; i++){
      *s *= v / static_cast<double>(i);// the initial value of v is rt/(2^m)
      MatrixTools::add(pijt_, *s, vPowExp_[i]);
      
    }

  }else{
    *s *= v / static_cast<double>(pow);
    MatrixTools::add(pijt_, *s, vPowExp_[pow]);
 
  }
  

}
/******************************************************************************/
const Matrix<double>& ChromosomeSubstitutionModel::getPij_t_func4(double d) const{
  RowMatrix<double> pijt_temp;
  MatrixTools::getId(size_, pijt_temp);
  double s = 1.0;
  double v = rate_ * d;
  size_t m = 0;
  bool converged = false;
  while (v > 0.5)    // exp(r*t*A)=(exp(r*t/(2^m) A))^(2^m)
  {
    m += 1;
    v /= 2;
  }
  for (size_t iternum = 2; iternum <  vPowExp_.size(); iternum++){
    calculateExp_Qt(iternum, &s, v);

    if (iternum > 2){
      converged = checkIfReachedConvergence(pijt_, pijt_temp);
      if (converged){
        break;
      }
    }
    MatrixTools::copy(pijt_, pijt_temp);
    if (iternum == vPowExp_.size()-1 && !converged){  //need to add more powers to the matrix
      RowMatrix<double> new_pow;
      //new_pow.resize(size_, size_);
      MatrixTools :: mult(vPowExp_[vPowExp_.size()-1], generator_, new_pow);
      vPowExp_.push_back(new_pow);

    }

  }
  while (m > 0){  // recover the 2^m
      
    MatrixTools::mult(pijt_, pijt_, tmpMat_);
    MatrixTools::copy(tmpMat_, pijt_);

    m--;
  }
  return pijt_;

}
/*********************************************************************************/
// double ChromosomeSubstitutionModel::getInitValue(size_t i, int state) const
// {
//   if (i >= size_)
//     throw IndexOutOfBoundsException("ChromosomeSubstitutionModel::getInitValue", i, 0, size_ - 1);
//   if (state < 0 || !alphabet_->isIntInAlphabet(state))
//     throw BadIntException(state, "ChromosomeSubstitutionModel::getInitValue. Character " + alphabet_->intToChar(state) + " is not allowed in model.");
//   vector<int> states = alphabet_->getAlias(state);
//   for (size_t j = 0; j < states.size(); j++)
//   {
//      if (getAlphabetStateAsInt(i) == states[j]){
//        if (dynamic_cast<const IntegerAlphabet*>(alphabet_)){
//          const IntegerAlphabet* alpha = dynamic_cast<const IntegerAlphabet*>(alphabet_);
//          // it is a composite state
//          if (state > alpha->getMax() + 1){
//            return alpha->getProbabilityForState(state, states[j]);

//          }else{
//            return 1.0;
//          }

//        }else{
//          return 1.;
//        }

//      }
//   }
//   return 0.;
// }

const Matrix<double>& ChromosomeSubstitutionModel::getPijt_test(double t) const {
  RowMatrix<double> pijt_temp;
  MatrixTools::getId(size_, pijt_temp);
  double s = 1.0;
  double v = rate_ * t;
  double norm = v * firstNormQ_;
  size_t m = 0;
  bool converged = false;
  //while (v > 0.5)    // exp(r*t*A)=(exp(r*t/(2^m) A))^(2^m)
  while (norm > 0.5)
  {
    m += 1;
    v /= 2;
    norm /= 2;
  }
  for (size_t iternum = 2; iternum <  vPowExp_.size(); iternum++){
    calculateExp_Qt(iternum, &s, v);

    if (iternum > 2){
      converged = checkIfReachedConvergence(pijt_, pijt_temp);
      if (converged){
        break;
      }
    }
    MatrixTools::copy(pijt_, pijt_temp);
    if (iternum > 250){
      //std :: cout << "ERROR: Pijt did not reach convergence for t = "<< t <<"!"<<endl;
      throw Exception("ChromosomeSubstitutionModel: Taylor series did not reach convergence!");
      break;
    }
    if (iternum == vPowExp_.size()-1 && !converged){  //need to add more powers to the matrix
      RowMatrix<double> new_pow;
      //new_pow.resize(size_, size_);
      MatrixTools :: mult(vPowExp_[vPowExp_.size()-1], generator_, new_pow);
      vPowExp_.push_back(new_pow);

    }

  }
  while (m > 0){  // recover the 2^m
      
    MatrixTools::mult(pijt_, pijt_, tmpMat_);
    MatrixTools::copy(tmpMat_, pijt_);

    m--;
  }
  //just for test/////////////////////
  // bool correct = true;
  if (!pijtCalledFromDeriv_){
    for (size_t i = 0; i < size_; i++){
      for (size_t j = 0; j < size_; j++){
        if (pijt_(i,j) < 0){
          pijt_(i,j) = NumConstants::VERY_TINY(); // trying to do it exactly as in ChromEvol. Maybe the "nan" problem will be solved
          //pijt_(i,j) = 0;
        }
        else if (pijt_(i, j) > 1){
          pijt_(i,j) = 1.0;
        }

      }
    }

  }
  return pijt_;

}

/******************************************************************************/