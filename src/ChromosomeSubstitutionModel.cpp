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
// CompositeParameter *** *** *** **** *** *** *** *** **** *** *** *** *** 
/******************************************************************************/
double compositeParameter::getRate(double state) const{
  return func_->getRate(params_, state);
} 
/******************************************************************************/
ChromosomeNumberDependencyFunction* compositeParameter::setDependencyFunction(ChromosomeNumberDependencyFunction::FunctionType funcType, bool continuous){
  switch (funcType)
  {
  case ChromosomeNumberDependencyFunction::CONSTANT:
    return new ConstantDependencyFunction(continuous);
  case ChromosomeNumberDependencyFunction::LINEAR:
    return new LinearDependencyFunction(continuous);
  case ChromosomeNumberDependencyFunction::EXP:
    return new ExponentailDependencyFunction(continuous);
  case ChromosomeNumberDependencyFunction::LINEAR_BD:
    return new LinearBDDependencyFunction(continuous);
  case ChromosomeNumberDependencyFunction::LOGNORMAL:
    return new LognormalDependencyFunction(continuous);
  case ChromosomeNumberDependencyFunction::POLYNOMIAL:
    return new PolynomialDependencyFunction(continuous);
  case ChromosomeNumberDependencyFunction::REVERSE_SIGMOID:
    return new RevSigmoidDependencyFunction(continuous);
  case ChromosomeNumberDependencyFunction::LOGITNORMAL:
    return new LogitnormalDependencyFunction(continuous);
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
    defineFunctionsNames(rateChangeType);
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
  defineFunctionsNames(rateChangeType);
  updateParameters(mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::GAIN)], 
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::LOSS)], 
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DUPL)],
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)],
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)]);
  computeFrequencies(false);
  isScalable_ = false;    //in ChromEvol the matrix should be not normalized
  updateMatrices();


}
/**************************************************************************** */
// the following two constructors are for derived models (classes)
/**************************************************************************** */
ChromosomeSubstitutionModel::ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
  std::map<int, vector<double>> mapOfParamValues,
  int baseNum,
  unsigned int chrRange, 
  rootFreqType freqType,
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
    
{}
/**************************************************************************** */
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
    
{}

/******************************************************************************/
void ChromosomeSubstitutionModel::defineFunctionsNames(vector<int> &rateChangeType){
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

}
/******************************************************************************/
void ChromosomeSubstitutionModel::setRandomModel(
  const ChromosomeAlphabet* alpha,
  int &baseNumber,
  map<int, vector<double>> initParams,
  unsigned int chrRange,
  rootFreqType rootFrequenciesType,
  vector<int> rateChangeType,
  vector<int>& fixedParams,
  bool demiOnlyForEven,
  double parsimonyBound,
  std::map<int, vector<double>> &mapRandomParams,
  int &newBaseNumber,
  bool brownianModel,
  double minDomain,
  double maxDomain)
{
  newBaseNumber = baseNumber;
  size_t startCompositeParams = getNumberOfNonCompositeParams();
  double minFuncDomain;
  double maxFuncDomain;
  if (brownianModel){
    minFuncDomain = minDomain;
    maxFuncDomain = maxDomain;
  }else{
    minFuncDomain = static_cast<double>(alpha->getMin());
    maxFuncDomain = static_cast<double>(alpha->getMax());
  }

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
        if (initParams[i].size() != 0){
          ChromosomeNumberDependencyFunction::FunctionType funcType = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startCompositeParams]);
          ChromosomeNumberDependencyFunction* functionOp = compositeParameter::setDependencyFunction(funcType, brownianModel);
          functionOp->setDomainsIfNeeded(minFuncDomain, maxFuncDomain);

          auto numOfParameters = functionOp->getNumOfParameters();
          for (size_t j = 0; j < numOfParameters; j++){
            functionOp->getBoundsForInitialParams(j, paramValues, &lowerBound, &upperBound, maxFuncDomain);
            double upperBoundCandidate = functionOp->getParsimonyBound(paramValues, parsimonyBound, j, minFuncDomain, maxFuncDomain);
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
  int newBaseNumber;
  setRandomModel(alpha, baseNumber, initParams, chrRange, rootFrequenciesType, rateChangeType, fixedParams, demiOnlyForEven, parsimonyBound, mapRandomParams, newBaseNumber, false);


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
void ChromosomeSubstitutionModel::updateParameters(vector<double> &gain, vector<double> &loss, vector<double> &dupl, vector<double> &demi, vector<double> &baseNumR, bool continuous){
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
    gain_ = new compositeParameter(gainFunc_, "gain", gainParams, continuous);

  }
  if (lossFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    loss_ = new compositeParameter(lossFunc_, "loss", lossParams, continuous);
  }
  if (duplFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    dupl_ = new compositeParameter(duplFunc_, "dupl", duplParams, continuous);
  }
  if (demiFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    demiploidy_ = new compositeParameter(demiFunc_, "demi", demiParams, continuous);
  }
  if (baseNumRFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    baseNumR_ = new compositeParameter(baseNumRFunc_, "baseNumR", baseNumParams, continuous);
  }
  setAllFunctionsDomains();


}
/******************************************************************************/
void ChromosomeSubstitutionModel::setAllFunctionsDomains(){
  if (gainFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    gain_->func_->setDomainsIfNeeded(ChrMinNum_, ChrMaxNum_);

  }
  if (lossFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    loss_->func_->setDomainsIfNeeded(ChrMinNum_, ChrMaxNum_);
  }
  if (duplFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    dupl_->func_->setDomainsIfNeeded(ChrMinNum_, ChrMaxNum_);

  }
  if (demiFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    demiploidy_->func_->setDomainsIfNeeded(ChrMinNum_, ChrMaxNum_);

  }
  if (baseNumRFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
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
      if (baseNumR_->getRate((double)i) < 0){
        throw Exception("ChromosomeSubstitutionModel::correctBaseNumForSimulation():Negative base number rate!");
      }
      for (size_t j = i + 1; j < maxChrNum + 1; j ++){
        if (j == maxChrNum){
          if (((j-i) <= maxChrRange_) && ((int)(j-i) > baseNum_)){
            generator_(i-minChrNum, maxChrNum-minChrNum) -= baseNumR_->getRate((double)i);
          }


        }else{
          if ((j-i) % baseNum_ == 0){
            if (i > (size_t)maxChrNumInferred){
              if (((j-i) <= maxChrRange_) && ((int)(j-i) > baseNum_)){
                generator_(i - minChrNum, j - minChrNum) -= baseNumR_->getRate((double)i);
              }
            }else{
              if (((int)j-maxChrNumInferred) >= baseNum_ ){
                if ((j-i) <= maxChrRange_){
                  if (((int)i != maxChrNumInferred) || ((int)j-maxChrNumInferred != baseNum_)){
                    generator_(i - minChrNum, j - minChrNum) -= baseNumR_->getRate((double)i);

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
    if (demiploidy_->getRate((double)i) < 0){
      throw Exception("ChromosomeSubstitutionModel::updateQWithDemiDupl(): Negative demiploidy rate!");
    }
    if (i % 2 == 0 && (double)i * 1.5 <= (double)maxChrNum){

      generator_(i-minChrNum, (size_t)((double)i * 1.5)-minChrNum) += demiploidy_->getRate((double)i);

                        
    }else if (i % 2 != 0 && (size_t)ceil((double)i*1.5) <= maxChrNum){
      if (!demiOnlyForEven_){
        if (i == 1){
          generator_(i-minChrNum, (size_t)ceil((double)i * 1.5)-minChrNum) += demiploidy_->getRate((double)i);
        }else{
          generator_(i-minChrNum, (size_t)ceil((double)i * 1.5)-minChrNum) += demiploidy_->getRate((double)i)/2;
          generator_(i-minChrNum, (size_t)floor((double)i * 1.5)-minChrNum) += demiploidy_->getRate((double)i)/2;

        }

      }


    }else{
      if (i != maxChrNum){
        generator_(i-minChrNum, maxChrNum-minChrNum) += demiploidy_->getRate((double)i);
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
  double gainRate = gain_->getRate((double)i);
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
  double lossRate = loss_->getRate((double)i);
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
  double duplRate = dupl_->getRate((double)i);
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
  if ( baseNumR_->getRate((double)currChrNum) < 0){
    throw Exception("ChromosomeSubstitutionModel::updateQWithBaseNumParameters():Negative base number rate!");
  }
  for (size_t j = currChrNum + 1; j < maxChrNum + 1; j ++){
    if (j == maxChrNum){
      if ((j-currChrNum) <= maxChrRange_){
        generator_(currChrNum-minChrNum, maxChrNum-minChrNum) += baseNumR_->getRate((double)currChrNum);
      }
    }else{
      if ((j-currChrNum) % baseNum_ == 0){
        if ((j-currChrNum) <= maxChrRange_){
          generator_(currChrNum - minChrNum, j - minChrNum) += baseNumR_->getRate((double)currChrNum);
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
#ifdef USE_VERSION_EIGEN_PIJT

const Matrix<double>& ChromosomeSubstitutionModel::getPij_t(double t) const
{
  size_t minTaylorIterations = 5;
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

      if (iternum > minTaylorIterations){
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
  return pijt_;
}

#else

const Matrix<double>& ChromosomeSubstitutionModel::getPij_t(double t) const
{
  size_t minTaylorIterations = 5;
  if (t == 0)
  {
    MatrixTools::getId(size_, pijt_);
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

      if (iternum > minTaylorIterations){
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
  return pijt_;
}
#endif
/******************************************************************************/
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
/************************************************************************************************/
// ChromosomeBMSubstitutionModel: For chromosome models with a continuous brownian motion trait
/************************************************************************************************/
ChromosomeBMSubstitutionModel::ChromosomeBMSubstitutionModel(double mu,
  double sigma,
  double state,
  double minTraitState,
  double maxTraitState,
  double minTraitStateInData,
  double maxTraitStateInData,
  const ChromosomeAlphabet* alpha, 
  vector<double> gain, 
  vector<double> loss, 
  vector<double> dupl, 
  vector<double> demi,
  int baseNum,
  vector<double> baseNumR,
  unsigned int maxChrRange, 
  rootFreqType freqType,
  vector<int> rateChangeType,
  bool demiOnlyForEven,
  bool simulated):
  AbstractParameterAliasable("Chromosome."),
  ChromosomeSubstitutionModel(alpha, gain, loss, dupl, demi, baseNum, baseNumR, maxChrRange, freqType, demiOnlyForEven, simulated),
  mu_(mu),
  sigma_(sigma),
  state_(state),
  minTraitState_(minTraitState),
  maxTraitState_(maxTraitState){
    defineFunctionsNames(rateChangeType);
    if (simulated){
      if (mu > maxTraitStateInData){
        maxTraitStateInData = mu; 
      }else if (mu < minTraitStateInData){
        minTraitStateInData = mu;
      }
    }
    updateBMParameters(minTraitStateInData, maxTraitStateInData);
    updateParameters(gain, loss, dupl, demi, baseNumR, true); 
    computeFrequencies(false);
    isScalable_ = false;    //in ChromEvol the matrix should be not normalized
    updateMatrices();
}


ChromosomeBMSubstitutionModel::ChromosomeBMSubstitutionModel(double mu,
  double sigma,
  double state,
  double minTraitState,
  double maxTraitState,
  double minTraitStateInData,
  double maxTraitStateInData,
  const ChromosomeAlphabet* alpha, 
  std::map<int, vector<double>> mapOfParamValues,
  int baseNum,
  unsigned int maxChrRange, 
  rootFreqType freqType,
  vector<int> rateChangeType,
  bool demiOnlyForEven,
  bool simulated):
  AbstractParameterAliasable("Chromosome."),
  ChromosomeSubstitutionModel(alpha, mapOfParamValues, baseNum, maxChrRange, freqType, demiOnlyForEven, simulated),
  mu_(mu),
  sigma_(sigma),
  state_(state),
  minTraitState_(minTraitState),
  maxTraitState_(maxTraitState){
    defineFunctionsNames(rateChangeType);
    if (simulated){
      if (mu > maxTraitStateInData){
        maxTraitStateInData = mu; 
      }else if (mu < minTraitStateInData){
        minTraitStateInData = mu;
      }
    }
    updateBMParameters(minTraitStateInData, maxTraitStateInData);
    updateParameters(mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::GAIN)], 
      mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::LOSS)], 
      mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DUPL)],
      mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)],
      mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)], true);
    computeFrequencies(false);
    isScalable_ = false;    //in ChromEvol the matrix should be not normalized
    updateMatrices();
}
void ChromosomeBMSubstitutionModel::setAllFunctionsDomains(){
  if (gainFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    gain_->func_->setDomainsIfNeeded(minTraitState_, maxTraitState_);
  }
  if (lossFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    loss_->func_->setDomainsIfNeeded(minTraitState_, maxTraitState_);
  }
  if (duplFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    dupl_->func_->setDomainsIfNeeded(minTraitState_, maxTraitState_);
  }
  if (demiFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    demiploidy_->func_->setDomainsIfNeeded(minTraitState_, maxTraitState_);
  }
  if(baseNumRFunc_ != ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
    baseNumR_->func_->setDomainsIfNeeded(minTraitState_, maxTraitState_);
  }

}
/******************************************************************************/
bool ChromosomeBMSubstitutionModel::areRatesNegative(double minTraitState, double maxTraitState, map<int, vector<double>> &params, vector<int> &rateChangeType){
  vector<double> paramValues;
  for (size_t i = 1; i < ChromosomeSubstitutionModel::paramType::NUM_OF_CHR_PARAMS; i++){
    switch (i)
    {
      case ChromosomeSubstitutionModel::BASENUM:
        break;
      case ChromosomeSubstitutionModel::GAIN:
        paramValues = params[ChromosomeSubstitutionModel::GAIN];
        break;
      case ChromosomeSubstitutionModel::LOSS:
        paramValues = params[ChromosomeSubstitutionModel::LOSS];
        break;
      case ChromosomeSubstitutionModel::DUPL:
        paramValues = params[ChromosomeSubstitutionModel::DUPL];
        break;
      case ChromosomeSubstitutionModel::DEMIDUPL:
        paramValues = params[ChromosomeSubstitutionModel::DEMIDUPL];
        break;
      case ChromosomeSubstitutionModel::BASENUMR:
        paramValues = params[ChromosomeSubstitutionModel::BASENUMR];
        break;
   
      default:
        throw Exception("ChromEvolOptions::setFunctions: parameter not found !!!");
    }
    size_t startNonComposite = getNumberOfNonCompositeParams();
    if (rateChangeType[i-startNonComposite] == ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
      continue;
    }
    auto funcType = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
    ChromosomeNumberDependencyFunction* func = compositeParameter::setDependencyFunction(funcType, true);
    func->setDomainsIfNeeded(minTraitState, maxTraitState);
    double rateMin;
    double rateMax;
    try {
      rateMin = func->getRate(paramValues, minTraitState+1);
    } catch (const std::exception& e) {
      delete func;
      return true;
    }
    try{
      rateMax = func->getRate(paramValues, maxTraitState+1);
    }catch (const std::exception& e) {
      delete func;
      return true;
    }
    delete func;
    

  }
  return false;

}
/**************************************************************************** */
void ChromosomeBMSubstitutionModel::getRandomState(std::map<int, vector<double>> &params, vector<int> &rateChangeType, double &state){
  vector<double> paramValues;
  state = 0;
  for (size_t i = 1; i < ChromosomeSubstitutionModel::paramType::NUM_OF_CHR_PARAMS; i++){
    switch (i)
    {
      case ChromosomeSubstitutionModel::BASENUM:
        break;
      case ChromosomeSubstitutionModel::GAIN:
        paramValues = params[ChromosomeSubstitutionModel::GAIN];
        break;
      case ChromosomeSubstitutionModel::LOSS:
        paramValues = params[ChromosomeSubstitutionModel::LOSS];
        break;
      case ChromosomeSubstitutionModel::DUPL:
        paramValues = params[ChromosomeSubstitutionModel::DUPL];
        break;
      case ChromosomeSubstitutionModel::DEMIDUPL:
        paramValues = params[ChromosomeSubstitutionModel::DEMIDUPL];
        break;
      case ChromosomeSubstitutionModel::BASENUMR:
        paramValues = params[ChromosomeSubstitutionModel::BASENUMR];
        break;
   
      default:
        throw Exception("ChromEvolOptions::setFunctions: parameter not found !!!");
    }
    size_t startNonComposite = getNumberOfNonCompositeParams();
    if (rateChangeType[i-startNonComposite] == ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
      continue;
    }
    auto funcType = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[i-startNonComposite]);
    if (funcType == ChromosomeNumberDependencyFunction::FunctionType::LINEAR){
      double rootPoint = -paramValues[0]/paramValues[1];
      if (rootPoint > state){
        state = rootPoint;
      }

    }
  }
 
}
/******************************************************************************/
void ChromosomeBMSubstitutionModel::fixInitialStates(std::map<string, vector<double>> &newParamValues){
  const compositeParameter* demiToUpdate = getDemiDupl();
  const compositeParameter* gainToUpdate = getGain();
  const compositeParameter* lossToUpdate = getLoss();
  const compositeParameter* duplToUpdate = getDupl();
  const compositeParameter* baseNumToUpdate = getBaseNumR();
  if (demiToUpdate != 0){
    fixInitialStates(demiToUpdate->getParameterValues(), newParamValues["demi"], demiToUpdate);
  }
  if (gainToUpdate != 0){
    fixInitialStates(gainToUpdate->getParameterValues(), newParamValues["gain"], gainToUpdate);
  }
  if (lossToUpdate != 0){
    fixInitialStates(lossToUpdate->getParameterValues(), newParamValues["loss"], lossToUpdate);
  }
  if(duplToUpdate != 0){
    fixInitialStates(duplToUpdate->getParameterValues(), newParamValues["dupl"], duplToUpdate);
  }
  if (baseNumToUpdate != 0){
    fixInitialStates(baseNumToUpdate->getParameterValues(), newParamValues["baseNumR"], baseNumToUpdate);
  }


}
/******************************************************************************/
void ChromosomeBMSubstitutionModel::fixInitialStates(const vector<double> &oldParamValues, vector<double> &newParamValues, const compositeParameter* paramToUpdate){
  if (paramToUpdate->getFunction()->getName() == ChromosomeNumberDependencyFunction::FunctionType::LINEAR){
    double rootPoint1 = 0;
    double rootPoint2 = 0;
    double rootPoint3 = 0;
    if (oldParamValues[1] != 0){
      rootPoint2 = -newParamValues[0]/oldParamValues[1];
      
    }
    if (newParamValues[1] != 0){
      rootPoint1= -oldParamValues[0]/newParamValues[1];
      rootPoint3 = -newParamValues[0]/newParamValues[1];
      
    }
    
    double rootPoint = std::max(state_, std::max(std::max(rootPoint1, rootPoint2), rootPoint3));
    if (rootPoint > state_){
      setParameterValue("state", rootPoint);
    }
  }

}
/******************************************************************************/
std::vector<Parameter*> ChromosomeBMSubstitutionModel::createCompositeParameter(ChromosomeNumberDependencyFunction::FunctionType &func, std::string paramName, vector<double> &vectorOfValues){
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
    ChromosomeNumberDependencyFunction* functionOp =  compositeParameter::setDependencyFunction(func, true);
    functionOp->setDomainsIfNeeded(minTraitState_, maxTraitState_);
    functionOp->getAbsoluteBounds(i, &lowerBound, &upperBound, maxTraitState_);
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
/*******************************************************************************/
void ChromosomeBMSubstitutionModel::updateQWithDemiDupl(size_t i, size_t minChrNum, size_t maxChrNum){
  //double demiploidy;
  
  if (demiploidy_ != 0){
    if (demiploidy_->getRate(state_+1) < 0){
      throw Exception("ChromosomeSubstitutionModel::updateQWithDemiDupl(): Negative demiploidy rate!");
    }
    if (i % 2 == 0 && (double)i * 1.5 <= (double)maxChrNum){

      generator_(i-minChrNum, (size_t)((double)i * 1.5)-minChrNum) += demiploidy_->getRate(state_+1);

                        
    }else if (i % 2 != 0 && (size_t)ceil((double)i*1.5) <= maxChrNum){
      if (!demiOnlyForEven_){
        if (i == 1){
          generator_(i-minChrNum, (size_t)ceil((double)i * 1.5)-minChrNum) += demiploidy_->getRate(state_+1);
        }else{
          generator_(i-minChrNum, (size_t)ceil((double)i * 1.5)-minChrNum) += demiploidy_->getRate(state_+1)/2;
          generator_(i-minChrNum, (size_t)floor((double)i * 1.5)-minChrNum) += demiploidy_->getRate(state_+1)/2;

        }

      }


    }else{
      if (i != maxChrNum){
        generator_(i-minChrNum, maxChrNum-minChrNum) += demiploidy_->getRate(state_+1);
      }

    }

  }
}
/*******************************************************************************/
void ChromosomeBMSubstitutionModel::updateQWithGain(size_t i, size_t minChrNum){
  if (gain_ == 0){
    return;
  }
  double gainRate = gain_->getRate(state_+1);
  if (gainRate < 0){
    throw Exception ("ChromosomeSubstitutionModel::updateQWithGain(): negative gain rate!");
  }
  generator_(i-minChrNum, i+1-minChrNum) += gainRate;


}
/*******************************************************************************/
void ChromosomeBMSubstitutionModel::updateQWithLoss(size_t i, size_t minChrNum){
  //generator_(i-minChrNum, i-1-minChrNum) = loss_ + (lossR_* i);
  if (loss_ == 0){
    return;
  }
  double lossRate = loss_->getRate(state_+1);
  if (lossRate < 0){
    throw Exception ("ChromosomeSubstitutionModel::updateQWithLoss(): negative loss rate!");
  }
  generator_(i-minChrNum, i-1-minChrNum) += lossRate;


}
/*******************************************************************************/
void ChromosomeBMSubstitutionModel::updateQWithDupl(size_t i, size_t minChrNum, size_t maxChrNum){
  if (dupl_ == 0){
    return;
  }
  // if the transition is not to maxChr
  double duplRate = dupl_->getRate(state_+1);
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
void ChromosomeBMSubstitutionModel::updateQWithBaseNumParameters(size_t currChrNum, size_t minChrNum, size_t maxChrNum){
  if ( baseNumR_->getRate(state_+1) < 0){
    throw Exception("ChromosomeSubstitutionModel::updateQWithBaseNumParameters():Negative base number rate!");
  }
  for (size_t j = currChrNum + 1; j < maxChrNum + 1; j ++){
    if (j == maxChrNum){
      if ((j-currChrNum) <= maxChrRange_){
        generator_(currChrNum-minChrNum, maxChrNum-minChrNum) += baseNumR_->getRate(state_+1);
      }
    }else{
      if ((j-currChrNum) % baseNum_ == 0){
        if ((j-currChrNum) <= maxChrRange_){
          generator_(currChrNum - minChrNum, j - minChrNum) += baseNumR_->getRate(state_+1);
        }
      }
    }
  }
}


/*******************************************************************************/
void ChromosomeBMSubstitutionModel::correctBaseNumForSimulation(int maxChrNumInferred){
    //update generator matrix
    size_t maxChrNum = (size_t)(getMax());
    size_t minChrNum = (size_t)(getMin()); 
 
    // updating Q matrix
    for (size_t i = minChrNum; i < maxChrNum; i++){
      if (baseNumR_->getRate((double)i) < 0){
        throw Exception("ChromosomeSubstitutionModel::correctBaseNumForSimulation():Negative base number rate!");
      }
      for (size_t j = i + 1; j < maxChrNum + 1; j ++){
        if (j == maxChrNum){
          if (((j-i) <= maxChrRange_) && ((int)(j-i) > baseNum_)){
            generator_(i-minChrNum, maxChrNum-minChrNum) -= baseNumR_->getRate(state_+1);
          }


        }else{
          if ((j-i) % baseNum_ == 0){
            if (i > (size_t)maxChrNumInferred){
              if (((j-i) <= maxChrRange_) && ((int)(j-i) > baseNum_)){
                generator_(i - minChrNum, j - minChrNum) -= baseNumR_->getRate(state_+1);
              }
            }else{
              if (((int)j-maxChrNumInferred) >= baseNum_ ){
                if ((j-i) <= maxChrRange_){
                  if (((int)i != maxChrNumInferred) || ((int)j-maxChrNumInferred != baseNum_)){
                    generator_(i - minChrNum, j - minChrNum) -= baseNumR_->getRate(state_+1);

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

void ChromosomeBMSubstitutionModel::getParametersValues(){
    sigma_ = getParameterValue("sigma");
    mu_ = getParameterValue("mu");
    state_ = getParameterValue("state");
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
void ChromosomeBMSubstitutionModel::updateBMParameters(double minTraitState, double maxTraitState){
  std::shared_ptr<IntervalConstraint> interval_sigma = make_shared<IntervalConstraint>(0.001, 10000, false, true);
  addParameter_(new Parameter("Chromosome.sigma", sigma_, interval_sigma));
  std::shared_ptr<IntervalConstraint> interval_mu = make_shared<IntervalConstraint>(minTraitState, maxTraitState, true, true); // don't really know what bound should I put here
  addParameter_(new Parameter("Chromosome.mu", mu_, interval_mu));
  std::shared_ptr<IntervalConstraint> interval_state = make_shared<IntervalConstraint>(-1000000, 1000000, false, true);
  addParameter_(new Parameter("Chromosome.state", state_, interval_state));

}

/******************************************************************************/
ChromosomeBMSubstitutionModel* ChromosomeBMSubstitutionModel::initBMRandomModel(
  const ChromosomeAlphabet* alpha,
  int &baseNumber,
  map<int, vector<double>> initParams,
  unsigned int chrRange,
  rootFreqType rootFrequenciesType,
  vector<int> rateChangeType,
  vector<int>& fixedParams,
  bool demiOnlyForEven,
  double parsimonyBound,
  double minTraitState,
  double maxTraitState,
  double minTraitStateInData,
  double maxTraitStateInData,
  double sigmaRoughEstimator)
{
  std::map<int, vector<double>> mapRandomParams;
  int newBaseNumber;
  double mu =  RandomTools::giveRandomNumberBetweenTwoPoints(minTraitStateInData + 0.001, maxTraitStateInData-0.001);
  double sigma = RandomTools::giveRandomNumberBetweenTwoPoints(0.001, sigmaRoughEstimator);
  setRandomModel(alpha, baseNumber, initParams, chrRange, rootFrequenciesType, rateChangeType, fixedParams, demiOnlyForEven, parsimonyBound, mapRandomParams, newBaseNumber, true, minTraitState, maxTraitState);


  ChromosomeBMSubstitutionModel* model = new ChromosomeBMSubstitutionModel(mu, sigma, 0, minTraitState, maxTraitState, minTraitStateInData, maxTraitStateInData, alpha, mapRandomParams, newBaseNumber, chrRange, rootFrequenciesType, rateChangeType, demiOnlyForEven);//, useExtendedFloat);
  return model;

}