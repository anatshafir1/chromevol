//
// File: CromosomeSubstitutionModel.h
// Created by: Anat Shafir
// Created on: 2020
//


#ifndef CHROMEVOL_CHROMOSOMESUBSTITUTIONMODEL_H
#define CHROMEVOL_CHROMOSOMESUBSTITUTIONMODEL_H

#include <Bpp/Phyl/Model/AbstractSubstitutionModel.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Exceptions.h>
#include <regex>
//#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloatTools.h>

#define lowerBoundOfRateParam 0.0
#define lowerBoundOfExpParam -3.0
#define upperBoundOfRateParam 100.0
#define upperBoundLinearRateParam 5.0
#define upperBoundExpParam 4.6
#define logNormalDomainFactor 4
#define revSigmoidExpRateParam 2
#define IgnoreParam -999
#define DemiEqualDupl -2
#define EPSILON 2.22045e-016
extern int lowerLimitBaseNumber;
using namespace std;
namespace bpp
{
class ChromosomeNumberDependencyFunction{
  protected:
    // the min and max chromosome counts
    int domainMin_;
    int domainMax_;

  public:
    enum FunctionType {CONSTANT, LINEAR, LINEAR_BD, EXP, POLYNOMIAL, LOGNORMAL, REVERSE_SIGMOID, IGNORE};
    ChromosomeNumberDependencyFunction():domainMin_(0), domainMax_(0){}
    virtual ~ChromosomeNumberDependencyFunction(){}

    virtual FunctionType getName() const = 0;
    virtual double getRate(std::vector<Parameter*> params, size_t state) const = 0;
    virtual size_t getNumOfParameters() const = 0;
    virtual void setDomainsIfNeeded(int minChrNum, int maxChrNum){}

    virtual void updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, int maxChrNum){
      std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(paramsNames[index]).getConstraint());
      *lowerBound = interval->getLowerBound();
      *upperBound = interval->getUpperBound();
    }
    virtual void updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound){return;};
    virtual void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
      *lowerBound = lowerBoundOfRateParam;
      *upperBound = upperBoundOfRateParam;
    }
    virtual void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
      *lowerBound = lowerBoundOfRateParam;
      *upperBound = upperBoundOfRateParam;
    }
    virtual double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, int minChrNum, int maxChrNum){
      return parsimonyBound;

    }

};
class ConstantDependencyFunction :
  public virtual ChromosomeNumberDependencyFunction
{
  public:
    ConstantDependencyFunction():ChromosomeNumberDependencyFunction(){}
    virtual ~ConstantDependencyFunction(){}

    FunctionType getName() const{return FunctionType::CONSTANT;}
    double getRate(std::vector<Parameter*> params, size_t state) const;
    size_t getNumOfParameters() const{return 1;}
    double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, int minChrNum, int maxChrNum){
      return parsimonyBound;
    }

};
class LinearDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  public:

    LinearDependencyFunction():ChromosomeNumberDependencyFunction(){}
    virtual ~LinearDependencyFunction(){}

    FunctionType getName() const{return FunctionType::LINEAR;}
    double getRate(std::vector<Parameter*> params, size_t state) const;
    size_t getNumOfParameters() const{return 2;}
    void updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, int maxChrNum);
    void updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound);
    void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
    void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);
    double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, int minChrNum, int maxChrNum);

};
class LinearBDDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  public:
    LinearBDDependencyFunction():ChromosomeNumberDependencyFunction(){}
    virtual ~LinearBDDependencyFunction(){}

    FunctionType getName() const{return FunctionType::LINEAR_BD;}
    double getRate(std::vector<Parameter*> params, size_t state) const;
    size_t getNumOfParameters() const{return 1;}
    //double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, int minChrNum, int maxChrNum);
    double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, int minChrNum, int maxChrNum);

};
class ExponentailDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  public:
    ExponentailDependencyFunction():ChromosomeNumberDependencyFunction(){}
    virtual ~ExponentailDependencyFunction(){}

    FunctionType getName() const {return FunctionType::EXP;}
    double getRate(std::vector<Parameter*> params, size_t state) const;
    size_t getNumOfParameters() const{return 2;}
    void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
    void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);
    double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, int minChrNum, int maxChrNum);

};
class PolynomialDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  public:
    PolynomialDependencyFunction():ChromosomeNumberDependencyFunction(){}
    virtual ~PolynomialDependencyFunction(){}

    FunctionType getName() const{return FunctionType::POLYNOMIAL;}
    double getRate(std::vector<Parameter*> params, size_t state) const;
    size_t getNumOfParameters() const{return 3;}
    void setDomainsIfNeeded(int minChrNum, int maxChrNum){
      domainMin_ = minChrNum;
      domainMax_ = maxChrNum;
    }
    //void updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, int maxChrNum);
    //void updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound);
    void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
    void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);

};
class LognormalDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  private:
  //int maxChrNum_;
  public:
    LognormalDependencyFunction():ChromosomeNumberDependencyFunction(){}
    virtual ~LognormalDependencyFunction(){}

    FunctionType getName() const {return FunctionType::LOGNORMAL;}
    void setDomainsIfNeeded(int minChrNum, int maxChrNum){
      domainMin_ = minChrNum;
      domainMax_ = maxChrNum;
    }

    double getRate(std::vector<Parameter*> params, size_t state) const;
    size_t getNumOfParameters() const{return 3;}
    void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
    void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);

};
class RevSigmoidDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  public:
    RevSigmoidDependencyFunction():ChromosomeNumberDependencyFunction(){}
    virtual ~RevSigmoidDependencyFunction(){}

    FunctionType getName() const {return FunctionType::REVERSE_SIGMOID;}
    double getRate(std::vector<Parameter*> params, size_t state) const;
    size_t getNumOfParameters() const{return 3;}
    void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
    void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);
    void setDomainsIfNeeded(int minChrNum, int maxChrNum){
      domainMin_ = minChrNum;
      domainMax_ = maxChrNum;
    }

};

class ChromosomeSubstitutionModel;


class compositeParameter{
  //public:
    //enum FunctionType {CONSTANT, LINEAR, LINEAR_BD, EXP, POLYNOMIAL, LOGNORMAL, REVERSE_SIGMOID, IGNORE};
    //enum ParamName {BASENUMR, LOSS, GAIN, DUPL, DEMI_DUPL, PARAMNAME_COUNT}; // 24_08 ->use only the substitution model 
    //typedef void (compositeParameter::*functionOp)(size_t, double*, double*);

  private:
    std::vector<Parameter*> params_; // a vector of ChromosomeSubstitutionModel parameters that correspond to a given parameter type
    ChromosomeNumberDependencyFunction* func_; // A function which is used to calculate the bounds and the rate of the composite parameter
    std::string name_;  // The name of the rate parameter. For exampe, "gain" for a paremter that contains gain0, gain1, etc.
    //functionOp updateParamFunc_;
    //size_t size_; // the number of parameters that the composite parameter contains 
    //int maxChrNumber_; // max chromosome number. Used to set the bounds.
    //vector<double> values_;
    

  public:
    // get the values of the composite parameter parameters
    std::vector<double> getParameterValues() const;
    const ChromosomeNumberDependencyFunction* getFunction() const{return func_;}

    // Returns true if the parameter is ignored
    static bool isIgnored(compositeParameter* param){return param == 0;}

    // get the function of the composite parameter
    const ChromosomeNumberDependencyFunction::FunctionType getFuncType() const {return func_->getName();}

    // get the name of the general name for the composite parameter
    // for example, if the parameter is gain, name_ = "gain", and contains actual parameters gain0, gain1, etc.
    const std::string getName() const {return name_;}

    // get the number of parameters
    const size_t getSize() const {return params_.size();}

    // given a function type returns the number of parameters
    //static size_t getNumOfParameters(ChromosomeNumberDependencyFunction* func){return func->getNumOfParameters();} 

    // Returns the value of the overall independent/dependent rate on the number of chromosomes (for any function)
    double getRate(size_t state) const;
    static ChromosomeNumberDependencyFunction* setDependencyFunction(ChromosomeNumberDependencyFunction::FunctionType funcType);
    


    // switches from the enum of ChromosomeSubstitutionModel parameter types to the matching ones in compositeParameter
    //static compositeParameter::ParamName getCompositeRateType(int param);

    // A general function that updates upperBound and lowerBound according to the function to prevent from negative values in the Q matrix.
    // params = the reference to the parameters which include the parameter that should be updated
    // paramNames = the names of the parameters whose values are need in order to set the bounds
    // index = the index of the parameter within the composite parameter that should be updated
    // lowerBound = a pointer towards the lower bound that should be updated
    // upperBound = a pointer towards the upper bound that should be updated
    // funcType = the type of function
    // maxChrNum =  the max chromosome number
    //static void updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, FunctionType funcType, int maxChrNum);

    // A general function that updates upperBound and lowerBound according to the function to prevent from negative values in the Q matrix.
    // Unlike the previous function, this function operates directly on the likelihood instance to ensure matching bounds
    //static void updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound, FunctionType funcType);

    // Functions used by the updateBounds() functions
    //static void updateBoundsLinear(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum, bool random = false);
    //static void updateBoundsExp(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum);
    // static void updateBoundsPolynomial(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum, bool random = false){
    //   throw Exception("Not impelemented yet!");
    // }
    // static void updateBoundsLogNormal(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum, bool random = false){
    //   throw Exception("Not implemented yet!");
    // }
    // static void updateBoundsReverseSigmoid(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum, bool random = false){
    //   throw Exception("Not implemented yet!");
    // }

    // get initial bounds for the random sampling of the parameter values
    //static void getBoundsForInitialParams(FunctionType func, size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber, bool random = false);

    // set the lowest and the highest possible bounds to the model parameters. This important,
    // because the parameters of the substitution model remain constant throughout optimization,
    // and there could be an issue with values set out of bounds in case of the likelihood parameters.
    //static void getAbsoluteBounds(FunctionType func, size_t index, double* lowerBound, double* upperBound, int maxChrNumber);

    static std::vector<std::string> getRelatedParameterNames(ParameterList &params, std::string pattern);
    //const Parameter& getParameter(size_t index){return *(params_[index]);}
    //void setNameStr(ParamName paramName);
    //void getParamUpdateFunction(FunctionType funcType);
    //vector <double>& getValues(size_t index);
    //static void updateBounds(Function* f, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, FunctionType funcType, int maxChrNum);
    //void getBounds(size_t index, double* lowerBound, double* upperBound){return (this->*updateParamFunc_)(index, lowerBound, upperBound);}



  


    compositeParameter(ChromosomeNumberDependencyFunction::FunctionType func, std::string paramName, vector<Parameter*> &params):
      params_(), func_(0), name_(paramName)
    {
      for (size_t i = 0; i < params.size(); i++){
        params_.push_back(params[i]);
      }
      //getParamUpdateFunction(func);
      func_ = setDependencyFunction(func);
    }


    virtual ~compositeParameter(){ delete func_;}



  protected:
    //Parameter* getParameter_(size_t index){return (params_[index]);}
    void setName(std::string name){name_ = name;}
    //void setFunction(FunctionType func){func_ = func;}
    std::vector<Parameter*>& getParams(){return params_;}
    void setParams(std::vector<Parameter*> params){params_ = params;}
    
    

    // functions to update the parameters and their respective intervals
    //void getConstBounds(size_t index, double* lowerBound, double* upperBound);
    //void getLinearBounds(size_t index, double* lowerBound, double* upperBound);
    //void getLinearBDBounds(size_t index, double* lowerBound, double* upperBound);
    //void getExpBounds(size_t index, double* lowerBound, double* upperBound);
   

    // void getPolynomialBounds(size_t index, double* lowerBound, double* upperBound){
    //   throw Exception("Not implemented yet!");
    // }
    // void getLogNormalBounds(size_t index, double* lowerBound, double* upperBound){
    //   throw Exception("Not implemented yet!");
    // }
    // void getReverseSigmoidBounds(size_t index, double* lowerBound, double* upperBound){
    //   throw Exception("Not implemented yet!");
    // }
    friend ChromosomeSubstitutionModel;  
};


class ChromosomeSubstitutionModel :
  public AbstractSubstitutionModel
{
public:
  enum rootFreqType {UNIFORM, ROOT_LL, STATIONARY, FIXED};
  enum rateChangeFunc {LINEAR = 0, EXP = 1};
  enum typeOfTransition {GAIN_T = 0, LOSS_T = 1, DUPL_T = 2, DEMIDUPL_T = 3, BASENUM_T = 4, MAXCHR_T = 5, NUMTYPES = 6, ILLEGAL = 7};
  // All the non-composite parameters should come before the composite ones!!! For example, BASENUM is the first one
  enum paramType {BASENUM = 0, BASENUMR = 1, DUPL = 2, LOSS = 3, GAIN = 4, DEMIDUPL = 5, NUM_OF_CHR_PARAMS = 6};

private:
  compositeParameter* gain_;
  compositeParameter* loss_;
  compositeParameter* dupl_;
  compositeParameter* demiploidy_;
  int baseNum_;
  compositeParameter* baseNumR_;
  unsigned int maxChrRange_;
  rootFreqType freqType_;
  int ChrMinNum_;
  int ChrMaxNum_;
  double firstNormQ_;
  mutable bool pijtCalledFromDeriv_;
  ChromosomeNumberDependencyFunction::FunctionType gainFunc_;
  ChromosomeNumberDependencyFunction::FunctionType lossFunc_;
  ChromosomeNumberDependencyFunction::FunctionType duplFunc_;
  ChromosomeNumberDependencyFunction::FunctionType demiFunc_;
  ChromosomeNumberDependencyFunction::FunctionType baseNumRFunc_;
  bool simulated_;
 


protected:
  mutable std::vector< RowMatrix<double> > vPowExp_;



public:
  ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
    vector<double> gain, 
    vector<double> loss, 
    vector<double> dupl, 
    vector<double> demi,
    int baseNum,
    vector<double> baseNumR,
    unsigned int maxChrRange, 
    rootFreqType freqType,
    vector<int> rateChangeType,
    bool simulated = false);

  ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
    std::map<int, vector<double>> mapOfParamValues,
    int baseNum,
    unsigned int maxChrRange, 
    rootFreqType freqType,
    vector<int> rateChangeType,
    bool simulated = false);


  virtual ~ChromosomeSubstitutionModel() {
    if (gain_ != 0){
      delete gain_;
    }
    if (loss_ != 0){
      delete loss_;
    }
    if (dupl_ != 0){
      delete dupl_;
    }
    if (demiploidy_ != 0){
      delete demiploidy_;
    }
    if (baseNumR_ != 0){
      delete baseNumR_;
    }
  }
  ChromosomeSubstitutionModel(const ChromosomeSubstitutionModel& model):
    AbstractParameterAliasable(model),
    AbstractSubstitutionModel(model),
    gain_(0),
    loss_(0),
    dupl_(0),
    demiploidy_(0),
    baseNum_(model.baseNum_),
    baseNumR_(0),
    maxChrRange_(model.maxChrRange_),
    freqType_(model.freqType_),
    ChrMinNum_(model.ChrMinNum_),
    ChrMaxNum_(model.ChrMaxNum_),
    firstNormQ_(model.firstNormQ_),
    pijtCalledFromDeriv_(model.pijtCalledFromDeriv_),
    gainFunc_(model.gainFunc_),
    lossFunc_(model.lossFunc_),
    duplFunc_(model.duplFunc_),
    demiFunc_(model.demiFunc_),
    baseNumRFunc_(model.baseNumRFunc_),
    simulated_(model.simulated_),
    vPowExp_(model.vPowExp_)
  {
    std::vector<compositeParameter**> newModelParams = {&gain_, &loss_, &dupl_, &demiploidy_, &baseNumR_};
    std::vector<compositeParameter*> originalModelParams = {model.gain_, model.loss_, model.dupl_, model.demiploidy_, model.baseNumR_};
    for (size_t i = 0; i < originalModelParams.size(); i++){
      if (originalModelParams[i] == 0){
        continue;
      }
      vector<Parameter*> params = originalModelParams[i]->getParams();
      vector<Parameter*> newParams;
      for (size_t j = 0; j < params.size(); j++){
        std::string name = params[j]->getName();
        string noPrefixName = std::regex_replace(name, std::regex("Chromosome."), "");
        newParams.push_back(&(getParameter_(noPrefixName)));
      }

      *(newModelParams[i]) = new compositeParameter(originalModelParams[i]->getFuncType(), originalModelParams[i]->getName(), newParams);
      (*(newModelParams[i]))->func_->setDomainsIfNeeded(ChrMinNum_, ChrMaxNum_);

    }

  }


  ChromosomeSubstitutionModel* clone() const { return new ChromosomeSubstitutionModel(*this);}
  

  
public:
  static ChromosomeSubstitutionModel* initRandomModel(
    const ChromosomeAlphabet* alpha,
    int &baseNumber,
    std::map<int, vector<double>> initParams,
    unsigned int chrRange,
    rootFreqType rootFrequenciesType,
    std::vector<int> rateChangeType,
    std::vector<int>& fixedParams,
    double parsimonyBound = 0);





  const Matrix<double>& getPij_t    (double d) const;
  const Matrix<double>& getdPij_dt  (double d) const;
  const Matrix<double>& getd2Pij_dt2(double d) const;
  //The following four functions are just for test.. Finally will bw removed.
  const Matrix<double>& getPij_t_func2(double d) const;
  const Matrix<double>& getPijt_test(double d) const;
  const Matrix<double>& getPij_t_func3(double d) const;
  const Matrix<double>& getPij_t_func4(double d) const;
  static size_t getNumberOfNonCompositeParams(){return 1;} // currently only base number


  std::string getName() const { return "Chromosome"; }
  void setFreq(std::map<int, double>& freqs);
  size_t getNumberOfStates() const { return size_; }
  int getMin() const {return ChrMinNum_;}
  int getMax() const {return ChrMaxNum_;}
  unsigned int getMaxChrRange() const {return maxChrRange_;}
  bool checkIfReachedConvergence(const Matrix<double>& pijt, const Matrix<double>& mt_prev) const;
  //double getInitValue(size_t i, int state) const;
  int getBaseNumber() const {return baseNum_;}
  bool isIgnoredGain() const {return gain_ == 0;}
  bool isIgnoredLoss() const {return loss_ == 0;}
  bool isIgnoredDupl() const {return dupl_ == 0;}
  bool isIgnoredDemiDupl() const {return demiploidy_ == 0;}
  bool isIgnoredBaseNumR() const {return baseNumR_ == 0;}
  const compositeParameter* getDemiDupl() const {return demiploidy_;}
  const compositeParameter* getGain() const {return gain_;}
  const compositeParameter* getLoss() const {return loss_;}
  const compositeParameter* getDupl() const {return dupl_;}
  const compositeParameter* getBaseNumR() const {return baseNumR_;}
  //double getConstDupl () const{return dupl_;}
  //double getChangeRateDupl() const {return duplR_;}
  //double getConstGain() const {return gain_;}
  //double getChangeRateGain() const {return gainR_;}
  //double getConstLoss() const {return loss_;}
  //double getChangeRateLoss() const {return lossR_;}
  //double getBaseNumR() const {return baseNumR_;}
  //double getRate (size_t state, double constRate, double changeRate) const;
  //static void getSetOfFixedParameters(vector<double>& initParams, vector<unsigned int>& fixedParams, map<int, double>& setOfFixedParams);
  //static void getSetOfFixedParameters(vector<double>& initGain, vector<double>& initLoss, vector<double>& initDupl, vector<double>& initDemiDupl, vector<double>& initBaseNumR, int &baseNumber, vector<unsigned int>& fixedParams, map<int, double>& setOfFixedParams);

  
  //These functions should be used from chromsome number optimizer
  void checkParametersBounds() const;
  // this function is needed, because a large range of base number tends to lead to very high chromosome numbers in the simulation
  void correctBaseNumForSimulation(int maxChrNum);


protected:
  void addCompositeParameter(std::vector<Parameter*> parameters);
  void getCompositeParametersValues(std::string paramName, compositeParameter* param);
  void calculatePijtUsingEigenValues(double t) const;
  //static void getRandomParameter(paramType type, double initParamValue, vector<double>& randomParams, double upperBound, double upperBoundLinear, double upperBoundExp, rateChangeFunc rateFunc, int maxChrNum, unsigned int chrRange, map<int, double>& setOfFixedParameters);
  //void updateParameters();
  void updateParameters(vector<double> &gain, vector<double> &loss, vector<double> &dupl, vector<double> &demi, vector<double> &baseNumR);
  //void updateLinearParameters();
  //void updateExpParameters();
  //void updateBaseNumParameters(std::shared_ptr<IntervalConstraint> interval);
  void updateMatrices();
  void updateQWithBaseNumParameters(size_t currChrNum, size_t minChrNum, size_t maxChrNum);
  void updateQWithGain(size_t currChrNum, size_t minChrNum);
  void updateQWithLoss(size_t currChrNum, size_t minChrNum);
  void updateQWithDupl(size_t currChrNum, size_t minChrNum, size_t maxChrNum = 0);
  void updateQWithDemiDupl(size_t currChrNum, size_t minChrNum, size_t maxChrNum);
  void updateEigenMatrices();
  void getParametersValues();
  void calculateExp_Qt(size_t pow, double s, size_t m, double v) const;
  void calculateExp_Qt(size_t pow, double* s, double v) const;
  double getFirstNorm() const;
  double get_epsilon() const{ return 0.0001;};
  //void updateConstRateParameter(double paramValueConst, double paramValueChange, string parameterName, std::shared_ptr<IntervalConstraint> interval);
  //void updateLinearChangeParameter(double paramValueConst, double paramValueChange, string parameterName);
  //void setNewBoundsForLinearParameters(double &constRate, double &changeRate, string paramNameConst, string paramNameLinear);
  std::vector<Parameter*> createCompositeParameter(ChromosomeNumberDependencyFunction::FunctionType &func, std::string paramName, vector<double> &vectorOfValues);
  friend class compositeParameter;
};
} // end of namespace bpp.

#endif  // CHROMEVOL_CHROMOSOMESUBSTITUTIONMODEL_H