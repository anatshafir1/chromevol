// File: CromosomeSubstitutionModel.h
// Created by: Anat Shafir
// Created on: 2020
//


#ifndef CHROMEVOL_CHROMOSOMENUMBERDEPENDENCYFUNCTION_H
#define CHROMEVOL_CHROMOSOMENUMBERDEPENDENCYFUNCTION_H

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
using namespace std;
namespace bpp
{
class ChromosomeNumberDependencyFunction{
  protected:
    // the min and max chromosome counts
    double domainMin_;
    double domainMax_;
    bool continuous_;

  public:
    enum FunctionType {CONSTANT, LINEAR, LINEAR_BD, EXP, POLYNOMIAL, LOGNORMAL, REVERSE_SIGMOID, LOGITNORMAL, IGNORE};
    ChromosomeNumberDependencyFunction():domainMin_(0), domainMax_(0), continuous_(false){}
    ChromosomeNumberDependencyFunction(bool continuous):domainMin_(0), domainMax_(0), continuous_(continuous){}
    virtual ~ChromosomeNumberDependencyFunction(){}

    virtual FunctionType getName() const = 0;
    virtual double getRate(std::vector<Parameter*> params, double state) const = 0;
    virtual size_t getNumOfParameters() const = 0;
    virtual void setDomainsIfNeeded(double minChrNum, double maxChrNum){}
    double getMinDomain(){return domainMin_;}
    double getMaxDomain(){return domainMax_;}
    

    virtual void updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, double maxChrNum){
      std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(paramsNames[index]).getConstraint());
      *lowerBound = interval->getLowerBound();
      *upperBound = interval->getUpperBound();
    }
    virtual void updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound){return;};
    virtual void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber){
      *lowerBound = lowerBoundOfRateParam;
      *upperBound = upperBoundOfRateParam;
    }
    virtual void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber){
      *lowerBound = lowerBoundOfRateParam;
      *upperBound = upperBoundOfRateParam;
    }
    virtual double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, double minChrNum, double maxChrNum){
      return parsimonyBound;

    }

};
class ConstantDependencyFunction :
  public virtual ChromosomeNumberDependencyFunction
{
  public:
    ConstantDependencyFunction():ChromosomeNumberDependencyFunction(){}
    ConstantDependencyFunction(bool continuous):ChromosomeNumberDependencyFunction(continuous){}
    virtual ~ConstantDependencyFunction(){}

    FunctionType getName() const{return FunctionType::CONSTANT;}
    double getRate(std::vector<Parameter*> params, double state) const;
    size_t getNumOfParameters() const{return 1;}
    double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, double minChrNum, double maxChrNum){
      return parsimonyBound;
    }

};
class LinearDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  public:

    LinearDependencyFunction():ChromosomeNumberDependencyFunction(){}
    LinearDependencyFunction(bool continuous):ChromosomeNumberDependencyFunction(continuous){}
    virtual ~LinearDependencyFunction(){}

    FunctionType getName() const{return FunctionType::LINEAR;}
    double getRate(std::vector<Parameter*> params, double state) const;
    size_t getNumOfParameters() const{return 2;}
    void setDomainsIfNeeded(double minChrNum, double maxChrNum){
      domainMin_ = minChrNum;
      domainMax_ = maxChrNum;
    }
    void updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, double maxChrNum);
    void updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound);
    void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber);
    void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber);
    double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, double minChrNum, double maxChrNum);

};
class LinearBDDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  public:
    LinearBDDependencyFunction():ChromosomeNumberDependencyFunction(){}
    LinearBDDependencyFunction(bool continuous):ChromosomeNumberDependencyFunction(continuous){}
    virtual ~LinearBDDependencyFunction(){}

    FunctionType getName() const{return FunctionType::LINEAR_BD;}
    double getRate(std::vector<Parameter*> params, double state) const;
    size_t getNumOfParameters() const{return 1;}
    //double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, int minChrNum, int maxChrNum);
    double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, double minChrNum, double maxChrNum);

};
class ExponentailDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  public:
    ExponentailDependencyFunction():ChromosomeNumberDependencyFunction(){}
    ExponentailDependencyFunction(bool continuous):ChromosomeNumberDependencyFunction(continuous){}
    virtual ~ExponentailDependencyFunction(){}

    FunctionType getName() const {return FunctionType::EXP;}
    double getRate(std::vector<Parameter*> params, double state) const;
    void setDomainsIfNeeded(double minChrNum, double maxChrNum){
      domainMin_ = minChrNum;
      domainMax_ = maxChrNum;
    }
    size_t getNumOfParameters() const{return 2;}
    void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber);
    void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber);
    double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, double minChrNum, double maxChrNum);

};
class PolynomialDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  public:
    PolynomialDependencyFunction():ChromosomeNumberDependencyFunction(){}
    PolynomialDependencyFunction(bool continuous):ChromosomeNumberDependencyFunction(continuous){}
    virtual ~PolynomialDependencyFunction(){}

    FunctionType getName() const{return FunctionType::POLYNOMIAL;}
    double getRate(std::vector<Parameter*> params, double state) const;
    size_t getNumOfParameters() const{return 3;}
    void setDomainsIfNeeded(double minChrNum, double maxChrNum){
      domainMin_ = minChrNum;
      domainMax_ = maxChrNum;
    }
    //void updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, int maxChrNum);
    //void updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound);
    void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber);
    void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber);

};
class LognormalDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  private:
  //int maxChrNum_;
  public:
    LognormalDependencyFunction():ChromosomeNumberDependencyFunction(){}
    LognormalDependencyFunction(bool continuous):ChromosomeNumberDependencyFunction(continuous){}
    virtual ~LognormalDependencyFunction(){}

    FunctionType getName() const {return FunctionType::LOGNORMAL;}
    void setDomainsIfNeeded(double minChrNum, double maxChrNum){
      domainMin_ = minChrNum;
      domainMax_ = maxChrNum;
    }

    double getRate(std::vector<Parameter*> params, double state) const;
    size_t getNumOfParameters() const{return 3;}
    void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber);
    void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber);

};
class RevSigmoidDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  public:
    RevSigmoidDependencyFunction():ChromosomeNumberDependencyFunction(){}
    RevSigmoidDependencyFunction(bool continuous):ChromosomeNumberDependencyFunction(continuous){}
    virtual ~RevSigmoidDependencyFunction(){}

    FunctionType getName() const {return FunctionType::REVERSE_SIGMOID;}
    double getRate(std::vector<Parameter*> params, double state) const;
    size_t getNumOfParameters() const{return 3;}
    void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber);
    void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber);
    void setDomainsIfNeeded(double minChrNum, double maxChrNum){
      domainMin_ = minChrNum;
      domainMax_ = maxChrNum;
    }

};
class LogitnormalDependencyFunction:
  public virtual ChromosomeNumberDependencyFunction
{
  private:
  //int maxChrNum_;
  public:
    LogitnormalDependencyFunction():ChromosomeNumberDependencyFunction(){}
    LogitnormalDependencyFunction(bool continuous):ChromosomeNumberDependencyFunction(continuous){}
    virtual ~LogitnormalDependencyFunction(){}

    FunctionType getName() const {return FunctionType::LOGITNORMAL;}
    void setDomainsIfNeeded(double minChrNum, double maxChrNum){
      domainMin_ = minChrNum;
      domainMax_ = maxChrNum;
    }

    double getRate(std::vector<Parameter*> params, double state) const;
    size_t getNumOfParameters() const{return 3;}
    void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, double maxChrNumber);
    void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, double maxChrNumber);

};


} // end of namespace bpp.

#endif  // CHROMEVOL_CHROMOSOMENUMBERDEPENDENCYFUNCTION_H