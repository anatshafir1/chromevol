#ifndef CHROMEVOL_CONTINUOUSTRAITMODEL_H
#define CHROMEVOL_CONTINUOUSTRAITMODEL_H


#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <random>
#include <Bpp/Numeric/Random/RandomTools.h>
#include "ChromosomeNumberDependencyFunction.h"
#include "ChromosomeSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The ContinuousTraitModel class: this class contains the various continuous trait models.
 *
 */ 

class ContinuousTraitModel:
    public virtual AbstractParameterAliasable
{
public:
    enum ParameterType {MU, SIGMA, ALPHA, THETA};
    
public:
    ContinuousTraitModel(string modelName) : AbstractParameterAliasable(modelName + "."){}
    virtual ~ContinuousTraitModel() {}

  /**
   * @brief update model parameters according to the new current state.
   *
   * @param state The current state of the character.
   */
    virtual double getParamValue(string &paramName, double state) const = 0;
    void createCompositeParameter(ChromosomeNumberDependencyFunction::FunctionType &func, 
        std::string paramName, 
        vector<double> &vectorOfValues, 
        string &prefix, 
        double minState, 
        double maxState, 
        compositeParameter** compositeParam);
      /**
   * @brief update model parameters according to the new current state.
   *
   * @param time The time duration for the simulation.
   * @param startState The current state.
   * @param dependencyState The state that the parameters depend on.
   * @param generator The random number generator.
   * @return The simulated state after time t.
   */
    virtual double simulate(double time, double startState, double dependencyState) = 0;//, std::default_random_engine &generator) = 0;
    virtual double sampleRootState() const = 0;
    static string getParameterName(ParameterType paramType){
        string paramName;
        switch(paramType){
            case MU:
                paramName = "mu";
                break;
            case SIGMA:
                paramName = "sigma";
                break;
            case ALPHA:
                paramName = "alpha";
                break;
            case THETA:
                paramName = "theta";
                break;
            default:
                throw Exception("ContinuousTraitModel::getParameterName: Unknown parameter type.");
        }
        return paramName;
    }
    static ParameterType getParameterType(const string &paramName){
        if(paramName == "mu")
            return MU;
        else if(paramName == "sigma")
            return SIGMA;
        else if(paramName == "alpha")
            return ALPHA;
        else if(paramName == "theta")
            return THETA;
        else
            throw Exception("ContinuousTraitModel::getParameterType: Unknown parameter name " + paramName + ".");
    }
    static vector<string> getModelSimpleParamNames(string modelName){
        vector<string> paramNames;
        if (modelName == "BrownianMotion"){
            paramNames = {"mu", "sigma"};
        }else{
            throw Exception("ContinuousTraitModel::getModelSimpleParamNames: Unknown model name " + modelName + ".");
            //paramNames = {"mu", "sigma", "alpha", "theta"};
        }
        return paramNames;
    }
    static vector<string> getModelNonCompositeParamNames(string modelName){
        vector<string> paramNames;
        if (modelName == "BrownianMotion"){
            paramNames = {"sigma"};
            
        }else{
            throw Exception("ContinuousTraitModel::getModelSimpleParamNames: Unknown model name " + modelName + ".");
        }
        return paramNames;
    }

};
class BrownianMotionModel:
  public virtual ContinuousTraitModel,
  public virtual AbstractParameterAliasable{
    protected:
        double mu_;
        compositeParameter* sigma_;
        ChromosomeNumberDependencyFunction::FunctionType sigmaFuncType_;
        //double state_;
    public:
        BrownianMotionModel(double mu, 
            vector<double> sigma, 
            double minDomain, 
            double maxDomain,
            ChromosomeNumberDependencyFunction::FunctionType sigmaFuncType):
            AbstractParameterAliasable("BrownianMotion."),
            ContinuousTraitModel("BrownianMotion."),
            mu_(mu),sigma_(0),sigmaFuncType_(sigmaFuncType) {
                string paramPrefix = "BrownianMotion.";
                createCompositeParameter(sigmaFuncType_, "sigma", sigma, paramPrefix, minDomain, maxDomain, &sigma_);
            }

        virtual ~BrownianMotionModel(){}

        BrownianMotionModel(const BrownianMotionModel& model):
            AbstractParameterAliasable(model),
            ContinuousTraitModel(model),
            mu_(model.mu_),
            sigma_(model.sigma_){}

        BrownianMotionModel* clone() const { return new BrownianMotionModel(*this);}

    
        double getParamValue(string &paramName, double state) const override;
        double simulate(double time, double startState, double dependencyState) override;
        double sampleRootState() const override{return mu_;};

  };

}
#endif // CHROMEVOL_CONTINUOUSTRAITMODEL_H