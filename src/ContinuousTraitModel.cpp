#include "ContinuousTraitModel.h"

using namespace std;
using namespace bpp;

/******************************************************************************/
void ContinuousTraitModel::createCompositeParameter(ChromosomeNumberDependencyFunction::FunctionType &func, 
    std::string paramName, 
    vector<double> &vectorOfValues, 
    string &prefix, 
    double minState, 
    double maxState, 
    compositeParameter** compositeParam)
{
    std::vector<Parameter*> params;
    if (func == ChromosomeNumberDependencyFunction::FunctionType::IGNORE){
        *compositeParam = 0;
    }
    for (size_t i = 0; i < vectorOfValues.size(); i++){
        auto paramValue = vectorOfValues[i];
        double lowerBound;
        double upperBound;
        ChromosomeNumberDependencyFunction* functionOp =  compositeParameter::setDependencyFunction(func);
        functionOp->setDomainsIfNeeded(minState, maxState);
        functionOp->getAbsoluteBounds(i, &lowerBound, &upperBound, maxState);
        delete functionOp;
        std::shared_ptr<IntervalConstraint> interval = make_shared<IntervalConstraint>(lowerBound, upperBound, true, true);
        Parameter* param = new Parameter(prefix + paramName + std::to_string(i), paramValue, interval);
        params.push_back(param);
        addParameter_(param);
    }
    *compositeParam = new compositeParameter(func, paramName, params, 0);

}
/******************************************************************************/

double BrownianMotionModel::getParamValue(string &paramName, double state) const{
    if (paramName == "BrownianMotion.sigma"){
        double sigmaValue = sigma_->getRate(state);
        return sigmaValue;
    }else if(paramName == "BrownianMotion.mu"){
        return mu_;
    }else{
        throw Exception("BrownianMotionModel::getParamValue(): Parameter name " + paramName + " not recognized!");
    }
}
/******************************************************************************/
double BrownianMotionModel::simulate(double time, double startState, double dependencyState){

    double sigmaValue = sigma_->getRate(dependencyState);
    if (sigmaValue < 0){
        throw Exception("BrownianMotionModel::simulate(): Negative sigma rate!");
    }
    double standard_deviation = std::sqrt(sigmaValue * time);
    std::normal_distribution<double> distribution(0.0, standard_deviation);
    double x = distribution(RandomTools::DEFAULT_GENERATOR);
    double newState = x + startState;
    return newState;
}

void GenomeSizeModel::updateState(double &currentTraitState, double &prevChromState, double &currentChromState){
    if (std::abs(currentChromState - prevChromState) > 1){
        currentTraitState = currentTraitState * (currentChromState/prevChromState);
    }
    return;
}