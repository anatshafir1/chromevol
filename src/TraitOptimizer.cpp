#include "TraitOptimizer.h"

using namespace bpp;
using namespace std;



    void TraitOptimizer::printLikelihoodParameters(JointPhyloLikelihood* lik, unsigned int optimized, vector<string> paramsNames) const{

        if (optimized == 0){
            std::cout << "Initial likelihood is : "<< lik->getValue() << std::endl;
            std::cout << "Initial trait likelihood is: " << lik->getAbstractPhyloLikelihood(lik->getNumbersOfPhyloLikelihoods()[0])->getValue() << std::endl;
            std::cout << "Initial chromosome likelihood is: " << lik->getPhylo2()->getValue() << std::endl;
            std::cout << "Parameters are:" << std::endl;
        }else{
            std::cout << "Optimized likelihood is : " << lik->getValue() << std::endl;
            std::cout << "Optimized trait likelihood is: " << lik->getAbstractPhyloLikelihood(lik->getNumbersOfPhyloLikelihoods()[0])->getValue() << std::endl;
            std::cout << "Optimized chromosome likelihood is: " << lik->getPhylo2()->getValue() << std::endl;
            std::cout << "Optimized parameters are:"<< std::endl;
        }   
        for (int i = 0; i < (int)(paramsNames.size()); i++){
            if (paramsNames[i].find("Chromosome.baseNum_") != std::string::npos){
                std::cout << paramsNames[i] << " = " << (int)(lik->getParameter(paramsNames[i]).getValue()) << std::endl;
            }else{
                std::cout << paramsNames[i] << " = " << lik->getParameter(paramsNames[i]).getValue() << std::endl;

            }    
        }
        std::cout <<  "***" << std::endl;

    }
    /***************************************************************************************************************************************/
    void TraitOptimizer::printLikelihoodParameters(SingleProcessPhyloLikelihood* lik, unsigned int optimized, vector<string> paramsNames) const{

        if (optimized == 0){
            std::cout << "Initial trait likelihood is: " << lik->getValue() << std::endl;
            std::cout << "Parameters are:" << std::endl;
        }else{
            std::cout << "Optimized trait likelihood is: " << lik->getValue() << std::endl;
            std::cout << "Optimized parameters are:"<< std::endl;
        }   
        for (int i = 0; i < (int)(paramsNames.size()); i++){

            std::cout << paramsNames[i] << " = " << lik->getParameter(paramsNames[i]).getValue() << std::endl;

             
        }
        std::cout <<  "***" << std::endl;

    }

/***************************************************************************************************************************************/
vector<string> TraitOptimizer::getTraitFixedParamFullNames(vector<string> &parameterNames, vector<string> &fixedTraitParams, bool fixedRootFreqs){
    vector<string> fixedParamFullNames;
    for (auto &fixedParam : fixedTraitParams){            
        for (auto &paramName : parameterNames){
            size_t found = paramName.find(fixedParam);
            if (found != std::string::npos){
                fixedParamFullNames.push_back(paramName);
            }
        }
    }
    if (fixedRootFreqs){
        std::regex rgx("theta\\d+");
        std::smatch match;
        for (auto &paramName : parameterNames){
            if (std::regex_search(paramName, match, rgx)){
                fixedParamFullNames.push_back(paramName);
            }

        }

    }
    return fixedParamFullNames;

}
