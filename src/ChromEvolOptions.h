#ifndef CHROMEVOL_CHROMEVOLOPTIONS_H
#define CHROMEVOL_CHROMEVOLOPTIONS_H

// From bpp-core:

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>


// From bpp-seq:
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include "ChromosomeSubstitutionModel.h"
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>



//standard libraries
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_map>
// the params file should contain vectors of parameters. For example,
//_gain = 2,1
//_gainFunc = LINEAR
// and so on...
// this implementation is meant to deal with more functions.
namespace bpp{

class ChromEvolOptions
{

public:
    static void initAllParameters(BppApplication& ChromEvol);
    static void getInitialValuesForComplexParams(std::map<uint, std::pair<int, std::map<int, vector<double>>>> &mapOfParams);
    static void getInitialValuesForComplexParamsForJointTraitModel(std::map<uint, std::pair<int, std::map<int, vector<double>>>> &mapOfParams, uint numOfModels);
    //static void initVectorOfChrNumParameters(vector<double>& paramVector);
    virtual ~ChromEvolOptions(){};
    
public:
    static string treeFilePath_;
    static string characterFilePath_;
    static int maxChrNum_;
    static int minChrNum_;
    static int numOfModels_;
    static double branchMul_;
    static std::vector <unsigned int> OptPointsNum_;
    static std::vector <unsigned int> OptIterNum_;
    static std::vector <unsigned int> OptPointsNumNextRounds_;
    static std::vector <unsigned int> OptIterNumNextRounds_;
    static std::map<uint, std::vector<double>> gain_;
    static std::map<uint, std::vector<double>> loss_;
    static std::map<uint, std::vector<double>> dupl_;
    static std::map<uint, std::vector<double>> demiDupl_;
    static std::map<uint, int> baseNum_;
    static std::map<uint, std::vector<double>> baseNumR_;
    static double tolerance_;
    static unsigned int maxIterations_;
    static bool maxParsimonyBound_;
    static bool standardOptimization_;
    static int BrentBracketing_;
    static string optimizationMethod_;
    static unsigned int maxAlpha_;
    static unsigned int minAlpha_;
    static int seed_;
    static std::vector <double> probsForMixedOptimization_;
    static string rootFreqs_;
    static string fixedFrequenciesFilePath_;
    static vector<int> rateChangeType_;
    //static bool optimizeBaseNumber_;
    static string baseNumOptimizationMethod_;
    static std::map<uint, std::vector<int>> fixedParams_; //1 if parameter should be fixed. The order corresponds to the one in the model definition.
    static int NumOfSimulations_;
    static int jumpTypeMethod_;
    static bool simulateData_;
    static int numOfDataToSimulate_;
    static string resultsPathDir_;
    static std::map<uint,uint> maxBaseNumTransition_; // needed for the simulator, since there is no data to infer it!
    static double treeLength_;
    static int maxNumOfTrials_; // to test the severity of the underflow problems
    // for heterogeneous model
    static int minCladeSize_;
    static int maxNumOfModels_;
    static std::map<uint, std::vector<uint>> mapOfNodeIdsPerModel_;
    static std::map<int, vector<std::pair<uint, int>>> sharedParameters_;
    static bool heterogeneousModel_;
    static double deltaAICcThreshold_;
    static std::map<uint, std::vector<uint>> mapModelNodesIds_;
    static string nodeIdsFilePath_;
    static std::vector<uint> initialModelNodes_;
    static std::vector<string> globalParams_;
    static bool parallelization_;
    static int maxChrInferred_;
    static bool backwardPhase_;
    static bool forwardPhase_;
    static bool runStochasticMapping_;
    static size_t numOfStochasticMappingTrials_;
    static size_t numOfFixingMappingIterations_;
    static bool useMaxBaseTransitonNumForOpt_;
    static string modelSelectionCriterion_;
    static size_t numOfSimulatedData_;
    static double fracAllowedFailedSimulations_;
    static bool correctBaseNumber_;
    static size_t numOfRequiredSimulatedData_;
    static string traitFilePath_;
    static string traitStateModel_;
    static vector<string> fixedTraitParams_;
    static std::map<string, double> traitParams_;
    static bool runOnlyJointModel_;
    static bool runOnlyIndependentModelWithTrait_;
    static int minBaseNumberBound_;
    static bool useMLReconstruction_;
    static int numberOfTraitStates_;
    static bool simulateTrait_;
    static bool heteroBootstrappingMode_;
    static std::unordered_map<std::string, string> sharedTraitParams_;
    static int numOfTraitConstraints_;
    static bool computeExpectations_;
    static bool fixedTraitRootFreqs_;

    // public functions
    static std::vector<int> translateStringParamsToInt(std::vector<string> &strParams);
    static std::shared_ptr<PhyloNode> getMRCA(PhyloTree* tree, std::vector<shared_ptr<PhyloNode>> nodes);
    //static std::string getParamName(int type);

private:
    static void initDefaultParameters();
    static void initParametersFromFile(BppApplication& ChromEvol);
    //static void setFixedParams(std::vector<unsigned int> fixedParams);
    
    static void setFunctions(std::string gainFunc, std::string lossFunc, std::string duplFunc, std::string demiDuplFunc, std::string baseNumRFunc);
    static int getFunctionFromString(string funcStr);
    static void setModelParameters(BppApplication& ChromEvol);
    //static void setSharedParametersInterModels();
    static void setSharedParametersPerModel(std::map<uint, std::map<int, std::pair<int, vector<double>>>> mapModelTypeValues, std::map<uint, std::pair<int, int>> mapModelBaseNumTypeAndVal, std::map<int, size_t> paramNumFreqs);
    static void setFixedParameters(BppApplication& ChromEvol);
    static void setFixedParametersTrait(BppApplication& ChromEvol);
    static void updateModelParameter(uint model, int type, vector<double> paramValues);
    //static string getParameterNameWithoutNamespace(const std::string& name);

};

}


#endif  // CHROMEVOL_CHROMEVOLOPTIONS_H