//standard libraries
#include <string>
#include <vector>
#include <iostream>
#include <time.h>

//from bpp-core
#include <Bpp/Version.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/BppApplication.h>




//from bpp-seq
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>




#include "ChromEvolOptions.h"
#include "ChromosomeNumberMng.h"




using namespace bpp;
using namespace std;




int main(int args, char **argv) {

    if (args == 1){
        std::cout << "No arguments provided"<<endl;
        return 0;
    }
    try{
        time_t t1;
        time(&t1);
        time_t t2;
        BppApplication ChromEvol(args, argv, "ChromEvol");
        ChromEvolOptions::initAllParameters(ChromEvol);
        ChromosomeNumberMng* mng = new ChromosomeNumberMng();
        if (!ChromEvolOptions::simulateData_){
            mng->getCharacterData(ChromEvolOptions::characterFilePath_);
        }
        mng->getTree(ChromEvolOptions::treeFilePath_, ChromEvolOptions::treeLength_);       
        //mng->runChromEvol();
        mng->runChromEvol();
        time(&t2);
        std::cout << "****** Max allowed chromosome number: "<< ChromEvolOptions::maxChrNum_ <<endl;
        std::cout <<"Total running time is: "<< static_cast<int>(t2-t1) <<endl;
        delete mng;

    }
    catch (exception& e)
    {
        cout << e.what() << endl;
        return 1;
    }

    return 0;
}
