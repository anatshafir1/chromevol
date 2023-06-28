#include "ChrFasta.h"
#include <fstream>
#include <Bpp/Seq/StringSequenceTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Io/FileTools.h>


using namespace bpp;
using namespace std;

// mainly code duplication from the fasta class. Not elegant- should think how to avoid it here!
VectorSequenceContainer* ChrFasta::readSequencesFromFile(const std::string& path , ChromosomeAlphabet* alpha){
    VectorSequenceContainer* vsc = new VectorSequenceContainer(alpha);
    std::ifstream input(path.c_str(), std::ios::in);
    if (!input)
        throw IOException("Fasta::appendFromStream: can't read from istream input");
    char c = '\n';
    char last_c;
    bool hasSeq = true;
    string line = "";
    while (!input.eof() && hasSeq)
    {
        last_c = c;
        input.get(c);
        if (c == '\n'){
            continue;
        }
        // Sequence detection
        if (c == '>' && last_c == '\n'){
            input.putback(c);
            c = last_c;
            BasicSequence tmpseq("", "", vsc->getAlphabet());
            hasSeq = nextSequence(input, tmpseq, alpha);
            vsc->addSequence(tmpseq, true);
        }
    }
    input.close();
    return vsc;
}
/*******************************************************************************/
bool ChrFasta::nextSequence(std::istream& input, Sequence& seq, ChromosomeAlphabet* alpha){
  if (!input)
    throw IOException("Fasta::nextSequence: can't read from istream input");
  string seqname = "";
  string content = "";
  short seqcpt = 0;
  string linebuffer = "";
  char c;
  while (!input.eof())
  {
    c = static_cast<char>(input.peek());
    if (input.eof())
      c = '\n';

    // Sequence begining detection
    if (c == '>')
    {
      // Stop if find a new sequence
      if (seqcpt++)
        break;
    }
    getline(input, linebuffer);
    if (c == '>')
    {
      // Get the sequence name line
      seqname = string(linebuffer.begin() + 1, linebuffer.end());
    }
    if (c != '>' && !TextTools::isWhiteSpaceCharacter(c)) {
      // Sequence content
      content += TextTools::toUpper(TextTools::removeWhiteSpaces(linebuffer));
    }
  }

  bool res = (!input.eof());
  seq.setName(seqname);
  //need to check if 
  alpha->setCompositeState(content);
  auto clean_content = alpha->extractChromosomeNumber(content);
  seq.setContent(clean_content);
  return res;

}