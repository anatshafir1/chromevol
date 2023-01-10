#ifndef CHROMEVOL_CHRFASTA_H
#define CHROMEVOL_CHRFASTA_H

#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Seq/Io/AbstractISequence.h>
#include <Bpp/Seq/Io/AbstractIAlignment.h>
#include <Bpp/Seq/Io/AbstractOSequence.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Io/ISequenceStream.h>
#include <Bpp/Seq/Io/OSequenceStream.h>
#include <Bpp/Seq/Io/SequenceFileIndex.h>


namespace bpp{


/**
 * @brief The fasta sequence file format.
 *
 * Read and write from/to Fasta files.
 */
class ChrFasta
{

  public:
    static VectorSequenceContainer* readSequencesFromFile(const std::string& path , ChromosomeAlphabet* alpha);
    static bool nextSequence(std::istream& input, Sequence& seq, ChromosomeAlphabet* alpha);
    virtual ~ChrFasta() {}


};

}

#endif // CHROMEVOL_CHRFASTA_H