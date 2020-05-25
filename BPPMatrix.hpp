//////////////////////////////////////////////////////////////
// BPPMatrix.hpp
//
// save the Base Pairing Probability Matrix for each sequences
//////////////////////////////////////////////////////////////

#ifndef __BBPMatrix_HPP__
#define __BBPMatrix_HPP__

#include <iostream>
#include <string>
#include "SparseMatrix.h"
#include "McCaskill.hpp"
#include "nrutil.h"

using namespace std;

class BPPMatrix {
private:

    int seqLength;       // sequence length;
    float cutOff;        // the threshold of probability

    BPPMatrix() {}
public:
    SparseMatrix bppMat; // base pairing probability matrix 1-origin(1..seqLength)
    BPPMatrix(const string &seq, int seqLength, float inCutOff = 0.0001) : seqLength (seqLength), cutOff(inCutOff) {
	setBPPMatrix(seq);
    }
    BPPMatrix(int seqLength, float inCutOff, const Trimat<float> & tmpBppMat) : seqLength(seqLength), cutOff(inCutOff) {
      bppMat.SetSparseMatrix(seqLength, seqLength, tmpBppMat, cutOff);
    }


    void setBPPMatrix(const string &seq) {
	const char *tmpSeq = seq.c_str();
	McCaskill mc(seqLength, &tmpSeq[1]);
	mc.calcPartitionFunction();
	Trimat<float> tmpBppMat(seqLength + 1);

	for (int i = 0; i < seqLength; i++) {
	  for(int j = i; j < seqLength; j++) {
	    tmpBppMat.ref(i+1, j+1) = mc.getProb(i,j);
	  }
	}
	bppMat.SetSparseMatrix(seqLength, seqLength, tmpBppMat, cutOff);
    }
    float GetProb(int i, int j) {
	return bppMat.GetValue(i,j);
    }

    float GetLength() const {
	return seqLength;
    }

    void updateBPPMatrix(const Trimat<float> &inbppMat) {
	bppMat.SetSparseMatrix(seqLength, seqLength, inbppMat, cutOff);
    }
};

#endif // __BPPMatrix_HPP__
