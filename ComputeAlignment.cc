/////////////////////////////////////////////////////////////////
// ComputeAlignment.cc
//
// Routines for (1) maximum weight trace alignment
//
/////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdio>
#include "SafeVector.h"
#include "SparseMatrix.h"
#include "MultiSequence.h"
#include<string>
#include<iostream>
#include<cfloat>

using namespace std;
const float LOG_ZERO = -2e20;
const float LOG_ONE = 0.0;

/////////////////////////////////////////////////////////////////
// ChooseBestOfThree()
//
// Store the largest of three values x1, x2, and x3 in *x.  Also
// if xi is the largest value, then store bi in *b.
/////////////////////////////////////////////////////////////////

inline void ChooseBestOfThree(float x1, float x2, float x3, char b1, char b2,
		char b3, float *x, char *b) {
	if (x1 >= x2) {
		if (x1 >= x3) {
			*x = x1;
			*b = b1;
			return;
		}
		*x = x3;
		*b = b3;
		return;
	}
	if (x2 >= x3) {
		*x = x2;
		*b = b2;
		return;
	}
	*x = x3;
	*b = b3;
}

/////////////////////////////////////////////////////////////////
// ComputeAlignment()
//
// Computes an alignment based on the given posterior matrix.
// This is done by finding the maximum summing path (or
// maximum weight trace) through the posterior matrix.  The
// final alignment is returned as a pair consisting of:
//    (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and
//        denote insertions in one of the two sequences and
//        B's denote that both sequences are present (i.e.
//        matches).
//    (2) a float indicating the sum achieved
/////////////////////////////////////////////////////////////////

pair < SafeVector < char>*, float> ComputeAlignment(int seq1Length,
		int seq2Length, const VF & posterior) {

	float *twoRows = new float[(seq2Length + 1) * 2];
	assert(twoRows);
	float *oldRow = twoRows;
	float *newRow = twoRows + seq2Length + 1;

	char *tracebackMatrix = new char[(seq1Length + 1) * (seq2Length + 1)];
	assert(tracebackMatrix);
	char *tracebackPtr = tracebackMatrix;

	VF::const_iterator posteriorPtr = posterior.begin() + seq2Length + 1;

	// initialization
	for (int i = 0; i <= seq2Length; i++) {
		oldRow[i] = 0;
		*(tracebackPtr++) = 'L';
	}

	// fill in matrix
	for (int i = 1; i <= seq1Length; i++) {

		// initialize left column
		newRow[0] = 0;
		posteriorPtr++;
		*(tracebackPtr++) = 'U';

		// fill in rest of row
		for (int j = 1; j <= seq2Length; j++) {
			ChooseBestOfThree(*(posteriorPtr++) + oldRow[j - 1], newRow[j - 1],
					oldRow[j], 'D', 'L', 'U', &newRow[j], tracebackPtr++);
		}

		// swap rows
		float *temp = oldRow;
		oldRow = newRow;
		newRow = temp;
	}

	// store best score
	float total = oldRow[seq2Length];
	delete[]twoRows;

	// compute traceback
	SafeVector < char>*alignment = new SafeVector < char>;
	assert(alignment);
	int r = seq1Length, c = seq2Length;
	while (r != 0 || c != 0) {
		char ch = tracebackMatrix[r * (seq2Length + 1) + c];
		switch (ch) {
		case 'L':
			c--;
			alignment->push_back('Y');
			break;
		case 'U':
			r--;
			alignment->push_back('X');
			break;
		case 'D':
			c--;
			r--;
			alignment->push_back('B');
			break;
		default:
			assert(false);
		}
	}

	delete[]tracebackMatrix;

	reverse(alignment->begin(), alignment->end());

	return make_pair(alignment, total);
}

/////////////////////////////////////////////////////////////////
// ComputeAlignmentWithGapPenalties()
//
// Similar to ComputeAlignment() except with gap penalties.
/////////////////////////////////////////////////////////////////

pair < SafeVector < char>*,
float> ComputeAlignmentWithGapPenalties(MultiSequence * align1,
		MultiSequence * align2, const VF & posterior, int numSeqs1,
		int numSeqs2, float gapOpenPenalty, float gapContinuePenalty) {
	int seq1Length = align1->GetSequence(0)->GetLength();
	int seq2Length = align2->GetSequence(0)->GetLength();
	SafeVector < SafeVector < char>::iterator>
			dataPtrs1(align1->GetNumSequences());
	SafeVector < SafeVector < char>::iterator>
			dataPtrs2(align2->GetNumSequences());

	// grab character data
	for (int i = 0; i < align1->GetNumSequences(); i++)
		dataPtrs1[i] = align1->GetSequence(i)->GetDataPtr();
	for (int i = 0; i < align2->GetNumSequences(); i++)
		dataPtrs2[i] = align2->GetSequence(i)->GetDataPtr();

	// the number of active sequences at any given column is defined to be the
	// number of non-gap characters in that column; the number of gap opens at
	// any given column is defined to be the number of gap characters in that
	// column where the previous character in the respective sequence was not
	// a gap
	SafeVector < int> numActive1(seq1Length + 1), numGapOpens1(seq1Length + 1);
	SafeVector < int> numActive2(seq2Length + 1), numGapOpens2(seq2Length + 1);

	// compute number of active sequences and gap opens for each group
	for (int i = 0; i < align1->GetNumSequences(); i++) {
		SafeVector < char>::iterator dataPtr = align1->GetSequence(i)->GetDataPtr();
		numActive1[0] = numGapOpens1[0] = 0;
		for (int j = 1; j <= seq1Length; j++) {
			if (dataPtr[j] != '-') {
				numActive1[j]++;
				numGapOpens1[j] += (j != 1 && dataPtr[j - 1] != '-');
			}
		}
	}
	for (int i = 0; i < align2->GetNumSequences(); i++) {
		SafeVector < char>::iterator dataPtr = align2->GetSequence(i)->GetDataPtr();
		numActive2[0] = numGapOpens2[0] = 0;
		for (int j = 1; j <= seq2Length; j++) {
			if (dataPtr[j] != '-') {
				numActive2[j]++;
				numGapOpens2[j] += (j != 1 && dataPtr[j - 1] != '-');
			}
		}
	}

	VVF openingPenalty1(numSeqs1 + 1, VF(numSeqs2 + 1));
	VF continuingPenalty1(numSeqs1 + 1);
	VVF openingPenalty2(numSeqs1 + 1, VF(numSeqs2 + 1));
	VF continuingPenalty2(numSeqs2 + 1);

	// precompute penalties
	for (int i = 0; i <= numSeqs1; i++)
		for (int j = 0; j <= numSeqs2; j++)
			openingPenalty1[i][j] = i * (gapOpenPenalty * j
					+ gapContinuePenalty * (numSeqs2 - j));
	for (int i = 0; i <= numSeqs1; i++)
		continuingPenalty1[i] = i * gapContinuePenalty * numSeqs2;
	for (int i = 0; i <= numSeqs2; i++)
		for (int j = 0; j <= numSeqs1; j++)
			openingPenalty2[i][j] = i * (gapOpenPenalty * j
					+ gapContinuePenalty * (numSeqs1 - j));
	for (int i = 0; i <= numSeqs2; i++)
		continuingPenalty2[i] = i * gapContinuePenalty * numSeqs1;

	float *twoRows = new float[6 * (seq2Length + 1)];
	assert(twoRows);
	float *oldRowMatch = twoRows;
	float *newRowMatch = twoRows + (seq2Length + 1);
	float *oldRowInsertX = twoRows + 2 * (seq2Length + 1);
	float *newRowInsertX = twoRows + 3 * (seq2Length + 1);
	float *oldRowInsertY = twoRows + 4 * (seq2Length + 1);
	float *newRowInsertY = twoRows + 5 * (seq2Length + 1);

	char *tracebackMatrix = new char[3 * (seq1Length + 1) * (seq2Length + 1)];
	assert(tracebackMatrix);
	char *tracebackPtr = tracebackMatrix;

	VF::const_iterator posteriorPtr = posterior.begin() + seq2Length + 1;

	// initialization
	for (int i = 0; i <= seq2Length; i++) {
		oldRowMatch[i] = oldRowInsertX[i] = (i == 0) ? 0 : LOG_ZERO;
		oldRowInsertY[i] = (i == 0) ? 0 : oldRowInsertY[i - 1]
				+ continuingPenalty2[numActive2[i]];
		*(tracebackPtr) = *(tracebackPtr + 1) = *(tracebackPtr + 2) = 'Y';
		tracebackPtr += 3;
	}

	// fill in matrix
	for (int i = 1; i <= seq1Length; i++) {

		// initialize left column
		newRowMatch[0] = newRowInsertY[0] = LOG_ZERO;
		newRowInsertX[0] = oldRowInsertX[0] + continuingPenalty1[numActive1[i]];
		posteriorPtr++;
		*(tracebackPtr) = *(tracebackPtr + 1) = *(tracebackPtr + 2) = 'X';
		tracebackPtr += 3;

		// fill in rest of row
		for (int j = 1; j <= seq2Length; j++) {

			// going to MATCH state
			ChooseBestOfThree(oldRowMatch[j - 1], oldRowInsertX[j - 1],
					oldRowInsertY[j - 1], 'M', 'X', 'Y', &newRowMatch[j],
					tracebackPtr++);
			newRowMatch[j] += *(posteriorPtr++);

			// going to INSERT X state
			ChooseBestOfThree(oldRowMatch[j]
					+ openingPenalty1[numActive1[i]][numGapOpens2
					[j]], oldRowInsertX[j] + continuingPenalty1[numActive1[i]],
					oldRowInsertY[j]
							+ openingPenalty1[numActive1[i]][numGapOpens2
							[j]], 'M', 'X', 'Y', &newRowInsertX[j],
					tracebackPtr++);

			// going to INSERT Y state
			ChooseBestOfThree(newRowMatch[j - 1]
					+ openingPenalty2[numActive2[j]][numGapOpens1
					[i]], newRowInsertX[j - 1]
					+ openingPenalty2[numActive2[j]][numGapOpens1
					[i]], newRowInsertY[j - 1]
					+ continuingPenalty2[numActive2[j]], 'M', 'X', 'Y',
					&newRowInsertY[j], tracebackPtr++);
		}

		// swap rows
		float *temp;
		temp = oldRowMatch;
		oldRowMatch = newRowMatch;
		newRowMatch = temp;
		temp = oldRowInsertX;
		oldRowInsertX = newRowInsertX;
		newRowInsertX = temp;
		temp = oldRowInsertY;
		oldRowInsertY = newRowInsertY;
		newRowInsertY = temp;
	}

	// store best score
	float total;
	char matrix;
	ChooseBestOfThree(oldRowMatch[seq2Length], oldRowInsertX[seq2Length],
			oldRowInsertY[seq2Length], 'M', 'X', 'Y', &total, &matrix);

	delete[]twoRows;

	// compute traceback
	SafeVector < char>*alignment = new SafeVector < char>;
	assert(alignment);
	int r = seq1Length, c = seq2Length;
	while (r != 0 || c != 0) {

		int offset = (matrix == 'M') ? 0 : (matrix == 'X') ? 1 : 2;
		char ch = tracebackMatrix[(r * (seq2Length + 1) + c) * 3 + offset];
		switch (matrix) {
		case 'Y':
			c--;
			alignment->push_back('Y');
			break;
		case 'X':
			r--;
			alignment->push_back('X');
			break;
		case 'M':
			c--;
			r--;
			alignment->push_back('B');
			break;
		default:
			assert(false);
		}
		matrix = ch;
	}

	delete[]tracebackMatrix;

	reverse(alignment->begin(), alignment->end());

	return make_pair(alignment, 1.0f);
}

/////////////////////////////////////////////////////////////////
// BuildPosterior()
//
// Builds a posterior probability matrix needed to align a pair
// of alignments.  Mathematically, the returned matrix M is
// defined as follows:
//    M[i,j] =     sum          sum      f(s,t,i,j)
//             s in align1  t in align2
// where
//                  [  P(s[i'] <--> t[j'])
//                  [       if s[i'] is a letter in the ith column of align1 and
//                  [          t[j'] it a letter in the jth column of align2
//    f(s,t,i,j) =  [
//                  [  0    otherwise
//
/////////////////////////////////////////////////////////////////


VF *BuildPosterior(MultiSequence * align1, MultiSequence * align2,
		const SafeVector < SafeVector <
		SparseMatrix *> >&sparseMatrices, float cutoff = 0.0f) {
	const int seq1Length = align1->GetSequence(0)->GetLength();
	const int seq2Length = align2->GetSequence(0)->GetLength();

	VF *posteriorPtr = new VF((seq1Length + 1) * (seq2Length + 1), 0);
	assert(posteriorPtr);
	VF & posterior = *posteriorPtr;

	// VF::iterator postPtr = posterior.begin();

	// for each s in align1
	for (int i = 0; i < align1->GetNumSequences(); i++) {
		int first = align1->GetSequence(i)->GetLabel();
		SafeVector < int>*mapping1 = align1->GetSequence(i)->GetMapping();

		// for each t in align2
		for (int j = 0; j < align2->GetNumSequences(); j++) {
			int second = align2->GetSequence(j)->GetLabel();
			SafeVector < int>*mapping2 = align2->GetSequence(j)->GetMapping();

			if (first < second) {

				// get the associated sparse matrix
				SparseMatrix *matrix = sparseMatrices[first][second];

				for (int ii = 1; ii <= matrix->GetSeq1Length(); ii++) {
					SafeVector < PIF>::iterator row = matrix->GetRowPtr(ii);
					int base = (*mapping1)[ii] * (seq2Length + 1);
					int rowSize = matrix->GetRowSize(ii);

					// add in all relevant values
					for (int jj = 0; jj < rowSize; jj++)
						posterior[base + (*mapping2)[row[jj].first]] += row[jj].second;

					// subtract cutoff
					for (int jj = 0; jj < matrix->GetSeq2Length(); jj++)
						posterior[base + (*mapping2)[jj]] -= cutoff;
				}

			} else {

				// get the associated sparse matrix
				SparseMatrix *matrix = sparseMatrices[second][first];

				for (int jj = 1; jj <= matrix->GetSeq1Length(); jj++) {
					SafeVector < PIF>::iterator row = matrix->GetRowPtr(jj);
					int base = (*mapping2)[jj];
					int rowSize = matrix->GetRowSize(jj);

					// add in all relevant values
					for (int ii = 0; ii < rowSize; ii++)
						posterior[base +
						(*mapping1)[row[ii].first] *
						(seq2Length + 1)] += row[ii].second;

					// subtract cutoff
					for (int ii = 0; ii < matrix->GetSeq2Length(); ii++)
						posterior[base +
						(*mapping1)[ii] * (seq2Length + 1)] -= cutoff;
				}

			}

			delete mapping2;
		}

		delete mapping1;
	}

	return posteriorPtr;
}
