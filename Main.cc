/////////////////////////////////////////////////////////////////
// Main.cc
//
// Main routines for PicXAA-R 1.0 program (Jan 2011).
//
/////////////////////////////////////////////////////////////////

#include "SafeVector.h"
#include "MultiSequence.h"
#include "SparseMatrix.h"
#include "ProbabilisticModel.h"
#include "Defaults.h"
#include "ScoreType.h"
#include "AlignGraph.h"
#include <string>
#include <iomanip>
#include <iostream>
#include <cerrno>
#include <sys/stat.h>
#include <sys/types.h>

#include "SparseMatrix.h"
#include "BPPMatrix.hpp"
#include "nrutil.h"
#include "AlifoldMEA.h"
#include <sstream>
#include <list>
#include <set>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <fstream>


//struct for column reliability

typedef struct {
	int columnNo;
	float probProduct;

} columnReliability;

columnReliability *column;

bool enableStockholmOutput = false;
bool enableMXSCARNAOutput = false;
bool enableVerbose = false;
bool enableAnnotation = false;
bool enableClustalWOutput = false;
bool enableAlignOrder = false;
bool NoSS = false;
float BaseProbThreshold= BASEPROBTHRESHOLD;
float BasePairConst= BASEPAIRCONST;
int BandWidth= BANDWIDTH;
SafeVector<string> sequenceNames;

string *ssCons= NULL;

const int MIN_ITERATIVE_REFINEMENT_REPS = 0;
const int MAX_ITERATIVE_REFINEMENT_REPS = 1000;

string parametersInputFilename = "";

int numInterConsistencyReps = 1;
int numIntraConsistencyReps = 1;
int numRefinementsReAligns = 100;

float cutoff = 0;
float gapOpenPenalty = 0;
float gapContinuePenalty = 0;

string infilename;

const int MIN_INTERCONSISTENCY_REPS = 0;
const int MAX_INTERCONSISTENCY_REPS = 5;
const int MIN_INTRACONSISTENCY_REPS = 0;
const int MAX_INTRACONSISTENCY_REPS = 5;
const int MIN_REFINEMENTSREALIGNS = 0;
const int MAX_REFINEMENTSREALIGNS = 1000;

float ALPHA=0.4;
float BETA=0.1;
float Tb=0.5;

///////////////////////////////
// global scoring matrix variables
//////////////////////////////
VF initDistrib(NumMatrixTypes);
VF gapOpen(2 * NumInsertStates);
VF gapExtend(2 * NumInsertStates);
VVF emitPairs(256, VF(256, 1e-10));
VF emitSingle(256, 1e-5);
string alphabet = alphabetDefault;

/////////////////////////////////////////////////////////////////
// Function prototypes
/////////////////////////////////////////////////////////////////

void PrintHeading();

void DoAlign(MultiSequence * sequences, SafeVector<string> sequenceNames);

pair<VVI, VF> ArrangePosteriorProbs(MultiSequence *sequences,
		const ProbabilisticModel &model,
		SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, VVF &distances,
		SafeVector<BPPMatrix*> &BPPMatrices, VVVF &BPPtob, int maxlength);

SafeVector<string> ParseParams(int argc, char **argv);

void ReadPHMMParameters();

void PrintPHMMParameters(const char *message, const VF &initDistrib,
		const VF &gapOpen, const VF &gapExtend, const VVF &emitPairs,
		const VF &emitSingle, const char *filename);

MultiSequence *AlignAlignments(MultiSequence * align1, MultiSequence * align2,
		const SafeVector<SafeVector<SparseMatrix *> >&sparseMatrices,
		const ProbabilisticModel &model);

SafeVector<SafeVector<SparseMatrix *> > DoRelaxation(MultiSequence * sequences,
		SafeVector<SafeVector<SparseMatrix *> >&sparseMatrices, VVF distances,
		bool nflag);

void Relax(SparseMatrix * matXZ, SparseMatrix * matZY, VF & posterior,
		float id_xz, float id_zy, bool nflag);

void Relax1(SparseMatrix * matXZ, SparseMatrix * matZY, VF & posterior,
		float id_xz, float id_zy, bool nflag);

SafeVector<BPPMatrix*> DoBasePairProbabilityRelaxation(MultiSequence *sequences,
		SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, SafeVector<
		BPPMatrix*> &BPPMatrices, VVF distances);

SafeVector<SafeVector<SparseMatrix *> > DoFourWayRelaxation(
		MultiSequence * sequences,
		SafeVector<SafeVector<SparseMatrix *> >&sparseMatrices, SafeVector<
		BPPMatrix*> &BPPMatrices, VVF distances, bool nflag, int RR);
void RelaxBP(SparseMatrix * matXZ, BPPMatrix* bpmatX, BPPMatrix* bpmatZ,
		VF & posterior, int Seqx, int Seqz, float id_zy, bool nflag, float mm);

void
		DoRefinement(MultiSequence* &alignment, SafeVector<SafeVector<
		SparseMatrix *> > &sparseMatrices, const ProbabilisticModel &model,
				VVF distances);
void FindSimilar(VVF distances, vector<set<int> > &SimSeqs);

float FindMaxPP(SafeVector<SafeVector<SparseMatrix *> > sparseMatrices);

bool GetInteger(char *data, int *val);

bool GetFloat(char *data, float *val);

//java gui related change
void WriteAnnotation(MultiSequence * alignment, const SafeVector<SafeVector<
SparseMatrix *> >&sparseMatrices);

//java gui related change
float ComputeScore(const SafeVector<pair<int, int> >&active, const SafeVector<
SafeVector<SparseMatrix *> >&sparseMatrices);

/////////////////////////////////////////////////////////////////
// main()
//
// Calls all initialization routines and runs the PicXAA-R aligner.
/////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

	// print PicXAA-R heading
	PrintHeading();

	// parse program parameters
	SafeVector<string> sequenceNames = ParseParams(argc, argv);
	infilename = sequenceNames[0];

	//Read Parameters for posterior alignment probabilities (PHMM)
	ReadPHMMParameters();
	PrintPHMMParameters("Using parameter set:", initDistrib, gapOpen,
			gapExtend, emitPairs, emitSingle, NULL);

	// now, we'll process all the files given as input.  If we are given
	// several filenames as input, then we'll load all of those sequences
	// simultaneously.


	// load all files together
	MultiSequence *sequences = new MultiSequence();
	assert(sequences);
	for (int i = 0; i < (int) sequenceNames.size(); i++) {
		cerr << "Loading sequence file: " << sequenceNames[i] << endl;
		sequences->LoadMFA(sequenceNames[i], true);
	}

	// now, we can perform the alignments using PicXAA-R algorithm
	DoAlign(sequences, sequenceNames);

	delete sequences;
}

/////////////////////////////////////////////////////////////////
// PrintHeading()
//
// Prints heading for PicXAA-R program.
/////////////////////////////////////////////////////////////////

void PrintHeading() {
	cerr << endl << "PicXAA-R Version 1.0 (Jan 2011) "
			<< "aligns multiple RNA sequences using structural information"
			<< endl
			<< "Written by Sayed Mohammad Ebrahim Sahraeian and Byung-Jun Yoon."
			<< endl << endl << endl;
}

/////////////////////////////////////////////////////////////////
// DoAlign()
// First computes all pairwise posterior probability matrices.
// Then apply the Consistency transformations
// Then construct the graph based on those probabilities
// Then, perform refinement step.
/////////////////////////////////////////////////////////////////

void DoAlign(MultiSequence * sequences, SafeVector<string> sequenceNames) {

	int numSeqs = sequences->GetNumSequences();

	SafeVector<SafeVector<SparseMatrix *> > sparseMatrices(numSeqs, SafeVector<
			SparseMatrix *>(numSeqs, NULL));

	VVF distances(numSeqs, VF(numSeqs, 0));

	// build new model for aligning, will be used for  PicXAA-R
	ProbabilisticModel model(initDistrib, gapOpen, gapExtend, emitPairs,
			emitSingle);

	// base-pairing probability matrices
	SafeVector<BPPMatrix*> BPPMatrices;

	// cumulative left and right base-pairing probability array
	VVVF BPProb;

	int maxlength = 0;
	for (int i = 0; i < sequences->GetNumSequences(); i++)
		if (sequences->GetSequence(i)->GetLength() > maxlength)
			maxlength = sequences->GetSequence(i)->GetLength();

	//compute the posterior pairwise alignment and base pairing probabilities
	pair<VVI, VF> alignsp = ArrangePosteriorProbs(sequences, model,
			sparseMatrices, distances, BPPMatrices, BPProb, maxlength);

	//Construct the alignment graph
	// Step1: Structural skeleton construction (usign base pairing probabilities)
	AlignGraph *MyGraph = new AlignGraph(sequences, sparseMatrices,
			BPPMatrices, BPProb, maxlength,Tb);
	// Step 2: Inserting highly probable local alignments (using base alignment probabilities)
	MyGraph->AlignGraph_PicXAA(alignsp.first, alignsp.second, sequences,
			sparseMatrices, maxlength);

	MyGraph->Graph2Align();
	MultiSequence *alignment = MyGraph->GetAlignment();

	//Perform the refinement step
	DoRefinement(alignment, sparseMatrices, model, distances);

	cerr<<BasePairConst<<endl;
	AlifoldMEA alifold(alignment, BPPMatrices, BasePairConst);
	alifold.Run();
	ssCons = alifold.getSScons();

	if (enableClustalWOutput) {
		if (NoSS)
			alignment->WriteALN(cout);
		else
			alignment->WriteMXSCARNA(cout, ssCons);
	} else if (enableStockholmOutput) {
		alignment->WriteSTOCKHOLM(cout, ssCons);
	} else {
		if (NoSS)
			alignment->WriteMFA2(cout);
		else
			alignment->WriteMFA(cout, ssCons);
	}

	if (enableAnnotation) {
		WriteAnnotation(alignment, sparseMatrices);
	}

	delete alignment;

}

/////////////////////////////////////////////////////////////////
// ArrangePosteriorProbs()
//
// Compute the posterior pairwise alignment and base pairing 
// probabilities
/////////////////////////////////////////////////////////////////

pair<VVI, VF> ArrangePosteriorProbs(MultiSequence *sequences,
		const ProbabilisticModel &model,
		SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, VVF &distances,
		SafeVector<BPPMatrix*> &BPPMatrices, VVVF &BPProb, int maxlength) {

	assert(sequences);

	const int numSeqs = sequences->GetNumSequences();

	int a, b;

	// do all pairwise alignments for posterior probability matrices
	for (a = 0; a < numSeqs - 1; a++) {
		for (b = a + 1; b < numSeqs; b++) {
			Sequence *seq1 = sequences->GetSequence(a);
			Sequence *seq2 = sequences->GetSequence(b);

			// if we are training, then we'll simply want to compute the
			// expected counts for each region within the matrix separately;
			// otherwise, we'll need to put all of the regions together and
			// assemble a posterior probability match matrix


			VF * posterior(NULL);

			VF *forward = model.ComputeForwardMatrix(seq1, seq2);
			assert (forward);
			VF *backward = model.ComputeBackwardMatrix(seq1, seq2);
			assert (backward);
			posterior = model.ComputePosteriorMatrix(seq1, seq2, *forward,
					*backward);
			assert (posterior);
			delete forward;
			delete backward;
			sparseMatrices[a][b] = new SparseMatrix(seq1->GetLength(),
					seq2->GetLength(), *posterior, 0.02);
			sparseMatrices[b][a] = NULL;

			pair<SafeVector<char> *, float> alignment;

			alignment = model.ComputeAlignment(seq1->GetLength(),
					seq2->GetLength(), *posterior);

			int alignlength = 0;
			for (SafeVector<char>::iterator iter = alignment.first->begin(); iter
					!= alignment.first->end(); ++iter)
				if (*iter == 'B')
					alignlength++;

			float distance = alignment.second / (float) alignlength;//  min(seq1->GetLength(),seq2->GetLength());
			distances[a][b] = distances[b][a] = distance;

			delete posterior;

		}
	}

	// compute the base pairing matrices for each sequences
	for (int i = 0; i < numSeqs; i++) {
		Sequence *tmpSeq = sequences->GetSequence(i);
		string seq = tmpSeq->GetString();
		int n_seq = tmpSeq->GetLength();
		BPPMatrix *bppmat = new BPPMatrix(seq, n_seq, BASEPROBTHRESHOLD);
		BPPMatrices.push_back(bppmat);
	}

	// compute the cumulative left and right base-pairing probability array for each sequences
	VF PZ(2, 0);
	for (int i = 0; i < numSeqs; i++) {
		VVF s;
		for (int j = 0; j < sequences->GetSequence(i)->GetLength(); j++)
			s.push_back(PZ);
		BPProb.push_back(s);
	}
	for (int i = 0; i < numSeqs; i++) {
		Sequence *tmpSeq = sequences->GetSequence(i);
		BPPMatrix *tmpBppMatrix = BPPMatrices[i];
		for (int j = 0; j < tmpSeq->GetLength(); j++)
			for (int k = j + 1; k < tmpSeq->GetLength(); k++) {
				BPProb[i][j][1] += tmpBppMatrix->GetProb(j + 1, k + 1);
				BPProb[i][k][0] += tmpBppMatrix->GetProb(j + 1, k + 1);
			}
	}

	// perform the inter-sequence consistency transformation  
	for (int r = 0; r < numInterConsistencyReps; r++) {
		SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices =
				DoRelaxation(sequences, sparseMatrices, distances, true);

		// now replace the old posterior matrices
		for (int i = 0; i < numSeqs; i++) {
			for (int j = 0; j < numSeqs; j++) {
				delete sparseMatrices[i][j];
				sparseMatrices[i][j] = newSparseMatrices[i][j];
			}
		}
	}

	// perform the intra-sequence consistency transformation  
	for (int r = 0; r < numIntraConsistencyReps; r++) {
		SafeVector<BPPMatrix *> newBPPMatrices =
				DoBasePairProbabilityRelaxation(sequences, sparseMatrices,
						BPPMatrices, distances);
		for (int i = 0; i < numSeqs; i++) {
			delete BPPMatrices[i];
			BPPMatrices[i] = newBPPMatrices[i];
		}
	}

	int numFWConsistencyReps = max((int) floor(5 - (float) numSeqs / 5), 1);
	// perform the Four-Way consistency transformation
	for (int r = 0; r < numFWConsistencyReps; r++) {
		SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices =
				DoFourWayRelaxation(sequences, sparseMatrices, BPPMatrices,
						distances, false, numFWConsistencyReps);

		// now replace the old posterior matrices
		for (int i = 0; i < numSeqs; i++) {
			for (int j = 0; j < numSeqs; j++) {
				delete sparseMatrices[i][j];
				sparseMatrices[i][j] = newSparseMatrices[i][j];
			}
		}

	}

	float maxpp = FindMaxPP(sparseMatrices);

	//Put the pairwise posterior probabilities in "alignp" and the pair residue positions in "aligns"
	VVI aligns;
	VF alignp;
	pair<VVI, VF> alignsp;
	int cnt = 0;
	for (int i = 0; i < numSeqs - 1; i++) {
		for (int j = i + 1; j < numSeqs; j++) {
			SparseMatrix *currSpMat = sparseMatrices[i][j];
			int numRows = currSpMat->GetSeq1Length();

			for (int k = 1; k <= numRows; k++) {
				for (int h = 0; h < currSpMat->GetRowSize(k); h++) {

					int ii = k - 1;
					int jj = (currSpMat->GetRowPtr(k)[h].first) - 1;
					if ((currSpMat->GetRowPtr(k)[h].second) >= (maxpp / 40)) {
						cnt++;
						VI newentry;
						newentry.push_back(i);
						newentry.push_back(k - 1);
						newentry.push_back(j);
						newentry.push_back((currSpMat->GetRowPtr(k)[h].first) - 1);
						aligns.push_back(newentry);
						alignp.push_back((currSpMat->GetRowPtr(k)[h].second));
					}
				}
			}
		}
	}

	alignsp.first = aligns;
	alignsp.second = alignp;
	return alignsp;
}

/////////////////////////////////////////////////////////////////
// ParseParams()
//
// Parse all command-line options.
/////////////////////////////////////////////////////////////////

SafeVector<string> ParseParams(int argc, char **argv) {

	if (argc < 2) {

		cerr << "PicXAA-R 1.0 comes with ABSOLUTELY NO WARRANTY." << endl
				<< "This is free software, and you are welcome to redistribute it under "
				<< endl
				<< "certain conditions. See the README file for details."
				<< endl << endl << "Usage:" << endl
				<< "       ./picxaa-r [options]  MFAFILE" << endl << "Example:"
				<< endl
				<< "       ./picxaa-r test/sample1.fasta > test/output.fasta"
				<< endl << endl << "Description:" << endl
				<< "       Align sequences in MFAFILE(s) and print result to standard output"
				<< endl << endl << "       -stockholm" << endl
				<< "              use STOCKHOLM output format instead of MFA"
				<< endl << endl << "       -clustalw" << endl
				<< "              use CLUSTALW output format instead of MFA"
				<< endl << endl << "       -noss" << endl
				<< "              do not show the consensus secondary structure"
				<< endl << endl << "       -columnscore" << endl
				<< "              write annotation for multiple alignment"
				<< endl << endl << "       -v, --verbose" << endl
				<< "              report progress while aligning (default: "
				<< (enableVerbose ? "on" : "off") << ")" << endl << endl
				<< "       -a, --alignment-order" << endl
				<< "              print sequences in alignment order rather than"
				<< endl << "              input order (default: "
				<< (enableAlignOrder ? "on" : "off") << ")" << endl << endl
				<< "       -c, --consistency REPS" << endl
				<< "              use " << MIN_INTERCONSISTENCY_REPS
				<< " <= REPS <= " << MAX_INTERCONSISTENCY_REPS << " (default: "
				<< numInterConsistencyReps << ") passes of" << endl
				<< "              inter-sequence consistency transformation"
				<< endl << endl << "       -ic, --intraconsistency REPS"
				<< endl << "              use " << MIN_INTERCONSISTENCY_REPS
				<< " <= REPS <= " << MAX_INTERCONSISTENCY_REPS << " (default: "
				<< numIntraConsistencyReps << ") passes of" << endl
				<< "              intra-sequence consistency transformation"
				<< endl << endl << "       -al, --alpha" << endl
				<< "              Sets the weight parameter (alpha) for "<<endl
				<< "              intra-sequence consistency transformation "
				<< endl << "              (0<alpha<1 , default: " << ALPHA
				<< ")" << endl << "       -bt, --beta" << endl
				<< "              Sets the weight parameter (beta) for "<<endl
				<< "              Four-Way consistency transformation " << endl
				<< "              (0<beta<1 , default: " << BETA << ")" << endl
				<< "       -Tb" << endl
				<< "              Sets the threshold value for identifying"
				<<endl
				<< "              the most probale base-pairing probabilities "
				<<endl << "              (0<Tb<1 , default: " << Tb << ")"
				<< endl << "       -r, --refinement REPS" << endl
				<< "              use " << MIN_REFINEMENTSREALIGNS
				<< " <= REPS <= " << MIN_REFINEMENTSREALIGNS << " (default: "
				<< numRefinementsReAligns << ") passes of " << endl
				<< "              refinement realignments" << endl << endl
				<< "       -p, --paramfile FILENAME" << endl
				<< "             read parameters for Pai-HMM from FILENAME "
				<< endl << "              (default: "
				<< parametersInputFilename << ")"<< endl << endl;
		exit(1);
	}

	SafeVector<string> sequenceNames;

	float tempFloat;
	int tempInt;

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {

			// clustalw output format
			if (!strcmp(argv[i], "-stockholm")) {
				enableStockholmOutput = true;
			}

			else if (!strcmp(argv[i], "-clustalw")) {
				enableClustalWOutput = true;
			}

			else if (!strcmp(argv[i], "-noss")) {
				NoSS = true;
			}

			// generate only column scores
			else if (!strcmp(argv[i], "-columnscore")) {
				enableAnnotation = true;
			}

			// verbose reporting
			else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
				enableVerbose = true;
			}

			// alignment order
			else if (!strcmp(argv[i], "-a") || !strcmp(argv[i],
					"--alignment-order")) {
				enableAlignOrder = true;
			}

			// number of inter-sequence consistency transformations
			else if (!strcmp(argv[i], "-c")
					|| !strcmp(argv[i], "--consistency")) {
				if (i < argc - 1) {
					if (!GetInteger(argv[++i], &tempInt)) {
						cerr << "ERROR: Invalid integer following option "
								<< argv[i - 1] << ": " << argv[i] << endl;
						exit(1);
					} else {
						if (tempInt < MIN_INTERCONSISTENCY_REPS || tempInt
>						MAX_INTERCONSISTENCY_REPS) {
							cerr << "ERROR: For option " << argv[i - 1]
							<< ", integer must be between "
							<< MIN_INTERCONSISTENCY_REPS << " and "
							<< MAX_INTERCONSISTENCY_REPS << "." << endl;
							exit(1);
						} else
						numInterConsistencyReps = tempInt;
					}
				} else {
					cerr << "ERROR: Integer expected for option " << argv[i]
					<< endl;
					exit(1);
				}
			}

			// number of inter-sequence consistency transformations
			else if (!strcmp(argv[i], "-ic")
			|| !strcmp(argv[i], "--intraconsistency")) {
				if (i < argc - 1) {
					if (!GetInteger(argv[++i], &tempInt)) {
						cerr << "ERROR: Invalid integer following option "
						<< argv[i - 1] << ": " << argv[i] << endl;
						exit(1);
					} else {
						if (tempInt < MIN_INTRACONSISTENCY_REPS || tempInt
>						MAX_INTRACONSISTENCY_REPS) {
							cerr << "ERROR: For option " << argv[i - 1]
							<< ", integer must be between "
							<< MIN_INTRACONSISTENCY_REPS << " and "
							<< MAX_INTRACONSISTENCY_REPS << "." << endl;
							exit(1);
						} else
						numIntraConsistencyReps = tempInt;
					}
				} else {
					cerr << "ERROR: Integer expected for option " << argv[i]
					<< endl;
					exit(1);
				}
			}

			else if (!strcmp(argv[i], "-al")
			|| !strcmp(argv[i], "--alpha")) {
				if (i < argc - 1) {
					if (!GetFloat(argv[++i], &tempFloat)) {
						cerr << "ERROR: Invalid floating-point value following option "
						<< argv[i - 1] << ": " << argv[i] << endl;
						exit(1);
					} else {
						if (tempFloat < 0 || tempFloat>1) {
							cerr << "ERROR: For option " << argv[i - 1]
							<< ", floating point value must be between 0 and 1."<< endl;
							exit(1);
						} else
						ALPHA = tempFloat;
					}
				} else {
					cerr << "ERROR: Floating-point value expected for option "
					<< argv[i] << endl;
					exit(1);
				}
			}

			else if (!strcmp(argv[i], "-bt")
			|| !strcmp(argv[i], "--beta")) {
				if (i < argc - 1) {
					if (!GetFloat(argv[++i], &tempFloat)) {
						cerr << "ERROR: Invalid floating-point value following option "
						<< argv[i - 1] << ": " << argv[i] << endl;
						exit(1);
					} else {
						if (tempFloat < 0 || tempFloat>1) {
							cerr << "ERROR: For option " << argv[i - 1]
							<< ", floating point value must be between 0 and 1." <<endl;
							exit(1);
						} else
						BETA = tempFloat;
					}
				} else {
					cerr << "ERROR: Floating-point value expected for option "
					<< argv[i] << endl;
					exit(1);
				}
			}

			else if (!strcmp(argv[i], "-Tb")) {
				if (i < argc - 1) {
					if (!GetFloat(argv[++i], &tempFloat)) {
						cerr << "ERROR: Invalid floating-point value following option "
						<< argv[i - 1] << ": " << argv[i] << endl;
						exit(1);
					} else {
						if (tempFloat < 0 || tempFloat>1) {
							cerr << "ERROR: For option " << argv[i - 1]
							<< ", floating point value must be between 0 and 1."<< endl;
							exit(1);
						} else
						Tb = tempFloat;
					}
				} else {
					cerr << "ERROR: Floating-point value expected for option "
					<< argv[i] << endl;
					exit(1);
				}
			}

			// number of realignments in refinement passes
			else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--refinement")) {
				if (i < argc - 1) {
					if (!GetInteger(argv[++i], &tempInt)) {
						cerr << "ERROR: Invalid integer following option "
						<< argv[i - 1] << ": " << argv[i] << endl;
						exit(1);
					} else {
						if (tempInt < MIN_REFINEMENTSREALIGNS || tempInt
>						MAX_REFINEMENTSREALIGNS) {
							cerr << "ERROR: For option " << argv[i - 1]
							<< ", integer must be between "
							<< MIN_REFINEMENTSREALIGNS << " and "
							<< MAX_REFINEMENTSREALIGNS << "." << endl;
							exit(1);
						} else
						numRefinementsReAligns = tempInt;
					}
				} else {
					cerr << "ERROR: Integer expected for option " << argv[i]
					<< endl;
					exit(1);
				}
			}
			else if (!strcmp(argv[i], "-p")
			|| !strcmp(argv[i], "--paramfile")) {
				if (i < argc - 1)
				parametersInputFilename = string(argv[++i]);
				else {
					cerr << "ERROR: Filename expected for option " << argv[i]
					<< endl;
					exit(1);
				}
			}

			// bad arguments
			else {
				cerr << "ERROR: Unrecognized option: " << argv[i] << endl;
				exit(1);
			}

		} else {
			sequenceNames.push_back(string(argv[i]));
		}
	}
	if (sequenceNames.size() == 0) {
		cerr << "ERROR: Please specify the input filename" << endl;
		exit(1);
	}
	return sequenceNames;

}

/////////////////////////////////////////////////////////////////
// ReadPHMMParameters()
//
// Read initial distribution, transition, and emission
// parameters for posterior alignment probabilities (PHMM)
/////////////////////////////////////////////////////////////////

void ReadPHMMParameters() {

	ifstream data;

	emitPairs = VVF(256, VF(256, 1e-10));
	emitSingle = VF(256, 1e-5);

	// read initial state distribution and transition parameters
	if (parametersInputFilename == string("")) {
		if (NumInsertStates == 1) {
			for (int i = 0; i < NumMatrixTypes; i++)
				initDistrib[i] = initDistrib1Default[i];
			for (int i = 0; i < 2 * NumInsertStates; i++)
				gapOpen[i] = gapOpen1Default[i];
			for (int i = 0; i < 2 * NumInsertStates; i++)
				gapExtend[i] = gapExtend1Default[i];
		} else if (NumInsertStates == 2) {
			for (int i = 0; i < NumMatrixTypes; i++)
				initDistrib[i] = initDistrib2Default[i];
			for (int i = 0; i < 2 * NumInsertStates; i++)
				gapOpen[i] = gapOpen2Default[i];
			for (int i = 0; i < 2 * NumInsertStates; i++)
				gapExtend[i] = gapExtend2Default[i];
		} else {
			cerr
					<< "ERROR: No default initial distribution/parameter settings exist"
					<< endl << "       for " << NumInsertStates
					<< " pairs of insert states.  Use --paramfile." << endl;
			exit(1);
		}

		alphabet = alphabetDefault;

		for (int i = 0; i < (int) alphabet.length(); i++) {
			emitSingle[(unsigned char) tolower(alphabet[i])] = emitSingleDefault[i];
			emitSingle[(unsigned char) toupper(alphabet[i])] = emitSingleDefault[i];
			for (int j = 0; j <= i; j++) {
				emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) tolower(
						alphabet[j])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) toupper(
						alphabet[j])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) tolower(
						alphabet[j])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) toupper(
						alphabet[j])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) tolower(
						alphabet[i])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) toupper(
						alphabet[i])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) tolower(
						alphabet[i])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) toupper(
						alphabet[i])] = emitPairsDefault[i][j];
			}
		}
	} else {
		data.open(parametersInputFilename.c_str());
		if (data.fail()) {
			cerr << "ERROR: Unable to read parameter file: "
					<< parametersInputFilename << endl;
			exit(1);
		}

		string line[3];
		for (int i = 0; i < 3; i++) {
			if (!getline(data, line[i])) {
				cerr
						<< "ERROR: Unable to read transition parameters from parameter file: "
						<< parametersInputFilename << endl;
				exit(1);
			}
		}
		istringstream data2;
		data2.clear();
		data2.str(line[0]);
		for (int i = 0; i < NumMatrixTypes; i++)
			data2 >> initDistrib[i];
		data2.clear();
		data2.str(line[1]);
		for (int i = 0; i < 2 * NumInsertStates; i++)
			data2 >> gapOpen[i];
		data2.clear();
		data2.str(line[2]);
		for (int i = 0; i < 2 * NumInsertStates; i++)
			data2 >> gapExtend[i];

		if (!getline(data, line[0])) {
			cerr << "ERROR: Unable to  alphabet from scoring matrix file: "
					<< parametersInputFilename << endl;
			exit(1);
		}

		// read alphabet as concatenation of all characters on alphabet line
		alphabet = "";
		string token;
		data2.clear();
		data2.str(line[0]);
		while (data2 >> token)
			alphabet += token;

		for (int i = 0; i < (int) alphabet.size(); i++) {
			for (int j = 0; j <= i; j++) {
				float val;
				data >> val;
				emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) tolower(
						alphabet[j])] = val;
				emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) toupper(
						alphabet[j])] = val;
				emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) tolower(
						alphabet[j])] = val;
				emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) toupper(
						alphabet[j])] = val;
				emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) tolower(
						alphabet[i])] = val;
				emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) toupper(
						alphabet[i])] = val;
				emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) tolower(
						alphabet[i])] = val;
				emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) toupper(
						alphabet[i])] = val;
			}
		}

		for (int i = 0; i < (int) alphabet.size(); i++) {
			float val;
			data >> val;
			emitSingle[(unsigned char) tolower(alphabet[i])] = val;
			emitSingle[(unsigned char) toupper(alphabet[i])] = val;
		}
		data.close();
	}
}

/////////////////////////////////////////////////////////////////
// PrintPHMMParameters()
//
// Print parameters for posterior alignment probabilities (PHMM) 
// to STDERR.  If a filename is specified, then the parameters 
// are also written to the file.
/////////////////////////////////////////////////////////////////

void PrintPHMMParameters(const char *message, const VF &initDistrib,
		const VF &gapOpen, const VF &gapExtend, const VVF &emitPairs,
		const VF &emitSingle, const char *filename) {

	// print parameters to the screen
	cerr << message << endl << "    initDistrib[] = { ";
	for (int i = 0; i < NumMatrixTypes; i++)
		cerr << setprecision(10) << initDistrib[i] << " ";
	cerr << "}" << endl << "        gapOpen[] = { ";
	for (int i = 0; i < NumInsertStates * 2; i++)
		cerr << setprecision(10) << gapOpen[i] << " ";
	cerr << "}" << endl << "      gapExtend[] = { ";
	for (int i = 0; i < NumInsertStates * 2; i++)
		cerr << setprecision(10) << gapExtend[i] << " ";
	cerr << "}" << endl << endl;

	// if a file name is specified
	if (filename) {

		// attempt to open the file for writing
		FILE *file = fopen(filename, "w");
		if (!file) {
			cerr << "ERROR: Unable to write parameter file: " << filename
					<< endl;
			exit(1);
		}

		// if successful, then write the parameters to the file
		for (int i = 0; i < NumMatrixTypes; i++)
			fprintf(file, "%.10f ", initDistrib[i]);
		fprintf(file, "\n");
		for (int i = 0; i < 2 * NumInsertStates; i++)
			fprintf(file, "%.10f ", gapOpen[i]);
		fprintf(file, "\n");
		for (int i = 0; i < 2 * NumInsertStates; i++)
			fprintf(file, "%.10f ", gapExtend[i]);
		fprintf(file, "\n");
		fprintf(file, "%s\n", alphabet.c_str());
		for (int i = 0; i < (int) alphabet.size(); i++) {
			for (int j = 0; j <= i; j++)
				fprintf(file, "%.10f ", emitPairs[(unsigned char) alphabet[i]][(unsigned char) alphabet[j]]);
			fprintf(file, "\n");
		}
		for (int i = 0; i < (int) alphabet.size(); i++)
			fprintf(file, "%.10f ", emitSingle[(unsigned char) alphabet[i]]);
		fprintf(file, "\n");
		fclose(file);
	}
}

/////////////////////////////////////////////////////////////////
// AlignAlignments()
//
// Returns the alignment of two MultiSequence objects.
/////////////////////////////////////////////////////////////////

MultiSequence *AlignAlignments(MultiSequence * align1, MultiSequence * align2,
		const SafeVector<SafeVector<SparseMatrix *> >&sparseMatrices,
		const ProbabilisticModel &model) {

	// print some info about the alignment
	if (enableVerbose) {
		for (int i = 0; i < align1->GetNumSequences(); i++)
			cerr << ((i == 0) ? "[" : ",") << align1->GetSequence(i)-> GetLabel();
		cerr << "] vs. ";
		for (int i = 0; i < align2->GetNumSequences(); i++)
			cerr << ((i == 0) ? "[" : ",") << align2->GetSequence(i)-> GetLabel();
		cerr << "]: ";
	}

	VF * posterior(NULL);

	posterior = model.BuildPosterior(align1, align2, sparseMatrices, cutoff);

	pair<SafeVector<char>*, float> alignment;
	if (gapOpenPenalty == 0 && gapContinuePenalty == 0)
		alignment = model.ComputeAlignment(align1->GetSequence(0)->GetLength(), align2->GetSequence(0)->GetLength(), *posterior);
	else
		alignment = model.ComputeAlignmentWithGapPenalties(align1, align2,
				*posterior, align1->GetNumSequences(),
				align2->GetNumSequences(), gapOpenPenalty, gapContinuePenalty);

	// choose the alignment routine depending on the "cosmetic" gap penalties used


	delete posterior;

	if (enableVerbose) {
		// compute total length of sequences
		int totLength = 0;
		for (int i = 0; i < align1->GetNumSequences(); i++)
			for (int j = 0; j < align2->GetNumSequences(); j++)
				totLength += min(align1->GetSequence(i)->GetLength(), align2->GetSequence(j)->GetLength());

		// give an "accuracy" measure for the alignment
		cerr << alignment.second / totLength << endl;
	}

	// now build final alignment
	MultiSequence *result = new MultiSequence();
	for (int i = 0; i < align1->GetNumSequences(); i++)
		result->AddSequence(align1->GetSequence(i)-> AddGaps(alignment.first, 'X'));
	for (int i = 0; i < align2->GetNumSequences(); i++)
		result->AddSequence(align2->GetSequence(i)-> AddGaps(alignment.first, 'Y'));
	if (!enableAlignOrder)
		result->SortByLabel();

	// free temporary alignment
	delete alignment.first;

	return result;
}

/////////////////////////////////////////////////////////////////
// DoRelaxation()
//
// Performs one round of the inter-sequence consistency 
// transformation.  
/////////////////////////////////////////////////////////////////

SafeVector<SafeVector<SparseMatrix *> > DoRelaxation(MultiSequence * sequences,
		SafeVector<SafeVector<SparseMatrix *> >&sparseMatrices, VVF distances,
		bool nflag) {
	const int numSeqs = sequences->GetNumSequences();

	SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices(numSeqs,
			SafeVector<SparseMatrix *>(numSeqs, NULL));

	// for every pair of sequences
	for (int i = 0; i < numSeqs; i++) {
		for (int j = i + 1; j < numSeqs; j++) {
			Sequence *seq1 = sequences->GetSequence(i);
			Sequence *seq2 = sequences->GetSequence(j);

			if (enableVerbose)
				cerr << "Relaxing (" << i + 1 << ") " << seq1->GetHeader()
						<< " vs. " << "(" << j + 1 << ") " << seq2->GetHeader()
						<< ": ";

			// get the original posterior matrix
			VF *posteriorPtr = sparseMatrices[i][j]->GetPosterior();
			assert(posteriorPtr);
			VF & posterior = *posteriorPtr;

			const int seq1Length = seq1->GetLength();
			const int seq2Length = seq2->GetLength();

			// contribution from the summation where z = x and z = y
			for (int k = 0; k < (seq1Length + 1) * (seq2Length + 1); k++) {
				posterior[k] += posterior[k];
				if (nflag)
					posterior[k] *= distances[j][i];
			}

			if (enableVerbose)
				cerr << sparseMatrices[i][j]->GetNumCells() << " --> ";

			float id_xyz = 0;
			// contribution from all other sequences
			for (int k = 0; k < numSeqs; k++)
				if (k != i && k != j) {
					if (k < i)
						Relax1(sparseMatrices[k][i], sparseMatrices[k][j],
								posterior, distances[k][i], distances[k][j],
								nflag);
					else if (k > i && k < j)
						Relax(sparseMatrices[i][k], sparseMatrices[k][j],
								posterior, distances[k][i], distances[k][j],
								nflag);
					else {
						SparseMatrix *temp =
								sparseMatrices[j][k]->ComputeTranspose();
						Relax(sparseMatrices[i][k], temp, posterior,
								distances[k][i], distances[k][j], nflag);
						delete temp;
					}
					id_xyz += distances[k][j] * distances[k][i];
				}
			// now renormalization
			id_xyz += 2 * distances[i][j];
			for (int k = 0; k < (seq1Length + 1) * (seq2Length + 1); k++) {
				if (nflag) {
					posterior[k] /= id_xyz;
				} else {
					posterior[k] /= (float) (numSeqs);

				}
			}

			// save the new posterior matrix
			newSparseMatrices[i][j] = new SparseMatrix(seq1->GetLength(),
					seq2->GetLength(), posterior, .005);
			newSparseMatrices[j][i] = NULL;

			if (enableVerbose)
				cerr << newSparseMatrices[i][j]->GetNumCells() << " -- ";

			delete posteriorPtr;

			if (enableVerbose)
				cerr << "done." << endl;
		}
	}

	return newSparseMatrices;
}

/////////////////////////////////////////////////////////////////
// Relax()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////

void Relax(SparseMatrix * matXZ, SparseMatrix * matZY, VF & posterior,
		float id_xz, float id_zy, bool nflag) {

	//use the conventional consistency transformation
	if (!nflag) {
		id_xz = 1;
		id_zy = 1;
	}
	float id_xyz = id_xz * id_zy;
	assert(matXZ);
	assert(matZY);

	int lengthX = matXZ->GetSeq1Length();
	int lengthY = matZY->GetSeq2Length();
	assert(matXZ->GetSeq2Length() == matZY->GetSeq1Length());

	// for every x[i]
	for (int i = 1; i <= lengthX; i++) {
		SafeVector<PIF>::iterator XZptr = matXZ->GetRowPtr(i);
		SafeVector<PIF>::iterator XZend = XZptr + matXZ->GetRowSize(i);

		VF::iterator base = posterior.begin() + i * (lengthY + 1);

		// iterate through all x[i]-z[k]
		while (XZptr != XZend) {
			SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
			SafeVector<PIF>::iterator ZYend = ZYptr
					+ matZY->GetRowSize(XZptr->first);
			const float XZval = XZptr->second;

			// iterate through all z[k]-y[j]
			while (ZYptr != ZYend) {
				base[ZYptr->first] += XZval * ZYptr->second * id_xyz;
				ZYptr++;
			}
			XZptr++;
		}
	}
}

/////////////////////////////////////////////////////////////////
// Relax1()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////

void Relax1(SparseMatrix * matZX, SparseMatrix * matZY, VF & posterior,
		float id_xz, float id_zy, bool nflag) {

	//use the conventional consistency transformation
	if (!nflag) {
		id_xz = 1;
		id_zy = 1;
	}
	float id_xyz = id_xz * id_zy;

	assert(matZX);
	assert(matZY);

	int lengthZ = matZX->GetSeq1Length();
	int lengthY = matZY->GetSeq2Length();

	// for every z[k]
	for (int k = 1; k <= lengthZ; k++) {
		SafeVector<PIF>::iterator ZXptr = matZX->GetRowPtr(k);
		SafeVector<PIF>::iterator ZXend = ZXptr + matZX->GetRowSize(k);

		// iterate through all z[k]-x[i]
		while (ZXptr != ZXend) {
			SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(k);
			SafeVector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(k);
			const float ZXval = ZXptr->second;
			VF::iterator base = posterior.begin() + ZXptr->first
					* (lengthY + 1);

			// iterate through all z[k]-y[j]
			while (ZYptr != ZYend) {
				base[ZYptr->first] += ZXval * ZYptr->second * id_xyz;
				ZYptr++;
			}
			ZXptr++;
		}
	}
}

///////////////////////////////////////////////////////////////
// DoBasePairProbabilityRelaxation()
//
// Performs one round of the intra-sequence consistency 
// transformation.  
/////////////////////////////////////////////////////////////////

SafeVector<BPPMatrix*> DoBasePairProbabilityRelaxation(MultiSequence *sequences,
		SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, SafeVector<
		BPPMatrix*> &BPPMatrices, VVF distances) {

	SafeVector<BPPMatrix*> newBPPMatrices;

	const int numSeqs = sequences->GetNumSequences();

	VF id_x(numSeqs, 0);
	for (int i = 0; i < numSeqs; i++) {
		for (int j = 0; j < numSeqs; j++) {
			if (i != j) {
				id_x[i] += distances[i][j];
			}
		}
	}

	// for each sequence
	for (int i = 0; i < numSeqs; i++) {

		Sequence *seq1 = sequences->GetSequence(i);
		BPPMatrix *seq1BppMatrix = BPPMatrices[seq1->GetLabel()];
		Trimat<float> consBppMat(seq1->GetLength() + 1);
		int seq1Length = seq1->GetLength();

		// contribution from all other sequences
		for (int j = 0; (j < numSeqs); j++) {
			if ((id_x[i] != 0) & (i != j)) {
				float dij = distances[i][j];
				Sequence *seq2 = sequences->GetSequence(j);
				BPPMatrix *seq2BppMatrix = BPPMatrices[seq2->GetLabel()];

				SparseMatrix *matchProb;
				if (i > j) {
					matchProb = sparseMatrices[j][i];
				} else {
					matchProb = sparseMatrices[i][j]->ComputeTranspose();
				}

				float M1 = dij / id_x[i] * (1-ALPHA);

				float sumProb = 0;
				vector<PIF2> &probs2 = seq2BppMatrix->bppMat.data2;
				for (int l = 0; l < (int) probs2.size(); l++) {
					float tmpProb2 = probs2[l].prob;
					int tmp2I = probs2[l].i;
					int tmp2J = probs2[l].j;

					float M2=tmpProb2 * M1;
					SafeVector<PIF>::iterator XYptri =
							matchProb->GetRowPtr(tmp2I);
					SafeVector<PIF>::iterator XYendi = XYptri
							+ matchProb->GetRowSize(tmp2I);
					while (XYptri != XYendi) {
						SafeVector<PIF>::iterator XYptrj =
								matchProb->GetRowPtr(tmp2J);
						SafeVector<PIF>::iterator XYendj = XYptrj
								+ matchProb->GetRowSize(tmp2J);
						while (XYptrj != XYendj) {
							if (XYptrj->first > XYptri->first) {
								consBppMat.ref(XYptri->first, XYptrj->first)
										+= XYptri->second * XYptrj->second * M2;
							}
							XYptrj++;
						}
						XYptri++;
					}
				}
			}
		}
		for (int k = 1; k <= seq1Length; k++) {

			SafeVector<PIF>::iterator XptrB = (seq1BppMatrix->bppMat).GetRowPtr(k);
			SafeVector<PIF>::iterator XptrBend = XptrB + (seq1BppMatrix->bppMat).GetRowSize(k);

			while (XptrB != XptrBend) {
				consBppMat.ref(k, XptrB->first) += ALPHA * XptrB->second;
				XptrB++;
			}
		}

		BPPMatrix *tmp_bppmat = new BPPMatrix(seq1Length, BASEPROBTHRESHOLD,
				consBppMat);
		newBPPMatrices.push_back(tmp_bppmat);
	}

	return newBPPMatrices;
}

///////////////////////////////////////////////////////////////
// DoFourWayRelaxation()
//
// Performs one round of the Four-Way consistency 
// transformation.  
/////////////////////////////////////////////////////////////////

SafeVector<SafeVector<SparseMatrix *> > DoFourWayRelaxation(
		MultiSequence * sequences,
		SafeVector<SafeVector<SparseMatrix *> >&sparseMatrices, SafeVector<
		BPPMatrix*> &BPPMatrices, VVF distances, bool nflag, int RR) {
	const int numSeqs = sequences->GetNumSequences();

	SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices(numSeqs,
			SafeVector<SparseMatrix *>(numSeqs, NULL));

	float TH=.005 / pow((double)3, (double)RR + 1);

	// for every pair of sequences
	for (int i = 0; i < numSeqs; i++) {
		for (int j = i + 1; j < numSeqs; j++) {
			Sequence *seq1 = sequences->GetSequence(i);
			Sequence *seq2 = sequences->GetSequence(j);

			// get the original posterior matrix
			VF *posteriorPtr = sparseMatrices[i][j]->GetPosterior();
			assert(posteriorPtr);
			VF & posterior = *posteriorPtr;

			const int seq1Length = seq1->GetLength();
			const int seq2Length = seq2->GetLength();

			for (int k = 0; k < (seq1Length + 1) * (seq2Length + 1); k++) {
				int ii = (int) floor(k / (seq2Length + 1));
				int jj = k % (seq2Length + 1);
				if (!((ii == 0) | (jj == 0)))
					posterior[k] *= BETA;

			}

			// contribution from all other sequences
			RelaxBP(sparseMatrices[i][j], BPPMatrices[i], BPPMatrices[j],
					posterior, i, j, distances[i][j], nflag, 1);

			//mask out positions not originally in the posterior matrix
			SparseMatrix *matXY = sparseMatrices[i][j];
			for (int y = 0; y <= seq2Length; y++)
				posterior[y] = 0;
			for (int x = 1; x <= seq1Length; x++) {
				SafeVector<PIF>::iterator XYptr = matXY->GetRowPtr(x);
				SafeVector<PIF>::iterator XYend = XYptr + matXY->GetRowSize(x);
				VF::iterator base = posterior.begin() + x * (seq2Length + 1);
				int curr = 0;
				while (XYptr != XYend) {

					// zero out all cells until the first filled column
					while (curr < XYptr->first) {
						base[curr] = 0;
						curr++;
					}
					// now, skip over this column
					curr++;
					++XYptr;
				}

				// zero out cells after last column
				while (curr <= seq2Length) {
					base[curr] = 0;
					curr++;
				}
			}

			// save the new posterior matrix
			newSparseMatrices[i][j] = new SparseMatrix(seq1->GetLength(),
					seq2->GetLength(), posterior, TH);
			newSparseMatrices[j][i] = NULL;

			if (enableVerbose)
				cerr << newSparseMatrices[i][j]->GetNumCells() << " -- ";

			delete posteriorPtr;

			if (enableVerbose)
				cerr << "done." << endl;
		}
	}

	return newSparseMatrices;
}

/////////////////////////////////////////////////////////////////
// RelaxBp()
//
// Computes the Four-Way consistency transformation for sequences 
// x and z.
/////////////////////////////////////////////////////////////////

void RelaxBP(SparseMatrix * matXZ, BPPMatrix* bpmatX, BPPMatrix* bpmatZ,
		VF & posterior, int Seqx, int Seqz, float id_zy, bool nflag, float mm) {

	VF *ppp = matXZ->GetPosterior();
	assert(ppp);
	VF & pp = *ppp;

	//use the conventional consistency transformation
	if (!nflag) {
		id_zy = 1;
	}
	assert(matXZ);

	int lengthX = matXZ->GetSeq1Length();
	int lengthZ = matXZ->GetSeq2Length();

	int aaa = 20000;


	vector<PIF2> &probs1 = bpmatX->bppMat.data2;
	vector<PIF2> &probs2 = bpmatZ->bppMat.data2;
	for (int i = 0; i < (int) probs1.size(); i++) {
		float tmpProb1 = probs1[i].prob;
		int tmp1I = probs1[i].i;
		int tmp1J = probs1[i].j;
		float M1=(1-BETA) * id_zy*tmpProb1;

		VF::iterator base = posterior.begin() + tmp1I * (lengthZ + 1);
		VF::iterator base2 = posterior.begin() + tmp1J * (lengthZ + 1);
		SafeVector<PIF>::iterator XZptri = matXZ->GetRowPtr(tmp1I);
		SafeVector<PIF>::iterator XZendi = XZptri + matXZ->GetRowSize(tmp1I);
		while (XZptri != XZendi) {

			SafeVector<PIF>::iterator XZptrj = matXZ->GetRowPtr(tmp1J);
			SafeVector<PIF>::iterator XZendj = XZptrj + matXZ->GetRowSize(tmp1J);
			while (XZptrj != XZendj) {
				base[XZptri->first] += bpmatZ->GetProb(XZptri->first,
						XZptrj->first) * XZptrj->second * M1;
				base2[XZptrj->first] += bpmatZ->GetProb(XZptri->first,
						XZptrj->first) * XZptri->second * M1;
				XZptrj++;
			}
			XZptri++;
		}
	}
}

//////////////////////////////////////////////////////////////////
// DoRefinement()
//
// Performs the refinement step
//////////////////////////////////////////////////////////////////

void DoRefinement(MultiSequence* &alignment, SafeVector<SafeVector<
SparseMatrix *> > &sparseMatrices, const ProbabilisticModel &model, VVF distances) {

	int NumSeqs = alignment->GetNumSequences();
	vector<set<int> > SimSeqs;

	FindSimilar(distances, SimSeqs);

	int cnt = 0;

	// Performe refinement for at least numRefinementsReAligns number of realignments
	while (cnt < (numRefinementsReAligns)) {
		srand(time(0));

		VI list(NumSeqs);

		for (int i = 0; i < NumSeqs; i++)
			list[i] = i;
		VI rnd_list;

		// obtain a random ordering of sequences
		while (list.size() > 0) {
			int index = rand() % (list.size());
			rnd_list.push_back(list[index]);
			list.erase(list.begin() + index);
		}

		//For each sequence, update the set of similar sequences (S_x) and then
		//align that to set of dissimilar sequences (N_x)
		for (int i = 0; i < alignment->GetNumSequences(); i++) {
			int si = rnd_list[i];
			set<int> groupOne, groupTwo;
			groupOne = SimSeqs[si];

			// project sequences to the two groups (S_x,N_x)
			for (int j = 0; j < alignment->GetNumSequences(); j++)
				if (groupOne.find(j) == groupOne.end())
					groupTwo.insert(j);
			cnt++;
			if ((groupOne.size() != 0) & (groupTwo.size() != 0)) {

				MultiSequence *groupOneSeqs = alignment->Project(groupOne);
				assert (groupOneSeqs);
				MultiSequence *groupTwoSeqs = alignment->Project(groupTwo);
				assert (groupTwoSeqs);

				//				//Find x
				set<int>::iterator it;
				int cnnt = 0;
				for (it = groupOne.begin(); it != groupOne.end(); ++it) {
					if (*it == si)
						break;
					cnnt++;
				}

				//Update S_x by aligning x with (S_x - x)
				if (groupOneSeqs->GetNumSequences() > 1) {
					set<int> groupOne2, groupTwo2;
					groupOne2.insert(cnnt);

					for (int k = 0; k < groupOneSeqs->GetNumSequences(); k++)
						if (k != cnnt)
							groupTwo2.insert(k);

					MultiSequence *groupOneSeqs2 =
							groupOneSeqs->Project(groupOne2);
					assert (groupOneSeqs2);
					MultiSequence *groupTwoSeqs2 =
							groupOneSeqs->Project(groupTwo2);
					assert (groupTwoSeqs2);
					delete groupOneSeqs;
					// realign
					groupOneSeqs = AlignAlignments(groupOneSeqs2,
							groupTwoSeqs2, sparseMatrices, model);
					cnt++;
				}
				delete alignment;

				// realign the updated similar set (S'_x) and N_x
				alignment = AlignAlignments(groupOneSeqs, groupTwoSeqs,
						sparseMatrices, model);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////
// Findsimilar()
//
// Find similar sequences for each sequence using kmeans algorithm
//////////////////////////////////////////////////////////////////

void FindSimilar(VVF distances, vector<set<int> > &SimSeqs) {

	int NumSeqs = distances.size();
	for (int i = 0; i < NumSeqs; i++)
		distances[i][i] = 1;

	// For each sequence x find the set of similar sequences S_x
	for (int i = 0; i < NumSeqs; i++) {
		set<int> c1, c2;

		float min_d = 1;
		float max_d = 0;
		int ii_min = 0;
		int ii_max = 0;
		for (int j = 0; j < NumSeqs; j++) {
			if ((distances[i][j] <= min_d)) {
				ii_min = j;
				min_d = distances[i][j];
			}
			if ((distances[i][j] >= max_d)) {
				ii_max = j;
				max_d = distances[i][j];
			}
		}

		c1.insert(ii_max); // The first cluster: similar sequences (S_x)
		c2.insert(ii_min); // The second cluster: dissimilar sequences (N_x)

		//initiate the clusters
		for (int j = 0; j < NumSeqs; j++) {
			if ((j != ii_min) & (j != ii_max)) {
				if (abs(distances[j][i] - max_d) < abs(distances[j][i] - min_d))
					c1.insert(j);
				else
					c2.insert(j);
			}
		}

		if (c1.find(i) == c1.end()) {
			c2.erase(i);
			c1.insert(i);
		}
		bool ch_flag = true;
		int cnt = 0;

		// Iterate 100 times to obtain the clusters using kmeans algorithm
		while ((cnt < 100) & (ch_flag)) {
			ch_flag = false;
			VI changes(NumSeqs, 0);

			//compute center of each cluster
			float m1 = 0;
			float m2 = 0;
			set<int>::iterator it;
			for (it = c1.begin(); it != c1.end(); ++it)
				m1 += distances[i][*it];
			for (it = c2.begin(); it != c2.end(); ++it)
				m2 += distances[i][*it];
			m1 /= c1.size();
			m2 /= c2.size();

			//update the clusters
			for (int j = 0; j < NumSeqs; j++) {
				if (j != i) {
					set<int>::iterator it;
					if (c1.find(j) != c1.end()) {
						if (abs(distances[j][i] - m1) > abs(distances[j][i]
								- m2)) {
							changes[j] = 1;
							ch_flag = true;
						}
					} else {
						if (abs(distances[j][i] - m2) > abs(distances[j][i]
								- m1)) {
							changes[j] = -1;
							ch_flag = true;
						}
					}
				}
			}
			if (ch_flag) {
				for (int j = 0; j < NumSeqs; j++) {
					if (changes[j] == 1) {
						c1.erase(j);
						c2.insert(j);
					} else if (changes[j] == -1) {
						c2.erase(j);
						c1.insert(j);
					}
				}
			}
			cnt++;
		}
		SimSeqs.push_back(c1);
	}

}
/////////////////////////////////////////////////////////////////
// FindMaxPP()
//
//Find the maximum posterior pairwise base alignment probability
/////////////////////////////////////////////////////////////////

float FindMaxPP(SafeVector<SafeVector<SparseMatrix *> > sparseMatrices) {

	float m = 0;
	int cnt = 0;
	for (int i = 0; i < sparseMatrices.size(); i++) {
		for (int j = i + 1; j < sparseMatrices.size(); j++) {
			for (int h = 0; h < sparseMatrices[i][j]->GetSeq1Length(); h++) {
				for (int k = 0; k < sparseMatrices[i][j]->GetRowSize(h); k++) {

					float tmpProb2 = sparseMatrices[i][j]->GetRowPtr(h)[k].second;
					if (tmpProb2 > m)
						m = tmpProb2;
				}
			}
		}
	}
	return m;

}

/////////////////////////////////////////////////////////////////
// GetInteger()
//
// Attempts to parse an integer from the character string given.
// Returns true only if no parsing error occurs.
/////////////////////////////////////////////////////////////////

bool GetInteger(char *data, int *val) {
	char *endPtr;
	long int retVal;

	assert(val);

	errno = 0;
	retVal = strtol(data, &endPtr, 0);
	if (retVal == 0 && (errno != 0 || data == endPtr))
		return false;
	if (errno != 0 && (retVal == LONG_MAX || retVal == LONG_MIN))
		return false;
if	(retVal < (long) INT_MIN || retVal> (long) INT_MAX)
	return false;
	*val = (int) retVal;
	return true;
}

/////////////////////////////////////////////////////////////////
// GetFloat()
//
// Attempts to parse a float from the character string given.
// Returns true only if no parsing error occurs.
/////////////////////////////////////////////////////////////////

bool GetFloat(char *data, float *val) {
	char *endPtr;
	double retVal;

	assert(val);

	errno = 0;
	retVal = strtod(data, &endPtr);
	if (retVal == 0 && (errno != 0 || data == endPtr))
		return false;
	if (errno != 0 && (retVal >= 1000000.0 || retVal <= -1000000.0))
		return false;
	*val = (float) retVal;
	return true;
}

/////////////////////////////////////////////////////////////////
// WriteAnnotation()
//
// Computes annotation for multiple alignment
/////////////////////////////////////////////////////////////////

void WriteAnnotation(MultiSequence * alignment, const SafeVector<SafeVector<
SparseMatrix *> >&sparseMatrices) {
	float probprodct = 0;
	const int alignLength = alignment->GetSequence(0)->GetLength();
	const int numSeqs = alignment->GetNumSequences();
	int i, j;

	SafeVector<int> position(numSeqs, 0);
	SafeVector<SafeVector<char>::iterator> seqs(numSeqs);
	for (i = 0; i < numSeqs; i++)
		seqs[i] = alignment->GetSequence(i)->GetDataPtr();
	SafeVector<pair<int, int> > active;
	active.reserve(numSeqs);

	column = new columnReliability[alignLength + 1];
	column[0].columnNo = 0;

	cerr << numSeqs << endl << alignLength << endl;
	// for every column
	for (i = 1; i <= alignLength; i++) {
		//initialize the column reliability structure
		column[i].columnNo = i;
		column[i].probProduct = 0;
		// find all aligned residues in this particular column

		active.clear();
		for (j = 0; j < numSeqs; j++) {
			if (seqs[j][i] != '-') {

				active.push_back(make_pair(j, ++position[j]));

				if (enableVerbose)
					printf("\nposition[j]=%d\n", position[j]);
			}
		}

		probprodct = ComputeScore(active, sparseMatrices);
		column[i].probProduct = probprodct;

		if (enableVerbose)
			printf("\ncolumn %d %f\n--\n", i, probprodct);

	}

	printf("probabilities ..(row column)\n");
	for (i = 1; i < alignLength; i++)
		printf("%d %f\n", column[i].columnNo, column[i].probProduct);

	delete[] column;

}

/////////////////////////////////////////////////////////////////
// ComputeScore()
//
// Computes the annotation score for a particular column.
/////////////////////////////////////////////////////////////////

float ComputeScore(const SafeVector<pair<int, int> >&active, const SafeVector<
SafeVector<SparseMatrix *> >&sparseMatrices) {

	if (active.size() <= 1)
		return 0;

	// ALTERNATIVE #1: Compute the average alignment score.


	float prob_product = 0;

	for (int i = 0; i < (int) active.size(); i++) {
		for (int j = i + 1; j < (int) active.size(); j++) {

			prob_product
					+= sparseMatrices[active[i].first][active[j].first]-> GetValue(
							active[i].second, active[j].second);

			if (enableVerbose)
				printf(
						"%d-%d %d-%d %1.3f %f\n",
						active[i].first,
						active[i].second,
						active[j].first,
						active[j].second,
						sparseMatrices[active[i].first][active[j].first]->GetValue(
								active[i].second, active[j].second),
						prob_product);
		}
	}

	if (enableVerbose)
		printf("active size= %d \n", (int) active.size());

	return 2 * prob_product / ((int) active.size() * ((int) active.size() - 1));
}
