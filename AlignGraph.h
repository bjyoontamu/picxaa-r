////////////////////////////////////////////////////////////////
// AlignGraph.h
//
// Utilities for constructing the alignment graph.
/////////////////////////////////////////////////////////////////

#ifndef ALIGNGRAPH_H
#define ALIGNGRAPH_H

#include "SafeVector.h"
#include "MultiSequence.h"
#include "Sequence.h"
#include "BPPMatrix.hpp"

/////////////////////////////////////////////////////////////////
// AlignGraph
//
// Class for multiple sequence alignment input/output.
/////////////////////////////////////////////////////////////////

class AlignGraph {

	//Sequences in the alignment
	MultiSequence *sequences;

	//Deduced alignment
	MultiSequence *alignment;

	SafeVector<char> *Profile_alignment;

	//The alignment graph
	// G[i][j] is the j'th child of the i'th node
	VVI G;

	//Columns(nodes) in Graph
	// cols[i]: vector of residues in node i
	VVVI cols;

	// Matrix of residues in each node:
	// IsPresent[i][j] is -1 if  j'th residue in i'th sequence
	// has not appeared in any node and is m (>=1) if it is
	// present in the  m'th node
	VVI IsPresent;

	// Matrix of bases appearance in stem components:
	// IsStem[i][j] is -1 if  j'th base in i'th sequence
	// has not appeared in any stem and is m (>=1) if it is
	// present in the  m'th stem component	
	VVI IsStem;

	// Matrix of bases position in stem components:
	// LRStem[i][j] is FALSE if (IsStem[i][j]==-1) or if
	// j'th base in i'th sequence is left stem component 
	// and is TRUE if it is right stem component 
	SafeVector<SafeVector<bool> > LRStem;

	// Number of stem components so far
	int StemSize;

	VVVF BPProb;
	SafeVector<BPPMatrix*> BPPMatrices;
	SafeVector<SafeVector<SparseMatrix *> > sparseMatrices;

	//Matrices of Ancestors and Descendents:
	// Ancs[i][j]=true if node j is ancestor of node i;  otherwise it is false
	// Descs[i][j]=true if node j is descendant of node i;  otherwise it is false
	SafeVector<SafeVector<bool> > Ancs, Descs;

	//Zero Vector
	SafeVector<bool> ZZ;

private:

	/////////////////////////////////////////////////////////////////
	// Partition()
	//
	// Partition the set for Quick Sort
	/////////////////////////////////////////////////////////////////

	int Partition(int low, int high, float arr[], int ind[]) {
		float high_vac, low_vac, pivot;
		int il, ih, ip;
		ip = ind[low];
		il = ind[low];
		ih = ind[high];
		ih = ind[high];
		pivot = arr[low];

		while (high > low) {
			high_vac = arr[high];
			ih = ind[high];
			while (pivot <= high_vac) {
				if (high <= low)
					break;
				high--;
				high_vac = arr[high];
				ih = ind[high];
			}
			arr[low] = high_vac;
			ind[low] = ih;
			low_vac = arr[low];
			il = ind[low];
			while (pivot >= low_vac) {
				if (high <= low)
					break;
				low++;
				low_vac = arr[low];
				il = ind[low];
			}
			arr[high] = low_vac;
			ind[high] = il;
		}
		arr[low] = pivot;
		ind[low] = ip;
		return low;
	}
	/////////////////////////////////////////////////////////////////
	// Quick_sort()
	//
	// Quick Sort Algorithm
	/////////////////////////////////////////////////////////////////

	void Quick_sort(int low, int high, float arr[], int ind[]) {
		int Piv_index;
		if (low < high) {
			Piv_index = Partition(low, high, arr, ind);
			Quick_sort(low, Piv_index - 1, arr, ind);
			Quick_sort(Piv_index + 1, high, arr, ind);
		}
	}

	/////////////////////////////////////////////////////////////////
	// FindCloseNodes()
	//
	// Find the preceding and succeeding nodes in the graph G that
	// contain residues in the sequence of residue x
	/////////////////////////////////////////////////////////////////

	VVI FindCloseNodes(VI x) {

		VVI Sx;

		// The preceding node
		VI Parent;

		// The succeeding node
		VI Child;
		int pi;
		int ci;
		int clx = -1;
		int crx = 10000; //inf
		for (int i = x[1] - 1; i >= (int) 0; i--)
			if (IsPresent[x[0]][i] != -1) {
				pi = IsPresent[x[0]][i];
				clx = i;
				break;
			}
		for (int i = x[1] + 1; i < (int) IsPresent[x[0]].size(); i++)
			if (IsPresent[x[0]][i] != -1) {
				ci = IsPresent[x[0]][i];
				crx = i;
				break;
			}
		if (clx != -1)
			Parent.push_back(pi);
		if (crx != 10000)
			Child.push_back(ci);

		Sx.push_back(Parent);
		Sx.push_back(Child);

		return Sx;
	}

	/////////////////////////////////////////////////////////////////
	// GiveNodes()
	//
	// Give the ancestors or descendents of a given node
	/////////////////////////////////////////////////////////////////

	VI GiveNodes(SafeVector<bool> A, int msize) {

		VI B;
		for (int i = 0; i < msize; i++)
			if (A[i])
				B.push_back(i);
		return B;

	}

	/////////////////////////////////////////////////////////////////
	// Union_VBVI()
	//
	// Perform the union of set a (in the boolean format)
	// with set b  (in integer format)
	/////////////////////////////////////////////////////////////////

	SafeVector<bool> Union_VBVI(SafeVector<bool> A, VI B) {
		SafeVector<bool> C(A);
		for (int i = 0; i < (int) B.size(); i++) {
			C[B[i]] = true;
		}
		return C;
	}

	/////////////////////////////////////////////////////////////////
	// Union()
	//
	// Perform the union of sets a and b ( both in the boolean format)
	/////////////////////////////////////////////////////////////////

	SafeVector<bool> Union(SafeVector<bool> A, SafeVector<bool> B, int msize) {
		SafeVector<bool> C(A);
		for (int i = 0; i < msize; i++) {
			C[i] = A[i] | B[i];
		}
		return C;
	}

	/////////////////////////////////////////////////////////////////
	// Union_VI()
	//
	// Perform the union of sets a and b ( both in the integer format)
	/////////////////////////////////////////////////////////////////

	VI Union_VI(VI A, VI B) {

		for (int j = 0; j < (int) B.size(); j++) {
			bool flag = false;
			for (int k = 0; k < (int) A.size(); k++)
				if (A[k] == B[j]) {
					flag = true;
					break;
				}
			if (!flag)
				A.push_back(B[j]);
		}
		return A;
	}

	/////////////////////////////////////////////////////////////////
	// Ismember()
	//
	// Check whether i is in set A (boolean format)
	/////////////////////////////////////////////////////////////////

	bool Ismember(SafeVector<bool> A, int i) {
		return (A[i]);
	}

	/////////////////////////////////////////////////////////////////
	// Ismember_VI()
	//
	// Check whether i is in set A (integer format)
	/////////////////////////////////////////////////////////////////

	bool Ismember_VI(VI A, int i) {
		for (int j = 0; j < (int) A.size(); j++)
			if (A[j] == i)
				return true;
		return false;
	}

	/////////////////////////////////////////////////////////////////
	// Remove()
	//
	// Remove i from set A (integer format)
	/////////////////////////////////////////////////////////////////

	VI Remove(VI A, int i) {
		VI B;
		for (int j = 0; j < (int) A.size(); j++)
			if (A[j] != i)
				B.push_back(A[j]);
		return B;
	}

	/////////////////////////////////////////////////////////////////
	// Update()
	//
	// Update the set A after we merges to components cx and cy
	/////////////////////////////////////////////////////////////////

	SafeVector<bool> Update(SafeVector<bool> A, int cy, int msize) {

		for (int j = cy; j < msize - 1; j++)
			A[j] = A[j + 1];
		A[msize - 1] = false;
		return A;
	}

	/////////////////////////////////////////////////////////////
	// Update()
	//
	// Update the value of i after we merges to components cx and cy
	/////////////////////////////////////////////////////////////////

	int Update(int i, int cx, int cy) {
		if (i < cy)
			return i;
		else if (i == cy)
			return cx;
		else
			return i - 1;
	}

	/////////////////////////////////////////////////////////////
	// GiveParent()
	//
	// Find the parents of the node i in graph G
	/////////////////////////////////////////////////////////////////

	VI GiveParent(VVI G, int i) {
		VI a;
		for (int j = 0; j < (int) G.size(); j++) {
			for (int k = 0; k < (int) G[j].size(); k++) {
				if (G[j][k] == i) {
					a.push_back(j);
					break;
				}
			}
		}
		return a;
	}

	/////////////////////////////////////////////////////////////
	// CheckAddNewNode()
	//
	// Procedure to insert the new node in the graph provided that
	// it does not introduce any cycle
	/////////////////////////////////////////////////////////////////

	void InsertSingleColumn(VI x) {

		//Develop the temporary Graph to check for cycles
		VVI Temp_G = G;

		//Find the closest nodes for residues x and y
		VVI Sx = FindCloseNodes(x);

		// Find the set of Parents and Children
		VI Parent = Sx[0];
		VI Child = Sx[1];

		//Add the Parents and Children to the graph
		Temp_G.push_back(Child);

		for (int j = 0; j < (int) Parent.size(); j++)
			Temp_G[Parent[j]].push_back(G.size());

		// Remove redundant edges
		for (int j = 0; j < (int) Parent.size(); j++)
			for (int k = 0; k < (int) Child.size(); k++)
				Temp_G[Parent[j]] = Remove(Temp_G[Parent[j]], Child[k]);

		//Update IsPresent
		IsPresent[x[0]][x[1]] = G.size();

		//Update G
		G = Temp_G;

		int Gsz = G.size();

		//New Node Ancestors
		SafeVector<bool> curr_Anc(ZZ);
		if (Parent.size() > 0)
			curr_Anc = Ancs[Parent[0]];
		if (Parent.size() == 2)
			curr_Anc = Union(curr_Anc, Ancs[Parent[1]], Gsz - 1);
		curr_Anc = Union_VBVI(curr_Anc, Parent);
		//Update Ancestors for the new node
		Ancs.push_back(curr_Anc);

		//New Node descendants
		SafeVector<bool> curr_Desc(ZZ);
		if (Child.size() > 0)
			curr_Desc = Descs[Child[0]];
		if (Child.size() == 2)
			curr_Desc = Union(curr_Desc, Descs[Child[1]], Gsz - 1);
		curr_Desc = Union_VBVI(curr_Desc, Child);

		//Update Descendants for the new node
		Descs.push_back(curr_Desc);

		//Update Ancestors and Descendants for other nodes
		VI DD, AA;
		for (int j = 0; j < Gsz; j++) {
			if (curr_Anc[j])
				AA.push_back(j);
			if (curr_Desc[j])
				DD.push_back(j);
		}
		for (int j = 0; j < (int) DD.size(); j++) {
			Ancs[DD[j]][Gsz - 1] = true;
			for (int k = 0; k < (int) AA.size(); k++) {
				Ancs[DD[j]][AA[k]] = true;
				Descs[AA[k]][DD[j]] = true;
			}
		}
		for (int k = 0; k < (int) AA.size(); k++)
			Descs[AA[k]][Gsz - 1] = true;
	}
	/////////////////////////////////////////////////////////////
	// CheckAddNewNode()
	//
	// Procedure to insert the new node in the graph provided that
	// it does not introduce any cycle
	/////////////////////////////////////////////////////////////////

	bool CheckAddNewNode(VI x, VI y) {

		//Develop the temporary Graph to check for cycles
		VVI Temp_G = G;

		//Find the closest nodes for residues x and y
		VVI Sx = FindCloseNodes(x);
		VVI Sy = FindCloseNodes(y);

		// Find the set of Parents and Children
		VI Parent = Union_VI(Sx[0], Sy[0]);
		VI Child = Union_VI(Sx[1], Sy[1]);

		//Add the Parents and Children to the graph
		Temp_G.push_back(Child);

		for (int j = 0; j < (int) Parent.size(); j++)
			Temp_G[Parent[j]].push_back(G.size());

		//Check for cycle
		bool Check_Cycle = true;
		//Check whether one residue's parent is the descendent of the other one
		if ((Sx[0].size() == 1) & (Sy[1].size() == 1))
			Check_Cycle = Check_Cycle & !Descs[Sy[1][0]][Sx[0][0]] & !(Sx[0][0]
					== Sy[1][0]);
		if ((Sy[0].size() == 1) & (Sx[1].size() == 1))
			Check_Cycle = Check_Cycle & !Descs[Sx[1][0]][Sy[0][0]] & !(Sy[0][0]
					== Sx[1][0]);

		// If no cycle appears update the graph
		if (Check_Cycle) {

			// Remove redundant edges
			if ((Sy[0].size() == 1) & (Sx[0].size() == 1)) {
				if (Ismember(Descs[Sx[0][0]], Sy[0][0]))
					Temp_G[Sx[0][0]] = Remove(Temp_G[Sx[0][0]], G.size());
				if (Ismember(Descs[Sy[0][0]], Sx[0][0]))
					Temp_G[Sy[0][0]] = Remove(Temp_G[Sy[0][0]], G.size());
			}
			if ((Sy[1].size()) & (Sx[1].size())) {
				if (Ismember(Descs[Sx[1][0]], Sy[1][0]))
					Temp_G[G.size()] = Remove(Temp_G[G.size()], Sy[1][0]);
				if (Ismember(Descs[Sy[1][0]], Sx[1][0]))
					Temp_G[G.size()] = Remove(Temp_G[G.size()], Sx[1][0]);
			}
			for (int j = 0; j < (int) Parent.size(); j++)
				for (int k = 0; k < (int) Child.size(); k++)
					Temp_G[Parent[j]] = Remove(Temp_G[Parent[j]], Child[k]);

			//Update IsPresent
			IsPresent[x[0]][x[1]] = G.size();
			IsPresent[y[0]][y[1]] = G.size();

			//Update G
			G = Temp_G;

			int Gsz = G.size();

			//New Node Ancestors
			SafeVector<bool> curr_Anc(ZZ);
			if (Parent.size() > 0)
				curr_Anc = Ancs[Parent[0]];
			if (Parent.size() == 2)
				curr_Anc = Union(curr_Anc, Ancs[Parent[1]], Gsz - 1);
			curr_Anc = Union_VBVI(curr_Anc, Parent);

			//Update Ancestors for the new node
			Ancs.push_back(curr_Anc);

			//New Node descendants
			SafeVector<bool> curr_Desc(ZZ);
			if (Child.size() > 0)
				curr_Desc = Descs[Child[0]];
			if (Child.size() == 2)
				curr_Desc = Union(curr_Desc, Descs[Child[1]], Gsz - 1);
			curr_Desc = Union_VBVI(curr_Desc, Child);

			//Update Descendants for the new node
			Descs.push_back(curr_Desc);

			//Update Ancestors and Descendants for other nodes
			VI DD, AA;
			for (int j = 0; j < Gsz; j++) {
				if (curr_Anc[j])
					AA.push_back(j);
				if (curr_Desc[j])
					DD.push_back(j);
			}
			for (int j = 0; j < (int) DD.size(); j++) {
				Ancs[DD[j]][Gsz - 1] = true;
				for (int k = 0; k < (int) AA.size(); k++) {
					Ancs[DD[j]][AA[k]] = true;
					Descs[AA[k]][DD[j]] = true;
				}
			}
			for (int k = 0; k < (int) AA.size(); k++)
				Descs[AA[k]][Gsz - 1] = true;
		}
		return Check_Cycle;
	}

	/////////////////////////////////////////////////////////////
	// CheckAddColumnEx()
	//
	// Procedure to extend an existing column in the graph provided that
	// it does not introduce any cycle
	/////////////////////////////////////////////////////////////////

	bool CheckAddColumnEx(VI y, int cx) {

		//Develop the temporary Graph to check for cycles
		VVI Temp_G = G;

		//Find the closest nodes for residue y
		VVI Sy = FindCloseNodes(y);

		// Find the set of Parents and Children
		VI Parent = Sy[0];
		VI Child = Sy[1];

		// If x is not already the child of parents of y insert
		// that as a child of them
		for (int j = 0; j < (int) Parent.size(); j++) {
			bool flag = false;
			for (int k = 0; k < (int) Temp_G[Parent[j]].size(); k++)
				if (Temp_G[Parent[j]][k] == cx) {
					flag = true;
					break;
				}
			if (!flag)
				Temp_G[Parent[j]].push_back(cx);
		}

		// If x is not already the parent of children of y insert
		// that as a parent of them
		for (int j = 0; j < (int) Child.size(); j++) {
			bool flag = false;
			for (int k = 0; k < (int) Temp_G[cx].size(); k++)
				if (Temp_G[cx][k] == Child[j]) {
					flag = true;
					break;
				}
			if (!flag)
				Temp_G[cx].push_back(Child[j]);
		}

		//Check for cycle
		bool Check_Cycle = true;
		//Check whether x is the descendant of any of y's children
		//or if any Descendants of x is a parent of y
		if (Child.size() > 0)
			Check_Cycle = !Descs[Child[0]][cx] & !(Child[0] == cx);
		if (Parent.size() > 0)
			Check_Cycle = Check_Cycle & !(Descs[cx][Parent[0]]) & !(Parent[0]
					== cx);

		// If no cycle appears update the graph
		if (Check_Cycle) {

			// Remove redundant edges
			if (Sy[0].size() == 1)
				if ((Ismember(Descs[Sy[0][0]], cx)) & (!Ismember_VI(
						G[Sy[0][0]], cx)))
					Temp_G[Sy[0][0]] = Remove(Temp_G[Sy[0][0]], cx);
			if (Sy[1].size() == 1)
				if ((Ismember(Descs[cx], Sy[1][0])) && (!Ismember_VI(G[cx],
						Sy[1][0])))
					Temp_G[cx] = Remove(Temp_G[cx], Sy[1][0]);
			if ((Sy[1].size() == 1) & (Sy[0].size() == 1))
				Temp_G[Sy[0][0]] = Remove(Temp_G[Sy[0][0]], Sy[1][0]);

			//Update IsPresent
			IsPresent[y[0]][y[1]] = cx;

			//Update G
			G = Temp_G;
			int Gsz = G.size();

			//cx Ancestors
			SafeVector<bool> curr_Anc(ZZ);
			if (Parent.size() > 0)
				curr_Anc = Ancs[Parent[0]];
			curr_Anc = Union_VBVI(curr_Anc, Parent);

			//Update cx's Ancestors
			Ancs[cx] = Union(Ancs[cx], curr_Anc, Gsz);

			//cx Descendants
			SafeVector<bool> curr_Desc(ZZ);
			if (Child.size() > 0)
				curr_Desc = Descs[Child[0]];
			curr_Desc = Union_VBVI(curr_Desc, Child);

			//Update cx's Descendants
			Descs[cx] = Union(Descs[cx], curr_Desc, Gsz);

			//Update Ancestors and Descendants for other nodes
			VI DD, AA;
			for (int j = 0; j < Gsz; j++) {
				if (Ancs[cx][j])
					AA.push_back(j);
				if (Descs[cx][j])
					DD.push_back(j);
			}
			for (int j = 0; j < (int) DD.size(); j++) {
				Ancs[DD[j]][cx] = true;
				for (int k = 0; k < (int) AA.size(); k++) {
					Ancs[DD[j]][AA[k]] = true;
					Descs[AA[k]][DD[j]] = true;
				}
			}
			for (int k = 0; k < (int) AA.size(); k++)
				Descs[AA[k]][cx] = true;

		}
		return Check_Cycle;
	}

	/////////////////////////////////////////////////////////////
	// CheckAddColumnMrg()
	//
	// Procedure to Merge two columns in the graph provided that
	// it does not introduce any cycle
	/////////////////////////////////////////////////////////////////

	bool CheckAddColumnMrg(int cx, int cy) {

		//Develop the temporary Graph to check for cycles
		VVI Temp_G;

		// Find the set of Children
		VI Child_x = G[cx];
		VI Child_y = G[cy];
		VI Child;
		Child_x = Union_VI(Child_x, Child_y);
		for (int j = 0; j < (int) Child_x.size(); j++)
			Child.push_back(Update(Child_x[j], cx, cy));

		//Construct the updated temporary Graph
		for (int j = 0; j < (int) G.size(); j++) {
			bool flag = false;
			VI nodes;
			for (int k = 0; k < (int) G[j].size(); k++)
				if ((G[j][k] == cx) | (G[j][k] == cy)) {
					if (flag == false) {
						nodes.push_back(cx);
						flag = true;
					}
				} else if (G[j][k] < cy)
					nodes.push_back(G[j][k]);
				else if (G[j][k] > cy)
					nodes.push_back(G[j][k] - 1);
			if (j == cx)
				Temp_G.push_back(Child);
			else if (j != cy)
				Temp_G.push_back(nodes);

		}

		//Check for cycle
		bool Check_Cycle;
		//Check whether any descendant of x is y or visevers
		Check_Cycle = !(Descs[cx][cy]) & !(Descs[cy][cx]);

		// If no cycle appears update the graph
		if (Check_Cycle) {

			// Remove redundant edges
			VI ax = GiveNodes(Ancs[cx], G.size());
			VI dy = GiveNodes(Descs[cy], G.size());
			for (int j = 0; j < (int) ax.size(); j++)
				for (int k = 0; k < (int) dy.size(); k++)
					if (Ismember_VI(G[ax[j]], dy[k]))
						Temp_G[Update(ax[j], cx, cy)] = Remove(Temp_G[Update(
								ax[j], cx, cy)], Update(dy[k], cx, cy));
			for (int j = 0; j < (int) ax.size(); j++)
				if (Ismember_VI(G[ax[j]], cy) & !Ismember_VI(G[ax[j]], cx))
					Temp_G[Update(ax[j], cx, cy)] = Remove(Temp_G[Update(ax[j],
							cx, cy)], cx);
			VI ay = GiveNodes(Ancs[cy], G.size());
			VI dx = GiveNodes(Descs[cx], G.size());
			for (int j = 0; j < (int) ay.size(); j++)
				for (int k = 0; k < (int) dx.size(); k++)
					if (Ismember_VI(G[ay[j]], dx[k])) {
						Temp_G[Update(ay[j], cx, cy)] = Remove(Temp_G[Update(
								ay[j], cx, cy)], Update(dx[k], cx, cy));
					}
			for (int j = 0; j < (int) ay.size(); j++)
				if (Ismember_VI(G[ay[j]], cx) & !Ismember_VI(G[ay[j]], cy)) {
					Temp_G[Update(ay[j], cx, cy)] = Remove(Temp_G[Update(ay[j],
							cx, cy)], cx);
				}
			VI pax = GiveParent(G, cx);
			VI chx = G[cx];
			VI pay = GiveParent(G, cy);
			VI chy = G[cy];
			for (int j = 0; j < (int) pax.size(); j++)
				if ((Ismember_VI(ay, pax[j])) & !(Ismember_VI(G[pax[j]], cy)))
					Temp_G[Update(pax[j], cx, cy)] = Remove(Temp_G[Update(
							pax[j], cx, cy)], cx);
			for (int j = 0; j < (int) pay.size(); j++)
				if ((Ismember_VI(ax, pay[j])) & !(Ismember_VI(G[pay[j]], cx)))
					Temp_G[Update(pay[j], cx, cy)] = Remove(Temp_G[Update(
							pay[j], cx, cy)], cx);
			for (int j = 0; j < (int) chx.size(); j++)
				if ((Ismember_VI(dy, chx[j])) & !(Ismember_VI(G[cy], chx[j])))
					Temp_G[cx] = Remove(Temp_G[cx], Update(chx[j], cx, cy));
			for (int j = 0; j < (int) chy.size(); j++)
				if ((Ismember_VI(dx, chy[j])) & !(Ismember_VI(G[cx], chy[j])))
					Temp_G[cx] = Remove(Temp_G[cx], Update(chy[j], cx, cy));

			//Update IsPresent
			for (int j = 0; j < (int) IsPresent.size(); j++)
				for (int k = 0; k < (int) IsPresent[j].size(); k++)
					IsPresent[j][k] = Update(IsPresent[j][k], cx, cy);

			//Update G
			G = Temp_G;
			int Gsz = G.size();

			//Merged node Ancestors and Descendants
			SafeVector<bool> curr_Anc(ZZ), curr_Desc(ZZ);
			curr_Anc = Union(Ancs[cx], Ancs[cy], Gsz + 1);
			curr_Desc = Union(Descs[cx], Descs[cy], Gsz + 1);

			//Update Ancestors and Descendants for other nodes
			SafeVector<SafeVector<bool> > new_Ancs, new_Descs;
			for (int j = 0; j < (int) Ancs.size(); j++) {
				if (j == cx) {
					curr_Anc = Update(curr_Anc, cy, Gsz + 1);
					new_Ancs.push_back(curr_Anc);
				} else if (j != cy) {
					SafeVector<bool> n_anc = Update(Ancs[j], cy, Gsz + 1);
					new_Ancs.push_back(n_anc);
				}
			}
			Ancs = new_Ancs;
			for (int j = 0; j < (int) Descs.size(); j++) {
				if (j == cx) {
					curr_Desc = Update(curr_Desc, cy, Gsz + 1);
					new_Descs.push_back(curr_Desc);
				} else if (j != cy) {
					SafeVector<bool> n_Desc = Update(Descs[j], cy, Gsz + 1);
					new_Descs.push_back(n_Desc);
				}
			}
			Descs = new_Descs;
			VI DD, AA;
			for (int j = 0; j < Gsz; j++) {
				if (curr_Anc[j])
					AA.push_back(j);
				if (curr_Desc[j])
					DD.push_back(j);
			}
			for (int j = 0; j < (int) DD.size(); j++) {
				Ancs[DD[j]][cx] = true;
				for (int k = 0; k < (int) AA.size(); k++) {
					Ancs[DD[j]][AA[k]] = true;
					Descs[AA[k]][DD[j]] = true;
				}
			}
			for (int k = 0; k < (int) AA.size(); k++)
				Descs[AA[k]][cx] = true;

		}
		return Check_Cycle;
	}

	/////////////////////////////////////////////////////////////
	// AddtoPath()
	//
	// Add node N2 after node N1 in the path
	/////////////////////////////////////////////////////////////////

	void AddtoPath(VI &Path, int N1, int N2) {

		int Ps = Path.size();
		int h;
		if (N1 == -1) {
			// If N2 is a root node
			h = -1;
		} else {
			for (h = 0; h < Ps; h++)
				if (Path[h] == N1)
					break;
		}
		if (h == (Ps - 1))
			Path.push_back(N2);
		else {
			Path.push_back(Path[Ps - 1]);
			for (int k = Ps - 1; k > (h + 1); k--)
				Path[k] = Path[k - 1];
			Path[h + 1] = N2;
		}
	}

	/////////////////////////////////////////////////////////////
	// FindPath()
	//
	// Update the path starting from the node N1 and considering the
	// passed nodes in "marked"
	/////////////////////////////////////////////////////////////////

	void FindPath(int N1, SafeVector<bool> &marked, VI &Path) {

		for (int j = 0; j < (int) G[N1].size(); j++) {
			if (!marked[G[N1][j]]) {
				marked[G[N1][j]] = true;
				AddtoPath(Path, N1, G[N1][j]);
				FindPath(G[N1][j], marked, Path);
			}
		}
	}

	/////////////////////////////////////////////////////////////
	// Path2Align()
	//
	// Convert path to a valid alignment
	/////////////////////////////////////////////////////////////////
	void Path2Align(VI Path, VVVI SRC, VVI ZeroPos, int SRCsize) {

		//Construct the alignment structure
		alignment = new MultiSequence();
		VI lengths(sequences->GetNumSequences());
		for (int i = 0; i < sequences->GetNumSequences(); i++) {

			Sequence *seq = sequences->GetSequence(i);
			Sequence *ret = new Sequence();
			assert (ret);

			ret->isValid = seq->isValid;
			ret->header = seq->header;
			ret->data = new SafeVector<char>;
			assert (ret->data);
			ret->length = (int) Path.size() + SRCsize;
			ret->sequenceLabel = seq->sequenceLabel;
			ret->inputLabel = seq->inputLabel;
			ret->data->push_back('@');
			alignment->AddSequence(ret);
		}

		VI temp_seqs(sequences->GetNumSequences());
		for (int i = 0; i < sequences->GetNumSequences(); i++)
			temp_seqs[i] = i;

		//insert single residue columns which should appear
		//at position zero
		for (int i = 0; i < (int) ZeroPos.size(); i++) {
			int seqnum = ZeroPos[i][0];
			int resnum = ZeroPos[i][1];
			alignment->GetSequence(seqnum)->data->push_back(sequences->GetSequence(seqnum)->GetPosition(resnum + 1));
			for (int k = 0; k < sequences->GetNumSequences(); k++)
				if (k != seqnum)
					alignment->GetSequence(k)->data->push_back('-');
		}

		//insert other columns
		for (int i = 0; i < (int) Path.size(); i++) {
			VI b(temp_seqs);
			//insert residues
			for (int j = 0; j < (int) cols[Path[i]].size(); j++) {
				int seqnum = cols[Path[i]][j][0];
				int resnum = cols[Path[i]][j][1];
				alignment->GetSequence(seqnum)->data->push_back(sequences->GetSequence(seqnum)->GetPosition(resnum + 1));
				b = Remove(b, seqnum);
			}

			//insert gaps
			for (int j = 0; j < (int) b.size(); j++) {
				alignment->GetSequence(b[j])->data->push_back('-');
			}

			//insert single residue columns which should
			// appear after i'th rgular column
			for (int j = 0; j < (int) SRC[i].size(); j++) {
				int seqnum = SRC[i][j][0];
				int resnum = SRC[i][j][1];
				alignment->GetSequence(seqnum)->data->push_back(sequences->GetSequence(seqnum)->GetPosition(resnum + 1));
				for (int k = 0; k < sequences->GetNumSequences(); k++)
					if (k != seqnum)
						alignment->GetSequence(k)->data->push_back('-');
			}
		}

	}

	/////////////////////////////////////////////////////////////
	// GetPrior()
	//
	// Find out whether base i in sequence x is more likely to 
	// base piar left or right
	/////////////////////////////////////////////////////////////////
	int GetPrior(int x, int i) {

		int max_x = 0;
		if ((BPProb[x][i][0]) < (BPProb[x][i][1]))
			max_x = 1;
		if ((1 - (BPProb[x][i][0]) - (BPProb[x][i][1])) > (BPProb[x][i][max_x]))
			max_x = 2;
		return max_x;
	}

	/////////////////////////////////////////////////////////////
	// FindBestAlign()
	//
	// Find best alignment pair (y,y') for base-pair (x[1],x[2])  
	// sequence x[0]
	/////////////////////////////////////////////////////////////////
	void FindBestAlign(VI x, VI &BestAlign, VVI BPairs, VF &AProb) {

		for (int i = 0; i < BPairs.size(); i++) {
			VI xx = x;
			if (BPairs[i][0] != x[0]) {
				VI y = BPairs[i];
				if (y[0] < xx[0]) {
					swap(xx, y);
				}
				int x0 = xx[0];
				int x1 = xx[1];
				int x2 = xx[2];
				if (x1 > x2)
					swap(x1, x2);
				int y0 = y[0];
				int y1 = y[1];
				int y2 = y[2];
				if (y1 > y2)
					swap(y1, y2);

				int yy0 = BPairs[i][0];
				int yy1 = BPairs[i][1];
				int yy2 = BPairs[i][2];
				if (yy1 > yy2)
					swap(yy1, yy2);
				AProb[i] = (sparseMatrices[x0][y0]->GetValue(x1 + 1, y1 + 1))
						* (sparseMatrices[x0][y0]->GetValue(x2 + 1, y2 + 1))
						* BPPMatrices[yy0]->GetProb(yy1 + 1, yy2 + 1);
			}
		}
		int n, low, high;
		float *a;
		int *ind;
		n = AProb.size();
		a = new float[n];
		ind = new int[n];
		for (int i = 0; i < n; i++) {
			a[i] = AProb[i];
			ind[i] = i;
		}
		high = n - 1;
		low = 0;
		Quick_sort(low, high, a, ind);
		for (int i = 0; i < AProb.size(); i++) {
			BestAlign.push_back(ind[AProb.size() - i - 1]);
		}
	}

	/////////////////////////////////////////////////////////////////
	// FindMaxBP()
	//
	//Find the maximum base pairing probability
	/////////////////////////////////////////////////////////////////
	float FindMaxBP(SafeVector<BPPMatrix*> BPPMatrices) {

		float m = 0;

		for (int i = 0; i < BPPMatrices.size(); i++) {
			vector<PIF2> &probs2 = BPPMatrices[i]->bppMat.data2;
			for (int j = 0; j < probs2.size(); j++) {
				float tmpProb2 = probs2[j].prob;
				if (tmpProb2 > m)
					m = tmpProb2;

			}

		}
		return m;

	}
	
	/////////////////////////////////////////////////////////////////
	// FindMaxPP()
	//
	//Find the maximum posterior pairwise base alignment probability
	/////////////////////////////////////////////////////////////////
	float FindMaxPP(SafeVector<SafeVector<SparseMatrix *> > sparseMatrices) {

		float m = 0;
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
public:

	/////////////////////////////////////////////////////////////////
	// AlignGraph::AlignGraph()
	//
	// Default constructor.
	/////////////////////////////////////////////////////////////////

	AlignGraph() {
	}

	/////////////////////////////////////////////////////////////////
	// AlignGraph::ConstGraph()
	//
	// Construct the alignment graph
	// Step1: Structural skeleton construction (usign base pairing 
	// probabilities)
	/////////////////////////////////////////////////////////////////
	AlignGraph(MultiSequence *sequences,
			SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
			SafeVector<BPPMatrix*> &BPPMatrices, VVVF BPProb, int maxlength, float Tb) :
		sequences(sequences), sparseMatrices(sparseMatrices),
				BPPMatrices(BPPMatrices), BPProb(BPProb), StemSize(0) {

		
		int numSeqs = sparseMatrices.size();

		// initialize the zero vector:
		// initial guess of number of columns in the alignment is 1.5*max_length
		for (int i = 0; i < 1.5 * maxlength; i++)
			ZZ.push_back(false);

		//initialize IsPresent (Matrix of residues in each node)
		for (int i = 0; i < numSeqs; i++) {
			VI s;
			for (int j = 0; j < maxlength; j++) {
				s.push_back(-1);
			}
			IsPresent.push_back(s);
		}

		//initialize IsStem (Matrix of Stem components)
		for (int i = 0; i < numSeqs; i++) {
			VI s;
			for (int j = 0; j < maxlength; j++) {
				s.push_back(-1);
			}
			IsStem.push_back(s);
		}
		//initialize LRStem (Matrix of Left-Right Stem components)
		for (int i = 0; i < numSeqs; i++) {
			SafeVector<bool> s;
			for (int j = 0; j < maxlength; j++) {
				s.push_back(false);
			}
			LRStem.push_back(s);
		}

		int cnt_NN = 0; //number of new nodes addition calls
		int cnt_CE = 0; //number of column extension calls
		int cnt_CM = 0; //number of column merging calls
		int cnt_tot = 0; // totla number of calls
		int cnt_NoCycle = 0;// number of calls with no cycle introduction

		float maxbp = FindMaxBP(BPPMatrices);
		float maxpp = FindMaxPP(sparseMatrices);

		//Start Graph construction from the largest base-pairing probabilities
		VVI BPairs;
		VF BPairp;
		for (int i = 0; i < numSeqs; i++) {
			SparseMatrix bppMat = BPPMatrices[i]->bppMat;
			int numRows = bppMat.GetSeq1Length();

			for (int k = 1; k <= numRows; k++) {
				for (int h = 0; h < bppMat.GetRowSize(k); h++) {

					int ii = k - 1;
					int jj = (bppMat.GetRowPtr(k)[h].first) - 1;
					float p = bppMat.GetRowPtr(k)[h].second;
					if (p >= (maxbp*Tb)) {
						VI newentry;
						newentry.push_back(i);
						newentry.push_back(ii);
						newentry.push_back(jj);
						BPairs.push_back(newentry);
						BPairp.push_back(p);
					}
				}
			}
		}

		//Sort the largest base-pairing probabilities
		float *a2;
		int *ind2;
		int n, low, high;
		n = BPairp.size();
		a2 = new float[n];
		ind2 = new int[n];
		for (int i = 0; i < n; i++) {
			a2[i] = BPairp[i];
			ind2[i] = i;
		}
		high = n - 1;
		low = 0;
		Quick_sort(low, high, a2, ind2);
		int CNNT = 0;

		for (int i = 0; i < (int) BPairs.size(); i++) {

			VI x = BPairs[ind2[BPairs.size() - i - 1]];
			VI x1(2), x2(2);
			x1[0] = x[0];
			x1[1] = min(x[1], x[2]);
			x2[0] = x[0];
			x2[1] = max(x[1], x[2]);

			VI BestAlign;
			VF AProb(BPairs.size(), 0);

			FindBestAlign(x, BestAlign, BPairs, AProb);

			for (int j = 0; j < BestAlign.size(); j++) {

				VI y = BPairs[BestAlign[j]];

				float p = AProb[BestAlign[j]];

				if ((CNNT>(maxlength/10*numSeqs))|(p < (maxbp / 2 * maxpp / 5 * maxpp
						/ 5)))
					break;
				else {
					VI y1(2), y2(2);
					y1[0] = y[0];
					y1[1] = min(y[1], y[2]);
					y2[0] = y[0];
					y2[1] = max(y[1], y[2]);

					//Check whether x1,x2,y1,y2 have appeared in any node
					int cx1 = -1;
					int cx2 = -1;
					int cy1 = -1;
					int cy2 = -1;
					int sx1 = -1;
					int sx2 = -1;
					int sy1 = -1;
					int sy2 = -1;
					bool fx1 = false;
					bool fx2 = false;
					bool fy1 = false;
					bool fy2 = false;
					if (IsPresent[x1[0]][x1[1]] != -1) {
						fx1 = true;
						cx1 = IsPresent[x1[0]][x1[1]];
						sx1 = IsStem[x1[0]][x1[1]];
					}
					if (IsPresent[x2[0]][x2[1]] != -1) {
						fx2 = true;
						cx2 = IsPresent[x2[0]][x2[1]];
						sx2 = IsStem[x2[0]][x2[1]];
					}
					if (IsPresent[y1[0]][y1[1]] != -1) {
						fy1 = true;
						cy1 = IsPresent[y1[0]][y1[1]];
						sy1 = IsStem[y1[0]][y1[1]];
					}
					if (IsPresent[y2[0]][y2[1]] != -1) {
						fy2 = true;
						cy2 = IsPresent[y2[0]][y2[1]];
						sy2 = IsStem[y2[0]][y2[1]];
					}
					if ((fx1 != fx2) | (fy1 != fy2) | (sy1 != sy2) | (sx1
							!= sx2))
						continue;

					bool Check_Cycle = false;

					bool fx = fx1;
					bool fy = fy1;
					int sx = sx1;
					int sy = sy1;
					//if niether x1,x2, nor y1,y2 have appeared:  Add a new node
					if ((!fx) & (!fy)) {
						VVI Temp_G = G;
						VVI Temp_IsPresent = IsPresent;
						SafeVector<SafeVector<bool> > Temp_Ancs = Ancs;
						SafeVector<SafeVector<bool> > Temp_Descs = Descs;
						Check_Cycle = CheckAddNewNode(x1, y1);
						if (Check_Cycle) {
							Check_Cycle = CheckAddNewNode(x2, y2);
						}
						if (Check_Cycle) {
							IsStem[x1[0]][x1[1]] = StemSize;
							IsStem[x2[0]][x2[1]] = StemSize;
							IsStem[y1[0]][y1[1]] = StemSize;
							IsStem[y2[0]][y2[1]] = StemSize;
							LRStem[x1[0]][x1[1]] = false;
							LRStem[x2[0]][x2[1]] = true;
							LRStem[y1[0]][y1[1]] = false;
							LRStem[y2[0]][y2[1]] = true;
							StemSize++;
							CNNT++;
						}


						if (!Check_Cycle) {
							G = Temp_G;
							IsPresent = Temp_IsPresent;
							Ancs = Temp_Ancs;
							Descs = Temp_Descs;
						} else {
							break;
						}
					}
					//if only one of x1,x2 or y1,y2 have appeared:  Extend the column
					else if (fx ^ fy) {

						//make x the residue which is already in the graph
						if (fy) {
							swap(x1, y1);
							swap(cx1, cy1);
							swap(sx1, sy1);
							swap(fx1, fy1);
							swap(x2, y2);
							swap(cx2, cy2);
							swap(sx2, sy2);
							swap(fx2, fy2);
							swap(sx, sy);
							swap(fx, fy);
						}

						//immediate check for cycle : check whether cx
						// has common residue in sequence of residue y
						bool iflag = false;
						for (int k = 0; k < (int) IsPresent[y1[0]].size(); k++)
							if ((IsPresent[y1[0]][k] == cx1)
									| (IsPresent[y1[0]][k] == cx2)) {
								iflag = true;
								break;
							}
						if (iflag == false) {
							VVI Temp_G = G;
							VVI Temp_IsPresent = IsPresent;
							SafeVector<SafeVector<bool> > Temp_Ancs = Ancs;
							SafeVector<SafeVector<bool> > Temp_Descs = Descs;
							Check_Cycle = CheckAddColumnEx(y1, cx1);
							if (Check_Cycle) {
								Check_Cycle = CheckAddColumnEx(y2, cx2);
							}
							if (Check_Cycle) {
								IsStem[y1[0]][y1[1]] = IsStem[x1[0]][x1[1]];
								IsStem[y2[0]][y2[1]] = IsStem[x2[0]][x2[1]];
								LRStem[y1[0]][y1[1]] = false;
								LRStem[y2[0]][y2[1]] = true;
								CNNT++;
							}
							if (!Check_Cycle) {
								G = Temp_G;
								IsPresent = Temp_IsPresent;
								Ancs = Temp_Ancs;
								Descs = Temp_Descs;
							} else {
								break;
							}
						} else
							continue;
					}
					//if both of x and y have appeared: Merge the columns
					else if ((fx & fy) & (sx != sy)) { //Column Merge

						//immediate check for cycle : check whether cx
						// has common residue in sequence of residue y or vice versa
						bool iflag = false;
						for (int k = 0; k < (int) IsPresent[y1[0]].size(); k++)
							if ((IsPresent[y1[0]][k] == cx1)
									| (IsPresent[y1[0]][k] == cx2)) {
								iflag = true;
								break;
							}
						if (iflag == false) {
							for (int k = 0; k < (int) IsPresent[x1[0]].size(); k++)
								if ((IsPresent[x1[0]][k] == cy1)
										| (IsPresent[x1[0]][k] == cy2)) {
									iflag = true;
									break;
								}
						}
						if (iflag == false) {
							VVI Temp_G = G;
							VVI Temp_IsPresent = IsPresent;
							SafeVector<SafeVector<bool> > Temp_Ancs = Ancs;
							SafeVector<SafeVector<bool> > Temp_Descs = Descs;
							Check_Cycle = CheckAddColumnMrg(min(cx1, cy1), max(
									cx1, cy1));
							if (Check_Cycle) {
								cx2 = Update(cx2, min(cx1, cy1), max(cx1, cy1));
								cy2 = Update(cy2, min(cx1, cy1), max(cx1, cy1));
								Check_Cycle = CheckAddColumnMrg(min(cx2, cy2),
										max(cx2, cy2));
							}
							if (Check_Cycle) {
								if (sx > sy)
									swap(sx, sy);
								//Update IsStem
								for (int k = 0; k < (int) IsStem.size(); k++)
									for (int h = 0; h < (int) IsStem[k].size(); h++)
										IsStem[k][h] = Update(IsStem[k][h], sx,
												sy);
							}

							if (!Check_Cycle) {
								G = Temp_G;
								IsPresent = Temp_IsPresent;
								Ancs = Temp_Ancs;
								Descs = Temp_Descs;
							} else {
								break;
							}
						}
					}

				}
			}
			if (G.size() > ZZ.size() - 10) {
				for (int j = 0; j < 100; j++)
					ZZ.push_back(false);
				for (int k = 0; k < (int) G.size(); k++)
					for (int j = 0; j < 100; j++) {
						Descs[k].push_back(false);
						Ancs[k].push_back(false);
					}
			}
		}
	}

	void AlignGraph_PicXAA(VVI aligns, VF alignp, MultiSequence *sequences,
			SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
			int maxlength) {

		int numSeqs = sequences->GetNumSequences();

		//sort the posterior probabilities using quick sort algorithm
		int n, low, high;
		float *a;
		int *ind;
		n = alignp.size();
		a = new float[n];
		ind = new int[n];
		for (int i = 0; i < n; i++) {
			a[i] = alignp[i];
			ind[i] = i;
		}
		high = n - 1;
		low = 0;
		Quick_sort(low, high, a, ind);

		// initialize counters
		int cnt_NN = 0; //number of new nodes addition calls
		int cnt_CE = 0; //number of column extension calls
		int cnt_CM = 0; //number of column merging calls
		int cnt_tot = 0; // totla number of calls
		int cnt_NoCycle = 0;// number of calls with no cycle introduction

		//Start Graph construction from the largest posteriro probabilities
		for (int i = 0; i < (int) aligns.size(); i++) {
			//residue pair to be aligned
			VI res_pair(aligns[ind[aligns.size() - i - 1]]);

			//residues x and y
			// x[0] sequence number of residue x
			// x[1] residue position of residue x
			VI x(2);
			x[0] = res_pair[0];
			x[1] = res_pair[1];
			VI y(2);
			y[0] = res_pair[2];
			y[1] = res_pair[3];

			//Check whether x and/or y have appeared in any node
			int cx = -1;
			int cy = -1;
			bool fx = false;
			bool fy = false;
			if (IsPresent[x[0]][x[1]] != -1) {
				fx = true;
				cx = IsPresent[x[0]][x[1]];
			}
			if (IsPresent[y[0]][y[1]] != -1) {
				fy = true;
				cy = IsPresent[y[0]][y[1]];
			}

			bool Check_Cycle = false;

			//if niether x nor y have appeared:  Add a new node
			if ((!fx) & (!fy)) {
				cnt_tot++;
				cnt_NN++;
				Check_Cycle = CheckAddNewNode(x, y);

				if (Check_Cycle)
					cnt_NoCycle++;
			}

			//if only one of x or y have appeared:  Extend the column
			else if (fx ^ fy) {
				//make x the residue which is already in the graph
				if (fy) {
					swap(x, y);
					swap(cx, cy);
					swap(fx, fy);
				}

				//immediate check for cycle : check whether cx
				// has common residue in sequence of residue y
				bool iflag = false;
				for (int j = 0; j < (int) IsPresent[y[0]].size(); j++)
					if (IsPresent[y[0]][j] == cx) {
						iflag = true;
						break;
					}
				if (iflag == false) {
					cnt_tot++;
					cnt_CE++;
					Check_Cycle = CheckAddColumnEx(y, cx);
					if (Check_Cycle)
						cnt_NoCycle++;
				}
			}

			//if both of x and y have appeared: Merge the columns
			else if ((fx & fy) & (cx != cy)) { //Column Merge
				//immediate check for cycle : check whether cx
				// has common residue in sequence of residue y or vice versa
				bool iflag = false;
				for (int j = 0; j < (int) IsPresent[y[0]].size(); j++)
					if (IsPresent[y[0]][j] == cx) {
						iflag = true;
						break;
					}
				for (int j = 0; j < (int) IsPresent[x[0]].size(); j++)
					if (IsPresent[x[0]][j] == cy) {
						iflag = true;
						break;
					}
				if (iflag == false) {
					if (cx > cy) {
						swap(x, y);
						swap(cx, cy);
						swap(fx, fy);
					}
					cnt_tot++;
					cnt_CM++;
					Check_Cycle = CheckAddColumnMrg(cx, cy);
					if (Check_Cycle)
						cnt_NoCycle++;
				}
			}

			// If the initial guess of the number of nodes
			// is close to reach extend number of possible nodes
			if (!((fx & fy) & (cx == cy))) {
				if (G.size() > ZZ.size() - 10) {
					for (int j = 0; j < 100; j++)
						ZZ.push_back(false);
					for (int k = 0; k < (int) G.size(); k++)
						for (int j = 0; j < 100; j++) {
							Descs[k].push_back(false);
							Ancs[k].push_back(false);
						}
				}
			}
		}

		//Construct the columns
		// cols[i]: vector of residues in node i

		VVI temp;
		for (int i = 0; i < (int) G.size(); i++)
			cols.push_back(temp);

		for (int i = 0; i < sequences->GetNumSequences(); i++)
			for (int j = 0; j < sequences->GetSequence(i)->GetLength(); j++)
				if (IsPresent[i][j] != -1) {
					VI x;
					x.push_back(i);
					x.push_back(j);
					cols[IsPresent[i][j]].push_back(x);
				}
	}

	/////////////////////////////////////////////////////////////////
	// AlignGraph::Graph2Align()
	//
	// Convert the graph to a valid alignment
	/////////////////////////////////////////////////////////////////
	void Graph2Align() {

		//vector of roots in G
		VI roots;
		for (int i = 0; i < (int) G.size(); i++)
			if (GiveParent(G, i).size() == 0)
				roots.push_back(i);

		//vector of passed nodes
		SafeVector<bool> marked(G.size());

		//Alignment Path deduced from G
		VI Path;
		for (int i = 0; i < (int) roots.size(); i++) {
			AddtoPath(Path, -1, roots[i]);
			FindPath(roots[i], marked, Path);
		}
		//Mapping the place of each node in Path
		VI PathMap(Path.size());
		for (int i = 0; i < (int) Path.size(); i++) {
			PathMap[Path[i]] = i;
		}
		// Single residue columns which should appear after
		// each node in Path
		VVVI SRC(Path.size());

		// single residue columns which should appear
		//at position zero
		VVI ZeroPos;

		//figure out the single residue columns
		int SRCsize = 0;
		for (int i = 0; i < sequences->GetNumSequences(); i++)
			for (int j = 0; j < sequences->GetSequence(i)->GetLength(); j++)
				if (IsPresent[i][j] == -1) {
					SRCsize++;
					VI x(2);
					x[0] = i;
					x[1] = j;
					int ct = j - 1;
					while (ct >= 0) {
						if (IsPresent[i][ct] != -1) {
							SRC[PathMap[IsPresent[i][ct]]].push_back(x);
							break;
						}
						ct--;
					}
					if (ct == -1)
						ZeroPos.push_back(x);
				}


		//Convert the path to valid alignment
		Path2Align(Path, SRC, ZeroPos, SRCsize);

	}

	/////////////////////////////////////////////////////////////////
	// AlignGraph::GetAlignment()
	//
	// Retrieve the alignment from the AlignGraph object.
	/////////////////////////////////////////////////////////////////

	MultiSequence* GetAlignment() {
		return alignment;
	}

	SafeVector<char> * GetProfileAlignment() {
		return Profile_alignment;
	}

	/////////////////////////////////////////////////////////////////
	// AlignGraph::~AlignGraph()
	//
	// Destructor.  Gets rid of objects contained in the graph.
	/////////////////////////////////////////////////////////////////

	~ AlignGraph() {

		if (sequences) {
			delete sequences;
			sequences = NULL;
		}
		if (alignment) {
			delete alignment;
			alignment = NULL;
		}
	}
};
#endif
