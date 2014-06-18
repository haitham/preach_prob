#ifndef GRAPH_H_
#define GRAPH_H_

#include <assert.h>

#include <vector>
#include <bitset>
#include <iostream>

#include "Stl.h"

#define MAX_NODES 256
#define MAX_EDGES 512

// Uncomment the next macro if detailed statistics of execution are needed
//#define STATISTICS


// macro to smoth up the use of bitsets
#define FOREACH_BS(v, vSet)	  \
	for (size_t v=vSet._Find_first(); v!=vSet.size(); v=vSet._Find_next(v))

using namespace std;
using namespace lemon;

typedef bitset<MAX_NODES> Nodes_T;
typedef bitset<MAX_EDGES> Edges_T;
typedef vector<Edges_T> EdgesSet;
typedef vector< pair<int,double> > EdgeList; // datastructure to get a set of edges and their weights

/* Encoding of a path */
class Path{
	Nodes_T nodes; // encoding of nodes touched by the path
	Edges_T edges; // encoding of edges
	int node; // the top node in the path

	double weight; // total weight of the path

public:

	// constructor to start a path from beginning
	Path(int source):node(source){
		nodes.set(source);
		weight=1.0;
	}

	// constructor to buid a path from existing path and step
	Path(Path& old, int _node, int _edge, double _weight):
		nodes(old.nodes), edges(old.edges), node(_node){
		nodes.set(_node), edges.set(_edge);
		weight=old.weight*_weight;
	}

	int GetNode(void){ return node; }
	bool ContainsNode(int _node){ return nodes[_node]; }
	bool ContainsEdge(int _edge){ return edges[_edge]; }
	void UnionNodes(Nodes_T& _nodes){ _nodes|=nodes; }
	Nodes_T& GetNodes(void){ return nodes; }

	Edges_T& GetEdges(void){ return edges; }

	void Print(void){ cout << "Node: " << node << "\tnodes=" << nodes << "\tedges=" << edges << endl; }
};

typedef vector<Path> PathSet;

/*  encoding of a cut */
class Cut {
	Nodes_T nodes; // set of nodes that specify the minimal cut (contains S implicitly)
	Edges_T edges; // set of edges forming the cut

	Nodes_T banned; // list of nodes that should not be included in future expansion
	Nodes_T choices; // set of nodes connected to current cut that are future choices
	// if chices = empty, this is a dead end

public:
	// consturctor for the special cut {S}
	Cut(int source, Nodes_T& _choices, Edges_T& _edges):
		edges(_edges), choices(_choices){
		nodes.set(source);
	}

	// constructor to build derived cut from existing cut
	Cut(Cut& prev, int _node, Edges_T& _edges, Nodes_T& _choices):
		edges(_edges), choices(_choices)
	{
		nodes = prev.nodes;
		nodes.set(_node);
		choices.reset(_node); // no circular choices
	}

	Cut(Edges_T& _edges){
		edges = _edges;
	}


	bool HasChoice(int _node){ return choices[_node]; }
	bool InsideCut(int _node){ return nodes[_node]; }
	Nodes_T& GetChoices(void){ return choices; }
	Edges_T& GetEdges(void){ return edges; }
	Nodes_T& GetBanned(void){ return banned; }
	Nodes_T& GetNodes(void){ return nodes; }
};

typedef vector<Cut> CutSet;

/* Helper class for Polynomial.

   Each Term encodes a term in the polynomial. We need to be able to
   create more terms from a term efficiently.

   x^* and y^* are special terms dealt with in the Polynomial.
*/
class Term {
	double coef; // coeficient of this term
	Edges_T x; // the x component of the term.
	// if x[i]=true, the term contains x_i
	Edges_T y; // same

public:

	// Constructor for a specific term
	// i: index
	// p: probability
	// The term added is one of p*x_i, (1-p)*y_i controlled by isX
	Term(){
		coef = 1.0;
	}

	// same signature as constructor
	void Multiply(int i, double p, bool isX){
		coef *= isX ? p : (1-p);
		if (isX) x.set(i); else y.set(i);
	}

	/* return -1.0 if no colapse, the coef otherwise
	   collapse if any path is contained in x part
	 */

	double XStarCollapse(PathSet& paths){
		FOREACH_STL(path, paths){
			Edges_T& pE = path.GetEdges();
			if ((~x & pE).none()) // path in x
				return coef;
		}END_FOREACH;
		return -1.0;  // not found
	}
	// version that uses Edges_T sets
	double XStarCollapse(EdgesSet& paths){
		FOREACH_STL(path, paths){
			if ((~x & path).none()) // path in x
				return coef;
		}END_FOREACH;
		return -1.0;  // not found
	}


	/* symetric of above */
	double YStarCollapse(CutSet& cuts){
		FOREACH_STL(cut, cuts){
			Edges_T& cE = cut.GetEdges();
			if ((~y & cE).none()) // cut in y
				return coef;
		}END_FOREACH;
		return -1.0; // not found
	}

	double YStarCollapse(EdgesSet& cuts){
		FOREACH_STL(cut, cuts){
			if ((~y & cut).none()) // cut in y
				return coef;
		}END_FOREACH;
		return -1.0; // not found
	}

};

/* Class to implement polynomial multiplicatio to carry out the
   probabilistic computation

   The Polynomial that needs to be produced is:
   \prod_{e\E} (p_e x_e+q_e y_e) with p_e is prob of edge e, q_e=1-p_e
*/

class Polynomial{
	typedef vector<Term> Terms;

	// paths and cuts. These do not change during computation
	PathSet& paths;
	CutSet& cuts;

	Edges_T covered; // set of covered edges so far
	bool collapse; // should we collapse?

	long double xStarCoef; // coefficient of xStar
	long double yStarCoef; // same for y
	Terms terms;
	int edgeCount;
 public:
	// main constructor. Requires paths and cuts
	// if collapse=true, we collapse as we go.
	Polynomial(PathSet& _paths, CutSet& _cuts, bool _collapse):
		paths(_paths), cuts(_cuts), collapse(_collapse),
		xStarCoef(0.0), yStarCoef(0.0){
		// Adding the empty term with 1.0 coefficient
		terms.push_back(Term());
		edgeCount = 0;
	}

	// Collapse using given sets of paths and edges
	// only these edges and cuts are needed to collapse
	void Collapse(EdgesSet& pathsC, EdgesSet& cutsC){
		Terms nTerms; // new terms
		FOREACH_STL(t, terms){
			// try an X* collapse
			double cCoef = t.XStarCollapse(pathsC);
			if (cCoef != -1.0){ // collapses
					xStarCoef+=cCoef;
#if 1 // sanity check. Remove latter
					assert(t.YStarCollapse(cutsC) == -1.0);
#endif
			} else { // try a Y* collapse
#if 1
				cCoef = t.YStarCollapse(cutsC);
#else
				cCoef == -1.0;
#endif
		if (cCoef != -1.0){ // collapses
					yStarCoef+=cCoef;
				} else { // both failed, keep term
					nTerms.push_back(t);
				}
			}
		}END_FOREACH;
		nTerms.swap(terms); // swap the new and old terms
		// nTerms gets distroyed with content of old terms
	}

	// version of AddEdge with no  collapsation
	void AddEdgeNoCollapse(int i, double p){
	Terms nTerms; // new terms
	edgeCount ++;
#ifdef STATISTICS
		cout << "Adding edge " << edgeCount << endl;
#endif
		FOREACH_STL(t, terms){
			Term xT = t;
			xT.Multiply(i, p, true);
			nTerms.push_back(xT);

			Term yT = t;
			yT.Multiply(i, p, false);
			nTerms.push_back(yT);

//			if (collapse){
//				xStarCoef = xStarCoef*p + xStarCoef*(1-p);
//				yStarCoef = yStarCoef*p + yStarCoef*(1-p);
//			}
		}END_FOREACH;
		nTerms.swap(terms); // swap the new and old terms
#ifdef STATISTICS
		cout << "Polynomial has " << terms.size() << " terms " << endl;
#endif
		// nTerms gets distroyed with content of old terms
	}


	/** Method to add an edge set.
	    Strategy:

	    1. Add each edge in the set without collapsation
	    2. Compute the set of paths and cuts that can be collapsed by this extra set.
	    3. Collapse using the sets above

	*/

	void AddSetEdges(EdgeList& edges, PathSet& paths, CutSet& cuts){
		//Step 1. Add each edge in the set without collapsation
		Edges_T newEdges;

		FOREACH_STL(edge, edges){
			if (covered[edge.first]){
				cout << "ERROR: Edges can provided only once" << endl;
				return;
			}
			newEdges.set(edge.first);
			AddEdgeNoCollapse(edge.first, edge.second);
		}END_FOREACH;

		if (collapse){
			Edges_T newCovered = covered|newEdges;

			//Step 2. Compute the set of paths and cuts that can be collapsed by this extra set.
			EdgesSet pathsC;
			FOREACH_STL(path, paths){
				Edges_T& ed = path.GetEdges();
				// a path is good only if included in newCovered but has at least one element from newEdges
				if ( (~newCovered & ed).none()/*in newCovered*/ &&
					 (newEdges & ed).any()/*one element from newEdges*/)
					pathsC.push_back(ed);
			}END_FOREACH;

			EdgesSet cutsC;
			FOREACH_STL(cut, cuts){
				Edges_T& ed = cut.GetEdges();
				// a path is good only if included in newCovered but has at least one element from newEdges
				if ( (~newCovered & ed).none()/*in newCovered*/ &&
					 (newEdges & ed).any()/*one element from newEdges*/)
					cutsC.push_back(ed);
			}END_FOREACH;

			//Step 3. Collapse using the sets above
			Collapse(pathsC, cutsC);

			covered=newCovered; // we not cover the new edges as well
		}
	}

	/* add an edge to the polynomial
	   i: the edge
	   p: the probability

	   Equations that have to be satisfied:
	   1. c x* (p x_i+q y_i) = cpx*+cqx*= c x* since p+q=1
	   2. c y* (p x_i+q y_i) = c y* since p+q=1
	   1-2 => new terms do not change x* and y* directly

	   3. If a term collapses on X, its coefficient is added to x*
	   c x* + term, if term => d x*, resutl (c+d)x*
	 */
	void AddEdge(int i, double p, PathSet& paths, CutSet& cuts){
		Terms nTerms; // new terms
		edgeCount ++;
#ifdef STATISTICS
		cout << "Adding edge " << edgeCount << endl;
#endif
		FOREACH_STL(t, terms){
			Term xT = t;
			xT.Multiply(i, p, true);
			Term yT = t;
			yT.Multiply(i, p, false);
			if (collapse){ // COLLAPSE CODE
				double cCoef = xT.XStarCollapse(paths);
				if (cCoef == -1.0) // no collapse
					nTerms.push_back(xT);
				else{
					xStarCoef+=cCoef;
#ifdef STATISTICS
					cout << "P";
#endif
				}
				// Symetric for Y
				cCoef = yT.YStarCollapse(cuts);
				if (cCoef == -1.0) // no collapse
					nTerms.push_back(yT);
				else {
					yStarCoef+=cCoef;
#ifdef STATISTICS
					cout << "C";
#endif
				}
			} else { // nonCOLLAPSE CODE
				nTerms.push_back(xT);
				nTerms.push_back(yT);
			}
		}END_FOREACH;
		nTerms.swap(terms); // swap the new and old terms
#ifdef STATISTICS
		cout << "Polynomial has " << terms.size() << " terms " << endl;
#endif
		// nTerms gets distroyed with content of old terms
	}


	// return probability that at least one path is active
	// If all edges in the graph have been corectly added, there should be
	// not Terms left and xStarCoef is the answer
	double ComputeProbability(void){
		if (collapse) { // smart algorithm
			if(terms.size()!=0){
//				cout << "WARNING: Result was extracted with an incomplete graph" << endl;
//				cout << "Number of terms left: " << terms.size() << endl;
			}

			double sum = xStarCoef+yStarCoef;
			const double tolerance = 10e-6;
			if ( (sum>1.0+tolerance) || (sum<1.0-tolerance) ){
//				cout << "WARNING: tolerance level for the solution exceeded. Sum=" << sum << endl;
			}
		} else {
			assert(xStarCoef==0.0); // we should not have collapsed anyting before
			// for each term, check if we can collapse the x part, that we contribute to the answer
			FOREACH_STL(term, terms){// each term is a possible world, all exclusive.
				double coef = term.XStarCollapse(paths); // coef is probability of the world
				if (coef != -1.0){ // successful collapse
					xStarCoef+=coef;
					assert(term.YStarCollapse(cuts) == -1.0);
				} else {
					coef = term.YStarCollapse(cuts); // coef is probability of the world
					if (coef != -1.0) // successful collapse
						yStarCoef+=coef;
				}
			}END_FOREACH;
		}
//		cout << "x*: " << xStarCoef << ", y*: " << yStarCoef << endl;
		return xStarCoef;
	}


};

#endif  // GRAPH_H_
