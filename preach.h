#ifndef GRAPH_H_
#define GRAPH_H_

#include <assert.h>

#include <vector>
#include <bitset>
#include <iostream>

#include <lemon/list_graph.h>
#include <lemon/core.h>
#include <lemon/bfs.h>
#include <lemon/dfs.h>

#include "Stl.h"

#define MAX_NODES 512
#define MAX_EDGES 1024

// Uncomment the next macro if detailed statistics of execution are needed
//#define STATISTICS


// macro to smoth up the use of bitsets
#define FOREACH_BS(v, vSet)	  \
	for (size_t v=vSet._Find_first(); v!=vSet.size(); v=vSet._Find_next(v))

using namespace std;

using lemon::ListDigraph;
using lemon::INVALID;
using lemon::Bfs;
using lemon::Dfs;

typedef bitset<MAX_NODES> Nodes_T;
typedef bitset<MAX_EDGES> Edges_T;
typedef vector<Edges_T> EdgesSet;
typedef vector< pair<int,double> > EdgeList; // datastructure to get a set of edges and their weights

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

	double XStarCollapse(int sourceNode, int targetNode, map<int, pair<int, int> > edgeTerminals){
        int edges[MAX_EDGES][2];
        int edgeCount = 0;
        FOREACH_BS(edge, x){
            pair<int, int> terminals = edgeTerminals[edge];
            edges[edgeCount][0] = terminals.first;
            edges[edgeCount][1] = terminals.second;
            edgeCount ++;
        }
        Nodes_T visited;
        visited.set(sourceNode);
        while(true){
            Nodes_T copy = visited;
            for (int i=0; i<edgeCount; i++){
                if (visited[edges[i][0]]){
                    visited.set(edges[i][1]);
                }
            }
            if (copy == visited) // no change
                break;
        }
        if (visited.test(targetNode)){
            return coef; //
        } else {
            return -1.0;  // not found
        }
	}

	/* symetric of above */
	double YStarCollapse(int sourceNode, int targetNode, map<int, pair<int, int> > edgeTerminals, Edges_T& allEdges){
		Edges_T yInverse = allEdges & ~y;
		int edges[MAX_EDGES][2];
        int edgeCount = 0;
        FOREACH_BS(edge, yInverse){
            pair<int, int> terminals = edgeTerminals[edge];
            edges[edgeCount][0] = terminals.first;
            edges[edgeCount][1] = terminals.second;
            edgeCount ++;
        }
        Nodes_T visited;
        visited.set(sourceNode);
        while(true){
            Nodes_T copy = visited;
            for (int i=0; i<edgeCount; i++){
                if (visited[edges[i][0]]){
                    visited.set(edges[i][1]);
                }
            }
            if (copy == visited) // no change
                break;
        }
        if (!visited.test(targetNode)){
            return coef; //
        } else {
            return -1.0;  // not found
        }
	}

	//get the term coeff
	double getCoef(){
	    return coef;
	}

};

/* Class to implement polynomial multiplicatio to carry out the
   probabilistic computation

   The Polynomial that needs to be produced is:
   \prod_{e\E} (p_e x_e+q_e y_e) with p_e is prob of edge e, q_e=1-p_e
*/

class Polynomial{
	typedef vector<Term> Terms;

    map<int, pair<int, int> > edgeTerminals;
    int sourceNode;
    int targetNode;
    Edges_T allEdges;
	bool collapse; // should we collapse?

	long double xStarCoef; // coefficient of xStar
	long double yStarCoef; // same for y
	Terms terms;
	int edgeCount;
 public:
	// main constructor. Requires paths and cuts
	// if collapse=true, we collapse as we go.
	Polynomial(map<int, pair<int, int> >& _edgeTerminals, int _source, int _target, Edges_T& _allEdges):
		sourceNode(_source), targetNode(_target), edgeTerminals(_edgeTerminals), allEdges(_allEdges),
		collapse(true), xStarCoef(0.0), yStarCoef(0.0){
		// Adding the empty term with 1.0 coefficient
		terms.push_back(Term());
		edgeCount = 0;
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
	void AddEdge(int i, double p){
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
				double cCoef = xT.XStarCollapse(sourceNode, targetNode, edgeTerminals);
				if (cCoef == -1.0) // no collapse
					nTerms.push_back(xT);
				else{
					xStarCoef+=cCoef;
#ifdef STATISTICS
					cout << "P";
#endif
				}
				// Symetric for Y
				cCoef = yT.YStarCollapse(sourceNode, targetNode, edgeTerminals, allEdges);
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

	void checkNoTerms(){
	    assert(terms.empty());
	}


	// get the coeff of xStar
	double xStarCoeff(){
	    //sanity check
	    assert(xStarCoef >= 0.0);
	    return xStarCoef;
	}

	// get the sum of uncollapsed coeffs
	double sumUncollapsed(){
	    double sum = 0.0;
	    FOREACH_STL(term, terms){
	        sum += term.getCoef();
	    }END_FOREACH;
	    //sanity check
	    assert(sum >= 0.0);
	    return sum;
	}

};

#endif  // GRAPH_H_
