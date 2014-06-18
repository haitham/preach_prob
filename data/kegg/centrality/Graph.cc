/** Compile with
    g++ -o graph Graph.cc -O3 -g

    Run with
    ./graph g1.txt a d
 */

#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>
#include <iostream>
#include <fstream>
#include <map>
#include <assert.h>

// uncomment the next line if a deep check of cuts is desired
#define CHECK_CUTS

using namespace std;

#include "Graph.h"

typedef ListDigraph::ArcMap<double> WeightMap;
typedef ListDigraph::NodeMap<string> NodeNames;
typedef map<string, ListDigraph::Node> NameToNode;


const string SOURCE = "SOURCE";
const string SINK = "SINK";
const string METHOD_IE = "ie";
const string METHOD_PM = "pm";
const string METHOD_PMC = "pmc";
const string PRE_YES = "pre";
const string PRE_NO = "nopre";


// map names to nodes
ListDigraph::Node FindNode(string name, ListDigraph& g,
                           NodeNames& nMap,
                           NameToNode& nodeMap
                           ){

	// look in map
	if (nodeMap.find(name) == nodeMap.end()){ // not found, add
		ListDigraph::Node node = g.addNode();
//		cout << "node: " << g.id(node) << endl;
		nodeMap[name]=node;
	}

	ListDigraph::Node node = nodeMap[name];
	nMap[node]=name;

	return node;
}

bool EdgeExists(ListDigraph& g, ListDigraph::Node& source, ListDigraph::Node& target){
	for (ListDigraph::OutArcIt arc(g, source); arc != INVALID; ++arc){
		if (g.id(g.target(arc)) == g.id(target))
			return true;
	}
	return false;
}

// Create the graph from the file
void CreateGraph(char* filename, ListDigraph& g,
                 ListDigraph::NodeMap<string>& nMap,
                 NameToNode& nodeMap,
                 WeightMap& wMap){
	fstream in(filename);

	while (!in.eof()){
		string start;
		string stop;
		double weight = -1.0;
		in >> start >> stop >> weight;

		if (weight==-1.0)
			continue;

//		cout << start << "\t" << stop << "\t" << weight << endl;

		ListDigraph::Node sN = FindNode(start, g, nMap, nodeMap );
		ListDigraph::Node tN = FindNode(stop, g, nMap, nodeMap);
		if (!EdgeExists(g, sN, tN)){
			ListDigraph::Arc a = g.addArc(sN, tN);
//			cout << "edge: " << g.id(a) << endl;
			wMap[a]=weight;
		}
	}
	in.close();
}


/** Function to check that cuts are correct w.r.t a set of paths
    Cut finding algorithm is more complicated, thus more error prone.

    Strategy: for each cut, make sure it intersects ALL paths
    (interupts all paths if removed from the graph).

    If we find a cut that is wrong we scream. This should never happen
    theoretically.
 */

bool checkCutCorrectness(Cut& cut, PathSet& paths){
	FOREACH_STL(path, paths){
		if((path.GetEdges() & cut.GetEdges()).none()){
			return false;
		}
	}END_FOREACH;
	return true;
}

void CheckCuts(PathSet& paths, CutSet& cuts){
	FOREACH_STL(cut, cuts){
		assert(checkCutCorrectness(cut, paths));
	}END_FOREACH;

}

/* Recursive function to find all paths.
   Call initially with an empty path starting at source.
   This is a DFS solution to the problem
*/
void FindAllPaths(Path& partialPath, PathSet& paths, ListDigraph& g,
                  WeightMap& wMap, ListDigraph::Node target){
	ListDigraph::Node curr = g.nodeFromId(partialPath.GetNode());

	// find all neighbors of curr
	for (ListDigraph::OutArcIt a(g, curr); a != INVALID; ++a){
		// for each neightbor
		ListDigraph::Node next = g.target(a);

		// 		assert(next!=curr);
		if (!partialPath.ContainsNode(g.id(next))){ //test if not visited
			Path newP(partialPath, g.id(next), g.id(a), wMap[a]);
			// if not, test if destination (create final path, add to paths)
			if (next==target){
				paths.push_back(newP);
			} else {// if not, form the path with this extra node and recursivelly call function
				FindAllPaths(newP, paths, g, wMap, target);
			}
		}
	}
}

/* Cut finding algorithm
   Initially FindAllCuts is called with the cut {S} alone

   Parameters:
     goodNodes: set of nodes that belong to at least one simple path

 */
void FindAllCuts(Cut& curr, CutSet& cuts,  ListDigraph& g,
                 ListDigraph::Node target, Nodes_T& goodNodes,
                 PathSet& paths /* required only for debugging */
                 ){

	Nodes_T banned = curr.GetBanned(); // list of banned nodes
	Nodes_T choices = curr.GetChoices(); // set of choices available

	// used to eliminate duplicates in cuts
	// go through all the choices

	for (size_t i=choices._Find_first();
	     i!=choices.size(); i=choices._Find_next(i)){
		// skip target choice
		if (g.id(target) == i){
			continue;
		}

		// add i to the cut and recurse
		choices.reset(i); // i not considered by other choices in future
		banned.set(i); // i banned from now on
		Edges_T edges = curr.GetEdges();
		// for each incoming edge into i, delete from edges
		for (ListDigraph::InArcIt a(g, g.nodeFromId(i)); a != INVALID; ++a){
			//				assert(edges[g.id(a) || g.id(a) == i]); // edge must have been on
			edges.reset(g.id(a)); // turn off edge
		}
		Nodes_T choicesN = choices; // new set of choices
		// for each outgoing edge into i, insert into edges
		for (ListDigraph::OutArcIt a(g, g.nodeFromId(i)); a != INVALID; ++a){
			ListDigraph::Node next = g.target(a);
			int _next = g.id(next);
			// for a node to be a future choice it must not be banned and must not be inside
			// or current node
			if (goodNodes[_next] && !banned[_next] && !curr.InsideCut(_next) && (_next)!=i){
				edges.set(g.id(a)); // turn on edge
				choicesN.set(_next);// mark as future choice
			}
		}

		Cut newC(curr, i, edges, choicesN);
		cuts.push_back(newC);
		FindAllCuts(newC, cuts, g, target, goodNodes, paths);
		banned.reset(i);
	}
}


// auxiliary function for InclusionExclusion
double fB(Edges_T& edges, ListDigraph& g, WeightMap& wMap){
	double res=1.0;
	for (size_t i=edges._Find_first();
	     i!=edges.size(); i=edges._Find_next(i)){
		double weight=wMap[g.arcFromId(i)];
		res*=weight;
	}
	return res;
}

/* function to compute the probability using inclusion-exclusion
   Probablility that at least one path is alive.
   A- set of paths
   \sum_{B\subseteq A} -1^{|B|+1} f(B)
   with
   S_B = \cup_{x\in B} E_x with x a path and E_x the set of edges
   f(B)=\prod_{e\in S_B} p(e) with p(e) the probability of edge e

   need to correct final result by substracting 1.0 for empty set
 */
void InclusionExclusion(PathSet& paths, int cPath, Edges_T edges,
                        ListDigraph& g, WeightMap& wMap,
                        bool positive, double& result, unsigned long long int& counter){
	if (cPath >= paths.size()){ // we reached the bottom
		result += positive ? fB(edges, g, wMap) : -fB(edges, g, wMap);
		const unsigned long long stepC = 1ULL<<24;
		counter++;
#ifdef STATISTICS
		if (counter % stepC == 0)
			cerr << "Analyzed " << counter << " path overlaps" << endl;
#endif
		return; // no recursion
	}

	// cPath not included
	InclusionExclusion(paths, cPath+1, edges, g, wMap, positive, result, counter);

	// cPath included
	edges|=paths[cPath].GetEdges(); // add edges to set
	InclusionExclusion(paths, cPath+1, edges, g, wMap, !positive /*flip sign*/, result, counter);
}


/*Function to merge edge sets from paths and cuts*/
void MergeEdgeSets(PathSet& paths, CutSet& cuts, EdgesSet& sets){
	Edges_T set;
	FOREACH_STL(path, paths){
		set = path.GetEdges();
		sets.push_back(set);
	}END_FOREACH;
	FOREACH_STL(cut, cuts){
		set = cut.GetEdges();
		sets.push_back(set);
	}END_FOREACH;
}

/*Function to get list of edge-weight pair from bitset of edges*/
void GetEdgeList(Edges_T& edges, ListDigraph& g, WeightMap& wMap, EdgeList& list){
	FOREACH_BS(e, edges){
	    ListDigraph::Arc arc = g.arcFromId(e);
	    double weight = wMap[arc];
		list.push_back(pair<int,double>(e, weight));
	}
}

/** Function to compute the probability that at
    least a path survives using PGFs

    Parameters:
      goodEdges: bitset of edges in at least one path
      paths: the paths in the graph
      cuts: the cuts
      collapse: if true, use the collapsing algorithm

    Return:
      probability computed
*/
double PGF(ListDigraph& g, Edges_T& goodEdges, WeightMap& wMap, PathSet& paths, CutSet& cuts, bool collapse){
	Polynomial poly(paths, cuts, collapse);

	// Get a unified set of edgesets from paths and cuts
	EdgesSet sets;
	MergeEdgeSets(paths, cuts, sets);

	//set of covered edges
	Edges_T covered;
	while((covered ^ goodEdges).any()){ //at least one good edge not yet covered
		// get the min-count edge set
	  int min = 100000; // -1; was buggy (haitham)
		Edges_T minSet;
		//		cout << "min before: " << min << endl;
		for (int i=0; i<sets.size(); i++){
			Edges_T set = sets.at(i);
			//mask the set with what's covered until now
			set &= ~covered;
			if (set.none()){ //empty set already covered, remove it and continue
				sets.erase(sets.begin()+i);
				i--;
				continue;
			}
			/*
			set = sets.at(i);  // recover the set if it
					   // remains (tamer)
					   */
			int foo = set.count();
			if (foo < min){
				min = set.count();
				minSet = set;
			}


		}
		//		cout << "min after: " << min << endl;
		// now we have a minSet.. cover it and add its edges
		covered |= minSet;
		EdgeList list;
		GetEdgeList(minSet, g, wMap, list);
		poly.AddSetEdges(list, paths, cuts);
	}
	return poly.ComputeProbability();
}



/* Printing of paths. For now the printing of the nodes is in order of
   ID not path format */

void PrintNodes(Nodes_T& nodes, ListDigraph& g, NodeNames& nNames){
	for (size_t i=nodes._Find_first();
	     i!=nodes.size(); i=nodes._Find_next(i)){
		cout << nNames[g.nodeFromId(i)]<< "\t";
	}
	cout << endl;
}

void  PrintPaths(PathSet& paths,  ListDigraph& g, NodeNames& nNames){
	FOREACH_STL(path, paths){
		PrintNodes(path.GetNodes(), g, nNames);
	}END_FOREACH;
}

void PrintCuts(CutSet& cuts, ListDigraph& g, NodeNames& nNames){
	FOREACH_STL(cut, cuts){
		PrintNodes(cut.GetNodes(), g, nNames);
	}END_FOREACH;
}

void ReadList(string filename, vector<string>& list) {
	fstream in(filename.data());
	string item;
//	cout << "reading " << filename <<endl;
	while (!in.eof()) {
		item = "";
		in >> item;

		if (item.empty())
			continue;
//		cout << item << endl;
		list.push_back(item);
	}
	in.close();
}

/*Reads sources and targets and adds a unified source and unified sink to the graph
 * HOW: adds a new SOURCE node to the graph and a 1.0-weight edge to all sources
 * same with SINK and all targets*/
void UnifyTerminals(ListDigraph& g,
		WeightMap& wMap,
		NodeNames& nMap,
		NameToNode& nodeMap,
		string sourcesFile,
		string targetsFile){
	vector<string> sources;
	vector<string> targets;

	//read sources and targets
	ReadList(sourcesFile, sources);
	ReadList(targetsFile, targets);
//	cout << sources.size() << " sources, " << targets.size() << " targets" << endl;

	//create unified source and sink nodes
	ListDigraph::Node source = FindNode(SOURCE, g, nMap, nodeMap);
	ListDigraph::Node sink = FindNode(SINK, g, nMap, nodeMap);

	// add an edge from the new source to all sources
	FOREACH_STL(nodeName, sources){
		ListDigraph::Node node = FindNode(nodeName, g, nMap, nodeMap);
		if (!EdgeExists(g, source, node)){
			ListDigraph::Arc arc = g.addArc(source, node);
			wMap[arc] = 1.0;
		}
	}END_FOREACH;

	// add an edge from all targets to the new sink
	FOREACH_STL(nodeName, targets){
		ListDigraph::Node node = FindNode(nodeName, g, nMap, nodeMap);
		if (!EdgeExists(g, node, sink)){
			ListDigraph::Arc arc = g.addArc(node, sink);
			wMap[arc] = 1.0;
		}
	}END_FOREACH;
}

/*Gets the In-degree of a node*/
int getNodeInDegree(ListDigraph& g, ListDigraph::Node& node){
	int count = 0;
	for (ListDigraph::InArcIt arc(g, node); arc != INVALID; ++arc)
		count++;
	return count;
}

/*Gets the Out-degree of a node*/
int getNodeOutDegree(ListDigraph& g, ListDigraph::Node& node){
	int count = 0;
	for (ListDigraph::OutArcIt arc(g, node); arc != INVALID; ++arc)
		count++;
	return count;
}

/*Removes "Isolated" nodes from the graph
 * An isolated node is the nodes that have no incoming or no outgoing edges
 * cannot participate in any path*/
void RemoveIsolatedNodes(ListDigraph& g, NodeNames& nMap){
	// repeat until nothing is there to remove
	bool changing = true;
	while(changing){
		changing = false;
		vector<ListDigraph::Node> isolatedNodes;
		for (ListDigraph::NodeIt node(g); node != INVALID; ++node){
			if (nMap[node] == SOURCE || nMap[node] == SINK)
				continue;
			if (getNodeInDegree(g, node) == 0){ // handle dead starts
				// mark the node for removal
				isolatedNodes.push_back(node);
			} else if (getNodeOutDegree(g, node) == 0){ // handle dead ends
				// mark the node for removal
				isolatedNodes.push_back(node);
			}
		}
		// remove the marked nodes
		FOREACH_STL(node, isolatedNodes){
			g.erase(node);
			changing = true;
		}END_FOREACH;
	}
}

/*Collapses all elementary paths
 * An elementary path: a --> x --> b , with x not connected to anything else
 * we delete x and create a new link a --> b with weight w = weight(a-->x)*weight(x-->b)
 * if an edge a --> b already exists before with weight w', we merge the old edge with the new one with
 * a weight = 1-(1-w)(1-w')
 * */
void CollapseELementaryPaths(ListDigraph& g, WeightMap& wMap, NodeNames& nMap){
	// repeat until nothing changes
	bool changing = true;
	while(changing){
		changing = false;
		vector<ListDigraph::Node> elementaryNodes;
		for (ListDigraph::NodeIt node(g); node != INVALID; ++node){
			if (nMap[node] == SOURCE || nMap[node] == SINK)
				continue;
			if (getNodeInDegree(g, node) == 1 && getNodeOutDegree(g, node) == 1){
				// elementary path, mark node to be removed
				elementaryNodes.push_back(node);
			}
		}
		// handle marked nodes: remove their edges and delete them
		FOREACH_STL(node, elementaryNodes){
			//link before with after
			for (ListDigraph::OutArcIt outArc(g, node); outArc != INVALID; ++outArc){
				for (ListDigraph::InArcIt inArc(g, node); inArc != INVALID; ++inArc){
					double newWeight;
					bool found = false;
					//Find existing arc between before and after
					for (ListDigraph::OutArcIt arc(g, g.source(inArc)); arc != INVALID; ++arc){
						if (g.target(arc) == g.target(outArc)){
							// a link already exists
							wMap[arc] = 1 - (1 - wMap[arc]) * (1 - wMap[inArc]*wMap[outArc]);
							found = true;
							break;
						}
					}
					if (!found){ // no existing link.. add one
						ListDigraph::Arc newArc = g.addArc(g.source(inArc), g.target(outArc));
						wMap[newArc] = wMap[inArc]*wMap[outArc];
					}
				}
			}
			g.erase(node);
			changing = true;
		}END_FOREACH;
	}
}


/*Does the needed preprocessing of the graph:
 * adding source & sink
 * remove isolated nodes
 * collapse elementary paths
 * see comments of each function for details
 * */
void Preprocess(ListDigraph& g,
		WeightMap& wMap,
		NodeNames& nMap,
		NameToNode& nodeMap,
		string sourcesFile,
		string targetsFile,
		string pre){

	UnifyTerminals(g, wMap, nMap, nodeMap, sourcesFile, targetsFile);
	if (pre == PRE_YES){
		RemoveIsolatedNodes(g, nMap);
		CollapseELementaryPaths(g, wMap, nMap);
	}
}

/*function to compare two <edge, score> pairs*/
bool compare(pair<int, double> a, pair<int, double> b){
	return (b.second < a.second);
}

/* Enumerate all the edges that participate in the paths
 * HOW: Union the edges of all paths
 * score each edge in terms of number of paths/cuts it matches
 * as well as size of the matching paths/cuts
 * Then enumerate those in a pair set of <edge, weight>
 * */
void EnumerateEdges(ListDigraph& g, PathSet& paths, CutSet& cuts, WeightMap& wMap, EdgeList& result){
	Edges_T goodEdges;
	// Union all path edges
	FOREACH_STL(path, paths){
		goodEdges |= path.GetEdges();
	}END_FOREACH;

	// enumerate all edges
	vector< pair<int, double> > edgeScores;
	for (size_t i=goodEdges._Find_first(); i!=goodEdges.size(); i=goodEdges._Find_next(i)){
		double score = 0.0;
//		FOREACH_STL(path, paths){
//			if (path.GetEdges()[i])
//				score += (1.0/path.GetEdges().count());
//		}END_FOREACH;
//		FOREACH_STL(cut, cuts){
//			if (cut.GetEdges()[i])
//				score += (1.0/cut.GetEdges().count());
//		}END_FOREACH;
		edgeScores.push_back(pair<int, double>(i, score));
	}
//	sort(edgeScores.begin(), edgeScores.end(), compare);

	FOREACH_STL(p, edgeScores){
		ListDigraph::Arc arc = g.arcFromId(p.first);
		result.push_back(pair<int, double>(p.first, wMap[arc]));
	}END_FOREACH;


}

/* Eliminates non-minimal cuts from the cut set
 * iterates over the edge representation of each cut and compares it with all other cuts
 * if the set of edges of a cut is a superset of another cut, delete it */
void RefineCuts(CutSet& cuts, PathSet& paths){
	// check first for correctness of all cuts
//	for (int i=0; i<cuts.size(); i++){
//			Cut currenti = cuts.at(i);
//			if (!checkCutCorrectness(currenti, paths)){
//				cuts.erase(cuts.begin() + i);
//				i--;
//			}
//	}
	// check for non-minimality: containment of cuts in other cuts
	for (int i=0; i<cuts.size(); i++){
		Cut currenti = cuts.at(i);
		for (int j=i+1; j<cuts.size(); j++){
			Cut currentj = cuts.at(j);
			// see if one of them is contained in the other
			if ((~currenti.GetEdges() & currentj.GetEdges()).none()){
				// j contained in i: delete i and break
				cuts.erase(cuts.begin() + i);
				i--;
				break;
			}
			if ((~currentj.GetEdges() & currenti.GetEdges()).none()){
				//i contained in j: delete j
				cuts.erase(cuts.begin() + j);
				j--;
			}
		}
	}
}

/*count number of nodes in graph!*/
void GetGraphSize(ListDigraph& g, int& nodes, int& edges){
	nodes = edges = 0;
	for (ListDigraph::NodeIt node(g); node != INVALID; ++node){
		nodes ++;
		for (ListDigraph::OutArcIt arc(g, node); arc != INVALID; ++arc){
			edges ++;
		}
	}
}

int main(int argc, char** argv) {
	if (argc < 6) {
		// argument 4 is the method you want to run
		// ie = inclusion exclusion
		// pm = polynomial multiplication (no collapse)
		// pmc = polynomial multiplication with collapse
		// argument 5 is whether to do preprocessing
		// pre = yes
		// nopre = no
		cout << "Usage: Graph graph-file sources-file targets-file "
				<< "[" << METHOD_IE << "|" << METHOD_PM << "|" << METHOD_PMC << "]"
				<< "[" << PRE_YES << "|" << PRE_NO << "]"
				<< endl;
		return -1;
	}

	ListDigraph g;
	WeightMap wMap(g); // keeps track of weights
	NodeNames nNames(g); // node names
	NameToNode nodeMap; // mapping from names to nodes in the graph

	double result = 1.0;
	string method = argv[4];

	CreateGraph(argv[1], g, nNames, nodeMap, wMap);
	int numNodes, numEdges;
	GetGraphSize(g, numNodes, numEdges);
//	cout << endl << "Original graph size: " << numNodes << " nodes, " << numEdges << " edges" << endl;


	// Read sources and targets and preprocess
	Preprocess(g, wMap, nNames, nodeMap, argv[2], argv[3], argv[5]);
	GetGraphSize(g, numNodes, numEdges);
	cout << endl << "Modified graph size: " << numNodes << " nodes, " << numEdges << " edges" << endl << endl;

	ListDigraph::Node source = FindNode(SOURCE, g, nNames, nodeMap);
	ListDigraph::Node target = FindNode(SINK, g, nNames, nodeMap);

	// print the graph in a file
	// graphToEps(g, "graph.eps");

	Path startP(g.id(source));
	PathSet paths; // collection of paths found by the algorithm
	FindAllPaths(startP, paths, g, wMap, target);

	cout << "Found " << paths.size() << " paths" << endl;
	/*	cout << "path sizes:" << endl;
	FOREACH_STL(path, paths){
		cout << "haitham " << path.GetEdges().count() << endl;
		}END_FOREACH; */
//	PrintPaths(paths, g, nNames);

	if (method == METHOD_IE){
//		cout << "Starting Inclusion-Exclusion" << endl;
		Edges_T empty;

		unsigned long long int counter = 0ULL;
		InclusionExclusion(paths, 0, empty, g, wMap, false, result, counter);
		cout << "Analyzed " << counter << " path combinations" << endl;
		cout << /*"Probability=" << */result << endl;
		return 0;
	}

	// compute the set of nodes that at part of at least one simple path
	Nodes_T goodNodes;
	FOREACH_STL(path, paths)
		{
			path.UnionNodes(goodNodes);
		}END_FOREACH;
//	cout << "Set of good nodes: ";
//	PrintNodes(goodNodes, g, nNames);

	CutSet cuts;
	Nodes_T choices; // steps from source
	Edges_T edges; // edges from source
	// find nodes reachable from source to form the first cut
	for (ListDigraph::OutArcIt a(g, source); a != INVALID; ++a) {
		ListDigraph::Node next = g.target(a);
		int _next = g.id(next);
		if (goodNodes[_next]) {
			choices.set(_next);
			edges.set(g.id(a));
		}

	}
	Cut startC(g.id(source), choices, edges);
	cuts.push_back(startC);
	FindAllCuts(startC, cuts, g, target, goodNodes, paths);

	/* The code above finds sometimes cuts that are not minimal. This
	 happens if the cut isolates groups of nodes from target. A
	 simple fix is to scan from target and try to reach tall the
	 nodes on target side of the cut without getting into the
	 cut. If this does not happen, the cut can be eliminated.
	 */

//	cerr << "Found " << cuts.size() << " non-minimal cuts:" << endl;
//	PrintCuts(cuts, g, nNames);

	// Remove non-minimal redundant cuts
	RefineCuts(cuts, paths);

//	cerr << "Found " << cuts.size() << " minimal cuts:" << endl;

	//cout << "cut sizes:" << endl;
	//	FOREACH_STL(cut, cuts){
	//	cout << "haitham " << cut.GetEdges().count() << endl;
	//}END_FOREACH;

//	PrintCuts(cuts, g, nNames);

	// now check if the cuts are correct
//	CheckCuts(paths, cuts);

//	 //Enumerate the set of edges
//	EdgeList edgeList;
//	EnumerateEdges(g, paths, cuts, wMap, edgeList);
	//enumerating good edges: edges appearing in at least one path
	Edges_T goodEdges;
	// Union all path edges
	FOREACH_STL(path, paths){
		goodEdges |= path.GetEdges();
	}END_FOREACH;
	bool collapse = (method == METHOD_PMC);
//	cout << "Starting polynomial multiplication with collapse = " << collapse << endl;
	result = PGF(g, goodEdges, wMap, paths, cuts, collapse);


	cout << /*"Probability=" << */ ">>" << result << endl;

	return 0;
}
