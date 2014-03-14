/** Compile with
    g++ -o graph Graph.cc -O3 -g

    Run with
    ./graph g1.txt a d
 */


#include <iostream>
#include <fstream>
#include <map>
#include <assert.h>
#include <ctime>
#include <sys/time.h>
#include <cmath>

#include <sys/time.h> // for gettimeofday
#include <stdlib.h> // drand48

// uncomment the next line if a deep check of cuts is desired
#define CHECK_CUTS

using namespace std;

#include "preach.h"

typedef ListDigraph::ArcMap<double> WeightMap;
typedef ListDigraph::NodeMap<string> NodeNames;
typedef map<string, ListDigraph::Node> NameToNode;

#ifdef WIN32
void srand48(int seed){srand(seed);}
double drand48(){return ((double)rand())/RAND_MAX;}
#endif

double getRefWeight(ListDigraph::Arc &arc, WeightMap &wMap, ListDigraph::ArcMap<ListDigraph::Arc> &arcRef){
    return wMap[arcRef[arc]];
}

/*This class encapsulates a variation of the graph, which is defined by a source and target pair,
  along with a current variable edge*/
class STVariation{
    public:
    ListDigraph g; // graph
    WeightMap wMap; // Edge weights of this graph variation,
                    //which will change in preprocessing without touching the original wMap
    ListDigraph::Node source; // source node here
    ListDigraph::Node target; // target node here
    ListDigraph::Arc varArc; // The variable arc here
    ListDigraph::Node sourceRef; // reference to source in original graph
    ListDigraph::Node targetRef; // reference to target in original graph
    ListDigraph::Arc varArcRef; // Reference to the variable arc in the original graph
    bool hasVarArc;
    double constant; // Result: constant part of the probability (a + bp)
    double coeff; // Result: coeff of p variable of the probability (a + bP)


    /*Gets the In-degree of a node*/
    int getNodeInDegree(ListDigraph::Node& node){
        int count = 0;
        for (ListDigraph::InArcIt arc(g, node); arc != INVALID; ++arc)
            count++;
        return count;
    }

    /*Gets the Out-degree of a node*/
    int getNodeOutDegree(ListDigraph::Node& node){
        int count = 0;
        for (ListDigraph::OutArcIt arc(g, node); arc != INVALID; ++arc)
            count++;
        return count;
    }

    /*Reverses the graph: replaces each edge by its reverse edge*/
    void reverseGraph(){
        // Collect a list of all edges
        vector<ListDigraph::Arc> arcs;
        for (ListDigraph::ArcIt arc(g); arc != INVALID; ++arc){
            arcs.push_back(arc);
        }
        FOREACH_STL(arc, arcs){
            g.reverseArc(arc);
        }END_FOREACH;
    }

    /*Removes "Isolated" nodes from the graph
     * An isolated node is the nodes that are not reachable from source or
     * can't reach to sink*/
    void RemoveIsolatedNodes(){
        // First: Make a forward traversal and mark the reachable nodes from source
        Nodes_T forward;
        Bfs<ListDigraph> bfs(g);
        bfs.run(source);
        for (ListDigraph::NodeIt node(g); node != INVALID; ++node){
            if (bfs.reached(node)){
                forward.set(g.id(node));
            }
        }
        // Second: reverse the graph and make a backward traversal
        // and mark the reachable nodes from the sink
        ListDigraph reverseG;
        Nodes_T backward;
        reverseGraph();
        bfs = Bfs<ListDigraph>(g);
        bfs.run(target);
        for (ListDigraph::NodeIt node(g); node != INVALID; ++node){
            if (bfs.reached(node)){
                backward.set(g.id(node));
            }
        }
        // reverse the graph again to return it to original state
        reverseGraph();

        //collect bad nodes
        vector<ListDigraph::Node> badNodes;
        for (ListDigraph::NodeIt node(g); node != INVALID; ++node){
            if (g.id(node) != g.id(source) && g.id(node) != g.id(target) && !(forward[g.id(node)] && backward[g.id(node)]))
                badNodes.push_back(node);
        }

        // Erase all bad nodes
        FOREACH_STL(node, badNodes){
            g.erase(node);
        }END_FOREACH;
    }


    /*
    Removes edges that are self cycles
    */
    void RemoveSelfCycles(){
        for (ListDigraph::ArcIt arc(g); arc != INVALID; ++arc){
            if (g.source(arc) == g.target(arc)){
                g.erase(arc);
            }
        }
    }


    /*Collapses all elementary paths
     * An elementary path: a --> x --> b , with x not connected to anything else
     * we delete x and create a new link a --> b with weight w = weight(a-->x)*weight(x-->b)
     * if an edge a --> b already exists before with weight w', we merge the old edge with the new one with
     * a weight = 1-(1-w)(1-w')
     * */
    void CollapseELementaryPaths(){
        RemoveSelfCycles();
        // repeat until nothing changes
        bool changing = true;
        while(changing){
            changing = false;
            vector<ListDigraph::Node> elementaryNodes;
            for (ListDigraph::NodeIt node(g); node != INVALID; ++node){
                if (node == source || node == target)
                    continue;
                if (hasVarArc && (node == g.source(varArc) || node == g.target(varArc)))
                    continue;
                if (getNodeInDegree(node) == 1 && getNodeOutDegree(node) == 1){
                    // elementary path, mark node to be removed
                    elementaryNodes.push_back(node);
                }
            }
            // handle marked nodes: remove their edges and delete them
            FOREACH_STL(node, elementaryNodes){
                //link before with after
                for (ListDigraph::OutArcIt outArc(g, node); outArc != INVALID; ++outArc){
                    for (ListDigraph::InArcIt inArc(g, node); inArc != INVALID; ++inArc){
                        if (g.source(inArc) == g.target(outArc)) // self loop, just remove the node later
                            continue;
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
    void Preprocess(){
        //As a sanity check: put a crazy value for the variable edge probability
        //to make sure it is never multiplied, EVER!
        if (hasVarArc)
            wMap[varArc] = -1000000.0;

        RemoveIsolatedNodes();
        CollapseELementaryPaths();
        RemoveSelfCycles();
    }


    /*Constructor version where there is a variable edge*/
    STVariation(ListDigraph& _g, ListDigraph::Node& _source, ListDigraph::Node& _target, ListDigraph::Arc& _varArc, WeightMap& _wMap):
                sourceRef(_source), targetRef(_target), varArcRef(_varArc), wMap(g), hasVarArc(true){
        //Build the new st variation
        digraphCopy(_g, g).node(_source, source).node(_target, target).arc(_varArc, varArc).arcMap(_wMap, wMap).run();
        //Preprocess based on the current pair
        Preprocess();
    }

    /*Constructor version where there is no variable edge*/
    STVariation(ListDigraph& _g, ListDigraph::Node& _source, ListDigraph::Node& _target, WeightMap& _wMap):
                sourceRef(_source), targetRef(_target), wMap(g), hasVarArc(false){
        //Build the new st variation
        digraphCopy(_g, g).node(_source, source).node(_target, target).arcMap(_wMap, wMap).run();
        //Preprocess based on the current pair
        Preprocess();
    }

    void Solve(){
        Edges_T allEdges;
        map<int, pair<int, int> > edgeTerminals;
        // Union all edges
        for (ListDigraph::ArcIt arc(g); arc != INVALID; ++arc){
            allEdges.set(g.id(arc));
            edgeTerminals[g.id(arc)] = pair<int, int>(g.id(g.source(arc)), g.id(g.target(arc)));
        }
        // The solver polynomial
        Polynomial poly(edgeTerminals, g.id(source), g.id(target), allEdges);

        // Traverse the graph DFS, add edges to the polynomial one at a time
        // except if it is the variable edge
        Dfs<ListDigraph> dfs(g);
        dfs.init();
        dfs.addSource(source);
        while(!dfs.emptyQueue()){
            ListDigraph::Arc arc = dfs.processNextArc();
            if (hasVarArc && arc == varArc)
                continue; // don't add the variable arc to the polynomial
            poly.AddEdge(g.id(arc), wMap[arc]);
        }
        constant = poly.xStarCoeff();
        coeff = poly.sumUncollapsed();
        //sanity check
        if (!hasVarArc){
            assert(coeff == 0.0);
            poly.checkNoTerms();
        }
    }
};

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

// Read expression file
// return expression map
void ReadExpression(ListDigraph& g, char* filename, ListDigraph::NodeMap<string>& nMap, NameToNode& nodeMap,
                    map<pair<ListDigraph::Node, ListDigraph::Node>, double>& exprMap){
    fstream in(filename);
    while (!in.eof()){
        string s;
        string t;
        double expr = -1.0;
        in >> s >> t >> expr;
        if (expr == -1.0)
            continue;
        ListDigraph::Node sN = FindNode(s, g, nMap, nodeMap );
		ListDigraph::Node tN = FindNode(t, g, nMap, nodeMap);
		exprMap[pair<ListDigraph::Node, ListDigraph::Node>(sN, tN)] = expr;
    }
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
		in >> start >> stop;
//		cout << start << "\t" << stop << "\t" << weight << endl;

		if (start.empty() || stop.empty())
            continue;
		ListDigraph::Node sN = FindNode(start, g, nMap, nodeMap );
		ListDigraph::Node tN = FindNode(stop, g, nMap, nodeMap);
		if (!EdgeExists(g, sN, tN)){
			ListDigraph::Arc a = g.addArc(sN, tN);
//			cout << "edge: " << g.id(a) << endl;
			wMap[a] = drand48();
		}
	}
	in.close();
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

#define POPULATION_SIZE 50
#define GENETIC_ROUNDS 100

class Species{

    /*Calculates the quality of this species based on the wMap values compared to the exprMap values*/
/*    void calculateQuality(ListDigraph& g, map<pair<ListDigraph::Node, ListDigraph::Node>, double>& exprMap){
        double numerator = 0.0;
        double denominator = 0.0;
        for (map<pair<ListDigraph::Node, ListDigraph::Node>, double>::iterator it=exprMap.begin(); it!=exprMap.end(); ++it){
            ListDigraph::Node source = it->first.first;
            ListDigraph::Node target = it->first.second;
            double expr = it->second;
            //create an STVariation and solve it
            STVariation variation(g, source, target, wMap);
            variation.Solve();
            //update numerator and denomenator
            numerator += (expr - variation.constant) * (expr - variation.constant);
            denominator ++;
        }
        quality = 1.0 - sqrt(numerator) / denominator;
    }
*/
    public:

    WeightMap wMap;
    double quality;

    /*Construct a species by random generation*/
    Species(ListDigraph& _g, map<pair<ListDigraph::Node, ListDigraph::Node>, double>& exprMap):wMap(_g){
        /*
        for (ListDigraph::ArcIt arc(g); arc != INVALID; ++arc){
            wMap[arc] = drand48();
        }
        calculateQuality(g, exprMap);
        */
    }

    double getQuality(){
        return quality;
    }
};

/*DOes the cross over step for the genetic algorithm*/
void CrossOver(vector<Species>& population){
    // Use this vector to model the preference of every species when it comes to crossover
    vector<double> preference;
    double sum = 0.0;
    for (vector<Species>::iterator s=population.begin(); s!=population.end(); ++s){
        sum += s->getQuality();
        preference.push_back(sum);
    }
}

/*Provides a good starting point for wMap, based on a genetic algorithm*/
void Geneticize(ListDigraph& g, map<pair<ListDigraph::Node, ListDigraph::Node>, double>& exprMap, WeightMap& wMap){
    // Build initial random population
    vector<Species> population;
    for (int i=0; i<POPULATION_SIZE; i++){
        Species s(g, exprMap);
        population.push_back(s);
    }
    for (int i=0; i<GENETIC_ROUNDS; i++){
    }
}


int main(int argc, char** argv) {
	if (argc < 4) {
		cout << "Usage: prob graph-file expression-file output-file" << endl;
		return -1;
	}

	// initialize srand
    timeval time;
    gettimeofday(&time,NULL);
    srand48((time.tv_sec * 1000) + (time.tv_usec / 1000));

    ListDigraph g;
	WeightMap wMap(g); // keeps track of weights
	NodeNames nNames(g); // node names
	NameToNode nodeMap; // mapping from names to nodes in the graph

	CreateGraph(argv[1], g, nNames, nodeMap, wMap);
	int numNodes, numEdges;
	GetGraphSize(g, numNodes, numEdges);
	cout << endl << "Original graph size: " << numNodes << " nodes, " << numEdges << " edges" << endl;

    //Read expression
    map<pair<ListDigraph::Node, ListDigraph::Node>, double> exprMap;
    ReadExpression(g, argv[2], nNames, nodeMap, exprMap);

	//This map will keep info about edges that don't relate to any s-t graph
	map<ListDigraph::Arc, bool> helplessEdges;

	// Do the following until nothing changes
	bool changing = true;
	while (changing){
	    changing = false;
	    double changeSum = 0.0;
        /*Here we loop over all edges, all sources and all targets, create a graph variant,
          and solve for the current edge as variable*/
        for (ListDigraph::ArcIt arc(g); arc != INVALID; ++arc){
            // if the edge is recorded as helpless, ignore
            if (helplessEdges[arc])
                continue;
            double numerator = 0.0;
            double denominator = 0.0;
            for (map<pair<ListDigraph::Node, ListDigraph::Node>, double>::iterator it=exprMap.begin(); it!=exprMap.end(); ++it){
                ListDigraph::Node source = it->first.first;
                ListDigraph::Node target = it->first.second;
                double expr = it->second;
                //create an STVariation and solve it
                STVariation variation(g, source, target, arc, wMap);
                variation.Solve();

                //update numerator and denomenator
                numerator += variation.coeff * (expr - variation.constant);
                denominator += variation.coeff * variation.coeff;
            }
            //If the edge turns out to be helpless (not related to any s-t, record it)
            if (denominator == 0.0){
                helplessEdges[arc] = true;
                continue;
            }
            // if not, update the arc probability if it's changed
            double newProb = numerator / denominator;
            if (newProb < 0.0)
                newProb = 0.0;
            else if (newProb > 1.0)
                newProb = 1.0;
            if (newProb > wMap[arc] + 0.00001 || newProb < wMap[arc] - 0.00001){
                changing = true;
                changeSum += abs(wMap[arc] - newProb);
                wMap[arc] = newProb;
            }
        }
        cout << "average change in prob: " << changeSum / numEdges << endl;
	}
	//output the results
	ofstream out(argv[3]);
	for (ListDigraph::ArcIt arc(g); arc != INVALID; ++arc){
	    if (helplessEdges[arc]){
	        out << nNames[g.source(arc)] << " " << nNames[g.target(arc)] << " " << "UNKNOWN" << endl;
	    } else {
	        out << nNames[g.source(arc)] << " " << nNames[g.target(arc)] << " " << wMap[arc] << endl;
	    }
	}
	out.close();
	return 0;
}
