// Project Identifier: 5949F553E20B650AB0FB2266D3C0822B13D248B0
#include <iostream>
#include <vector>
#include <getopt.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <queue>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>

using namespace std;

enum Terrain {
    land,
    coastline,
    sea
};

struct node
{
    pair<double, double> placement;
    bool kv = false;
    double dv = numeric_limits<double>::infinity();
    size_t pv;
    Terrain terrain;
    node(double x, double y, bool kvVal = false, 
         double dvVal = numeric_limits<double>::infinity())
        : placement(make_pair(x, y)), kv(kvVal), dv(dvVal){}
};

static option long_options[] = {
    {"mode", required_argument, nullptr, 'm'},
    {"help",       no_argument, nullptr, 'h'},
    {nullptr,                0, nullptr,  0}
};

string get_optarg_argument_as_string() {
    if (optarg == nullptr) {
        throw "Error: No mode specified";
    }
    std::string str(optarg); // convert from const char* to std::string.

    // assume the command line will otherwise be correct
    if (str == "MST" || str == "FASTTSP" || str == "OPTTSP") {
        return str;
    } else {
        // no this case will get error
        throw "Error: Invalid command line option";
    }
}

bool MST(vector<node>& graph, bool flight);
double FASTTSP(vector<node>& graph, bool print, vector<size_t>& path);
void OPTTSP(vector<node>& graph);
double MST_distance(const node& n1, const node& n2, bool flight);
double calculateDistance(const node& n1, const node& n2);
void terrain(node& n);
void two_opt(vector<size_t> &tour, const vector<node>& graph);
void furthestInsertion(vector<size_t>& tour, vector<node>& graph, size_t num_vertices);
void randomInsertion(vector<size_t>& tour, vector<node>& graph, size_t num_vertices);

class opttsp
{
public:
    vector<node> graph;
    vector<size_t> path;
    vector<size_t> best_tour;
    vector<vector<double>> dist;
    double best_cost = numeric_limits<double>::infinity();

    // Constructor
    opttsp(const vector<node>& inputGraph) {
        graph = inputGraph;
        size_t num_vertices = graph.size();
        
        // Initialize path with indices {0, 1, 2, ..., num_vertices-1}
        path.resize(num_vertices);
        iota(path.begin(), path.end(), 0);
        
        // Initialize distance matrix
        dist.resize(num_vertices, vector<double>(num_vertices, 0.0));
        for (size_t i = 0; i < num_vertices; ++i) {
            for (size_t j = 0; j < num_vertices; ++j) {
                dist[i][j] = calculateDistance(graph[i], graph[j]);
            }
        }
    }

    void genPerms(size_t permLength, double current_cost);
    bool promising(size_t permLength, double current_cost);
};

int main(int argc, char** argv) {
    // std::ios_base::sync_with_stdio(false); // Turning off Synchronized I/O to speed up
    cout << std::setprecision(2); //Always show 2 decimal places
    cout << std::fixed; //Disable scientific notation for large numbers

    try
    {
        string mode;

        // we remember whether we've already set them in these variables:
        bool mode_was_set = false;


        int choice = 0;
        while ((choice = getopt_long(argc, argv, "m:h", long_options, nullptr)) != -1) {
            switch (choice) {
            case 'm':
                if (mode_was_set) {
                    throw "Multiple modes specified";
                }
                mode_was_set = true;
                mode = get_optarg_argument_as_string();
                break;
            case 'h':
                std::cout <<
                " --mode (MST|FASTTSP|OPTTSP) or -m (MST|FASTTSP|OPTTSP)\n"
                "         MST: Find the minimum spanning tree (MST) of the map\n"
                "         FASTTSP: Find a fast, but not necessarily optimal, solution to the TSP\n"
                "         OPTTSP: Find the optimal solution to the TSP\n"
                " --help or -h\n"
                "         Prints this input specification.\n";
                return 0; // return from main with success
            default:
                // unrecognized option
                throw "Error: Invalid command line option";
            }
        }

        size_t num_vertices;
        vector<node> graph;
        double x, y;

        cin >> num_vertices;
        for (size_t i = 0; i < num_vertices; ++i) {
            cin >> x >> y;
            if (cin) {
                graph.emplace_back(x, y);
            }
        }
        
        if (mode == "MST")
        {
            if (!MST(graph, false))
            {
                throw "Cannot construct MST";
            }
            
            double total_weight = 0.0; 
            for (size_t i = 1; i < num_vertices; i++) {
                total_weight += sqrt(graph[i].dv);
            }
            // print the weight
            cout << total_weight << endl;

            for (size_t i = 1; i < num_vertices; i++)
            {   
                cout << min(i, graph[i].pv) << " " << max(i, graph[i].pv) << endl;
            }  
        }

        if (mode == "FASTTSP")
        {
            vector<size_t> tour(graph.size(), 0);
            FASTTSP(graph, true, tour);
        }
        
        if (mode == "OPTTSP")
        {
            OPTTSP(graph);
        }
    }
    catch(const char* err)
    {
        std::cerr << err << std::endl; // print the error message to std::cerr
        return 1; // returning 1 from main indicates that an error occurred
    }
}

bool MST(vector<node>& graph, bool flight) {
    size_t num_vertices = graph.size();
    graph[0].kv = true;
    graph[0].dv = 0;

    for (size_t i = 0; i < num_vertices; i++)
    {
        terrain(graph[i]);
    }
    

    for (size_t i = 1; i < num_vertices; i++)
    {
        graph[i].dv = MST_distance(graph[0], graph[i], flight);
        graph[i].pv = 0;
    }

    for (size_t i = 0; i < num_vertices; i++)
    {   
        // Find the closest outer between each inner
        double min_dist = numeric_limits<double>::infinity();
        size_t next = 0;
        for (size_t j = 0; j < num_vertices; j++)
        {
            if (graph[j].kv) continue;
            if (graph[j].dv < min_dist)
            {
                min_dist = graph[j].dv;
                next = j;
            }
        }

        // Update the each outer point's closest inner point
        graph[next].kv = true;
        for (size_t j = 0; j < num_vertices; j++)
        {
            if (graph[j].kv) continue;
            double temp_dist = MST_distance(graph[next], graph[j], flight);
            if (temp_dist < graph[j].dv)
            {
                graph[j].dv = temp_dist;
                graph[j].pv = next;
            }
        }
    }

    for (size_t i = 1; i < num_vertices; i++) {
        if (!graph[i].kv)
        {
            return false;
        }
    }

    return true;
}

double MST_distance(const node& n1, const node& n2, bool flight) {
    if (!flight)
    {
        if ((n1.terrain == Terrain::sea && n2.terrain == Terrain::land) ||
            (n1.terrain == Terrain::land && n2.terrain == Terrain::sea)) {
            return numeric_limits<double>::infinity();
        }
    }

    double dx = n1.placement.first - n2.placement.first;
    double dy = n1.placement.second - n2.placement.second;
    return dx * dx + dy * dy;
}

double calculateDistance(const node& n1, const node& n2) {
    double dx = n1.placement.first - n2.placement.first;
    double dy = n1.placement.second - n2.placement.second;
    return sqrt(dx * dx + dy * dy);
}


void terrain(node& n) {
    if (n.placement.first < 0 && n.placement.second < 0) {
        n.terrain = Terrain::sea;
    } else if ((n.placement.first <= 0 && n.placement.second == 0) || 
               (n.placement.first == 0 && n.placement.second <= 0)) {
        n.terrain = Terrain::coastline;
    } else {
        n.terrain = Terrain::land;
    }
}

double FASTTSP(vector<node>& graph, bool print, vector<size_t>& tour) {
    size_t num_vertices = graph.size();

    // furthestInsertion(tour, graph, num_vertices);
    randomInsertion(tour, graph, num_vertices);
    // tour.resize(num_vertices);
    // iota(tour.begin(), tour.end(), 0);

    // two_opt(tour, graph);

    double cost = 0;
    for (size_t i = 0; i < num_vertices; i++)
    {
        cost += calculateDistance(graph[tour[i % num_vertices]], graph[tour[(i + 1) % num_vertices]]);
    }
    if (print)
    {
        cout << cost << "\n";
        for (size_t i = 0; i < num_vertices; i++)
        {
            cout << tour[i] << " ";
        }
        cout << "\n"; 
    } 
    return cost;
}

void two_opt(std::vector<size_t> &tour, const std::vector<node> &nodes) {
    for (size_t i = 0; i < tour.size() - 2; i++)
    {   
        size_t j = i + 2;
        for (;j < tour.size() - 1; j++)
        {
            double change = calculateDistance(nodes[tour[i]], nodes[tour[j]]) + calculateDistance(nodes[tour[i + 1]], nodes[tour[j + 1]])
                           -calculateDistance(nodes[tour[i]], nodes[tour[i + 1]]) - calculateDistance(nodes[tour[j]], nodes[tour[j + 1]]);
            if (change < 0)
            {
                reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                j = i + 2;
            }   
        }  
    }
}

void furthestInsertion(vector<size_t>& tour, vector<node>& graph, size_t num_vertices) {
    tour[0] = 0;
    graph[0].kv = true;
    for (size_t i = 1; i < num_vertices; i++)
    {
        graph[i].dv = calculateDistance(graph[i], graph[0]);
        graph[i].pv = 0;
    }


    // Loop for every vertile O(n)
    for (size_t i = 1; i < num_vertices; i++)
    {   
        // Find the max distance outer point O(n)
        size_t max_idx = 0;
        double max_distance = 0;
        for (size_t j = 0; j < num_vertices; j++)
        {
            if (graph[j].kv) continue;
            double temp_distance = graph[j].dv;
            if (temp_distance > max_distance)
            {
                max_distance = temp_distance;
                max_idx = j;
            }
        }

        graph[max_idx].kv = true;

        // Find the optimal connect place O(n)
        double min_cost = numeric_limits<double>::infinity();
        size_t min_idx = 0;
        for (size_t j = 0; j < i; j++)
        {
            size_t v1 = tour[j];
            size_t v2 = tour[(j + 1) % i];
            double old_edge = calculateDistance(graph[v1], graph[v2]);
            double new_edges = calculateDistance(graph[v1], graph[max_idx])
                            + calculateDistance(graph[v2], graph[max_idx]);
            double temp_cost = new_edges - old_edge;
            if (temp_cost < min_cost) {
                min_cost = temp_cost;
                min_idx = j;
            }
        }
        tour.insert(tour.begin() + min_idx + 1, max_idx);
        // Insert the node at the best position
        for (size_t j = i; j > min_idx; --j) {
            tour[j] = tour[j - 1];
        }

        // Insert the new node at the best position
        tour[min_idx] = max_idx;
        
        // Updata the max distance between inner and outer O(n)
        for (size_t j = 0; j < num_vertices; j++)
        {
            if (graph[j].kv) continue;
            double temp_max_distance = calculateDistance(graph[j], graph[max_idx]);
            if (temp_max_distance > graph[j].dv)
            {
                graph[j].dv = temp_max_distance;
                graph[j].pv = max_idx;
            }
        }
    }
}

void randomInsertion(vector<size_t>& tour, vector<node>& graph, size_t num_vertices) {
    std::vector<bool> visited(num_vertices, false);

    // Random number generator setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, int(num_vertices) - 1);

    // Start with two random nodes
    size_t first = 0;
    size_t second;
    do {
        second = dist(gen);
    } while (second == first);

    tour[0] = first;
    tour[1] = second;
    visited[first] = true;
    visited[second] = true;

    // Loop for every vertile O(n)
    for (size_t i = 2; i < num_vertices; i++)
    {   
        size_t next_node;
        do {
            next_node = dist(gen);
        } while (visited[next_node]);

        visited[next_node] = true;

        // Find the best place to insert it in the tour
        double min_increase = std::numeric_limits<double>::infinity();
        size_t best_position = 0;
        for (size_t i = 0; i < tour.size(); ++i) {
            size_t a = tour[i];
            size_t b = tour[(i + 1) % tour.size()];
            double increase = calculateDistance(graph[a], graph[next_node]) +
                              calculateDistance(graph[next_node], graph[b]) -
                              calculateDistance(graph[a], graph[b]);
            if (increase < min_increase) {
                min_increase = increase;
                best_position = i + 1;
            }
        }

        // Insert the node at the best position
        for (size_t j = i; j > best_position; --j) {
            tour[j] = tour[j - 1];
        }

        // Insert the new node at the best position
        tour[best_position] = next_node;
    }
}

void OPTTSP(vector<node>& graph) {
    vector<size_t> tour(graph.size(), 0);
    double lowerbound = FASTTSP(graph, false, tour);
    opttsp optTSP(graph);
    optTSP.best_cost = lowerbound;
    optTSP.best_tour = tour;
    optTSP.genPerms(1, 0.0);
    cout << optTSP.best_cost << "\n";
    for (size_t i = 0; i < graph.size(); i++)
    {
        cout << optTSP.best_tour[i] << " ";
    }
    cout << "\n";
} 

void opttsp::genPerms(size_t permLength, double current_cost) {
  if (permLength == path.size()) {
    current_cost += dist[path.back()][path[0]];
    if (current_cost < best_cost)
    {
        best_cost = current_cost;
        best_tour = path;
    }
    return;
  }

  if (!promising(permLength, current_cost)) return;

  for (size_t i = permLength; i < path.size(); ++i) {
    swap(path[permLength], path[i]);
    genPerms(permLength + 1, current_cost + dist[path[permLength - 1]][path[permLength]]); // it should be dist[permLength - 1][i] before swap;
    swap(path[permLength], path[i]);
  }
}

bool opttsp::promising(size_t permLength, double current_cost) {
    if (current_cost >= best_cost) return false;
    if (path.size() - permLength < 4) return true;
    

    vector<size_t> unvisited;
    vector<node> mst_graph;
    for (size_t i = permLength; i < path.size(); ++i) {
        unvisited.push_back(path[i]);
        mst_graph.push_back(graph[path[i]]);
    } 

    if (!MST(mst_graph, true)) // will not happen
    {
        throw "Cannot construct MST";
    }

    double mst_cost = 0;
    for (size_t i = 0; i < mst_graph.size(); i++)
    {
        mst_cost += sqrt(mst_graph[i].dv);
    }

    double connect_to_mst = numeric_limits<double>::infinity();
    double connect_to_start = numeric_limits<double>::infinity();
    for (size_t u : unvisited) {
        connect_to_mst = min(connect_to_mst, dist[path[permLength - 1]][u]);
        connect_to_start = min(connect_to_start, dist[u][path[0]]);
    }

    double total_lower_bound = current_cost + mst_cost + connect_to_mst + connect_to_start;

    return total_lower_bound < best_cost;
}