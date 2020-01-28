#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <list>

#define INF std::numeric_limits<int>::max()

class Graph{
public:
    unsigned int num_edges;
    unsigned int num_vertices;
    std::string file_name;
    std::vector<std::vector<double>> matrix;
    std::vector<std::set<std::pair<unsigned int,double>>> adj_list;
    //konstruktor, destruktor
    Graph() : num_edges(0), num_vertices(0), file_name(""), matrix(), adj_list() {}
    ~Graph(){}

    //ostale funkcije
    bool load(std::string);     //funkcija za ucitavanje grafa iz datoteke (true ako uspjesno, false inace)
    void print_in_file(std::string);    //funkcija za printanje grafa u datoteku
    void print_graph() const;       //funkcija za printanje (samo za provjeru)
    double adjacent(unsigned int, unsigned int);        //funkcija za provjeru susjedstva, vraca tezinu ako su cvorovi susjedi, 0 inace
    std::pair<std::vector<int>, double> shorthest_path(unsigned int,unsigned int);       //najkraci put od jednog do drugog vrha
    unsigned int find_degree(unsigned int);     //stupanj vrha
};

bool Graph::load(std::string file){
    std::ifstream myfile;
    std::string line;
    myfile.open(file);
    if(!myfile.is_open()) return false;
    file_name=file;
    getline(myfile, line);      //broj cvorova i broj bridova
    std::stringstream ss1(line);
    ss1 >> num_vertices >> num_edges;
    int a,b;
    double c;

    std::vector<std::vector<double>> v(num_vertices, std::vector<double>(num_vertices,static_cast<double>(0)));
    matrix=v;
    std::vector<std::set<std::pair<unsigned int,double>>> ad(num_vertices);
    adj_list=ad;
    for(unsigned int i=0; i<num_edges; ++i){
        getline(myfile, line);
        std::stringstream ss2(line);
        ss2 >> a >> b >> c;
        unsigned int u_a=static_cast<unsigned int>(a);
        unsigned int u_b=static_cast<unsigned int>(b);
        matrix[u_a][u_b]=c;
        matrix[u_b][u_a]=c;
        adj_list[u_b].insert(std::make_pair(u_a,c));
        adj_list[u_a].insert(std::make_pair(u_b,c));
    }
    myfile.close();
    std::cout<<num_vertices<<" vrhova."<<std::endl;
    std::cout<<num_edges<<" bridova."<<std::endl;
    return true;
}

void Graph::print_in_file(std::string outname){
    //printanje cvorova u datoteku
    std::ofstream file_n;
    file_n.open(outname+"_vertices.csv");
    file_n<<"Id;Label\n";
    for(unsigned int i=0; i<num_vertices; ++i)
        file_n<<i<<";"<<i<<"\n";
    file_n.close();
    //printanje bridova u datoteku, bez tezine (nije nam potrebna)
    std::ofstream file_v;
    file_v.open(outname+"_edges.csv");
    file_v<<"Source;Target;Type\n";
    for(unsigned int i=0; i<num_vertices; ++i)
        for(std::set<std::pair<unsigned int,double>>::iterator it=adj_list[i].begin(); it!=adj_list[i].end(); ++it)
            if(i<it->first) file_v<<i<<";"<<it->first<<";"<<"Undirected\n";
    file_v.close();
}

void Graph::print_graph() const{
    std::cout<<"Datoteka: "<<file_name<<std::endl;
    for(unsigned int i=0; i<num_vertices; ++i){
        std::cout<<i<<": ";
        for (std::set<std::pair<unsigned int, double>>::iterator it=adj_list[i].begin(); it!=adj_list[i].end();++it)
            std::cout<<"("<<it->first<<", "<<it->second<<")";
        std::cout<<std::endl;
        }
}

double Graph::adjacent(unsigned int i, unsigned int j){
    return matrix[i][j];
}

//vraca tezinu puta i vektor parent koji omogucava rekonstrukciju puta po potrebi
std::pair<std::vector<int>, double> Graph::shorthest_path(unsigned int i, unsigned int j){
    std::vector<int> parent(num_vertices, -2);  //parent[n] sadrzi cvor iz kojeg smo dosli u cvor n
    parent[i]=-1;
    std::vector<double> dist(num_vertices, INF); // dist[n] sadrzi duljinu najkraceg puta od i do n
    std::vector<bool> sptSet(num_vertices, false); // sptSet[n]=true ako se cvor n nalazi u najkracem putu do j, false inace

    // udaljenost od i do i je uvijek 0
    dist[i]=0;
    for(unsigned int count=0; count < num_vertices; ++count){
        unsigned int u, min_index=0;
        double min=INF;
        for(u=0; u<num_vertices; ++u)
            if (sptSet[u] == false && dist[u] <= min){ min = dist[u]; min_index = u; }
        u=min_index;
        if(u==j) break;
        sptSet[u]=true;
        for(std::set<std::pair<unsigned int,double>>::iterator it=adj_list[u].begin(); it!=adj_list[u].end(); ++it){
            if( !sptSet[it->first] && dist[u] < INF && dist[u]+matrix[u][it->first] < dist[it->first]){
                dist[it->first]=dist[u]+matrix[u][it->first];
                parent[it->first]=static_cast<int>(u);
            }
        }
    }
    return std::make_pair(parent, dist[j]);
}

unsigned int Graph::find_degree(unsigned int n){
    return static_cast<unsigned int>(adj_list[n].size());
}

#endif // GRAPH_H
