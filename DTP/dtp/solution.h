#ifndef SOLUTION_H
#define SOLUTION_H

#include "graph.h"
#include <iostream>
#include <set>
#include <random> //mersenne-twister
#include <functional>
#include <vector>
#include <limits>
#include <chrono>
#include <algorithm>
#include <queue>
#include <map>

extern long long seed;          //seed i random_number_generator se pozove samo jednom na pocetku izvodenja programa
extern std::mt19937 random_number_generator;

#define INF std::numeric_limits<int>::max()

bool double_equals(double a, double b, double epsilon = 0.001){
    return std::abs(a-b)<epsilon;
}

unsigned int random_num(unsigned int i, unsigned int j){
    std::uniform_int_distribution<unsigned int> number_distribution(i,j);
    return number_distribution(random_number_generator);
}

unsigned int random_from_set(std::set<unsigned int> S){
    int br=static_cast<int>(random_num(0,static_cast<unsigned int>(S.size()-1)));
    return *(std::next(S.begin(), br));
}

class Solution{
public:
    double fitness;
    unsigned int num_edges;
    unsigned int num_vertices;
    std::vector<std::set<std::pair<unsigned int,double>>> adj_list;
    std::set<unsigned int> vertices_in_solution;

    //konstruktor:
    Solution(Graph&, int);     //int=1 iz ABC, int=2 iz SSGA, inace prazno dijete
    //konstruktor za testiranje, graf kopira u rjesenje
    Solution(Graph& G) : fitness(0), num_edges(G.num_edges), num_vertices(G.num_vertices), adj_list(G.adj_list), vertices_in_solution() {
        for(unsigned int i=0; i<num_vertices; ++i) vertices_in_solution.insert(i);
        for(unsigned int i=0; i<G.num_vertices; ++i)
            for(std::set<std::pair<unsigned int,double>>::iterator it=G.adj_list[i].begin(); it!=G.adj_list[i].end(); ++it)
                if(i<it->first) fitness+=it->second;
    }
    //copy konstruktor:
    Solution(const Solution& s) : fitness(s.fitness), num_edges(s.num_edges), num_vertices(s.num_vertices), adj_list(s.adj_list), vertices_in_solution(s.vertices_in_solution){}
    //destruktor:
    ~Solution(){}
    //operator pridruzivanja:
    Solution& operator=(const Solution& s){
        fitness=s.fitness;
        num_edges=s.num_edges;
        num_vertices=s.num_vertices;
        adj_list=s.adj_list;
        vertices_in_solution=s.vertices_in_solution;
        return *this;
    }

    //operatori usporedivanja:
    friend bool operator==(const Solution& first, const Solution& second){
        if(!double_equals(first.fitness, second.fitness) || first.num_edges!=second.num_edges) return false;
        if(first.adj_list!=second.adj_list) return false;
        return true;
    }
    friend bool operator!=(const Solution& first, const Solution& second){
        return !(first==second);
    }
    friend bool operator<(const Solution& first,const Solution& second){
        if(!double_equals(first.fitness, second.fitness) && first.fitness<second.fitness) return true;         //rjesenje je bolje ako je fitness manji
        if(double_equals(first.fitness, second.fitness) && first.num_edges<second.num_edges) return true;      //ili ako su fitness jednaki, ali ima manje bridova
        return false;
    }

    //ostale funkcije:
    bool feasible(Graph&) const;    //dopustivost
    void print_in_file(std::string,unsigned int) const;      //funkcija za printanje grafa u datoteku
    void print_solution(Graph&) const;        //funkcija za printanje rjesenja, samo za provjeru
    double adjacent(unsigned int, unsigned int);        //provjera jesu li dva cvora susjedna i ako jesu vraca tezinu brida (mislim da ne treba)
    unsigned int find_degree(unsigned int);         //funkcija za racunanje stupnja cvora u rjesenju
    bool cut(Graph&);       //funkcija za rezanje listova (vraca true ako je odrezao nesto, false inace)
    void MST(Graph&);       //funkcija za odredivanje minimalnog razapinjuceg stabla na podgrafu induciranom skupom S (Primov algoritam)
    double get_fitness() const { return fitness; }        //funkcija koja vraca fitness
    friend Solution cross(const Solution&,const Solution&, Graph&);     //funkcija koja kriza 2 rjesenja, vraca dijete
    friend Solution mutate(const Solution&m, Graph&);        //funkcija koja mutira rjesenje, vraca dijete
};

Solution::Solution(Graph& G, int f) : fitness(0), num_edges(0), num_vertices(0), adj_list(G.num_vertices), vertices_in_solution() {
    //vertices_in_solution=S, edges=DT, problem=GV (u radu)
    if(f==1){
        std::vector<unsigned int> visited(G.num_vertices, 0);
        std::vector<unsigned int> D;
        std::set<unsigned int> V;
        for(unsigned int i=0; i<G.num_vertices; ++i) V.insert(i);

        //slucajno biramo vrh v
        unsigned int v=random_num(0,G.num_vertices-1);

        //dodajemo ga u skup S i posjecujemo i izabcujemo iz V=V/S
        vertices_in_solution.insert(v); visited[v]=1; V.erase(v);

        //posjecujemo sve njegove susjede
        for(std::set<std::pair<unsigned int, double>>::iterator it=G.adj_list[v].begin();it!=G.adj_list[v].end(); ++it) visited[it->first]=1;

        //update D
        D.clear();
        for(std::set<unsigned int>::iterator it=V.begin(); it!=V.end(); ++it){
            int d=0;
            for(std::set<std::pair<unsigned int, double>>::iterator it1=G.adj_list[*it].begin();it1!=G.adj_list[*it].end(); ++it1)
                if(visited[it1->first]==0) ++d;
            if(d) D.push_back(*it);
        }

        int stop;
        while(1){
            stop=1;
            unsigned int ind_u=random_num(0,static_cast<unsigned int>(D.size()-1));
            unsigned int u=D[ind_u];
            std::vector<int> path=G.shorthest_path(v,u).first;
            unsigned int a=u;
            unsigned int b=static_cast<unsigned int>(path[u]);
            while(path[a]!=-1){
                vertices_in_solution.insert(a);
                V.erase(a);
                visited[a]=1;
                for(std::set<std::pair<unsigned int, double>>::iterator it=G.adj_list[a].begin(); it!=G.adj_list[a].end(); ++it) visited[it->first]=1;

                std::pair<unsigned int, double> p1(a, G.adjacent(a,b));
                std::pair<unsigned int, double> p2(b, G.adjacent(a,b));

                if(adj_list[b].insert(p1).second) { ++num_edges; fitness+=p1.second; };
                adj_list[a].insert(p2);

                a=static_cast<unsigned int>(path[a]);
                b=static_cast<unsigned int>(path[a]);
            }
            vertices_in_solution.insert(a);
            V.erase(a);
            visited[a]=1;
            for(std::set<std::pair<unsigned int, double>>::iterator it=G.adj_list[a].begin(); it!=G.adj_list[a].end(); ++it) visited[it->first]=1;

            D.clear();
            for(std::set<unsigned int>::iterator it=V.begin(); it!=V.end(); ++it){
                int d=0;
                for(std::set<std::pair<unsigned int, double>>::iterator it1=G.adj_list[*it].begin();it1!=G.adj_list[*it].end(); ++it1)
                    if(visited[it1->first]==0) ++d;
                if(d) D.push_back(*it);
            }

            for(unsigned int i=0; i<G.num_vertices; ++i) if(visited[i]==0) stop=0;
            if(stop) break;
        }
        num_vertices=static_cast<unsigned int>(vertices_in_solution.size());
        //rezemo listove iz DT
        while(1){
            if(!cut(G)) break;
        }
        MST(G);
        while(1){
            if(!cut(G)) break;
        }
    }

    if(f==2){
        std::set<unsigned int> U;   //posjeceni cvorovi koji nisu u rjesenju
        std::vector<unsigned int> in_solution(G.num_vertices,0);
        //odaberemo random cvor iz grafa, njega dodamo u rjesenje, a njegove susjede u U
        unsigned int rand=random_num(0,G.num_vertices-1);
        vertices_in_solution.insert(rand);
        in_solution[rand]=1;
        for(std::set<std::pair<unsigned int, double>>::iterator it=G.adj_list[rand].begin(); it!=G.adj_list[rand].end(); ++it) U.insert(it->first);
        //sve dok svi cvorovi ne budu ili u rjesenju ili u U
        while(U.size()+vertices_in_solution.size()!=G.num_vertices){
            //pronadi brid izmedu jednog cvora iz U i jednog iz S=vertices_in_solution i posjeti susjede
            unsigned int vx=random_from_set(vertices_in_solution);
            unsigned int vy=random_from_set(U);
            if(G.adjacent(vx,vy)>0){
                U.erase(vy);
                vertices_in_solution.insert(vy);
                in_solution[vy]=true;
                adj_list[vx].insert(std::make_pair(vy,G.adjacent(vx,vy)));
                adj_list[vy].insert(std::make_pair(vx,G.adjacent(vx,vy)));
                fitness+=G.adjacent(vx,vy);
                ++num_edges;
                for(std::set<std::pair<unsigned int, double>>::iterator it=G.adj_list[vx].begin(); it!=G.adj_list[vx].end(); ++it)
                    if(!in_solution[it->first]) U.insert(it->first);
                for(std::set<std::pair<unsigned int, double>>::iterator it=G.adj_list[vy].begin(); it!=G.adj_list[vy].end(); ++it)
                    if(!in_solution[it->first]) U.insert(it->first);
            }

        }
        num_vertices=static_cast<unsigned int>(vertices_in_solution.size());
        while(1){
          if(!cut(G)) break;
        }
        MST(G);
    }
}


void Solution::print_in_file(std::string name, unsigned int num) const{
    //printanje bridova u datoteku, bez tezine(za sad)
    std::cout<<"printa: "<<name<<std::endl;
    std::ofstream file_v;
    file_v.open(name+"_edges.csv");
    file_v<<"Source;Target;Type\n";
    for(unsigned int i=0; i<num; ++i)
        for(std::set<std::pair<unsigned int,double>>::iterator it=adj_list[i].begin(); it!=adj_list[i].end(); ++it)
            if(i<it->first) file_v<<i<<";"<<it->first<<";"<<"Undirected\n";
    file_v.close();
}

void Solution::print_solution(Graph &G) const{
    for(unsigned int i=0; i<G.num_vertices; ++i)
        for(std::set<std::pair<unsigned int,double>>::iterator it=adj_list[i].begin(); it!=adj_list[i].end(); ++it)
            if(i<it->first) std::cout<<"("<<i<<","<<it->first<<","<<it->second<<")"<<std::endl;
}

bool Solution::feasible(Graph& G) const{
    for(int unsigned i=0; i<G.num_vertices; ++i){
        int flag=0;
        //ako cvor nije u rjesenju pogledaj je li neki od njegovih susjeda u rjesenju
        if(vertices_in_solution.find(i)==vertices_in_solution.end()){
            for(std::set<std::pair<unsigned int, double>>::iterator it=G.adj_list[i].begin(); it!=G.adj_list[i].end(); ++it){
                if(vertices_in_solution.find(it->first)!=vertices_in_solution.end()) { flag=1; break; }
            }
            if(flag==0) return false;
        }
    }
    return true;
}

bool Solution::cut(Graph& G){
    std::set<unsigned int> vertices_to_delete;
    for(std::set<unsigned int>::iterator it=vertices_in_solution.begin(); it!=vertices_in_solution.end(); ){
        int flag=0;
        if(find_degree(*it)==1){
            flag=1;
            for(std::set<std::pair<unsigned int, double>>::iterator it1=G.adj_list[*it].begin(); it1!=G.adj_list[*it].end(); ++it1){
                int flag2=0;
                //mislim da je ova verzija brza jer je manje vrhova u rjesenju nego sto je susjeda u G (u vecim grafovima)
                for(std::set<unsigned int>::iterator it2=vertices_in_solution.begin(); it2!=vertices_in_solution.end(); ++it2){
                    if(*it==*it2) continue;
                    if(G.adjacent(it1->first,*it2)>0) { flag2=1; break; }
                }
                if(!flag2) { flag=0; break; }
            }
        }
        if(flag){ vertices_to_delete.insert(*it);  vertices_in_solution.erase(*(it++)); --num_vertices; }
        else ++it;
    }
    for(std::set<unsigned int>::iterator it=vertices_to_delete.begin(); it!=vertices_to_delete.end(); ++it){
        for(std::set<std::pair<unsigned int,double>>::iterator it1=adj_list[*it].begin(); it1!=adj_list[*it].end(); ++it1){
            adj_list[it1->first].erase(std::make_pair(*it,it1->second));
            fitness-=it1->second;
            --num_edges;
        }
        adj_list[*it].clear();
    }
    if(vertices_to_delete.size()>0) return true;
    return false;
}


unsigned int Solution::find_degree(unsigned int n){
    return static_cast<unsigned int>(adj_list[n].size());
}

void Solution::MST(Graph &G){
    std::vector<int> parent(num_vertices,-2);   //za rekonstrukciju mst
    std::vector<unsigned int> mstSet(num_vertices, 0);   //da znamo koji jos vrhovi nisu u mst
    std::vector<double> key(num_vertices, INF);

    unsigned int rand=random_num(0,num_vertices-1);
    key[rand]=0;
    parent[rand]=-1;   //odabrani cvor je korijen mst

    for(unsigned int i=0; i<num_vertices-1; ++i) {
        //odabiremo cvor s najmanjom vrijednost key od onih koji vec nisu u mst
        double min_key=INF;
        unsigned int u=INF;
        for(unsigned int v=0; v<num_vertices; ++v)
            if(mstSet[v]==0 && key[v]<min_key) { min_key=key[v]; u=v; }

        // dodajemo taj cvor u mst
        mstSet[u]=1;

        // update key vrijednosti i postavi roditelje
        unsigned int vertex_u=*(std::next(vertices_in_solution.begin(), static_cast<int>(u)));
        //stavila G.adj_list umjesto adj_list
        for(std::set<std::pair<unsigned int,double>>::iterator it=G.adj_list[vertex_u].begin(); it!=G.adj_list[vertex_u].end(); ++it){
            unsigned int vertex_v=it->first;
            if(vertices_in_solution.find(vertex_v)!=vertices_in_solution.end()){
                unsigned int v=0;
                for(std::set<unsigned int>::iterator it1=vertices_in_solution.begin(); it1!=vertices_in_solution.end(); ++it1){
                    if(*it1==vertex_v) break;
                    ++v;
                }
                if(mstSet[v]==0 && it->second < key[v]) { parent[v]=static_cast<int>(u); key[v]=it->second; }
            }
        }
    }
    fitness=0;
    num_edges=0;
    for(unsigned int i=0; i<G.num_vertices; ++i) adj_list[i].clear();
    for(unsigned int ver=0; ver<num_vertices; ++ver){
        if(ver==rand) continue;
        unsigned int vertex=*(std::next(vertices_in_solution.begin(), static_cast<int>(ver)));
        unsigned int par=static_cast<unsigned int> (parent[ver]);
        unsigned int vertex_parent=*(std::next(vertices_in_solution.begin(), static_cast<int>(par)));
        double weight=G.adjacent(vertex,vertex_parent);
        adj_list[vertex].insert(std::make_pair(vertex_parent,weight));
        adj_list[vertex_parent].insert(std::make_pair(vertex,weight));
        fitness+=weight;
        ++num_edges;
    }
}

double Solution::adjacent(unsigned int i, unsigned int j){
    unsigned int smaller=std::min(i,j);
    unsigned int bigger=std::max(i,j);
    for(std::set<std::pair<unsigned int,double>>::iterator it=adj_list[bigger].begin(); it!=adj_list[bigger].end(); ++it)
        if(it->first==smaller) return it->second;
    return 0;
}

Solution cross(const Solution& parent1,const Solution& parent2,Graph& problem){
    Solution child(problem,0);  //prazno dijete
    std::vector<unsigned int> visited(problem.num_vertices, 0); //svi vrhovi u grafu neoznaceni
    unsigned int rand=random_num(0,1); //random briramo izmedu 1 i 2 roditelja
    std::set<unsigned int> vertices1=parent1.vertices_in_solution;
    std::set<unsigned int> vertices2=parent2.vertices_in_solution;

    //uzimamo prvi vrh iz odabranog roditelja
    unsigned int v;
    if(!rand){
        v=*(vertices1.begin());
        vertices1.erase(v);
    }
    else{
        v=*(vertices2.begin());
        vertices2.erase(v);
    }

    child.vertices_in_solution.insert(v);
    //posjecujemo v i sve njegove susjede
    visited[v]=1;
    for(std::set<std::pair<unsigned int,double>>::iterator it=problem.adj_list[v].begin(); it!=problem.adj_list[v].end(); ++it) visited[it->first]=1;

    //sve dok ne posjetimo sve vrhove, odabiremo random roditelja i uzimamo njegov prvi vrh
    int stop=0;
    while(1){
        stop=1;
        if(vertices1.size() && vertices2.size()){
            rand=random_num(0,1);
            if(!rand){
                v=*(vertices1.begin());
                vertices1.erase(v);
                if(child.vertices_in_solution.find(v)!=child.vertices_in_solution.end()) continue;  //ako je vrh vec u rjesenju idemo dalje
            }
            else {
                v=*(vertices2.begin());
                vertices2.erase(v);
                if(child.vertices_in_solution.find(v)!=child.vertices_in_solution.end()) continue;  //ako je vrh vec u rjesenju idemo dalje
            }
        }
        if(vertices1.size() && !vertices2.size()){
            v=*(vertices1.begin());
            vertices1.erase(v);
            if(child.vertices_in_solution.find(v)!=child.vertices_in_solution.end()) continue;  //ako je vrh vec u rjesenju idemo dalje
        }

        if(!vertices1.size() && vertices2.size()){
            v=*(vertices2.begin());
            vertices2.erase(v);
            if(child.vertices_in_solution.find(v)!=child.vertices_in_solution.end()) continue;  //ako je vrh vec u rjesenju idemo dalje
        }
        //ako postoji brid izmedu v i nekog vrha u rjesenju povezi ih i posjeti susjede
        int flag=0;
        for(std::set<unsigned int>::iterator it=child.vertices_in_solution.begin(); it!=child.vertices_in_solution.end(); ++it)
            if(problem.adjacent(v,*it)>0){
                child.vertices_in_solution.insert(*it);
                if(child.adj_list[v].insert(std::make_pair(*it, problem.adjacent(v,*it))).second) { ++child.num_edges; child.fitness+=problem.adjacent(v,*it); };
                child.adj_list[*it].insert(std::make_pair(v, problem.adjacent(v,*it)));
                visited[*it]=1;
                visited[v]=1;
                for(std::set<std::pair<unsigned int, double>>::iterator it2=problem.adj_list[v].begin(); it2!=problem.adj_list[v].end(); ++it2) visited[it2->first]=1;
                for(std::set<std::pair<unsigned int, double>>::iterator it3=problem.adj_list[*it].begin(); it3!=problem.adj_list[*it].end(); ++it3) visited[it3->first]=1;
                flag=1;
                break;
            }
        //ako ne postoji brid, pronadi put s najvecim potencijalom od v do nekog vrha u rjesenju
        //potencijal=(broj neposjecenih vrhova koji su u putu)/(tezina puta)
        if(flag==0){
            std::pair<std::vector<int>,double> sp, sp_best;
            unsigned int sp_best_target=problem.num_vertices; //postavimo na nesto da se kompajler ne buni
            double potential=-1;
            for(std::set<unsigned int>::iterator it=child.vertices_in_solution.begin(); it!=child.vertices_in_solution.end(); ++it){
                sp=problem.shorthest_path(v,*it);
                int cnt_unvisited=0;
                unsigned int a=*it;
                //brojimo neposjecene cvorove u putu
                while(sp.first[a]!=-1){
                    if(!visited[a]) ++cnt_unvisited;
                    a=static_cast<unsigned int>(sp.first[a]);
                }
                if(!visited[a]) ++cnt_unvisited;

                if(static_cast<double>(cnt_unvisited)/sp.second > potential){      //ako je potencijal bolji, spremi ga u sp_best
                    potential=static_cast<double>(cnt_unvisited)/sp.second;
                    sp_best=sp;
                    sp_best_target=*it;     //moramo zapamtiti krajnji cvor da mozemo reproducirati najbolji put
                }
            }
            //cijeli put i vrhove dodaj u child
            unsigned int a=sp_best_target;
            unsigned int b=static_cast<unsigned int>(sp_best.first[a]);
            while(sp_best.first[a]!=-1){
                child.vertices_in_solution.insert(a);
                visited[a]=1;
                for(std::set<std::pair<unsigned int, double>>::iterator it=problem.adj_list[a].begin(); it!=problem.adj_list[a].end(); ++it) visited[it->first]=1;

                std::pair<unsigned int, double> p1(a, problem.adjacent(a,b));
                std::pair<unsigned int, double> p2(b, problem.adjacent(a,b));

                if(child.adj_list[b].insert(p1).second) { ++child.num_edges; child.fitness+=problem.adjacent(a,b); };
                child.adj_list[a].insert(p2);

                a=static_cast<unsigned int>(sp_best.first[a]);
                b=static_cast<unsigned int>(sp_best.first[a]);
            }
            child.vertices_in_solution.insert(a);
            visited[a]=1;
            for(std::set<std::pair<unsigned int, double>>::iterator it=problem.adj_list[a].begin(); it!=problem.adj_list[a].end(); ++it) visited[it->first]=1;
        }
        child.vertices_in_solution.insert(v);

        for(unsigned int i=0; i<problem.num_vertices; ++i) if(visited[i]==0) stop=0;
        if(stop) break;
    }
    child.num_vertices=static_cast<unsigned int>(child.vertices_in_solution.size());
    return child;
}

Solution mutate(const Solution& parent,Graph& problem,double pm){
    Solution C(parent);
    std::set<unsigned int> V; //svi vrhovi grafa
    for(unsigned int i=0; i<problem.num_vertices; ++i) V.insert(i);

    //V/C
    for(std::set<unsigned int>::iterator it=C.vertices_in_solution.begin(); it!=C.vertices_in_solution.end(); ++it) V.erase(*it);

    //iz skupa V=V/C slucajno izaberemo vrhove Vm i to njih Pm*min(|C|, |V\C|), Pm je parametar odreden empirijski
    unsigned int n=static_cast<unsigned int>(pm*std::min(C.num_vertices, static_cast<unsigned int>(V.size())));

    //slucajno biramo n vrhova iz Vm
    std::set<unsigned int> Vm;
    while(Vm.size()!=n){
        unsigned int vs=random_from_set(V);
        while(problem.find_degree(vs)<=1){      //ne uzimamo vrhove koji su stupnja 1 u grafu jer cemo njih ionako brisati s cut
            vs=random_from_set(V);
        }
        Vm.insert(vs);
    }
    while(!Vm.empty()){
        //trazimo vs koji ima najkraci put do nekog random odabranog cvora iz vertices_in_solution
        std::pair<std::vector<int>, double> min_sp;
        double min_cost=INF;
        unsigned int vs_min;
        unsigned int vd=random_from_set(C.vertices_in_solution);
        for(std::set<unsigned int>::iterator it=Vm.begin(); it!=Vm.end(); ++it){
            //trazimo vrh iz Vm koji je "najblizi" vrhu vd
           std::pair<std::vector<int>, double> sp=problem.shorthest_path(*it, vd);
                if(sp.second<min_cost){
                    min_sp=sp;
                    min_cost=sp.second;
                    vs_min=*it;
                }
        }
        //dodajemo vrhove iz puta u C
        unsigned int a=vd;
        unsigned int b=static_cast<unsigned int>(min_sp.first[vd]);
        while(min_sp.first[a]!=-1){
            std::pair<unsigned int, double> p1(a, problem.adjacent(a,b));
            std::pair<unsigned int, double> p2(b, problem.adjacent(a,b));

            C.vertices_in_solution.insert(a);

            if(C.adj_list[b].insert(p1).second==true) { C.fitness+=p1.second; ++C.num_edges; }
            C.adj_list[a].insert(p2);

            a=static_cast<unsigned int>(min_sp.first[a]);
            b=static_cast<unsigned int>(min_sp.first[a]);
        }
        C.vertices_in_solution.insert(a);

        //micemo vs iz Vm;
        Vm.erase(vs_min);
    }

    C.num_vertices=static_cast<unsigned int>((C.vertices_in_solution).size());
    //kaze da se nakon toga primjeni prim
    C.MST(problem);
    return C;
}

#endif // SOLUTION_H
