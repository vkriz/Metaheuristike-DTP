#ifndef POPULATION_H
#define POPULATION_H

#include "solution.h"
#include <set>
#include <string>

extern unsigned int cross_cnt;
extern unsigned int mutate_cnt;

class Population{
public:
    unsigned int pop_size;
    std::set<Solution> solutions;
    Graph problem;

    //konstruktori:
    Population(Graph& G) : pop_size(0), solutions(), problem(G){}
    Population(unsigned int n, Graph& G, int f) : pop_size(n), solutions(), problem(G){        //stvara populaciju zadane velicine pomocu heuristickog algoritma
        while(solutions.size()!=n){
            Solution sol(problem, f);
            solutions.insert(sol);
            //u populaciji ce rjesenja biti poredana od najboljeg do najlosijeg
        }
    }

    //copy konstruktor:
    Population(const Population& p) : pop_size(p.pop_size), solutions(p.solutions){}

    //destruktor:
    ~Population(){}

    //operatori usporedivanja:
    friend bool operator==(const Population& first, const Population& second){
        if(first.pop_size!=second.pop_size) return false;
        if(first.solutions!=second.solutions) return false;
        return true;
    }

    friend bool operator!=(const Population& first, const Population& second){
        return !(first==second);
    }

    //operator pridruzivanja:
    Population& operator=(const Population& p){
        pop_size=p.pop_size;
        solutions=p.solutions;
        return *this;
    }

    //printanje u datoteku
    void print_population_in_file(std::string) const;

    //vraca najbolje rjesenje u populaciji
    Solution best(){ return *(solutions.begin()); }

    //bira 2 random rjesenja i vraca bolje s vjerojatnoscu p
    Solution select(double);

    //funkcija koja nareduje krizanje
    void crossover(double);

    //funkcija koja nareduje mutaciju
    void mutation(double, double);
};

void Population::print_population_in_file(std::string name) const{
    int i=1;
    for(std::set<Solution>::iterator it=solutions.begin(); it!=solutions.end(); ++it) {
       std::string str=std::to_string(i);
       (*it).print_in_file(name+"_solution_"+str,problem.num_vertices);
       ++i;
    }
}

Solution Population::select(double p){
    unsigned int rand1=random_num(0,pop_size-1);
    unsigned int rand2=random_num(0,pop_size-1);
    while(rand1==rand2) rand2=random_num(0,pop_size-1);
    unsigned int r=random_num(0,100);
    unsigned int selected;
    if(double(r)/100 <= p) selected=std::min(rand1,rand2);
    else selected=std::max(rand1,rand2);
    return *(std::next(solutions.begin(), static_cast<int>(selected)));
}

void Population::crossover(double p){
    Solution parent1=select(p);
    Solution parent2=select(p);
    while(parent1==parent2) parent2=select(p);     //biramo novog parent2 sve dok nisu razliciti
    Solution child=cross(parent1,parent2,problem);
    while(1){
        if(!child.cut(problem)) break;
    }
    child.MST(problem);
    while(1){
        if(!child.cut(problem)) break;
    }
    //ako je child razlicit od svih rjesenja, izbaci najgore rjesenje i ubaci child
    if(solutions.insert(child).second==true) solutions.erase(--solutions.end());
}

void Population::mutation(double p, double pm){
    Solution parent=select(p);
    Solution child=mutate(parent,problem, pm);
    while(1){
        if(!child.cut(problem)) break;
    }
    child.MST(problem);
    while(1){
        if(!child.cut(problem)) break;
    }
    //ako je child razlicit od svih rjesenja, izbaci najgore rjesenje i ubaci child
    if(solutions.insert(child).second==true) solutions.erase(--solutions.end());
}

#endif // POPULATION_H
