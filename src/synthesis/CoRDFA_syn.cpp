#include "CoRDFA_syn.h"

using std::make_unique;
using std::move;
using std::set;
using std::shared_ptr;
using std::string;
using std::unique_ptr;
using std::vector;

CoRDFA_syn::CoRDFA_syn(shared_ptr<Cudd> m, string filename, string filenonbackup, string partfile)
{
    std::cout << " A " << filename << std::endl;
    unique_ptr<SSNFA> ssnfa = make_unique<SSNFA>(m);
    ssnfa->initialize_mona(filename, partfile);
    mgr = m;
    std::cout << " B" << std::endl;
    ssnfa->project_unobservables();
    ssnfa->complement();
    std::cout << " C" << std::endl;
    initializer(ssnfa);
    std::cout << " D" << std::endl;
    unique_ptr<DFA> dfa = make_unique<DFA>(m);
    dfa->initialize(filenonbackup, partfile, false);

    unique_ptr<DFA> finaldfa = make_unique<DFA>(m);
    finaldfa->init_from_cross_product(ssnfa.get(), dfa.get()); 
    std::cout << "Cross has " << finaldfa->nstates << " states" << std::endl;
    bdd = move(ssnfa);
    
}

CoRDFA_syn::~CoRDFA_syn() {}

void CoRDFA_syn::initializer(unique_ptr<SSNFA>& ssnfa){
  for(int i = 0; i < ssnfa->nbits; i++){
    BDD b = mgr->bddVar();
    ssnfa->bddvars.push_back(b);
  }
  BDD t = ssnfa->finalstatesBDD;
  W.push_back(ssnfa->finalstatesBDD);
  Wprime.push_back(ssnfa->finalstatesBDD);
  cur = 0;

  for (int i = 0; i < ssnfa->nstates; ++i) {
    ssnfa->res.push_back(mgr->bddZero());
    for (int j = 0; j < ssnfa->nstates; ++j) {
          std::cout << "i: " << i << " " << ssnfa->nstates << std::endl;
          std::cout << "j: " << j << " " << ssnfa->nstates << std::endl;

      BDD trans = ssnfa->bddvars[j] & ssnfa->labels[i][j];
      ssnfa->res[i] |= trans;
    }
  }
}
















