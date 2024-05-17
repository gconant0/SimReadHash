#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "gen_dna_funcs.cpp"
#include "read_seq.cpp"
#include "write_seq.cpp"
#include "libranlib.h"
#include "random.cpp"
#include <unordered_map>
#include <upcxx/upcxx.hpp>
#include <sys/time.h>
#include <sys/resource.h>

unsigned long RunTime()
{
    unsigned long retval;
    struct timespec t1;
    
    clock_gettime(CLOCK_MONOTONIC,  &t1);
    
    retval = (t1.tv_sec*1e9 + t1.tv_nsec)/1000000;
    return(retval);
}



using namespace::std;

#define K_LEN 31

//Functions to extract kmer sequences
string get_kmer_string(long kmer);
char get_last_base(long kmer);

//Holds a single kmer
class Kmer {
public:
    long my_kmer, rev_kmer, left, right;
    Kmer(Molecule_Sequence *data, int pos);
    
};


//Holds graph node for kmer walk
class Kmer_node {
public:
    Kmer_node();
    int have_kmer, have_next_kmers[4], has_next, has_pred;
    long get_kmer() {return(my_kmer);};
    long get_rev_kmer() {return(rev_kmer);};
    long get_next_kmer() {return(next_kmer);};
    long get_pred_kmer() {return(pred_kmer);};
    long get_next_kmer(long i);                     //Tests all 4 possible next kmers
    long get_next_rev_kmer(long i);                 //Tests all 4 possible next reverse kmers
    void set_kmer(long kmer);
    void set_next(long n_kmer) {has_next=1; next_kmer=n_kmer;};   //Found a valid next kmer
    void set_prev(long p_kmer) {has_pred=1; pred_kmer=p_kmer;};   //Found a valid previous kmer
private:
    long two_mask, kmask;
    long my_kmer, rev_kmer, next_kmer, pred_kmer;
    
};



class DistrMap {
private:
  // store the local unordered map in a distributed object to access from RPCs (Remote Procedure Calls)
    using dobj_map_t = upcxx::dist_object<std::unordered_map<long, long> >;

  
public:
    // map the key to a target process
    int get_target_rank(const long &key) {
      return std::hash<long>{}(key) % upcxx::rank_n();
    }
    
    dobj_map_t local_map;
  // initialize the local map
  DistrMap() : local_map({}) {}
    
  // insert a key-value pair into the hash table
  upcxx::future<> insert(const long &key, const long &val) {
    // the RPC returns an empty upcxx::future by default
    return upcxx::rpc(get_target_rank(key),
                      // lambda to insert the key-value pair
                      [](dobj_map_t &lmap, const long &key, const long &val) {
                        // insert into the local map at the target
                        lmap->insert({key, val});
                      }, local_map, key, val);
  }
    upcxx::future<> increment(const long &key) {
      // the RPC returns an empty upcxx::future by default
      return upcxx::rpc(get_target_rank(key),
                        // lambda to insert the key-value pair
                        [](dobj_map_t &lmap, const long &key) {
                          // insert into the local map at the target
          auto elem = lmap->find(key);
          long new_val=elem->second+1;
          lmap->erase(key);
          lmap->insert({key, new_val});
                        }, local_map, key);
    }
    
  // find a key and return associated value in a future
  upcxx::future<long> find(const long &key) {
    return upcxx::rpc(get_target_rank(key),
                      // lambda to find the key in the local map
                      [](dobj_map_t &lmap, const long &key) -> long {
                        auto elem = lmap->find(key);
                        if (elem == lmap->end()) return -1; // not found
                        else return elem->second; // key found: return value
                      }, local_map, key);
  }
};

//Make the local nodes and the map to them global for rpc function access
Kmer_node *access_nodes;
std::map<long, int> node_ids;

void rpc_set_prev (const long my_kmer, const long p_kmer)
{
    auto index=node_ids.find(my_kmer);
    
    if (index == node_ids.end()) {
        cout<<"ERROR: No kmer: "<<my_kmer<<" on node "<<upcxx::rank_me()<<endl;
        return;
    }
    
    //cout<<"RPC: "<<upcxx::rank_me()<<" Kmer: "<<my_kmer<<" to "<<p_kmer<<" Kmer index="<<index->second<<endl;
    access_nodes[index->second].set_prev(p_kmer);
}
int rpc_next_valid(const long my_kmer)
{
     auto index=node_ids.find(my_kmer);
     if (index == node_ids.end()) {
         cout<<"ERROR: No kmer: "<<my_kmer<<" on node "<<upcxx::rank_me()<<endl;
         return(0);
     }
     return(access_nodes[index->second].has_next);
     
 }
long rpc_next_kmer(const long my_kmer)
{
     auto index=node_ids.find(my_kmer);
     if (index == node_ids.end()) {
         cout<<"ERROR: No kmer: "<<my_kmer<<" on node "<<upcxx::rank_me()<<endl;
         return(0);
     }
     return(access_nodes[index->second].get_next_kmer());
     
 }




int main (int argc, char **argv)
{
    upcxx::init();
    
    int i,j, cnt, genome_size, num_reads, read_len, r_start, kmers_per_read, next_valid, ntaxa, nchars, valid_local_kmers=0, num_found, found[8];
    long next_kmer, my_kmer, val, b, t_start, t_end;
    string genome_file, outfile, name, new_contig;
    Random_Gen *mygen;
    DATATYPE cdata=NUCLEIC;
	Read_Sequence *read_seq;
	Write_Sequence *write_seq=0;
	Sequence_dataset *current_data, *dataset;
    Kmer ***the_kmers;
    DistrMap kmer_count;
   
    
	if (argc<4) {
		cerr<<"Usage: sim_read_hash <in genome> <out file> num_reads read_len\n";
		return(-1);
	}
	else	{
        genome_file=argv[1];
        std::stringstream ss;
        ss<<argv[2]<<" "<<argv[3];
        ss>>num_reads>>read_len;
        
        
        mygen=new Random_Gen;
        read_seq=new Read_FASTA;
        t_start=RunTime();
        current_data=read_seq->get_dataset(ntaxa, nchars, genome_file.c_str(), FALSE);
        t_end=RunTime();
        cout<<upcxx::rank_me()<<": IO time: "<<t_end-t_start<<endl;
        
        cout<<"Simulating "<<num_reads<<" reads of length: "<<read_len<<" from sequence of size "<<nchars<<"\n";
        
   
        t_start=RunTime();
        dataset=new Sequence_dataset(num_reads/upcxx::rank_n(), read_len, cdata);
        cnt=0;
        
		for(i=0; i<num_reads; i++) {
             name = "SimRead_Num";
            
            stringstream ss;
            ss <<upcxx::rank_me()<<"_"<< cnt;
            name = name + ss.str();
    
            r_start=ignuin(0,(*current_data)[0].Sequence_size()-read_len);
           
            if ((i%upcxx::rank_n()) == upcxx::rank_me()) {
                (*dataset)[cnt].Assign_name(name.c_str());
                for(j=0; j<read_len; j++)
                    (*dataset)[cnt].Assign_site(j, (*current_data)[0][r_start+j]);
                cnt++;
            }
           
        }
        upcxx::barrier();
        t_end=RunTime();
        cout<<upcxx::rank_me()<<": Simulation time: "<<t_end-t_start<<endl;
        
        t_start=RunTime();
        kmers_per_read=read_len-K_LEN+1;
        
        the_kmers=new Kmer**[num_reads];
        
        for(i=0; i<num_reads/upcxx::rank_n(); i++) {
            the_kmers[i]=new Kmer*[kmers_per_read];
            
            for(j=0; j<read_len-K_LEN; j++)
                the_kmers[i][j]=new Kmer(&(*dataset)[i], j);
        }
        upcxx::barrier();
        t_end=RunTime();
        cout<<upcxx::rank_me()<<": Kmer build time: "<<t_end-t_start<<endl;
        
        
        //Construct distributed hash table on the kmers
        t_start=RunTime();
        for(i=0; i<num_reads/upcxx::rank_n(); i++) {
            for(j=0; j<read_len-K_LEN; j++) {
                val =kmer_count.find(the_kmers[i][j]->my_kmer).wait();
                if (val != -1 )
                    kmer_count.increment(the_kmers[i][j]->my_kmer).wait();
                else
                    kmer_count.insert(the_kmers[i][j]->my_kmer, 1).wait();
            }
        }
        
        upcxx::barrier();
        t_end=RunTime();
        cout<<upcxx::rank_me()<<": Kmer hash table time: "<<t_end-t_start<<endl;
        
        //Walk the local hash table looking for kmers seen more than once
        cout <<upcxx::rank_me()<<" : "<<kmer_count.local_map->size()<<endl;
        t_start=RunTime();
        for (auto it=kmer_count.local_map->begin();  it != kmer_count.local_map->end(); ++it) {
            if (it->second >1) valid_local_kmers++;
        }
        
        //Local graph for holding the kmers and their next/last kmers
        access_nodes=new Kmer_node[valid_local_kmers];
        
        cnt=0;
        for (auto it=kmer_count.local_map->begin();  it != kmer_count.local_map->end(); ++it) {
            if (it->second >1) {
                access_nodes[cnt].set_kmer(it->first);
                for(i=0; i<8; i++) found[i]=0;
                
                //Check all four possible next kmers for validity in the hash table
                for(b=0; b<4; b++) {
                    next_kmer=access_nodes[cnt].get_next_kmer(b);
                    
                    const auto nextit = kmer_count.find(next_kmer).wait();
                    if (nextit >1) found[(int)b] =1;
                }
                
                //REVERSE CODE CURRENTLY OFF
#ifdef USE_REV
                for(b=0; b<4; b++) {
                    next_kmer=access_nodes[cnt].get_next_rev_kmer(b);
                    
                    const auto nextit = kmer_count.find(next_kmer).wait();
                    if (nextit >1) found[(int)b+4] =1;
                }
#endif
                //Check if exactly one next kmer is found validly in the hash
                num_found=0;
                for(i=0; i<8; i++) {
                    if (found[i] ==1) num_found++;
                }
                
                if (num_found==1) {
                    for(b=0; b<4; b++) {
                        if (found[(int)b] ==1) access_nodes[cnt].set_next(access_nodes[cnt].get_next_kmer(b));
                    }
#ifdef USE_REV
                    for(b=0; b<4; b++) {
                        if (found[(int)b+4] ==1)  access_nodes[cnt].set_next(access_nodes[cnt].get_next_rev_kmer(b));
                    }
#endif
                }
                
                //Using c++ map to relate local kmers to node array
                node_ids[access_nodes[cnt].get_kmer()]=cnt;
                cnt++;
            }
        }
        
        upcxx::barrier();
        t_end=RunTime();
        cout<<upcxx::rank_me()<<": Kmer next computation: "<<t_end-t_start<<endl;
        
        
        
        num_found=0;
        for(i=0; i<valid_local_kmers; i++) {
            if (access_nodes[i].has_next==1) num_found++;
        }
        cout<<upcxx::rank_me()<<": "<<num_found<<" of "<<valid_local_kmers<<" have a next"<<endl;
        
        //Reverse walk the graph looking for the precessor of each kmer
        t_start=RunTime();
        for(i=0; i<valid_local_kmers; i++) {
            if (access_nodes[i].has_next==1) {
                next_kmer=access_nodes[i].get_next_kmer();
                my_kmer=access_nodes[i].get_kmer();
                
                if (kmer_count.get_target_rank(next_kmer) ==upcxx::rank_me())
                    rpc_set_prev(next_kmer, my_kmer);
                else
                    upcxx::rpc(kmer_count.get_target_rank(next_kmer), rpc_set_prev, next_kmer, my_kmer).wait();
            }
        }

        upcxx::barrier();
        if (upcxx::rank_me()==0) {
            cout<<"Assigned previous kmers"<<endl;
        }
        upcxx::barrier();
        t_end=RunTime();
        cout<<upcxx::rank_me()<<": Kmer previous computation time: "<<t_end-t_start<<endl;
        
        //Walk the graph looking for kmers with no previous kmers: start contruction
        t_start=RunTime();
        for(i=0; i<valid_local_kmers; i++) {
            if ((access_nodes[i].has_pred == 0) && (access_nodes[i].has_next == 1)) {
                new_contig=get_kmer_string(access_nodes[i].get_kmer());
                
                next_kmer=access_nodes[i].get_next_kmer();
                next_valid=access_nodes[i].has_next;
                while(next_valid !=0) {
                    new_contig=new_contig+get_last_base(next_kmer);
                    
                    if(kmer_count.get_target_rank(next_kmer) == upcxx::rank_me()) {
                        next_valid=access_nodes[node_ids[next_kmer]].has_next;
                        next_kmer=access_nodes[node_ids[next_kmer]].get_next_kmer();
                    }
                    else {
                        upcxx::future<int> has_result = upcxx::rpc(kmer_count.get_target_rank(next_kmer), rpc_next_valid, next_kmer);
                        upcxx::future<long> result_is = upcxx::rpc(kmer_count.get_target_rank(next_kmer), rpc_next_kmer, next_kmer);
                        
                        next_valid = has_result.wait();
                        next_kmer = result_is.wait();
                    }
                    
                }
                if (new_contig.length()>40)
                    cout<<"Contig: ("<<new_contig.length()<<"): "<<new_contig<<endl;
            }
        }
        t_end=RunTime();
        cout<<upcxx::rank_me()<<": Kmer walk time: "<<t_end-t_start<<endl;
        
        //write_seq=new Write_FASTA(outfile.c_str(), cdata);
			
        //write_seq->write_dataset(2*num_reads, read_len, dataset);
			
	
		
        delete mygen;
		
		delete current_data;
        delete dataset;
        delete write_seq;
        
        upcxx::barrier(); // wait for finds to complete globally
        if (!upcxx::rank_me()) cout << "SUCCESS" << endl;
        upcxx::finalize();
        
		return(0);
    }
}//end main






Kmer::Kmer(Molecule_Sequence *data, int pos)
//Extracts the kmer at position pos from the read
//Left and right are computed but not used
//Reverse is computed but not used
{
    int i, two_mask=3, new_base;
    long temp, rev_left, rev_right, kmask;
    
    rev_kmer=0;
    my_kmer=0;
    
    
    for(i=pos; i<pos+K_LEN; i++) {
        my_kmer = my_kmer | ((*data)[i] & two_mask);
        my_kmer = my_kmer <<2;
        new_base=(*data)[i] ^ two_mask;
        rev_kmer= rev_kmer | (new_base << ((2*i)-pos));
        
    }
    my_kmer = my_kmer >>2;

    
    kmask=0;
    for(i=0; i<K_LEN; i++) {
        kmask=kmask | two_mask;
        kmask =kmask<<2;
    }
    kmask =kmask>>2;
    
    if (pos !=0) {
        left=my_kmer;
        left=left>>2;
        left = left | (((*data)[pos-1] & two_mask) << (2*(K_LEN-1)));
        
        rev_left=rev_kmer;
        rev_left = (rev_left <<2) & kmask;
        
        new_base=(*data)[pos-1] ^ two_mask;
        rev_left = rev_left | new_base;
       
        if (left >rev_left) {
            temp=left;
            left=rev_left;
            rev_left=temp;
            
        }
        
    }
    if (pos < data->Sequence_size()-1) {
        right=my_kmer;
        right = (right << 2) & kmask;
        right = right | ((*data)[pos+K_LEN] & two_mask);
        
        rev_right=rev_kmer;
        rev_right =rev_right>>2;
        new_base=(*data)[pos+K_LEN] ^ two_mask;
        rev_right = rev_right | (new_base << (2*(K_LEN-1)));
        
       
        if (rev_right < right) {
            temp=right;
            right=rev_right;
            rev_right=temp;
        }
        
    }
#ifdef USE_REV
    if (rev_kmer<my_kmer) {
        temp=my_kmer;
        my_kmer=rev_kmer;
        rev_kmer=temp;
    }
#endif
}


string get_kmer_string(long kmer)
{
    int i, two_mask=3, val;
    
    string ret_val;
    
    ret_val="";
    for(i=0; i<K_LEN; i++) {
        val=kmer & two_mask;
        ret_val = num_to_base(val) + ret_val;
        kmer =kmer >> 2;
    }
    return(ret_val);
}

char get_last_base(long kmer)
{
    int two_mask=3;
    return(num_to_base(kmer & two_mask));
}

Kmer_node::Kmer_node ()
{
    int i;
    my_kmer=-1;
    have_kmer=0;
    pred_kmer=0;
    has_pred=0;
    has_next=0;
    next_kmer=0;
    for(i=0; i<4; i++) {
        have_next_kmers[i]=0;
    }
    two_mask=3;
    kmask=0;
    for(i=0; i<K_LEN; i++) {
        kmask=kmask | two_mask;
        kmask =kmask<<2;
    }
    kmask =kmask>>2;
}

void Kmer_node::set_kmer(long kmer)
{
    int i, j,  val1, val2;
    long part_kmer, next_base, orig_base;
    my_kmer=kmer;
    rev_kmer=0;
    
    part_kmer=my_kmer;
    
    for(i=0; i<K_LEN; i++) {
        orig_base = part_kmer & two_mask;
        next_base = orig_base ^ two_mask;
        next_base = (next_base << (2*(K_LEN-i-1)));
        part_kmer =part_kmer >>2;
        rev_kmer= rev_kmer | next_base;

        
    }
    
}

long Kmer_node::get_next_kmer(long i)
{
    long ret_val;
    ret_val=my_kmer;
    ret_val = (ret_val << 2) & kmask;
    ret_val = ret_val | (i & two_mask);
    return(ret_val);
}


long Kmer_node::get_next_rev_kmer(long i)
{
    long ret_val, next_base;
     
    ret_val=rev_kmer;
    
    ret_val = ret_val >> 2;
    
    next_base=(i ^ two_mask) << (2*(K_LEN-1));
    ret_val = ret_val | next_base;
    
    return(ret_val);
}
