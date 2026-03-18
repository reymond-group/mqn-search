#include "reader.h"

#include <algorithm>
#include <math.h>
#include <set>
#include <unordered_map>
#include <iostream>

// Atom type/covalency/group information
#define s_2_val -22
#define o_2_val -19
#define n_3_val_pos -16
#define n_3_val -13
#define n_2_val -12
#define c_123_val -1
#define C_1234_val 1
#define C_2_bnd_amide_hba 2
#define C_2_bnd_amide_hbd 3
#define C_2_bnd_guanidine 4
#define C_2_bnd_guanidine_hba 5
#define C_2_bnd_guanidine_hbd 6
#define C_2_bnd_arylamidine 7
#define C_2_bnd_carboxy 8
#define C_2_bnd_ester 9
#define C_2_bnd_keto 10
#define C_3_bnd 11
#define C_2_bnd_thiamide 12
#define N_12_val 13
#define N_3_val 14
#define N_2_bnd 15
#define N_4_val 17
#define N_3_val_pos 18
#define N_3_bnd 19
#define O_12_val 20
#define O_2_bnd 21
#define O_1_val_neg 22
#define S_12_val 23
#define S_2_bnd 24
#define S_2x2_bnd_6_val 25
#define S_6_val 26 // e.g. pentafluorosulfanyl-R there are some in ChEMBL
#define P_3_val 27
#define P_1x2_bnd_5_val 28
#define P_5_val 29

typedef unsigned short ushort;
typedef unsigned int uint;

const std::map<char, std::unordered_set<char>> reader::cnos =
{
    //! Removed pair 'c' - 'S' (there are compounds like
    //! Vortioxetine where S is followed by an aromatic carbon...
    //! update: also removed labels 's' - 'C' and 'o' - 'C' for the same reason.
    //! and since it apparently IS possible 'n' - 'S'
    //! -> This would be relevant for elements Sc, Sn, Co
    {'c', std::unordered_set<char>({'A','T'})},
    {'n', std::unordered_set<char>({'M', 'Z', 'I', 'R'})},
    {'o', std::unordered_set<char>({'M', 'N', 'P', 'H'})},
    {'s', std::unordered_set<char>({'A', 'O', 'E'})},
    {'p', std::unordered_set<char>()}
};

const std::map<std::string, uint> reader::amass =
{
    {"H", 1},
    {"B", 11},
    {"C", 12},
    {"c", 12},
    {"N", 14},
    {"n", 14},
    {"O", 16},
    {"o", 16},
    {"F", 19},
    {"Mg", 24},
    {"Al", 27},
    {"Si", 28},
    {"P", 31},
    {"p", 31},
    {"S", 32},
    {"s", 32},
    {"Cl", 35},
    {"Ca", 40},
    {"Fe", 56},
    {"Cu", 64},
    {"Zn", 65},
    {"As", 75},
    {"Se", 79},
    {"Br", 80},
    {"I", 127}
};
const std::map<std::string, uint> reader::anum =
{
    {"H", 1},
    {"B", 5},
    {"C", 6},
    {"c", 6},
    {"N", 7},
    {"n", 7},
    {"O", 8},
    {"o", 8},
    {"F", 9},
    {"Mg", 12},
    {"Al", 13},
    {"Si", 14},
    {"P", 15},
    {"p", 15},
    {"S", 16},
    {"s", 16},
    {"Cl", 17},
    {"Ca", 20},
    {"Fe", 26},
    {"Cu", 29},
    {"Zn", 30},
    {"As", 33},
    {"Se", 34},
    {"Br", 35},
    {"I", 53}
};


reader::reader() :
m_cyclic_edge_count(0),
m_edges(0),
m_bonds_in_min_2_rings(0),
m_atoms_in_min_2_rings(0),
m_cyclic_single_bonds(0),
m_acyclic_single_bonds(0),
m_cyclic_double_bonds(0),
m_acyclic_double_bonds(0),
m_cyclic_triple_bonds(0),
m_acyclic_triple_bonds(0),
m_monovalent(0),
m_cyclic_divalent(0),
m_acyclic_divalent(0),
m_cyclic_trivalent(0),
m_acyclic_trivalent(0),
m_cyclic_tetravalent(0),
m_acyclic_tetravalent(0),
m_pos_charges(0),
m_neg_charges(0),
m_hba(0),
m_hba_sites(0),
m_hbd(0),
m_hbd_sites(0),
m_triple_shared(0)
{

}

reader::~reader(){

}

// ### PUBLIC METHODS ###
void reader::read( const std::string SMILES ){
    m_SMILES = SMILES;

    gen_adj_list();
    count_implicit_H();
    correct_valence();
    assign_mass_number();

    find_rings();

    gen_bulk_ring_adj_lst();
    ring_systems();
    sssr();
    gen_ring_adj_lst();

    spheroid_system_detection();

    single_bonds();
    double_bonds();
    triple_bonds();

    m_atom_types.resize( m_atom_labels.size(), 0 );

    // concretize atomic environment
    init_pattern();

    find_fix_amides();
    find_fix_carboxy();
    find_fix_guanidines_amidines(); // The returned bool is for distribution m_accepted_patterns
    find_carbonyl();
    find_fix_S_P();
/*
    if( has_NCOO() ||
        has_terminal_carbonate_carbamate() ||
        has_anhydride_like() ||
        bad_planarization() )
    {
        uint test = 0;
    }
*/
    node_order();
    charge();
    h_bond();
}

const std::vector<uint>& reader::get_mqn(){
    if( m_mqn.empty() ){
        m_mqn.reserve( 43 );
        m_mqn.emplace_back( 0 ); // sum; temporary
        m_mqn.emplace_back( m_cyclic_single_bonds );
        m_mqn.emplace_back( m_acyclic_single_bonds );
        m_mqn.emplace_back( m_cyclic_divalent );
        m_mqn.emplace_back( m_atoms_in_min_2_rings );
        m_mqn.emplace_back( m_hbd_sites );
        m_mqn.emplace_back( m_acyclic_divalent );
        m_mqn.emplace_back( m_cyclic_trivalent );
        m_mqn.emplace_back( m_hba_sites );
        m_mqn.emplace_back( rotatable_bonds() );
        m_mqn.emplace_back( m_bonds_in_min_2_rings );
        m_mqn.emplace_back( carbon_count() );
        m_mqn.emplace_back( m_hba );
        m_mqn.emplace_back( m_monovalent );
        m_mqn.emplace_back( m_hbd );
        m_mqn.emplace_back( cyclic_n_count() );
        m_mqn.emplace_back( m_acyclic_double_bonds );
        m_mqn.emplace_back( acyclic_n_count() );
        m_mqn.emplace_back( m_cyclic_double_bonds );
        m_mqn.emplace_back( acyclic_o_count() );
        m_mqn.emplace_back( count_size_5_rings() );
        m_mqn.emplace_back( m_acyclic_trivalent );
        m_mqn.emplace_back( m_cyclic_tetravalent );
        m_mqn.emplace_back( cyclic_o_count() );
        m_mqn.emplace_back( m_pos_charges );
        m_mqn.emplace_back( count_size_3_rings() );
        m_mqn.emplace_back( count_size_6_rings() );
        m_mqn.emplace_back( m_acyclic_triple_bonds );
        m_mqn.emplace_back( count_size_4_rings() );
        m_mqn.emplace_back( heavy_atom_count() );
        m_mqn.emplace_back( count_size_7_rings() );
        m_mqn.emplace_back( sulfur_count() );
        m_mqn.emplace_back( count_size_8_rings() );
        m_mqn.emplace_back( m_acyclic_tetravalent );
        m_mqn.emplace_back( m_neg_charges );
        m_mqn.emplace_back( count_size_9_rings() );
        m_mqn.emplace_back( m_cyclic_triple_bonds );
        m_mqn.emplace_back( count_macrocycles() );
        m_mqn.emplace_back( chlorine_count() );
        m_mqn.emplace_back( iodine_count() );
        m_mqn.emplace_back( bromine_count() );
        m_mqn.emplace_back( fluorine_count() );
        m_mqn.emplace_back( phosphorous_count() );

        // update sum
        for(uint i = 1; i < m_mqn.size(); ++ i){
            m_mqn[0] += m_mqn[i];
        }
    }
    return m_mqn;
}

uint reader::count_size_3_rings(){
    uint count = 0;
    for(const auto* ring: m_topol_SSSR ){
        if( ring->size() == 3 ){
            ++count;
        }
    }
    return count;
}
uint reader::count_size_4_rings(){
    uint count = 0;
    for(const auto* ring: m_topol_SSSR ){
        if( ring->size() == 4 ){
            ++count;
        }
    }
    return count;
}
uint reader::count_size_5_rings(){
    uint count = 0;
    for(const auto* ring: m_topol_SSSR ){
        if( ring->size() == 5 ){
            ++count;
        }
    }
    return count;
}
uint reader::count_size_6_rings(){
    uint count = 0;
    for(const auto* ring: m_topol_SSSR ){
        if( ring->size() == 6 ){
            ++count;
        }
    }
    return count;
}
uint reader::count_size_7_rings(){
    uint count = 0;
    for(const auto* ring: m_topol_SSSR ){
        if( ring->size() == 7 ){
            ++count;
        }
    }
    return count;
}
uint reader::count_size_8_rings(){
    uint count = 0;
    for(const auto* ring: m_topol_SSSR ){
        if( ring->size() == 8 ){
            ++count;
        }
    }
    return count;
}
uint reader::count_size_9_rings(){
    uint count = 0;
    for(const auto* ring: m_topol_SSSR ){
        if( ring->size() == 9 ){
            ++count;
        }
    }
    return count;
}
uint reader::count_macrocycles(){
    uint count = 0;
    for(const auto* ring: m_topol_SSSR ){
        if( ring->size() >= 10 ){
            ++count;
        }
    }
    return count;
}

uint reader::ring_system_count(){
    return m_ring_systems.size();
}
uint reader::get_sssr_size(){
    return m_topol_SSSR.size();
}

uint reader::count_monovalent(){
    uint count = 0;
    for(const auto& nbrs : m_adj_list ){
        if( nbrs.size() == 1 ){
            ++count;
        }
    }
    return count;
}
uint reader::count_divalent(){
    uint count = 0;
    for(const auto& nbrs : m_adj_list ){
        if( nbrs.size() == 2 ){
            ++count;
        }
    }
    return count;
}
uint reader::count_trivalent(){
    uint count = 0;
    for(const auto& nbrs : m_adj_list ){
        if( nbrs.size() == 3 ){
            ++count;
        }
    }
    return count;
}
uint reader::count_tetravalent(){
    uint count = 0;
    for(const auto& nbrs : m_adj_list ){
        if( nbrs.size() == 4 ){
            ++count;
        }
    }
    return count;
}

uint reader::cyclic_node_count(){
    return m_cyclic_nodes.size();
}
uint reader::cyclic_edge_count(){
    return m_cyclic_edge_count;
}


float reader::fraction_cyclic_edges(){
    float total_edge_count = m_edges;
    return m_cyclic_edge_count / total_edge_count;
}
uint reader::bonds_in_min_2_rings(){
    return m_bonds_in_min_2_rings;
}
uint reader::atoms_in_min_2_rings(){
    return m_atoms_in_min_2_rings;
}

float reader::molecular_weight(){
    uint i = 0;
    float mwt = 0;
    for(uint mass : m_atomic_masses){
        /*
        if( mass == 12 ){
            mwt += 12.01;
        }
        */
        mwt += mass;

        mwt += m_attached_H[i] * 1.008;
        ++i;
    }
    return mwt;
}
uint reader::heavy_atom_count(){
    return m_adj_list.size();
}
uint reader::rotatable_bonds(){
    uint bonds = 0;
    for(const auto& nbrs: m_adj_list){
        bonds += nbrs.size();
    }
    bonds /= 2;
    m_edges = bonds;
    bonds -= m_cyclic_edge_count;

    for( const auto& nbrs: m_adj_list ){
        if( nbrs.size() == 1 ){
            -- bonds;
        }
    }

    return bonds;
}

uint reader::heteroatom_count(){
    uint hetatm = 0;
    for(const std::string& label : m_atom_labels){
        if( label != "C" && label != "c"){
            ++hetatm;
        }
    }
    return hetatm;
}

float reader::fraction_aromatic(){
    float f_ar = 0.0;
    for(const std::string& label : m_atom_labels){
        if(islower(label.front())){
            ++f_ar;
        }
    }
    f_ar /= m_adj_list.size();
    return f_ar;
}

uint reader::carbon_count(){
    uint count_c = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 6 ){
            ++ count_c;
        }
    }
    return count_c;
}
uint reader::nitrogen_count(){
    uint count_n = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 7 ){
            ++ count_n;
        }
    }
    return count_n;
}
uint reader::oxygen_count(){
    uint count_o = 8;
    for( uint anum : m_atomic_numbers ){
        if( anum == 16 ){
            ++ count_o;
        }
    }
    return count_o;
}
uint reader::sulfur_count(){
    uint count_s = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 16 ){
            ++ count_s;
        }
    }
    return count_s;
}
uint reader::fluorine_count(){
    uint count_f = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 9 ){
            ++ count_f;
        }
    }
    return count_f;
}
uint reader::chlorine_count(){
    uint count_cl = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 17 ){
            ++ count_cl;
        }
    }
    return count_cl;
}
uint reader::bromine_count(){
    uint count_br = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 35 ){
            ++ count_br;
        }
    }
    return count_br;
}
uint reader::phosphorous_count(){
    uint count_p = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 15 ){
            ++ count_p;
        }
    }
    return count_p;
}
uint reader::iodine_count(){
    uint count_i = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 53 ){
            ++ count_i;
        }
    }
    return count_i;
}
uint reader::acyclic_n_count(){
    uint count_acn = 0;
    uint i = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 7 && m_cyclic_nodes.find(i) == m_cyclic_nodes.end() /*m_ring_membership[i] == 0*/ ){
            ++ count_acn;
        }
        ++ i;
    }
    return count_acn;
}
uint reader::cyclic_n_count(){
    uint count_cn = 0;
    uint i = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 7 && m_cyclic_nodes.find(i) != m_cyclic_nodes.end() /*m_ring_membership[i] == 1*/ ){
            ++ count_cn;
        }
        ++ i;
    }
    return count_cn;
}
uint reader::acyclic_o_count(){
    uint count_aco = 0;
    uint i = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 8 && m_cyclic_nodes.find(i) == m_cyclic_nodes.end() /*m_ring_membership[i] == 0*/ ){
            ++ count_aco;
        }
        ++ i;
    }
    return count_aco;
}
uint reader::cyclic_o_count(){
    uint count_co = 0;
    uint i = 0;
    for( uint anum : m_atomic_numbers ){
        if( anum == 8 && m_cyclic_nodes.find(i) != m_cyclic_nodes.end() /*m_ring_membership[i] == 1*/ ){
            ++ count_co;
        }
        ++ i;
    }
    return count_co;
}

void reader::cleanup(){
    m_SMILES.clear();
    m_adj_list.clear();
    m_bond_orders.clear();
    m_atom_labels.clear();
    m_atom_types.clear();
    m_valence.clear();
    m_charge.clear();
    m_attached_H.clear();
    m_atomic_masses.clear();
    m_atomic_numbers.clear();
    m_ring_membership.clear();
    m_nestedness.clear();
    m_aromatic_H.clear();

    m_orig = 0;
    m_temp_path.clear();
    m_rings.clear();
    m_bulk_ring_adj_lst.clear();
    m_ring_sys_membership.clear();
    m_sys_size_ring.clear();
    m_sys_size_ring_23.clear();
    m_topol_SSSR.clear();
    m_symm_SSSR.clear();
    m_symm_SSSR_23.clear();
    m_ring_adj_lst.clear();
    m_shared_nodes_lst.clear();
    m_bridging_nodes.clear();
    m_ring_systems.clear();

    //! Under construction
    m_largest_rings.clear();
    //! End under construction

    m_planarizable.clear();
    m_planarizable_nodes.clear();
    m_node_ring_membership.clear();

    m_cyclic_nodes.clear();
    m_cyclic_edge_count = 0;
    m_bonds_in_min_2_rings = 0;
    m_atoms_in_min_2_rings = 0;

    m_cyclic_single_bonds = 0;
    m_acyclic_single_bonds = 0;
    m_cyclic_double_bonds = 0;
    m_acyclic_double_bonds = 0;
    m_cyclic_triple_bonds = 0;
    m_acyclic_triple_bonds = 0;

    m_monovalent = 0;
    m_cyclic_divalent = 0;
    m_acyclic_divalent = 0;
    m_cyclic_trivalent = 0;
    m_acyclic_trivalent = 0;
    m_cyclic_tetravalent = 0;
    m_acyclic_tetravalent = 0;

    m_pos_charges = 0;
    m_neg_charges = 0;

    m_hba = 0;
    m_hba_sites = 0;
    m_hbd = 0;
    m_hbd_sites = 0;

    m_triple_shared = 0;
    m_mqn.clear();
}

// ### PRIVATE METHODS ###
bool reader::is_valid_atom_id(const uint i){
    char current = m_SMILES.at(i);
    char last = m_SMILES.at(i-1);

    auto forbidden_chars = cnos.find(current);
    if( forbidden_chars == cnos.end() ){
        return false;
    }
    else if( forbidden_chars->second.find( last ) != forbidden_chars->second.end() ){
        return false;
    }
    return true;
}

void reader::get_nestedness(){

    m_nestedness.resize(m_SMILES.size());

    uint nestedness = 0;
    int max_nestedness = 0;
    for(uint i = 0; i < m_nestedness.size(); ++i ){
        if( m_SMILES.at(i) == '(' ){

            ++ nestedness;

            m_nestedness.at(i) = nestedness;

            if(nestedness > max_nestedness){
                max_nestedness = nestedness;
            }
        }
        else if( m_SMILES.at(i) == ')' ){

            m_nestedness.at(i) = nestedness;

            -- nestedness;
        }
        else{
            m_nestedness.at(i) = nestedness;
        }
    }
}

bool reader::look_ahead(const uint i, std::vector<std::vector<uint>>& connectivity, std::vector<uint>& rid_count ){

    uint ori = i - 1;
    uint last_atom_id = ori; // was set to i which in case of e.g. '(' led to wrongly setting last atom id to that position, hence had no obvious effect
    bool is_two_letter_id = false;
    bool sense_ring_start_ids = true;
    bool sense_ring_stop_ids = false;
    int parentheses = 0;
    uint neighbors_detected = 0;
    uint ring_closures = 0;
    std::set<char> ring_ids;
    bool search_direct_nbr = true;
    bool next_valid_is_nbr = false;
    bool got_in_line_nbr = false;
    bool got_dbond = false;
    bool got_tbond = false;
    bool got_pos_chrg = false;
    bool got_neg_chrg = false;
    bool poss_pos_end = false;
    bool poss_neg_end = false;
    bool in_brackets = false;
    bool search_stop_before_bracket = false;
    bool skip_multi_bond = false;

    bool left_sidechain = false;


    for( uint j = i; j < m_SMILES.size(); ++ j ){

        char chr = m_SMILES.at(j);
        if( chr == '[' ){
            if( j == i || ( j==1 && is_two_letter_id) ){
                got_in_line_nbr = true;
            }
            if(sense_ring_stop_ids){
                search_stop_before_bracket = true;
            }
            sense_ring_stop_ids = false;

            in_brackets = true;
        }
        //! 05.01.2023: commenting out - below seems to be justified.
        //! 16.1.2023: absolutely. the only possible problem is it crashes when not done.
        //! I guess this is for the in brackes + multiple explicit H case where you temporarly ignore numbers
        //! to not confuse them with ring ids.
        else if( chr == 'H' ){
            //sense_ring_start_ids = false;
            sense_ring_stop_ids = false;
            m_aromatic_H.insert(last_atom_id);
        }
        else if( chr == ']' ){
            if(search_stop_before_bracket){
                sense_ring_stop_ids = true;
                search_stop_before_bracket = false;
            }
            //! 05.09.2023: Moved below conditional up here. Looks a bit counter-intuitive
            //! to to it when in_brackets == false. It makes sense though:
            //! Brackets are sensed strating from every atom to their left X...[...]
            //! + & - are always sensed (also if - connects rings). However, the charge is only set
            //! when
            //! - charge was detected,
            //! - hitting a ] before [ which only happens when starting from an atom inside brackets.
            //! inside_brackets is only set true upon meeting a [ first and if true nothing happens
            //! and it's again turned off below.
            if( !in_brackets ){
                if( got_pos_chrg ){
                    ++m_charge[ori];
                    got_pos_chrg = false;
                }
                else if( got_neg_chrg ){
                    --m_charge[ori];
                    got_neg_chrg = false;
                }
            }
            in_brackets = false;
        }
        else if( isdigit(chr) ){
            if( sense_ring_start_ids ){
                if(rid_count.at(j) % 2 == 0){
                    ring_ids.insert( chr );
                    sense_ring_stop_ids = true;
                }
                ++ rid_count.at(j);
                if( got_dbond || got_tbond ){
                    skip_multi_bond = true;
                }
            }
            else if( sense_ring_stop_ids && ring_ids.find( chr ) != ring_ids.end() ){
                connectivity[ ori ].push_back( last_atom_id );
                connectivity[ last_atom_id ].push_back( ori );

                //! 13.09.2023: Added below check to account for ugly but valid SMILES like:
                //! C=1NCCCC1. Also assured that such double bond is not set twice
                //! (between first C and N, see downstream)
                if( skip_multi_bond ){
                    if( got_dbond ){
                        m_bond_orders[ std::pair<uint, uint>( ori, last_atom_id ) ] = 2;
                        got_dbond = false;
                    }
                    else if( got_tbond ){
                        m_bond_orders[ std::pair<uint, uint>( ori, last_atom_id ) ] = 3;
                        got_tbond = false;
                    }
                    skip_multi_bond = false;
                }
                else{
                    m_bond_orders[ std::pair<uint, uint>( ori, last_atom_id ) ] = 1;
                }

                m_ring_membership[ ori ] = 1;
                m_ring_membership[ last_atom_id ] = 1;

                // unset all ring membership flags of "higher nestedness" than at closure (i.e. count open parentheses)
                // that may have been set during previous ring closure search
                // e.g. unset C4 (and the rest in parentheses) in C1CC(Cc2ccccc2)CCC1 that was set
                // cyclic while looking for first closure (at the and).
                // The rest in parentheses will be reset when searching closure for 2.
                // Note: it must be degree of closure as it can be higher than at opening.
                // - Only doing it when start and closure are nested equally detects (cc2) in:
                // c2cc(ccc2)CCCCCCCCC
                // yet flags (C) cyclic in:
                // c2cc(c(C)cc2)CCCCCCCCC
                for( uint k = ori+1; k < last_atom_id; ++k ){

                    if( m_nestedness[k] >= m_nestedness[last_atom_id]){
                        m_ring_membership[k] = 0;
                    }
                }

                ring_ids.erase(chr);
                if(ring_ids.empty()){
                    sense_ring_stop_ids = false;
                }
                ++ rid_count[j];
                ++ neighbors_detected;
            }
        }
        else if( parentheses == 0 && chr == '(' ){
            ++ parentheses;
            if( search_direct_nbr ){
                next_valid_is_nbr = true;
            }
            sense_ring_start_ids = false;
        }
        else if( chr == '(' && parentheses < 0 ){
            ++ parentheses;
            next_valid_is_nbr = false;
            sense_ring_start_ids = false;
        }
        else if( chr == '(' ){
            ++ parentheses;
            sense_ring_start_ids = false;
        }
        else if( parentheses == 0 && chr == ')' ){
            if( ring_ids.empty() ){
                return is_two_letter_id;
            }
            left_sidechain = true;
            -- parentheses;
        }
        else if( chr == ')' ){
            -- parentheses;
            sense_ring_start_ids = false;
            if( parentheses == 0 && !got_in_line_nbr ){
                next_valid_is_nbr = true;
            }
        }
        else if( isupper(chr) || ( isalpha(chr) && is_valid_atom_id( j )) ){

            if( !left_sidechain && (((parentheses == 0 && search_direct_nbr ) ||
                (next_valid_is_nbr && (parentheses==0 || parentheses==1 )))) )
            {
                if( parentheses == 0 ){
                    got_in_line_nbr = true;
                }
                connectivity[ ori ].push_back( j );
                connectivity[ j ].push_back( ori );

                ++ neighbors_detected;
                search_direct_nbr = false;
                next_valid_is_nbr = false;

                if( j == last_atom_id || (j == 1 && is_two_letter_id ) ){
                    got_in_line_nbr = true;
                }
                if( !skip_multi_bond ){
                    if(got_dbond){
                        ++m_valence.at(ori);
                        ++m_valence.at(j);
                        m_bond_orders[ std::pair<uint, uint>(ori, j) ] = 2;
                        got_dbond = false;
                    }
                    else if(got_tbond){
                        m_valence.at(ori) += 2;
                        m_valence.at(j) += 2;
                        m_bond_orders[ std::pair<uint, uint>(ori, j) ] = 3;
                        got_tbond = false;
                    }
                    else{
                        m_bond_orders[ std::pair<uint, uint>(ori, j) ] = 1;
                    }
                }
                else{
                    m_bond_orders[ std::pair<uint, uint>(ori, j) ] = 1;
                }
            }
            if(sense_ring_stop_ids){
                m_ring_membership.at(j) = 1;
            }
            sense_ring_start_ids = false;
            last_atom_id = j;
        }
        else if( j == i && islower(chr) ){
            is_two_letter_id = true;
        }

        else if( chr == '=' && ((next_valid_is_nbr && parentheses >= 1 ) || (!got_in_line_nbr && parentheses == 0 ))){
            got_dbond = true;
        }
        //! 05.01.2023: added the same second block as for double bonds, not 100% clear if this holds.
        else if( chr == '#' && ((next_valid_is_nbr && parentheses >= 1 ) || (!got_in_line_nbr && parentheses == 0 ))){
            got_tbond = true;
        }
        else if( chr == '+' ){
            got_pos_chrg = true;
        }
        else if( chr == '-' /*&& in_brackets*/ ){
            got_neg_chrg = true;
        }
    }
    return is_two_letter_id;
}

void reader::gen_adj_list(){

    m_valence.clear();
    m_valence.resize( m_SMILES.size() );

    m_charge.clear();
    m_charge.resize( m_SMILES.size(), 0 );

    m_ring_membership.clear();
    m_ring_membership.resize( m_SMILES.size(), 0 );

    m_atom_labels.clear();


    get_nestedness();

    uint i=0;
    std::vector<uint> rid_count;
    std::vector<std::vector<uint>> connectivity;
    connectivity.resize( m_SMILES.size() );
    rid_count.resize( m_SMILES.size(), 0 );

    for( char chr : m_SMILES ){
        if( chr != 'H' && ( isupper(chr) || (islower(chr) && ( i== 0 || is_valid_atom_id( i )))) ){

            bool is_two_letter_atom_id = look_ahead( i+1, connectivity, rid_count );

            m_atom_labels.push_back( std::string() );
            m_atom_labels.back() = chr;

            if( is_two_letter_atom_id ){
                m_atom_labels.back() += m_SMILES.at(i+1);
            }
        }
        ++ i;
    }

    // convert the likely sparse connectivity list with SMILES
    // indices to a proper adjacency list (contiguous interval)
    m_adj_list.clear();
    m_adj_list.reserve( m_atom_labels.size() );

    // the same for a possibly sparse valence list (at this stage only holds bond order correctors)
    std::vector<uint> dense_valence;
    dense_valence.reserve( m_atom_labels.size() );

    std::vector<int> dense_charge;
    dense_charge.reserve( m_atom_labels.size() );

    std::vector<uint> dense_ring_membership;
    dense_ring_membership.reserve( m_atom_labels.size() );


    uint lowest_id = 0;
    std::map<uint, uint> old_new_id;
    for( uint sparse_id= 0; sparse_id < connectivity.size(); ++ sparse_id ){
        if( ! connectivity.at( sparse_id ).empty() ){
            old_new_id[sparse_id] = lowest_id;
            ++ lowest_id;
        }
    }

    for( uint sparse_id= 0; sparse_id < connectivity.size(); ++ sparse_id ){
        if( ! connectivity.at( sparse_id ).empty() ){
            m_adj_list.push_back( std::vector<uint>() );
            for( uint id : connectivity.at(sparse_id) ){
                m_adj_list.back().push_back( old_new_id[id] );
            }
            dense_valence.push_back( m_valence.at(sparse_id) );
            dense_charge.push_back( m_charge.at(sparse_id) );
            dense_ring_membership.push_back( m_ring_membership.at(sparse_id) );

            // C++17
            /*
            auto node_handle = m_aromatic_H.extract(sparse_id);
            if( node_handle ){
                m_aromatic_H.insert( old_new_id[sparse_id] );
            }
            */
            std::size_t erased = m_aromatic_H.erase(sparse_id);
            if( erased == 1 ){
                m_aromatic_H.insert( old_new_id[sparse_id] );
            }

        }
    }
    m_valence.swap( dense_valence );
    m_charge.swap( dense_charge );
    m_ring_membership.swap( dense_ring_membership );

    // efficient...
    std::map<std::pair<uint,uint>, uint> dense_bond_order;
    for( const std::pair<std::pair<uint,uint>, uint>& bond_pair : m_bond_orders ){
        dense_bond_order[std::pair<uint,uint>( old_new_id[ bond_pair.first.first], old_new_id[ bond_pair.first.second ])] = bond_pair.second;
    }
    std::swap( m_bond_orders, dense_bond_order );
}

void reader::count_implicit_H(){

    m_attached_H.clear();
    m_attached_H.reserve( m_atom_labels.size() );

    uint atom = 0;
    for(const std::vector<uint>& neighbors : m_adj_list){

        std::string label = m_atom_labels[ atom ];
        uint valence = m_valence[ atom ];
        uint nbr_count = neighbors.size();
        int charge = m_charge[ atom ];

        if( label == "C" ){
            if( valence == 0 ){
                if( nbr_count == 3 ){
                    m_attached_H.push_back( 1 );
                }
                else if( nbr_count == 2 ){
                    m_attached_H.push_back( 2 );
                }
                else if( nbr_count == 1 ){
                    m_attached_H.push_back( 3 );
                }
                else{ // 4 or 0 neighbors
                    m_attached_H.push_back( 0 );
                }
            }
            else if( valence == 1 ){

                if( nbr_count == 2 ){
                    m_attached_H.push_back( 1 );
                }
                //! 11.09.2023: added else, otherwise double push_back
                //! for single atoms is possible
                else if( nbr_count == 3 ){
                    m_attached_H.push_back( 0 );
                }
                else{
                    m_attached_H.push_back( 2 );
                }
            }
            //! 05.09.2023: C#C-H
            else if( valence == 2 && nbr_count == 1){
                m_attached_H.push_back( 1 );
            }
            else{
                m_attached_H.push_back(0);
            }
        }
        else if( label == "N" ){
            bool sp2_nbr = false;
            bool sp1_nbr = false;
            for( uint nbr : neighbors ){
                if( m_valence[nbr] == 1 ){
                    sp2_nbr = true;
                }
                else if( m_valence[nbr] == 2 && m_adj_list[nbr].size() == 2 ){
                    sp1_nbr = true;
                }
            }
            if( sp2_nbr ){
                if( nbr_count == 1 ){
                    if( valence == 0 ){
                        m_attached_H.push_back( 2 );
                    }
                    else{
                        m_attached_H.push_back( 1 );
                    }
                }
                else if( nbr_count == 2 ){
                    m_attached_H.push_back( 1 );
                }
                else{
                    m_attached_H.push_back( 0 );
                }
            }
            else if( sp1_nbr ){
                m_attached_H.push_back( 0 );
            }
            else{
                // assuming we have unprotonated species
                if( nbr_count == 2 ){
                    m_attached_H.push_back( 1 );
                }
                else if( nbr_count == 1 ){
                    m_attached_H.push_back( 2 );
                }
                else if( nbr_count > 0 ){ // 3 or 4 neighbors
                    m_attached_H.push_back( 0 );
                }
                else{ // ammonia
                    m_attached_H.push_back( 4 );
                }
            }
        }
        else if( label == "O" ){

            if( valence == 1 || nbr_count == 2 || charge <= -1 ){
                m_attached_H.push_back( 0 );
            }
            else{
                m_attached_H.push_back( 1 );
            }
        }
        else if( label == "c" ){

            if( nbr_count == 2 ){
                m_attached_H.push_back( 1 );
            }
            else{
                m_attached_H.push_back( 0 );
            }
        }
        else if( label == "n" ){
            if( m_aromatic_H.find(atom) != m_aromatic_H.end() ){
                m_attached_H.push_back(1);
            }
            else{
                m_attached_H.push_back( 0 );
            }
        }
        else if( label == "o" ){

            m_attached_H.push_back( 0 );
        }
        else if( label == "S" ){
            // valence 1 at this point means double bonded S
            // neighbor count >= 2 is for anything other than sulfhydryl
            if( valence == 1 || nbr_count >= 2 ){
                m_attached_H.push_back( 0 );
            }
            // R-SH
            else if( nbr_count == 1 ){
                m_attached_H.push_back( 1 );
            }
            // H-SH special case
            else{
                m_attached_H.push_back( 2 );
            }
        }
        else if( label == "s" ){

            m_attached_H.push_back( 0 );
        }
        /*
        else if( label == "P" ){

        }
        */
        else{
            m_attached_H.push_back( 0 );
        }

        ++ atom;
    }
}

void reader::correct_valence(){

    uint i = 0;
    for( const std::vector<uint>& adjacency : m_adj_list){

        m_valence.at(i) += adjacency.size();
        ++ i;
    }
}

void reader::assign_mass_number(){

    m_atomic_masses.clear();
    m_atomic_masses.reserve( m_atom_labels.size() );
    m_atomic_numbers.clear();
    m_atomic_numbers.reserve( m_atom_labels.size() );

    for( const std::string& label : m_atom_labels ){

        m_atomic_masses.push_back( amass.at( label ) );
        m_atomic_numbers.push_back( anum.at( label ) );
    }
}

void reader::single_bonds(){

    m_acyclic_single_bonds = 0;
    m_cyclic_single_bonds = 0;

    for( std::pair<const std::pair<uint, uint>, uint>& bond_order : m_bond_orders ){

        if( bond_order.second == 1 ){
            if( m_cyclic_nodes.find( bond_order.first.first ) == m_cyclic_nodes.end() ||
                m_cyclic_nodes.find( bond_order.first.second ) == m_cyclic_nodes.end() )
            {
                ++ m_acyclic_single_bonds;
            }
            else{
                bool both_same_system = false;
                // const std::pair<const std::vector<uint> *, std::vector<const std::vector<uint>*>*>& members_ring
                // first is a vector* of unique system member nodes:
                //! 03.02.2025: Added extended csb detection to correctly assign "C1CCCCC1-C1CCCCC1"-type bonds as asb.
                //! Note: The nalogous is applied to reader::double_bonds(), seems unnecessary for reader::triple_bonds().
                //! Note: This is primarily meant to be a straight-forward solution.
                //! Maybe there is minor potential for optimization.
                for(auto& members_ring : m_ring_sys_membership){

                    bool got_bond_1st = false;
                    bool got_bond_2nd = false;
                    // test if both cyclic nodes are present in the same ring(-system)
                    for( const uint member : *members_ring.first ){
                        if(
                            member == bond_order.first.first )
                        {
                            got_bond_1st = true;
                        }
                        else if(
                           member == bond_order.first.second )
                        {
                            got_bond_2nd = true;
                        }
                    }
                    if( got_bond_1st && got_bond_2nd ){
                        both_same_system = true;
                        break;
                    }
                }
                if( both_same_system ){
                    ++ m_cyclic_single_bonds;
                }
                else{
                    ++ m_acyclic_single_bonds;
                }
            }
        }
    }
}

void reader::double_bonds(){

    m_acyclic_double_bonds = 0;
    m_cyclic_double_bonds = 0;

    for( std::pair<const std::pair<uint, uint>, uint>& bond_order : m_bond_orders ){
        if( bond_order.second == 2 ){

            if( m_cyclic_nodes.find( bond_order.first.first ) == m_cyclic_nodes.end() ||
                m_cyclic_nodes.find( bond_order.first.second ) == m_cyclic_nodes.end() )
            {
                ++ m_acyclic_double_bonds;
            }
            else{
                bool both_same_system = false;
                for(auto& members_ring : m_ring_sys_membership){

                    bool got_bond_1st = false;
                    bool got_bond_2nd = false;
                    for( const uint member : *members_ring.first ){
                        if(
                            member == bond_order.first.first )
                        {
                            got_bond_1st = true;
                        }
                        else if(
                            member == bond_order.first.second )
                        {
                            got_bond_2nd = true;
                        }
                    }
                    if( got_bond_1st && got_bond_2nd ){
                        both_same_system = true;
                        break;
                    }
                }
                if( both_same_system ){
                    ++ m_cyclic_double_bonds;
                }
                else{
                    ++ m_acyclic_double_bonds;
                }
            }
        }
    }
    // correct cyclic double bonds in case we have aromatic SMILES input.
    for( const std::vector<const std::vector<uint> *>& system : m_ring_systems){
        uint aromatized_nodes = 0;
        std::unordered_set<uint> processed;

        uint aromatic_amide_cn = 0;

        for( const std::vector<uint> * ring : system ){
            for( uint node : *ring ){

                if( islower(m_atom_labels[node][0]) &&
                    processed.emplace( node ).second )
                {                    
                    ++ aromatized_nodes;

                    bool aromatic_amide = false;
                    if( m_atom_labels[node][0] == 'c' ){
                        bool aromatic_cn = false;
                        bool aromatic_carbonyl = false;
                        for( uint nbr : m_adj_list[node] ){
                            if( m_atom_labels[nbr] == "n" )
                            {
                                aromatic_cn = true;
                            }
                            else if( m_atom_labels[nbr] == "O" || m_atom_labels[nbr] == "S" ){
                                // bond partners are always ordered lower index first
                                if( node < nbr &&
                                    m_bond_orders[std::pair<uint, uint>(node,nbr)] == 2 )
                                {
                                    aromatic_carbonyl = true;
                                }
                                else if( m_bond_orders[std::pair<uint, uint>(nbr,node)] == 2)
                                {
                                    aromatic_carbonyl = true;
                                }
                            }
                        }
                        if( aromatic_cn && aromatic_carbonyl ){
                            aromatic_amide = true;
                        }
                    }
                    if( aromatic_amide ){
                        ++ aromatic_amide_cn;
                    }
                }
            }
        }
        //! 11.09.2023: added below. condition assures it's only done when
        //! aromatic SMILES representation is used.
        if( aromatized_nodes >= 1 ){

            uint delta = floor(aromatized_nodes / 2.0f);

            m_cyclic_double_bonds += delta;
            m_cyclic_double_bonds -= aromatic_amide_cn;

            m_cyclic_single_bonds -= delta;
            m_cyclic_single_bonds += aromatic_amide_cn;
        }
    }
}

void reader::triple_bonds(){

    m_acyclic_triple_bonds = 0;
    m_cyclic_triple_bonds = 0;

    for( std::pair<const std::pair<uint, uint>, uint>& bond_order : m_bond_orders ){
        if( bond_order.second == 3 ){

            if( m_cyclic_nodes.find( bond_order.first.first ) == m_cyclic_nodes.end() ||
                m_cyclic_nodes.find( bond_order.first.second ) == m_cyclic_nodes.end() )
            {
                ++ m_acyclic_triple_bonds;
            }
            else{
                ++ m_cyclic_triple_bonds;
            }
        }
    }
}

void reader::node_order(){

    uint node = 0;
    for( const std::vector<uint>& nbr : m_adj_list){
        if( nbr.size() == 1 ){
            ++ m_monovalent;
        }
        else if( nbr.size() == 2){
            if( m_cyclic_nodes.find( node ) == m_cyclic_nodes.end() ){
                ++ m_acyclic_divalent;
            }
            else{
                ++ m_cyclic_divalent;
            }
        }
        else if( nbr.size() == 3){
            if( m_cyclic_nodes.find( node ) == m_cyclic_nodes.end() ){
                ++ m_acyclic_trivalent;
            }
            else{
                ++ m_cyclic_trivalent;
            }
        }
        else if( nbr.size() == 4){
            if( m_cyclic_nodes.find( node ) == m_cyclic_nodes.end() ){
                ++ m_acyclic_tetravalent;
            }
            else{
                ++ m_cyclic_tetravalent;
            }
        }
        ++node;
    }
}

void reader::charge(){
    for(int charge : m_charge){
        if( charge < 0){
            m_neg_charges += abs(charge);
        }
        else if( charge > 0 ){
            m_pos_charges += charge;
        }
    }
}

void reader::h_bond(){

    uint node = 0;
    for( const std::string& label : m_atom_labels ){
        char first = label[0];
        uint attached_H = m_attached_H[node];
        if( first == 'N' ){
            // hbd
            if( attached_H > 0 ){
                m_hbd_sites += attached_H;
                ++ m_hbd;
            }
            // hba
            bool pair_undeloc = true;
            bool guanidine_N = false;
            for( uint nbr : m_adj_list[node] ){
                char neighbor_type = m_atom_types[nbr];
                if( neighbor_type == C_2_bnd_amide_hba || // note: the a here implies that its N is not a donor, only the O an acceptor.
                    neighbor_type == C_2_bnd_amide_hbd ||
                    neighbor_type == S_2x2_bnd_6_val ||
                    neighbor_type == P_1x2_bnd_5_val ||
                    neighbor_type == C_2_bnd_thiamide
                ){
                    pair_undeloc = false;
                }
                // captures guanidine N aswell as arylamidines; other amidines are not detected
                else if( neighbor_type >= 4 && neighbor_type <= 7 ){
                    guanidine_N = true;
                }
            }
            if( pair_undeloc ){
                // just as marvin does it (pH independent):
                // every N of a guanidine that has less than 3 heavy neighbors
                // gets 1 hba-site assigned
                if( guanidine_N && m_adj_list[node].size() <= 2 ){
                    ++ m_hba;
                    ++ m_hba_sites;
                }
                else if( m_valence[ node ] <= 3){
                    ++ m_hba;
                    m_hba_sites += 4 - m_valence[node] - attached_H;
                }
            }
        }
        else if( first == 'n' ){
            // hbd
            if( attached_H > 0 ){ // assuming assuming valid input, will be == 1
                m_hbd_sites += attached_H;
                ++ m_hbd;
            }
            //! Initially as out-commented below, then switches to <= 2 and back again
            //! to have it like the RDKit. However, ChemAxon does not count such.
            // hba
            //else if( m_valence[node] <= 3 ){
            else if( m_valence[node] <= 2 ){
                ++ m_hba_sites;
                ++ m_hba;
            }
        }
        else if( first == 'O' ){
            // hbd
            if( attached_H > 0 ){
                m_hbd_sites += attached_H;
                ++ m_hbd;
            }
            // hba
            bool ester_sp3_O = false;
            ushort aryl_neighbors = 0;
            for( uint nbr : m_adj_list[node] ){
                char neighbor_type = m_atom_types[nbr];
                if( (neighbor_type == C_2_bnd_ester ||
                     neighbor_type == C_2_bnd_amide_hba ||
                     neighbor_type == C_2_bnd_amide_hbd ) &&
                     m_adj_list[node].size() >= 2 )
                {
                    ester_sp3_O = true;
                    //break;
                }
                else if( neighbor_type < 0){
                    ++ aryl_neighbors;
                }
            }
            if( !ester_sp3_O ){
                if( aryl_neighbors <= 1 ){
                    m_hba_sites += 2;
                    m_hba_sites += abs(m_charge[ node ]);
                    ++ m_hba;
                }
            }

        }
        /*
        else if( first == 'o' ){ // omitting test for positive charge
            ++ m_hba_sites;
            ++ hba
        }
        */
        ++ node;
    }
}


uint reader::neighbor_count_of(uint node){
    if( node < m_adj_list.size() ){
        return m_adj_list[ node ].size();
    }
    return 0;
}

bool reader::append_if_new_ring(std::vector<std::vector<uint> >& size_container, const std::vector<uint> qry_ring ){
    for(const std::vector<uint>& db_ring : size_container){
        uint similars = 0;
        for(uint db_node : db_ring){
            for(uint qry_node : qry_ring){
                if(qry_node == db_node){
                    ++similars;
                }
            }
            if(similars == qry_ring.size()){
                return false;
            }
        }
    }
    size_container.emplace_back(qry_ring);

    return true;
}

bool reader::node_in_range(std::vector<uint>::const_iterator begin, std::vector<uint>::const_iterator end, uint lookup_node){
    while( begin <= end ){
        if( *begin == lookup_node ){
            return true;
        }
        ++begin;
    }
    return false;
}

uint reader::ring_traverse(uint current, uint depth){

    if( depth >= 3){
        if( current == m_orig ){
            return depth;
        }
    }

    uint node_coval = this->neighbor_count_of( current );
    if( node_coval <= 1 || depth >= 17 ){
        return 0;
    }

    m_temp_path.at( depth ) = current;
    uint closure_depth = 0;

    for( const uint neighbor : m_adj_list[ current ] ){

        if((depth < 2 && neighbor != m_orig) || depth >= 2 ){

            if( depth < 2 || !node_in_range(m_temp_path.begin()+1, m_temp_path.begin()+depth-1, neighbor) ){

                closure_depth = ring_traverse( neighbor, depth+1 );

                if( closure_depth ){
                    std::vector<uint>::iterator begin = m_temp_path.begin();
                    std::vector<uint>::iterator end = begin+closure_depth;

                    if(closure_depth == 17){
                        --closure_depth;
                    }
                    append_if_new_ring( m_rings[ closure_depth-3 ], std::vector<uint>(begin, end) );
                    closure_depth = 0;
                }
            }
        }
    }
    return closure_depth;
}

bool reader::find_rings(){

    m_temp_path.resize(17,0);
    m_rings.assign(14,std::vector<std::vector<uint>>());

    for( uint node = 0; node < m_adj_list.size(); ++node ){
        m_orig = node;
        ring_traverse(node, 0);
    }
    if( !m_rings.empty() ){
        return true;
    }
    return false;
}


void reader::gen_ring_adj_lst(){

    std::vector<const std::vector<uint>*> all_smallest_rings;

    for( const std::vector<uint>* ring : m_symm_SSSR ){
        all_smallest_rings.emplace_back( ring );
    }

    for( const std::vector<uint>* ring : m_symm_SSSR_23 ){
        bool add = true;
        for( const std::vector<uint>* added : all_smallest_rings ){
            if( ring == added ){
                add = false;
                break;
            }
        }
        if( add ){
            all_smallest_rings.emplace_back( ring );
        }
    }

    for( uint l = 0; l < all_smallest_rings.size(); ++l ){
        if( m_bulk_ring_adj_lst.at( all_smallest_rings[ l ]).empty() ){
            m_ring_adj_lst.insert(std::pair<const std::vector<uint>*,std::vector<const std::vector<uint>*>>( all_smallest_rings[ l ],std::vector<const std::vector<uint>*>()) );
            m_shared_nodes_lst.insert(std::pair<const std::vector<uint>*, std::vector<std::vector<uint>>>( all_smallest_rings[ l ], std::vector<std::vector<uint>>()) );

        }
        else{

            for( const std::vector<uint>* adjacent : m_bulk_ring_adj_lst.at( all_smallest_rings[ l ] ) ){
                bool set_member = false;
                for( uint r = l+1; r < all_smallest_rings.size(); ++r ){
                    if( adjacent == all_smallest_rings[ r ] ){
                        set_member = true;
                        break;
                    }
                }
                if( set_member ){
                    m_ring_adj_lst[ all_smallest_rings[ l ] ].emplace_back( adjacent );
                    m_ring_adj_lst[ adjacent ].emplace_back( all_smallest_rings[ l ] );

                    std::vector<uint> shared_nodes;
                    for( uint l_node : *all_smallest_rings[ l ] ){
                        for( uint r_node : *adjacent ){
                            if( l_node == r_node ){
                                shared_nodes.emplace_back( l_node );
                                break;
                            }
                        }
                    }

                    m_shared_nodes_lst[ all_smallest_rings[ l ] ].emplace_back( shared_nodes );
                    m_shared_nodes_lst[ adjacent ].emplace_back( shared_nodes );
                }
            }
        }
    }
}

void reader::gen_bulk_ring_adj_lst(){

    uint count_rings = 0;
    for(auto& container : m_rings){
        count_rings += container.size();
    }
    std::vector<const std::vector<uint>*> bulk_rings;
    bulk_rings.reserve( count_rings );
    for(auto& container : m_rings){
        for(auto& ring : container){
            bulk_rings.emplace_back(&ring);
        }
    }

    for( const std::vector<std::vector<uint>>& container : m_rings ){
        for( const std::vector<uint>& ring : container ){
            m_bulk_ring_adj_lst[ &ring ] = std::vector<const std::vector<uint>*>();
        }
    }
    for(uint ref_id = 0; ref_id < bulk_rings.size(); ++ref_id){
        const std::vector<uint>* ref_ring = bulk_rings[ ref_id ];
        for(uint cmp_id = ref_id+1; cmp_id <  bulk_rings.size(); ++cmp_id){
            const std::vector<uint>* cmp_ring = bulk_rings[ cmp_id ];
            bool done_compare = false;
            for( uint ref_node : *ref_ring ){
                for( uint cmp_node : *cmp_ring ){
                    if( ref_node == cmp_node ){
                        m_bulk_ring_adj_lst[ ref_ring ].emplace_back( cmp_ring );
                        m_bulk_ring_adj_lst[ cmp_ring ].emplace_back( ref_ring );
                        done_compare = true;
                        break;
                    }
                }
                if(done_compare){
                    break;
                }
            }
        }
    }
}

const std::vector<uint>* reader::bulk_ring_system_traversal(const std::vector<uint>* current, std::vector<const std::vector<uint> *>* const current_system, const std::vector<uint>* last){

    const std::vector<uint>* last_ring_added_per_instance = nullptr;

    for( const std::vector<uint>* adj_ring : m_bulk_ring_adj_lst.at( current ) ){

        if(adj_ring != last){

            auto ring_assoc_system = m_ring_sys_membership.find(adj_ring);
            if(ring_assoc_system == m_ring_sys_membership.end()){

                current_system->emplace_back(adj_ring);
                m_ring_sys_membership[ adj_ring ] = current_system;
                bulk_ring_system_traversal( adj_ring,current_system, current );
                last_ring_added_per_instance = adj_ring;
            }
        }
    }
    return last_ring_added_per_instance;
}

void reader::ring_systems(){

    bool all_rings_assigned = false;
    const std::vector<uint>* last_ring_added_dfs = nullptr;

    std::map<const std::vector<uint>*,std::vector<const std::vector<uint>*> >::reverse_iterator iter_ring_adj_other = m_bulk_ring_adj_lst.rbegin();
    for( ;iter_ring_adj_other != m_bulk_ring_adj_lst.rend(); iter_ring_adj_other++ ){

        auto& ring_adj_other = *iter_ring_adj_other;

        if( last_ring_added_dfs != ring_adj_other.first ){
            auto ring_assoc_system = m_ring_sys_membership.find(ring_adj_other.first);

            if( ring_assoc_system == m_ring_sys_membership.end() ){

                m_sys_size_ring.emplace_back(std::vector<const std::vector<uint>*>({ring_adj_other.first}));
                m_ring_sys_membership[ ring_adj_other.first]  = &(m_sys_size_ring.back());
                for( const std::vector<uint>* nb_ring : m_bulk_ring_adj_lst.at(ring_adj_other.first) ){

                    if(last_ring_added_dfs != nb_ring){
                        m_sys_size_ring.back().emplace_back(nb_ring);
                        m_ring_sys_membership[nb_ring] = &(m_sys_size_ring.back());
                    }
                }

                for( const std::vector<uint>* nb_ring : m_bulk_ring_adj_lst.at(ring_adj_other.first) ){

                    if(m_ring_sys_membership.size() >= m_bulk_ring_adj_lst.size()){

                        all_rings_assigned = true;
                        break;
                    }

                    else if(last_ring_added_dfs != nb_ring){
                        last_ring_added_dfs = bulk_ring_system_traversal(nb_ring, &m_sys_size_ring.back(), ring_adj_other.first);
                    }
                }
            }
        }
        if( all_rings_assigned ){
            break;
        }
    }

    for( std::vector<const std::vector<uint>*>& system : m_sys_size_ring){
        std::sort(system.begin(), system.end(), [](const std::vector<uint>* a, const std::vector<uint>* b){ return a->size() < b->size(); });
    }
}

void reader::ring_system_traversal_23( const std::vector<uint>* current, std::vector<const std::vector<uint> *>& current_system ){
    current_system.emplace_back(current);
    for( const std::vector<uint>* adj_ring : m_bulk_ring_adj_lst.at( current ) ){

        bool possible_new_member = true;
        for(int idx = current_system.size()-1; idx >= 0; --idx){

            if( adj_ring == current_system[ idx ] ){
                possible_new_member = false;
                break;
            }
        }
        if( possible_new_member ){
            bool quat_free_neighbor = true;
            for(uint node : *adj_ring){
                if( neighbor_count_of( node ) >= 4 ){
                    quat_free_neighbor = false;
                }
            }
            if(quat_free_neighbor){
                ring_system_traversal_23( adj_ring, current_system );
            }
        }
    }
}

void reader::sssr(){

    std::vector<std::vector<std::map<uint,std::set<uint>> >> rings_as_edges;
    std::vector<std::map<uint, std::set<uint>>> unique_edges;
    std::vector<std::set<uint>> unique_nodes;

    std::vector<std::map<uint,std::set<uint>> > edges_sssr;

    std::vector<std::vector<std::map<uint,std::set<uint>> >> rings_as_edges_23;
    std::vector<std::map<uint, std::set<uint>*>> unique_edges_23;
    std::vector<std::set<uint>> unique_nodes_23;

    for( std::vector<const std::vector<uint>*>& system : m_sys_size_ring ){

        rings_as_edges.emplace_back( std::vector<std::map<uint,std::set<uint>> >() );
        unique_edges.emplace_back( std::map<uint, std::set<uint>>() );
        unique_nodes.emplace_back( std::set<uint>() );

        for(const std::vector<uint>* ring : system){

            bool quat_free = true;
            rings_as_edges.back().emplace_back( std::map<uint,std::set<uint> >() );

            std::vector<uint>::const_iterator last = ring->end()-1;
            std::vector<uint>::const_iterator current = ring->begin();

            do{
                if( *current < *last ){
                    rings_as_edges.back().back()[*current].insert(*last);
                    unique_edges.back()[*current].insert(*last);
                }
                else{
                    rings_as_edges.back().back()[*last].insert(*current);
                    unique_edges.back()[*last].insert(*current);
                }

                unique_nodes.back().insert(*current);

                if( this->neighbor_count_of(*current) >= 4 ){
                    quat_free = false;
                }

                last = current;
                current ++;

            }while( current != ring->end() );

            if( quat_free ){

                bool new_quat_free_sys = true;
                if( !m_sys_size_ring_23.empty() ){
                    uint sys_idx = 0;
                    for( std::vector<const std::vector<uint>*>& system : m_sys_size_ring_23 ){
                        for( const std::vector<uint>* const ring_in_latest_23_sys : system){
                            if( ring_in_latest_23_sys == ring ){
                                new_quat_free_sys = false;

                                break;
                            }
                        }
                        ++ sys_idx;
                    }
                }
                if( new_quat_free_sys ){

                    m_sys_size_ring_23.emplace_back( std::vector<const std::vector<uint>*>() );
                    rings_as_edges_23.emplace_back( std::vector<std::map<uint,std::set<uint>> >() );
                    rings_as_edges_23.back().reserve( rings_as_edges.back().size() );
                    unique_edges_23.emplace_back( std::map<uint, std::set<uint>*>() );
                    unique_nodes_23.emplace_back( std::set<uint>() );

                    ring_system_traversal_23( ring, m_sys_size_ring_23.back() );
                }
            }
        }
    }

    uint s = 0;
    for( const std::vector<const std::vector<uint>*>& sys_23 : m_sys_size_ring_23 ){
        for( const std::vector<uint>* ring_23 : sys_23 ){
            uint g = 0;
            for(const std::vector<const std::vector<uint>*>& sys : m_sys_size_ring){

                for( uint i = 0; i < sys.size(); ++i ){

                    if( sys[ i ] == ring_23 ){

                        rings_as_edges_23.at(s).emplace_back( rings_as_edges[ g ].at(i) );

                        for( const auto& adj : rings_as_edges_23[ s ].back() ){
                            unique_nodes_23[ s ].insert(adj.first);

                            for(uint neighbor : adj.second){
                                unique_nodes_23[ s ].insert(neighbor);
                            }

                            unique_edges_23[ s ].insert(std::pair<uint, std::set<uint>*>(adj.first, &unique_edges[ g ].at( adj.first )));
                        }
                        break;
                    }
                }
                ++g;
            }
        }
        ++s;
    }
    uint sys = 0;
    for( std::vector<const std::vector<uint>*>& system : m_sys_size_ring ){

        uint edge_count = 0;
        for(const std::pair<uint, std::set<uint>>& edges : unique_edges[ sys ]){
            edge_count += edges.second.size();
        }
        uint rings_expected = edge_count - unique_nodes[ sys ].size() + 1;
        uint rings_expected_symm = edge_count - unique_nodes[ sys ].size() + 2;

        m_cyclic_edge_count += edge_count;

        const std::vector<uint>* last_discarded = nullptr;
        const std::vector<uint>* before_last_discarded = nullptr;
        const std::vector<uint>* p_ring = nullptr;

        uint subsys_size = 0;
        uint last_discarded_shares_nodes_with_n = 0;

        for( int ref_id = rings_as_edges[ sys ].size()-1; ref_id >= 0; --ref_id ){
            std::map<uint,std::set<uint>>& ref_ring = rings_as_edges[ sys ].at(ref_id);
            p_ring = system[ ref_id ];
            uint edges_in_smaller = 0;

            for(auto& edge : ref_ring){
                uint e0 = edge.first;

                for(uint e1 : edge.second){

                    for( int cmp_id = 0; cmp_id < ref_id; ++cmp_id ){
                        std::map<uint,std::set<uint>>& cmp_ring = rings_as_edges[ sys ].at(cmp_id);

                        if( cmp_ring.find(e0) != cmp_ring.end() &&
                            cmp_ring.at(e0).find(e1) != cmp_ring.at(e0).end() )
                        {
                            ++edges_in_smaller;
                            break;
                        }
                    }
                }
            }

            if( edges_in_smaller < p_ring->size() ){
                edges_sssr.emplace_back(ref_ring);
                m_topol_SSSR.emplace_back( p_ring );
                m_symm_SSSR.emplace_back( p_ring );
                ++ subsys_size;
            }
            else if( ref_id < rings_expected_symm ){

                std::vector<uint> node_in_n_rings = std::vector<uint>( m_adj_list.size(), 0 );
                std::vector<const std::vector<uint>*>::iterator smaller_ring_it = system.begin();
                while( *smaller_ring_it != p_ring ){

                    for( uint node_in_smaller : **smaller_ring_it ){
                        ++node_in_n_rings[ node_in_smaller ];
                    }
                    smaller_ring_it++;
                }

                uint shared_among_smaller = 0;
                uint shared_with_current = 0;
                uint node_sys = 0;
                for( uint n : node_in_n_rings ){
                    if( n >= 2 ){
                        ++shared_among_smaller;
                        bool node_shared_w_current = false;
                        for( uint node_in_current : *p_ring ){
                            if( node_in_current == node_sys ){
                                node_shared_w_current = true;
                                break;
                            }
                        }
                        if( node_shared_w_current ){
                            ++shared_with_current;
                        }
                    }
                    ++node_sys;
                }
                if( shared_among_smaller >= 3 &&
                    shared_with_current < shared_among_smaller )
                {
                    before_last_discarded = last_discarded;
                    last_discarded = p_ring;
                }
            }
        }

        if(subsys_size < rings_expected && last_discarded != nullptr ){

            m_topol_SSSR.emplace_back(last_discarded);
        }

        if( subsys_size < rings_expected_symm ){
            if( last_discarded != nullptr && last_discarded_shares_nodes_with_n < subsys_size ){

                m_symm_SSSR.emplace_back( last_discarded );
            }
            if( before_last_discarded != nullptr && last_discarded_shares_nodes_with_n < subsys_size ){

                m_symm_SSSR.emplace_back( before_last_discarded );
            }
        }
        ++sys;
    }

    sys = 0;
    for( std::vector<const std::vector<uint>*>& system : m_sys_size_ring_23 ){

        if( system.size() == 1 ){
            m_symm_SSSR_23.emplace_front( system.front() );
        }
        else{

            uint edge_count_23 = 0;
            for(const std::pair<uint, std::set<uint>*>& edges : unique_edges_23[ sys ]){
                edge_count_23 += edges.second->size();
            }
            uint rings_expected_23 = edge_count_23 - unique_nodes_23[ sys ].size() + 1;

            const std::vector<uint>* last_discarded = nullptr;
            const std::vector<uint>* before_last_discarded = nullptr;
            const std::vector<uint>* p_ring = nullptr;
            uint subsys_size = 0;
            for( int ref_id = rings_as_edges_23[ sys ].size()-1; ref_id >= 0; --ref_id ){
                const std::map<uint,std::set<uint>>& ref_ring = rings_as_edges_23[ sys ].at(ref_id);
                p_ring = system[ ref_id ];
                uint edges_in_smaller = 0;

                for(auto& edge : ref_ring){
                    uint e0 = edge.first;

                    for(uint e1 : edge.second){

                        for( int cmp_id = 0; cmp_id < ref_id; ++cmp_id ){
                            const std::map<uint,std::set<uint>>& cmp_ring = rings_as_edges_23[ sys ].at(cmp_id);

                            if( cmp_ring.find(e0) != cmp_ring.end() && cmp_ring.at(e0).find(e1) != cmp_ring.at(e0).end() ){

                                ++edges_in_smaller;
                                break;
                            }
                        }
                    }
                }

                if( edges_in_smaller < ref_ring.size() ){
                    m_symm_SSSR_23.emplace_front( p_ring );
                    subsys_size++;
                }
                else{
                    before_last_discarded = last_discarded;
                    last_discarded = p_ring;
                }
            }

            if(subsys_size < rings_expected_23){

                if( last_discarded != nullptr && unique_nodes_23[ sys ].size() > last_discarded->size() ){
                    m_symm_SSSR_23.emplace_back(last_discarded);
                }
                if( before_last_discarded != nullptr && unique_nodes_23[ sys ].size() > last_discarded->size() ){
                    m_symm_SSSR_23.emplace_back(before_last_discarded);
                }
            }
        }
        ++sys;
    }

    uint in_min_2 = 0;
    for(auto& system : unique_edges){
        for(auto& edges : system){
            for(auto& end : edges.second){
                uint occurrences = 0;
                for(auto& ring : edges_sssr){

                    if( ring.find(edges.first) != ring.end() ){

                        if(ring.at(edges.first).find(end) != ring.at(edges.first).end()){
                            ++occurrences;
                        }
                    }
                }
                if (occurrences >= 2){
                    ++ in_min_2;
                }
            }
        }
    }
    m_bonds_in_min_2_rings = in_min_2;

    for(const std::vector<uint>* ring : m_topol_SSSR){
        for(uint cyclic_node : *ring){
            m_cyclic_nodes.insert(cyclic_node);
        }
    }

    in_min_2 = 0;
    for( uint l_node : m_cyclic_nodes ){
        uint occurrences = 0;
        for( auto& ring : m_topol_SSSR ){
            for( uint r_node : *ring){
                if( l_node == r_node ){
                    ++ occurrences;
                    break;
                }
            }
        }
        if (occurrences >= 2){
            ++ in_min_2;
        }
    }
    m_atoms_in_min_2_rings = in_min_2;

    std::set<const std::vector<uint>*> processed_rings;

    while(processed_rings.size() < m_topol_SSSR.size()){
        bool added_rings = true;

        m_ring_systems.emplace_back( std::vector<const std::vector<uint>*>() );
        std::vector<const std::vector<uint>*>& current_system = m_ring_systems.back();
        while(added_rings){
            added_rings = false;
            for(const std::vector<uint>* ring : m_topol_SSSR){
                if( processed_rings.find(ring) == processed_rings.end() ){
                    if(current_system.empty()){
                        current_system.emplace_back(ring);
                        processed_rings.insert(ring);
                        added_rings = true;
                    }
                    else{
                        for(uint node : *ring){
                            for(auto& present_ring: m_ring_systems.back()){
                                for(uint present_node : *present_ring){
                                    if(node == present_node){
                                        m_ring_systems.back().emplace_back(ring);
                                        processed_rings.insert(ring);
                                        added_rings = true;
                                        break;
                                    }
                                }
                                if(added_rings){
                                    break;
                                }
                            }
                            if(added_rings){
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

// Below come methods related to planarizability check
bool reader::in_di_tri_val_rings(const std::vector<uint>* lookup_ring){

    for( const std::vector<uint>* di_tri_valent : m_symm_SSSR_23){
        if( lookup_ring == di_tri_valent ){
            return true;
        }
    }
    return false;
}

bool reader::ring_in_set(std::vector<const std::vector<uint>*> ring_set, const std::vector<uint>* lookup_ring){

    for( const std::vector<uint>* member : ring_set){
            if( lookup_ring == member ){
                return true;
            }
        }
        return false;
}

uint reader::shared_with_two_rings(const std::vector<uint> * first, const std::vector<uint> * second, uint idx_2nd){

    std::vector<const std::vector<uint> *> shared_rings;
    for(const std::vector<uint>* nb_first : m_ring_adj_lst.at(first) ){
        for(const std::vector<uint>* nb_second : m_ring_adj_lst.at(second) ){
            if( nb_first == nb_second ){
                uint shared_3rd_1st = 0;
                for(uint node_1 : *first){
                    for( uint node_3 : *nb_first ){
                        if(node_1 == node_3){
                            ++shared_3rd_1st;
                        }
                    }
                }
                if( shared_3rd_1st >= first->size() -1 ){
                    break;
                }
                else{
                    uint shared_3rd_2nd = 0;
                    for(uint node_2 : *second){
                        for( uint node_3 : *nb_first ){
                            if(node_2 == node_3){
                                ++shared_3rd_2nd;
                            }
                        }
                    }
                    if( shared_3rd_2nd >= second->size() -1 ){
                        break;
                    }
                }

                shared_rings.emplace_back( nb_first );
                break;
            }
        }
    }
    uint max_shared = 0;

    std::vector<uint> unique_12;
    unique_12.reserve( first->size() + second->size() );

    for( uint node_1st : *first){
        bool unique = true;
        for( uint node_shared : m_shared_nodes_lst.at(first)[ idx_2nd ] ){
            if( node_1st == node_shared ){
                unique = false;
                break;
            }
        }
        if( unique ){
            unique_12.emplace_back( node_1st );
        }
    }
    for( uint node_2nd : *second){
        bool unique = true;
        for( uint node_shared : m_shared_nodes_lst.at(first)[ idx_2nd ] ){
            if( node_2nd == node_shared ){
                unique = false;
                break;
            }
        }
        if( unique ){
            unique_12.emplace_back( node_2nd );
        }
    }

    for(const std::vector<uint>* shared : shared_rings){
        m_triple_shared = 0;

        uint shared_3rd = 0;
        for( uint qry_node : *shared ){
            for( uint us12_node : unique_12 ){
                if( qry_node == us12_node ){
                    ++shared_3rd;
                    break;
                }
            }
            for( uint s12_node : m_shared_nodes_lst.at(first)[ idx_2nd ] ){
                if( qry_node == s12_node ){
                    shared_3rd++;
                    ++m_triple_shared;
                    break;
                }
            }
        }

        if( shared_3rd > max_shared ){
            max_shared = shared_3rd;
        }
    }

    return max_shared;
}

void reader::extract_planarizable_subgraph( const std::vector<uint> * const smallest_ring ){

    std::vector<const std::vector<uint>*>* current_part = nullptr;

    for( std::vector<const std::vector<uint>*>& part : m_planarizable ){
        for( const std::vector<uint>* present : part ){
            for( const std::vector<uint>* neighbor : m_ring_adj_lst.at(present) ){
                if( smallest_ring == neighbor ){
                    current_part = &part;
                    break;
                }
            }
            if( current_part ){
                break;
            }
        }
        if(current_part){
            break;
        }
    }
    bool isplanarizable = true;

    for( int nb_idx = 0; nb_idx < m_shared_nodes_lst.at( smallest_ring ).size(); ++nb_idx ){

        uint ref_nb_shared = m_shared_nodes_lst.at( smallest_ring )[ nb_idx ].size();

        if(ref_nb_shared >= 4 && ref_nb_shared < smallest_ring->size()-1){
            isplanarizable = false;
            break;
        }
        else if(ref_nb_shared == 3){
            uint ref_frag = smallest_ring->size() - ref_nb_shared;
            uint cmp_frag = m_ring_adj_lst.at( smallest_ring )[ nb_idx ]->size() - ref_nb_shared;

            uint bridge_len = ref_nb_shared;
            uint large_frag_1 = ref_frag;
            uint large_frag_2 = cmp_frag;

            if( ref_frag <= ref_nb_shared ){
                bridge_len = ref_frag;
                large_frag_1 = ref_nb_shared;
            }
            if( cmp_frag <= bridge_len ){
                large_frag_2 = bridge_len;
                bridge_len = cmp_frag;
            }
            bridge_len -= 2;

            if( bridge_len > 1 || ( large_frag_1 + large_frag_2 < 8 ) ){

                isplanarizable = false;
                break;
            }
        }
        if( isplanarizable && m_ring_adj_lst.at(smallest_ring).size() >= 3 ){
            uint two_shared_with_third = shared_with_two_rings(smallest_ring, m_ring_adj_lst.at( smallest_ring )[ nb_idx ], nb_idx);

            uint unq_nodes_count_12 = smallest_ring->size() + m_ring_adj_lst.at( smallest_ring )[ nb_idx ]->size() - 2*ref_nb_shared;

            if( ref_nb_shared == 2 ){
                if( two_shared_with_third-m_triple_shared >= 4 && two_shared_with_third < unq_nodes_count_12-1 ){

                    isplanarizable = false;
                    break;
                }
            }
            else{
                if( two_shared_with_third >= 4 && two_shared_with_third < unq_nodes_count_12-1 ){

                    isplanarizable = false;
                    break;
                }

            }
        }
    }

    if( isplanarizable ){

        if( current_part != nullptr ){
            current_part->emplace_back(smallest_ring);
        }
        else{
            m_planarizable.emplace_back(std::vector<const std::vector<uint>*>( {smallest_ring} ));
        }
        for(uint planarizable_node : *smallest_ring){
            m_planarizable_nodes.insert(planarizable_node);
        }
    }
}

uint reader::ring_system_traversal(const std::vector<uint>* current, std::vector<const std::vector<uint> *>& current_system, uint depth){

    current_system.emplace_back( current );
    for( const std::vector<uint>* adj_ring : m_ring_adj_lst.at( current ) ){
        bool new_member = true;
        for(const std::vector<uint>* member : current_system){
            if( adj_ring == member ){
                new_member = false;
                break;
            }
        }
        if( new_member ){
            depth += ring_system_traversal(adj_ring, current_system, depth+1);
        }
    }
    return depth;
}

void reader::spheroid_system_detection(){

    for( uint ref_pos = 0; ref_pos < m_symm_SSSR_23.size(); ++ref_pos ){

        bool unknown = true;

        for(const std::vector<const std::vector<uint> *>& system : m_planarizable){
            for(const std::vector<uint>* member_ring : system){
                if(m_symm_SSSR_23[ ref_pos ] == member_ring){
                    unknown = false;
                    break;
                }
            }
            if( !unknown ){
                break;
            }
        }

        if( unknown ){

            m_ring_systems.emplace_back( std::vector<const std::vector<uint> *>() );
            extract_planarizable_subgraph( m_symm_SSSR_23[ ref_pos ] );
        }
    }
}

void reader::node_ring_membership(){

    m_node_ring_membership.resize( m_adj_list.size() );
    for( uint node : m_cyclic_nodes ){
        for( const std::vector<uint>* const ring : m_symm_SSSR ){
            for( uint ring_node : *ring ){
                if( ring_node == node ){
                    m_node_ring_membership[ node ].emplace_back(ring);
                    break;
                }
            }
        }
    }
}

bool reader::is_planarizable_node(uint query_node){

    if(m_cyclic_nodes.find(query_node) == m_cyclic_nodes.end()){
        return true;
    }
    if(m_planarizable_nodes.find(query_node) != m_planarizable_nodes.end()){
        return true;
    }
    return false;
}


void reader::init_pattern(){
    uint i = 0;
    for( uint label : m_atomic_numbers ){
        if( label == 6 ){
            if( m_atom_labels[i] == "c" ){
                m_atom_types[i] = c_123_val;
            }
            else{
                m_atom_types[i] = C_1234_val;
            }
        }
        else if( label == 7 ){
            if( m_atom_labels[i] == "n" ){
                if( m_adj_list[i].size() == 2 ){
                    m_atom_types[i] = n_2_val;
                }
                else if( m_charge[i] == 0 ){
                    m_atom_types[i] = n_3_val;
                }
                else{
                    m_atom_types[i] = n_3_val_pos;
                }
            }
            else{
                if( m_adj_list[i].size() <= 2 ){
                    m_atom_types[i] = N_12_val;
                }
                else if( m_adj_list[i].size() == 3 ){
                    if( m_charge[i] == 0 ){
                        m_atom_types[i] = N_3_val;
                    }
                    else{
                        m_atom_types[i] = N_3_val_pos;
                    }
                }
                else{
                    m_atom_types[i] = N_4_val;
                }
            }
        }
        else if( label == 8 ){
            if( m_atom_labels[i] == "o" ){
                m_atom_types[i] = o_2_val;
            }
            else{
                if( m_charge[i] == 0 ){
                    m_atom_types[i] = O_12_val;
                }
                else{
                    m_atom_types[i] = O_1_val_neg;
                }
            }
        }
        ++ i;
    }
}


//! Under construction

void reader::largest_ring_of_group(
    const std::map<uint, std::set<uint> >& node_adj_lst,
    std::map<uint, std::pair<const std::vector<uint>*,const std::vector<uint>*>>& nodes_parent_rings/*,
    const std::set<uint>& bridging_nodes*/)
{
    m_largest_rings.emplace_back(std::vector<uint>());
    std::vector<uint>* largest_ring = &m_largest_rings.back();
    largest_ring->reserve( node_adj_lst.size() );

    const std::vector<uint>* current_ring = nullptr;

    uint orig = 0;
    uint current = 0;

    for( const std::pair<const uint, std::set<uint>>& adjacency : node_adj_lst ){
        orig = adjacency.first;

        // get a node that's shared by two rings only
        std::map<uint, std::pair< const std::vector<uint>*, const std::vector<uint>* >>::const_iterator result = nodes_parent_rings.find( orig );
        if( result != nodes_parent_rings.end() ){

            // watch out for a neighbor node that's just in the current ring
            //! is this always working? the other ring?
            //! 25.02.2021: Above is a good point: the answer is no. Therefore,
            //! all rings must be checked.
            current_ring = result->second.first;
            bool neighbor_in_same_ring = false;
            for(uint neighbor : node_adj_lst.at( orig )){

                if( nodes_parent_rings.find( neighbor ) == nodes_parent_rings.end() ){
                    for(uint node : *current_ring){
                        if(node == neighbor){
                            neighbor_in_same_ring = true;
                            break;
                        }
                    }
                    if(neighbor_in_same_ring){
                        current = neighbor;
                        break;
                    }
                }
            }
            if( neighbor_in_same_ring ){

                largest_ring->emplace_back( orig );
                largest_ring->emplace_back( current );
                break;
            }
            else{
                current_ring = result->second.second;
                bool neighbor_in_same_ring = false;
                for(uint neighbor : node_adj_lst.at( orig )){

                    if( nodes_parent_rings.find( neighbor ) == nodes_parent_rings.end() ){
                        for(uint node : *current_ring){
                            if(node == neighbor){
                                neighbor_in_same_ring = true;
                                break;
                            }
                        }
                        if(neighbor_in_same_ring){
                            current = neighbor;
                            break;
                        }
                    }
                }
                if( neighbor_in_same_ring ){

                    largest_ring->push_back( orig );
                    largest_ring->push_back( current );
                    break;
                }
            }
        }
    }

    bool find_back_to_orig = true;
    uint last = orig;
    uint depth = 2;
    bool shared_node = true;
    // iterate neighbors of current, picking the right one as next current and repeat until
    // finding back to origin at the right depth
    do{
        bool added = false;
        for( uint neighbor : node_adj_lst.at( current ) ){

            if( neighbor == orig ){
                // This extra indentation is made on purpose
                // to treat the rare exception of a largest ring being smaller than the adjacency list
                // that can occur e.g. with coronene: permutation of all outer excluding the central ring
                // has the same largest ring as the one including the central one.
                // Explanation:
                // Never do anything when orig occurs, so it's not added when there's no other possibility
                // Also avoid orig downstream when !added. This way, by the end of do-while and if still !added, can return.
                if( depth == node_adj_lst.size() ){
                    find_back_to_orig = false;
                    break;
                }
            }
            // prevent return to last
            else if( neighbor != last ){
                // if reach a node that's member of two rings, exit current and move to the other ring
                // -> e.g. in naphthalene, avoid traversing second fused after the first, which would lead to an 8...
                std::map<uint, std::pair< const std::vector<uint>*, const std::vector<uint>* >>::const_iterator result = nodes_parent_rings.find( neighbor );
                if( shared_node && result != nodes_parent_rings.end() ){

                    largest_ring->emplace_back( neighbor );
                    last = current;
                    current = neighbor;

                    added = true;
                    ++depth;
                    shared_node = false;
                    // set current ring to the the other ring (the one that's not current)
                    if( nodes_parent_rings.at( neighbor ).first != current_ring ){
                        current_ring = nodes_parent_rings.at( neighbor ).second;
                    }
                    else{
                        current_ring = nodes_parent_rings.at( neighbor ).first;
                    }

                    break;
                }
                else if( !shared_node && result == nodes_parent_rings.end() ){

                    largest_ring->emplace_back( neighbor );
                    last = current;
                    current = neighbor;

                    // ... but try getting out of a ring as soon as possible after entering
                    shared_node = true;
                    ++depth;
                    added = true;
                    break;
                }
            }
        }
        // This happens when no node turned out to have the required constellation
        // i.e. looking for a shared node and not finding any or unsuccessfully looking for an unshared node
        // in case of first, can add whichever unshared node that's not last is found first
        // (orig at the righ depth would have been caught upstream)
        // in the second case, it must be assured that
        // 1) the shared node is not a member of current ring
        // 2) Say we have phenanthrene rings A-B-C, started at the 'outer bend' node being in A^B, traversed A and reach the 'inner bend' A^B.
        // (1) assures we don't return to orig. If however current ring isn't properly set to C but B instead, we might decide wrongly and continue with B
        // (at least in an even more complex case)
        if( find_back_to_orig && !added ){
            if( !shared_node ){

                const std::vector<uint>* other_ring = nullptr;
                if(nodes_parent_rings.find(current) != nodes_parent_rings.end() ){

                    if(nodes_parent_rings.at( current ).first == current_ring){
                        other_ring = nodes_parent_rings.at( current ).second;
                    }
                    else{
                        other_ring = nodes_parent_rings.at( current ).first;
                    }
                }

                for( uint neighbor : node_adj_lst.at( current ) ){
                    if( neighbor != last && neighbor != orig){
                        bool found_right_neighbor = false;
                        if( nodes_parent_rings.at( neighbor ).first != current_ring &&
                            nodes_parent_rings.at( neighbor ).first != other_ring )
                        {
                            current_ring = nodes_parent_rings.at( neighbor ).first;
                            found_right_neighbor = true;
                        }
                        else if( nodes_parent_rings.at(neighbor).second != current_ring &&
                                 nodes_parent_rings.at(neighbor).second != other_ring )
                        {
                            current_ring = nodes_parent_rings.at( neighbor ).first;
                            found_right_neighbor = true;
                        }
                        if( found_right_neighbor ){
                            largest_ring->emplace_back( neighbor );
                            added = true;
                            last = current;
                            current = neighbor;
                            ++depth;
                            break;
                        }
                    }
                }
            }
            else{
                for( uint neighbor : node_adj_lst.at( current ) ){
                    if( neighbor != last && neighbor != orig ){
                        largest_ring->emplace_back( neighbor );
                        added = true;
                        last = current;
                        current = neighbor;
                        ++depth;
                        break;
                    }
                }
            }
        }
        // If absolutely impossible to add a further node,
        // stop even if largest_ring < node_adj_list
        if(!added){
            find_back_to_orig = false;
        }
    }while( find_back_to_orig );
}

void reader::aromaticity_detection(){

    for( auto& group : m_planarizable ){
        std::vector<const std::vector<uint>*> subgroup;
        subgroup.reserve( group.size() ); // is this really faster for such small vectors?

        for( auto l_ring : group ){

            subgroup.emplace_back( l_ring );

        }
    }
}
//! End under construction



bool reader::has_NCOO(){
    uint node = 0;
    char neighbor_type;
    ushort sp3_het_nbr_count = 0;
    ushort arom_nbr_count = 0;
    for( char atom_type : m_atom_types ){
        if( atom_type == C_1234_val ){
            sp3_het_nbr_count = 0;

            for( uint neighbor : m_adj_list[ node ] ){

                neighbor_type = m_atom_types.at( neighbor );
                if( neighbor_type >= 12 ){
                    ++ sp3_het_nbr_count;
                }
                else if( neighbor_type <= -12 ){
                    ++ arom_nbr_count;
                }
            }
        }
        if( sp3_het_nbr_count >= 2 || (sp3_het_nbr_count > 0 && arom_nbr_count > 0) ){
            return true;
        }
        ++node;
    }
    return false;
}

bool reader::has_terminal_carbonate_carbamate(){
    for(uint node = 0; node < m_atom_types.size(); ++node){
        char node_type = m_atom_types[ node ];
        if( node_type == O_12_val &&
            m_adj_list[node].size() == 1 )
        {
            for( uint nbr : m_adj_list[node] ){

                char nbr_type = m_atom_types[ nbr ];
                if( nbr_type == C_2_bnd_amide_hba ||
                    nbr_type == C_2_bnd_amide_hbd ||
                    nbr_type == C_2_bnd_ester || nbr_type == C_2_bnd_carboxy )
                {
                    char count_O = 0;
                    char count_N = 0;
                    for( uint nbr_of_nbr : m_adj_list[nbr] ){

                        char nbr_of_nbr_type = m_atom_types[ nbr_of_nbr ];

                        if( nbr_of_nbr_type == N_12_val || nbr_of_nbr_type == N_3_val || nbr_of_nbr_type == N_3_val_pos){
                            ++ count_N;
                        }
                        else if( nbr_of_nbr_type == O_2_bnd || nbr_of_nbr_type == O_12_val || nbr_of_nbr_type == O_1_val_neg ){
                            ++ count_O;
                        }
                    }
                    if( count_O >= 3 || ( count_O >= 2 && count_N >= 1 ) ){
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

bool reader::has_anhydride_like(){
    uint count_crit_sp2_nbr = 0;
    for( uint i = 0; i < m_atom_types.size(); ++i ){
        if( m_atom_types[i] == N_12_val ||
            m_atom_types[i] == N_3_val ||
            m_atom_types[i] == N_4_val ||
            m_atom_types[i] == O_12_val )
        {
            for( uint nbr : m_adj_list[i] ){
                if( m_atom_types[nbr] == C_2_bnd_amide_hba ||
                    m_atom_types[nbr] == C_2_bnd_amide_hbd ||
                    m_atom_types[nbr] == C_2_bnd_ester ||
                    m_atom_types[nbr] == C_2_bnd_carboxy )
                {
                    ++ count_crit_sp2_nbr;
                }
            }
            if( count_crit_sp2_nbr >= 2){
                return true;
            }
            count_crit_sp2_nbr = 0;
        }
    }
    return false;
}

bool reader::bad_planarization(){
    uint i = 0;
    for( char type : m_atom_types ){
        if( type < 0 && m_planarizable_nodes.find( i ) == m_planarizable_nodes.end() ){
            return true;
        }
        ++ i;
    }
    return false;
}

void reader::find_fix_amides(){

    for(uint node = 0; node < m_adj_list.size(); ++node){
        bool potential_hbd = false;

        if( m_atom_types[ node ] == 1 && m_adj_list[node].size() == 3 ){
            bool adj_terminal_oxygen = false;
            uint double_bond_id = 0;
            ushort adj_nitrogen = 0;
            for( uint neighbor : m_adj_list[node] ){

                if( m_atom_types[ neighbor ] == O_12_val && m_adj_list[neighbor].size() == 1 ){
                    adj_terminal_oxygen = true;
                    double_bond_id = neighbor;
                }

                else if( m_atom_types[ neighbor ] == N_12_val ){

                    ++adj_nitrogen;
                    potential_hbd = true;
                }
                else if( m_atom_types[ neighbor ] == N_3_val  ){
                    ++adj_nitrogen;
                }
            }

            if( adj_terminal_oxygen && adj_nitrogen > 0 ){
                m_atom_types[ double_bond_id ] = O_2_bnd;
                if( potential_hbd ){
                    m_atom_types[ node ] = C_2_bnd_amide_hbd;
                }
                else{
                    m_atom_types.at( node ) = C_2_bnd_amide_hba;
                }
            }
        }
    }
}

void reader::find_fix_carboxy(){

    if( m_atom_types[4] == N_3_val && m_atom_types[7] == O_12_val && m_atom_types[5] == O_12_val ){
        uint test = 0;
    }

    for(uint node = 0; node < m_adj_list.size(); ++node){

        if( m_atom_types[ node ] == C_1234_val && m_adj_list[node].size() == 3 ){

            bool adj_terminal_oxygen = false;
            uint double_bond_id = 0;
            ushort adj_oxygen = 0;
            ushort adj_O_hac_nbrs = 1;

            for( uint neighbor : m_adj_list[node] ){
                if( m_atom_types[ neighbor ] == O_12_val ){

                    if( m_adj_list[neighbor].size() == 1 ){
                        adj_terminal_oxygen = true;
                        double_bond_id = neighbor;
                    }
                    else{
                        ++adj_O_hac_nbrs;
                    }
                    ++adj_oxygen;
                }
                else if(m_atom_types[ neighbor ] >=12 && !is_planarizable_node(neighbor)){
                    adj_terminal_oxygen = false;
                    break;
                }
            }

            if( adj_terminal_oxygen && adj_oxygen >= 2 ){
                if( adj_O_hac_nbrs == 1 ){
                    m_atom_types[ double_bond_id ] = O_2_bnd;
                    m_atom_types[ node ] = C_2_bnd_carboxy;
                }
                else{
                    m_atom_types[ double_bond_id ] = O_2_bnd;
                    m_atom_types[ node ] = C_2_bnd_ester;
                }
            }
        }
    }
}

bool reader::find_fix_guanidines_amidines(){
    bool new_guanidine_amidine = false;
    bool neighbors_planarizable;
    for(uint node = 0; node < m_adj_list.size(); ++node){

        if( m_atom_types[ node ] == 1 && m_adj_list[node].size() == 3 && is_planarizable_node(node) ){
            neighbors_planarizable = true;
            ushort lowest_coval_nbr = 0;
            ushort lowest_coval = 4;
            ushort highest_coval = 0;
            ushort adj_123_val_nitrogen = 0;
            ushort adj_arom_nitrogen = 0;
            ushort adj_arom_carbon = 0;

            for( uint neighbor : m_adj_list[node] ){

                if(!is_planarizable_node(neighbor)){
                    neighbors_planarizable = false;
                    break;
                }
                if( m_atom_types.at(neighbor) == N_12_val || m_atom_types[ neighbor ] == N_3_val ){
                    ++adj_123_val_nitrogen;
                    uint current_coval = m_adj_list[neighbor].size();
                    if(current_coval < lowest_coval){
                        lowest_coval_nbr = neighbor;
                        lowest_coval = current_coval;
                    }
                    if( current_coval > highest_coval ){
                        highest_coval = current_coval;
                    }
                    if(current_coval == 3 && !is_planarizable_node(neighbor)){
                        break;
                    }
                }

                else if( m_atom_types[ neighbor ] == N_2_bnd ){
                    ++adj_123_val_nitrogen;
                }

                else if( m_atom_types[ neighbor ] == n_2_val || m_atom_types[ neighbor ] == n_3_val || m_atom_types[ neighbor ] == n_3_val_pos ){
                    ++adj_123_val_nitrogen;
                    ++adj_arom_nitrogen;
                }
                else if( m_atom_types[ neighbor ] < 0 && m_atom_types[ neighbor ] >= -11 ){
                    ++adj_arom_carbon;
                }
                else{
                    break;
                }
            }

            if( adj_123_val_nitrogen == 3 && neighbors_planarizable ){

                if( lowest_coval < 3 && adj_arom_nitrogen <= 1 ){
                    m_atom_types[ lowest_coval_nbr ] = N_2_bnd;
                    if( lowest_coval <= 2 ){
                        m_atom_types[ node ] = C_2_bnd_guanidine_hbd;
                    }
                    else{
                        m_atom_types[ node ] = C_2_bnd_guanidine_hba;
                    }
                    new_guanidine_amidine = true;
                }
                else{
                    m_atom_types[ lowest_coval_nbr ] = N_3_val_pos;
                    m_atom_types[ node ] = C_2_bnd_guanidine;

                    new_guanidine_amidine = true;
                }
            }
            else if( adj_123_val_nitrogen == 2 && adj_arom_carbon == 1 && lowest_coval < 3 && neighbors_planarizable ){
                m_atom_types[ lowest_coval_nbr ] = N_2_bnd;
                m_atom_types[ node ] = C_2_bnd_arylamidine;

                new_guanidine_amidine = true;
            }
        }
    }
    return new_guanidine_amidine;
}

void reader::find_carbonyl(){
    for( uint i = 0; i < m_adj_list.size(); ++i ){
        uint bond_order = 1;

        bond_order += m_valence[i] - m_adj_list[i].size();
        uint o_count = 0;
        uint n_count = 0;
        uint o_pos = 0;

        if( bond_order == 2 && m_atomic_numbers[i] == 6 ){
            uint j = 0;
            for( uint neighbor : m_adj_list[i] ){
                if( m_atomic_numbers[neighbor] == 7 ){
                    ++ n_count;
                }
                else if( m_atomic_numbers[neighbor] == 8 ){
                    ++ o_count;
                    o_pos = j;
                }
                ++ j;
            }
        }
        if( n_count == 0 && o_count == 1 ){
            m_atom_types[i] = C_2_bnd_keto;
            m_atom_types[m_adj_list[i][o_pos]] = O_2_bnd;
        }
    }
}

void reader::find_fix_S_P(){
    for( uint i = 0; i < m_adj_list.size(); ++i ){
        if( m_atomic_numbers[i] == 16 ){
            uint neighbor_count = m_adj_list[i].size();
            if( neighbor_count == 1 && m_valence[i] > 1 ){
                m_atom_types[i] = S_2_bnd;
                m_atom_types[ m_adj_list[i][0] ] =  C_2_bnd_thiamide;
            }
            else if( neighbor_count <= 2 ){
                m_atom_types[i] = S_12_val;
            }
            else if( neighbor_count == 3 ){
                m_atom_types[i] = S_2_bnd;
            }
            else if( neighbor_count == 4 ){
                m_atom_types[i] = S_2x2_bnd_6_val;
            }
            else if( neighbor_count == 6 ){
                m_atom_types[i] = S_6_val;
            }
        }
        else if( m_atomic_numbers[i] == 15 ){
            uint neighbor_count = m_adj_list[i].size();
            if( neighbor_count == 4 ){
                m_atom_types[i] = P_1x2_bnd_5_val;
            }
            else if( neighbor_count == 5 ){
                m_atom_types[i] = P_5_val;
            }
            else if( neighbor_count == 3 ){
                m_atom_types[i] = P_3_val;
            }
        }
    }
}
