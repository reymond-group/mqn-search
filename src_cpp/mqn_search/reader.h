#ifndef READER_H
#define READER_H

#include <vector>
#include <deque>
#include <string>
#include <map>
#include <set>
#include <unordered_set>

class reader
{
    typedef unsigned int uint;

public:
    reader();
    ~reader();

    void read( const std::string SMILES );
    const std::vector<uint>& get_mqn();


    uint count_size_3_rings();
    uint count_size_4_rings();
    uint count_size_5_rings();
    uint count_size_6_rings();
    uint count_size_7_rings();
    uint count_size_8_rings();
    uint count_size_9_rings();
    uint count_macrocycles();

    uint ring_system_count();
    uint get_sssr_size();

    uint count_monovalent();
    uint count_divalent();
    uint count_trivalent();
    uint count_tetravalent();

    uint cyclic_node_count();
    uint cyclic_edge_count();

    float fraction_cyclic_edges();
    uint bonds_in_min_2_rings();
    uint atoms_in_min_2_rings();

    float molecular_weight();
    uint heavy_atom_count();
    uint rotatable_bonds();

    uint heteroatom_count();
    float fraction_aromatic();

    uint carbon_count();
    uint nitrogen_count();
    uint oxygen_count();
    uint sulfur_count();
    uint fluorine_count();
    uint chlorine_count();
    uint bromine_count();
    uint phosphorous_count();
    uint iodine_count();
    uint acyclic_n_count();
    uint cyclic_n_count();
    uint acyclic_o_count();
    uint cyclic_o_count();

    void cleanup();

private:

    static const std::map<char, std::unordered_set<char>> cnos;
    static const std::map<std::string, uint> amass;
    static const std::map<std::string, uint> anum;

    std::string m_SMILES;
    std::vector<std::vector<uint>> m_adj_list;
    std::vector<std::string> m_atom_labels;
    std::vector<char> m_atom_types;
    std::vector<uint> m_valence;
    std::map<std::pair<uint,uint>, uint> m_bond_orders;
    std::vector<int> m_charge;
    std::vector<uint> m_attached_H;
    std::vector<uint> m_atomic_masses;
    std::vector<uint> m_atomic_numbers;
    std::vector<uint> m_ring_membership;
    std::vector<uint> m_nestedness;
    std::set<uint> m_aromatic_H;

    bool is_valid_atom_id(const uint i);
    void get_nestedness();
    bool look_ahead(const uint i, std::vector<std::vector<uint>>& connectivity, std::vector<uint>& rid_count );
    void gen_adj_list();
    void count_implicit_H();
    void correct_valence();
    void assign_mass_number();
    void get_bond_orders();

    void single_bonds();
    void double_bonds();
    void triple_bonds();

    void node_order();
    void charge();
    void h_bond();

    struct ring_comparator {
        bool operator()(const std::vector<uint>* const v1, const std::vector<uint>* const v2) const {
            return *v1 < *v2;
        }
    };

    uint m_orig; // origin node of ring traversal
    std::vector<uint> m_temp_path; // traversal path
    std::vector<std::vector<std::vector<uint> > > m_rings; // size containers

    std::map<const std::vector<uint>*,std::vector<const std::vector<uint>*>, ring_comparator > m_bulk_ring_adj_lst;
    std::map<const std::vector<uint>*, std::vector<const std::vector<uint>*>*> m_ring_sys_membership;
    std::vector<std::vector<const std::vector<uint>*>> m_sys_size_ring;
    std::vector<std::vector<const std::vector<uint>*>> m_sys_size_ring_23;
    std::deque<const std::vector<uint>* > m_topol_SSSR;
    std::deque<const std::vector<uint>* > m_symm_SSSR;
    std::deque<const std::vector<uint>* > m_symm_SSSR_23;
    std::map<const std::vector<uint>*,std::vector<const std::vector<uint>*> > m_ring_adj_lst;
    std::map<const std::vector<uint>*, std::vector<std::vector<uint> > > m_shared_nodes_lst;
    std::map<const std::vector<uint>*, std::set<uint> > m_bridging_nodes;
    std::vector<std::vector<const std::vector<uint> *> > m_ring_systems;

    //! Under construction
    std::vector<std::vector<uint>> m_largest_rings;
    //! End under construction

    std::set<uint> m_cyclic_nodes;
    uint m_cyclic_edge_count;
    uint m_edges;
    uint m_bonds_in_min_2_rings;
    uint m_atoms_in_min_2_rings;

    int m_cyclic_single_bonds;
    int m_acyclic_single_bonds;
    int m_cyclic_double_bonds;
    int m_acyclic_double_bonds;
    int m_cyclic_triple_bonds;
    int m_acyclic_triple_bonds;

    int m_monovalent;
    int m_cyclic_divalent;
    int m_acyclic_divalent;
    int m_cyclic_trivalent;
    int m_acyclic_trivalent;
    int m_cyclic_tetravalent;
    int m_acyclic_tetravalent;

    int m_pos_charges;
    int m_neg_charges;

    int m_hba;
    int m_hba_sites;
    int m_hbd;
    int m_hbd_sites;

    uint neighbor_count_of(uint node);
    bool append_if_new_ring(std::vector<std::vector<uint> >& size_container, const std::vector<uint> qry_ring );
    bool node_in_range(std::vector<uint>::const_iterator begin, std::vector<uint>::const_iterator end, uint lookup_node);
    uint ring_traverse(uint current, uint depth);
    bool find_rings();
    void gen_ring_adj_lst();

    void gen_bulk_ring_adj_lst();
    const std::vector<uint>* bulk_ring_system_traversal(const std::vector<uint>* current, std::vector<const std::vector<uint> *>* const current_system, const std::vector<uint>* last);
    void ring_systems();
    void ring_system_traversal_23(const std::vector<uint>* current, std::vector<const std::vector<uint> *>& current_system);
    void sssr();

    bool in_di_tri_val_rings(const std::vector<uint>* lookup_ring);
    bool ring_in_set(std::vector<const std::vector<uint>*> ring_set, const std::vector<uint>* lookup_ring);
    const std::vector<uint>* spheroid_traversal(const std::vector<uint>* current, std::vector<const std::vector<uint> *>& forbidden, uint depth);

    uint shared_with_two_rings(const std::vector<uint> * first, const std::vector<uint> * second, uint idx_2nd);
    void extract_planarizable_subgraph( const std::vector<uint> * const smallest_ring );
    void spheroid_system_detection();
    uint ring_system_traversal(const std::vector<uint>* current, std::vector<const std::vector<uint> *>& current_system, uint depth);
    void node_ring_membership();
    bool is_planarizable_node(uint query_node);

    std::vector<std::vector<const std::vector<uint>*>> m_planarizable;
    uint m_triple_shared;
    std::set<uint> m_planarizable_nodes;
    std::vector<std::vector<const std::vector<uint>*> > m_node_ring_membership;

    //! Under construction
    void largest_ring_of_group(const std::map<uint, std::set<uint> >& node_adj_lst, std::map<uint, std::pair<const std::vector<uint>*, const std::vector<uint>*>>& nodes_parent_rings/*, const std::set<uint>& bridging_nodes*/);
    void aromaticity_detection();
    //! End under construction

    void init_pattern();

    bool has_NCOO();
    bool has_terminal_carbonate_carbamate();
    bool has_anhydride_like();

    void find_fix_amides();
    void find_fix_carboxy();
    bool find_fix_guanidines_amidines();
    void find_carbonyl();

    void find_fix_S_P();

    bool bad_planarization();

    std::vector<uint> m_mqn;
};

#endif // READER_H
