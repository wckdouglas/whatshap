#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>
#include <unordered_set>

#include "pedigreecolumncostcomputer.h"
#include "pedigreedptable.h"

using namespace std;

PedigreeDPTable::PedigreeDPTable(ReadSet* read_set, const vector<unsigned int>& recombcost, const Pedigree* pedigree, unsigned int ploidy, bool distrust_genotypes, const vector<unsigned int>* allele_counts, const vector<unsigned int>* positions, const vector<unsigned int>* precomputed_partitioning) :
	read_set(read_set),
	recombcost(recombcost),
	pedigree(pedigree),
	ploidy(ploidy),
	distrust_genotypes(distrust_genotypes),
	optimal_score(0u),
	optimal_score_index(0u),
	input_column_iterator(*read_set, positions),
	allele_counts(allele_counts)
{
	read_set->reassignReadIds();

	// create all pedigree partitions
	for (size_t i=0; i<std::pow(4, pedigree->triple_count()); ++i) {
		pedigree_partitions.push_back(new PedigreePartitions(*pedigree, i, ploidy));
	}

	// translate all individual ids to individual indices
	for (size_t i=0; i<read_set->size(); ++i) {
		read_sources.push_back(pedigree->id_to_index(read_set->get(i)->getSampleID()));
	}

	if (precomputed_partitioning == nullptr){
		compute_table();
	} else {
		// map read_id -> assigned partition
		for (unsigned int i = 0; i < read_set->size(); i++) {
			read_to_partition[ read_set->get(i)->getID() ] = (*precomputed_partitioning)[i];
		}
		set_index_path();
	}
}


PedigreeDPTable::~PedigreeDPTable() {
	init(projection_column_table, 0);
	init(index_backtrace_table, 0);
	init(transmission_backtrace_table, 0);
	init(indexers, 0);
	init(pedigree_partitions, 0);
	delete allele_counts;
}


unique_ptr<vector<unsigned int> > PedigreeDPTable::extract_read_ids(const vector<const Entry *>& entries) {
	unique_ptr<vector<unsigned int> > read_ids(new vector<unsigned int>());
	for (size_t i=0; i<entries.size(); ++i) {
		read_ids->push_back(entries[i]->get_read_id());
	}
	return read_ids;
}


size_t PedigreeDPTable::popcount(size_t x) {
	unsigned int count = 0;
	for (;x; x >>= 1) {
		count += x & 1;
	}
	return count;
}


void PedigreeDPTable::clear_table() {
	size_t column_count = input_column_iterator.get_column_count();

	init(projection_column_table, column_count);
	init(index_backtrace_table, column_count);
	init(transmission_backtrace_table, column_count);
	init(indexers, column_count);

	index_path.clear();

	optimal_score = numeric_limits<unsigned int>::max();
	optimal_score_index = 0;
	optimal_transmission_value = 0;
	previous_transmission_value = 0;
}


void PedigreeDPTable::set_index_path() {
	clear_table();

	optimal_score = 0;
	optimal_score_index = 0;
	// empty read-set, nothing to phase, MEC is 0
	if (input_column_iterator.get_column_count() == 0) {
		return;
	}

	// store partitioning for each column
	input_column_iterator.jump_to_column(0);
	unique_ptr<vector<const Entry *> > current_input_column = input_column_iterator.get_next();
	unique_ptr<vector<unsigned int> > current_read_ids = extract_read_ids(*current_input_column);

cout << "HERE 1" << endl;

	size_t column_count = input_column_iterator.get_column_count();
	index_path.assign(column_count, index_and_inheritance_t());

cout << "HERE 2" << endl;

	for (size_t column_index = 0; column_index < column_count; ++column_index) {
		// construct partitioning of reads in this column
		unsigned int partitioning = 0;
		unsigned int offset = 1;
		for (size_t i = 0; i < current_read_ids->size(); i++) {
			unsigned int id = (*current_read_ids)[i];
			partitioning += read_to_partition[id]*offset;
			offset *= ploidy;
		}

		// store in index path 
		index_and_inheritance_t v;
		v.index = partitioning;
		v.inheritance_value = 0;

		index_path[column_index] = v;

		// compute optimal cost
		PedigreeColumnCostComputer cost_computer(*current_input_column, column_index, read_sources, pedigree, *pedigree_partitions[v.inheritance_value], distrust_genotypes, allele_counts->at(column_index));
		cost_computer.set_partitioning(v.index);
		optimal_score += cost_computer.get_cost();

		if (input_column_iterator.has_next()){
			current_input_column = input_column_iterator.get_next();
			current_read_ids = extract_read_ids(*current_input_column);
		}
	}
}


void PedigreeDPTable::compute_table() {
	clear_table();

	// empty read-set, nothing to phase, so MEC score is 0
	if (input_column_iterator.get_column_count() == 0) {
		optimal_score = 0;
		optimal_score_index = 0;
		return;
	}

	input_column_iterator.jump_to_column(0);
	unique_ptr<vector<const Entry *> > current_input_column;
	unique_ptr<vector<const Entry *> > next_input_column;
	// get the next column ahead of time
	next_input_column = input_column_iterator.get_next();
	unique_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_input_column);
	ColumnIndexingScheme* next_indexer = new ColumnIndexingScheme(0, *next_read_ids, ploidy);
	indexers[0] = next_indexer;

	// forward pass: create a sparse table, storing values at every sqrt(#columns)-th position,
	size_t k = (size_t)sqrt(input_column_iterator.get_column_count());
	for (size_t column_index=0; column_index<input_column_iterator.get_column_count(); ++column_index) {
		// make former next column the current one
		current_input_column = std::move(next_input_column);
		unique_ptr<vector<unsigned int> > current_read_ids = std::move(next_read_ids);

		ColumnIndexingScheme* current_indexer = next_indexer;
		// peek ahead and get the next column
		if (input_column_iterator.has_next()) {
			next_input_column = input_column_iterator.get_next();
			next_read_ids = extract_read_ids(*next_input_column);
			next_indexer = new ColumnIndexingScheme(current_indexer,*next_read_ids,ploidy);
			current_indexer->set_next_column(next_indexer);
			indexers[column_index + 1] = next_indexer;
		} else {
			assert(next_input_column.get() == 0);
			assert(next_read_ids.get() == 0);
			next_indexer = 0;
		}

		compute_column(column_index, std::move(current_input_column));

		// determine whether to delete previous column (to save space)
		if ((k>1) && (column_index > 0) && (((column_index-1)%k) != 0)) {
			delete index_backtrace_table[column_index-1];
			delete transmission_backtrace_table[column_index-1];
			delete projection_column_table[column_index-1];
			index_backtrace_table[column_index-1] = nullptr;
			transmission_backtrace_table[column_index-1] = nullptr;
			projection_column_table[column_index-1] = nullptr;
		}
	}

	// perform a backtrace to get optimal path
	index_path.assign(indexers.size(), index_and_inheritance_t());
	index_and_inheritance_t v;
	unsigned int prev_inheritance_value = previous_transmission_value;
	v.index = optimal_score_index;
	v.inheritance_value = optimal_transmission_value;
	index_path[indexers.size()-1] = v;

//	// remember when to break
//	bool set_break = false;

	for(size_t i = indexers.size()-1; i > 0; --i) { // backtrack through table
		// ensure that index_backtrace_table[i-1] and transmission_backtrace_table[i-1] exist
		if (projection_column_table[i-1] == nullptr) {
			// compute index of last previous column that has been stored
			size_t j = (i-1) / k * k;
			assert(projection_column_table[j] != nullptr);
			for (j=j+1; j<i; ++j) {
				compute_column(j);
			}
		}
		// compute index and transmission value for the current column
		unique_ptr<ColumnIndexingIterator> iterator = indexers[i]->get_iterator();
		unsigned int backtrace_index = iterator->index_backward_projection(v.index);

/**		// determine how many sets are empty in each column	
		unsigned int nonempty_backward = empty_partitions(v.index, indexers[i]->get_read_ids()->size());
		unsigned int nonempty_projection = empty_partitions(backtrace_index, indexers[i]->get_backward_projection_width());
		unsigned int nonempty_forward = empty_partitions(index_backtrace_table[i-1]->at(backtrace_index, prev_inheritance_value), indexers[i-1]->get_read_ids()->size());

		// if start of a connected component, skip
		if (nonempty_projection == ploidy){
			set_break = false;
			continue;
		}

		// determine whether to break haplotype apart
		if (nonempty_projection > 1){
			if (nonempty_backward < nonempty_projection){
				set_break = true;
			}
			if (nonempty_forward < nonempty_projection){
				if(set_break) block_boundaries.push_back(i);
			}
		}
**/
		v.index = index_backtrace_table[i-1]->at(backtrace_index, prev_inheritance_value);
		v.inheritance_value = prev_inheritance_value;
		prev_inheritance_value = transmission_backtrace_table[i-1]->at(backtrace_index, v.inheritance_value);
		index_path[i-1] = v;
		// free parts of the DP table no longer needed
		if (i%k == 0) {
			for (size_t j=i; (j<i+k) && (j<input_column_iterator.get_column_count()-1); ++j) {
				assert(projection_column_table[j] != nullptr);
				delete index_backtrace_table[j];
				delete transmission_backtrace_table[j];
				delete projection_column_table[j];
				index_backtrace_table[j] = nullptr;
				transmission_backtrace_table[j] = nullptr;
				projection_column_table[j] = nullptr;
			}
		}
	}
}

/**
unsigned int PedigreeDPTable::empty_partitions(unsigned int index, unsigned int size){
	unordered_set<unsigned int> elements;
	while(size > 0) {
		elements.insert(index % ploidy);
		index /= ploidy;
		size -= 1;
	}

	return ploidy-elements.size();
}
**/

void PedigreeDPTable::compute_column(size_t column_index, unique_ptr<vector<const Entry*>> current_input_column) {
	assert(column_index < input_column_iterator.get_column_count());

	// check whether requested column is already there
	if (projection_column_table[column_index] != nullptr) {
		assert(index_backtrace_table[column_index] != nullptr);
		assert(transmission_backtrace_table[column_index] != nullptr);
		return;
	}

	ColumnIndexingScheme* current_indexer = indexers[column_index];
	assert(current_indexer != nullptr);

	// compute the number of different transmission vectors
	unsigned int transmission_configurations = std::pow(4, pedigree->triple_count());

	// if current input column was not provided, then create it
	if (current_input_column.get() == nullptr) {
		input_column_iterator.jump_to_column(column_index);
		current_input_column = input_column_iterator.get_next();
	}

	// reserve memory for the current DP column
	Vector2D<unsigned int> dp_column(current_indexer->column_size(), transmission_configurations, 0);

	// obtain previous projection column (which is assumed to have been already computed)
	Vector2D<unsigned int>* previous_projection_column = nullptr;
	if (column_index > 0) {
		previous_projection_column = projection_column_table[column_index - 1];
	}

	// initialize forward projection column and associated backtrace columns,
	// if existing (i.e. if not last column)
	Vector2D<unsigned int>* current_projection_column = nullptr;
	Vector2D<unsigned int>* transmission_backtrace_column = nullptr;
	Vector2D<unsigned int>* index_backtrace_column = nullptr;
	if (column_index + 1 < input_column_iterator.get_column_count()) {
		current_projection_column = new Vector2D<unsigned int>(
			current_indexer->forward_projection_size(),
			transmission_configurations,
			numeric_limits<unsigned int>::max()
		);
		transmission_backtrace_column = new Vector2D<unsigned int>(
			current_indexer->forward_projection_size(),
			transmission_configurations,
			numeric_limits<unsigned int>::max()
		);
		index_backtrace_column = new Vector2D<unsigned int>(
			current_indexer->forward_projection_size(),
			transmission_configurations,
			numeric_limits<unsigned int>::max()
		);
	}

	// create column cost computers
	vector<PedigreeColumnCostComputer> cost_computers;
	cost_computers.reserve(transmission_configurations);

	for(unsigned int i = 0; i < transmission_configurations; ++i) {
		cost_computers.emplace_back(*current_input_column, column_index, read_sources, pedigree, *pedigree_partitions[i], distrust_genotypes, allele_counts->at(column_index));
	}

	// iterate over all bipartitions
	unique_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator();
	while (iterator->has_next()) {
		int bit_changed = -1;
		int partition_changed = -1;
		iterator->advance(&bit_changed, &partition_changed);
		if (bit_changed >= 0) {
			for(auto& cost_computer : cost_computers) {
				cost_computer.update_partitioning(bit_changed, partition_changed);
			}
		} else {
			for(auto& cost_computer : cost_computers) {
				cost_computer.set_partitioning(iterator->get_partition());
			}
		}

		// Determine index in backward projection column from where to fetch the previous cost
		size_t backward_projection_index = 0;
		if (column_index > 0) {
			backward_projection_index = iterator->get_backward_projection();
		}
		// Determine index in the current DP column to be written
		size_t current_index = iterator->get_index();

		// Compute aggregate cost based on cost in previous and cost in current column
		vector<unsigned int> min_recomb_index(transmission_configurations);
		bool found_valid_transmission_vector = false;
		for (size_t i = 0; i < transmission_configurations; ++i) {
			// Compute cost incurred by current cell of DP table
			unsigned int current_cost = cost_computers[i].get_cost();
			unsigned int min = numeric_limits<unsigned int>::max();
			size_t min_index = 0;
			if (current_cost < numeric_limits<unsigned int>::max()) {
				found_valid_transmission_vector = true;
			}
			for (size_t j = 0; j < transmission_configurations; ++j) {
				// Step 1: add up cost from current_cost column and previous columns
				unsigned int val;
				unsigned int previous_cost = 0;
				if (column_index > 0) {
					previous_cost = previous_projection_column->at(backward_projection_index,j);
				}
				if ((current_cost < numeric_limits<unsigned int>::max()) && (previous_cost < numeric_limits<unsigned int>::max())) {
					val = current_cost + previous_cost;
				} else {
					val = numeric_limits<unsigned int>::max();
				}
				// Step 2: add further cost incurred by recombination
				// change in bit 0 --> recombination in mother
				size_t x = i ^ j; // count the number of bits set in x

				if (val < numeric_limits<unsigned int>::max()) {
					val += popcount(x) * recombcost[column_index];
				}

				// check for new minimum
				if (val < min) {
					min = val;
					min_index = j;
				}
			}
			dp_column.set(current_index, i, min);
			min_recomb_index[i] = min_index;
		}
		if (!found_valid_transmission_vector) {
			throw std::runtime_error("Error: Mendelian conflict");
		}

		// if last DP column, then check for new optimal score, otherwise update forward projection and backtrace columns
		if (current_projection_column == 0) {
			// update running optimal score index
			for (size_t i = 0; i < transmission_configurations; ++i) {
				if (dp_column.at(current_index, i) < optimal_score) {
					optimal_score = dp_column.at(current_index, i);
					optimal_score_index = iterator->get_index();
					optimal_transmission_value = i;
					previous_transmission_value = min_recomb_index[i];
				}
			}
		} else {
			unsigned int forward_index = iterator->get_forward_projection();
			unsigned int it_idx = iterator->get_index();
			for (unsigned int i = 0; i < transmission_configurations; ++i) {
				if (dp_column.at(current_index, i) < current_projection_column->at(forward_index,i)) {
					current_projection_column->set(forward_index, i, dp_column.at(current_index, i));
					index_backtrace_column->set(forward_index, i, it_idx);
					transmission_backtrace_column->set(forward_index,i, min_recomb_index[i]);
				}
			}
		}
	}

	// if not last column, then store computed tables
	if (current_projection_column != 0) {
		index_backtrace_table[column_index] = index_backtrace_column;
		transmission_backtrace_table[column_index] = transmission_backtrace_column;
		projection_column_table[column_index] = current_projection_column;
	}
}


unsigned int PedigreeDPTable::get_optimal_score() {
	//if (backtrace_table.empty()) throw runtime_error("Empty backtrace table");
	return optimal_score;
}

/**
vector<unsigned int>* PedigreeDPTable::get_block_boundaries(){
	cout <<  "BLOCK BOUNDARIES: " << block_boundaries.size() << endl;
	vector<unsigned int>* output = new vector<unsigned int>(block_boundaries.rbegin(), block_boundaries.rend());
	return output;
}
**/

void PedigreeDPTable::get_super_reads(std::vector<ReadSet*>* output_read_set, vector<unsigned int>* transmission_vector) {
	assert(output_read_set != nullptr);
	assert(output_read_set->size() == pedigree->size());
	assert(transmission_vector != nullptr);
	transmission_vector->clear();

	input_column_iterator.jump_to_column(0);
	const vector<unsigned int>* positions = input_column_iterator.get_positions();

	std::vector<std::vector<Read*>> superreads;
	for (unsigned int i = 0; i < pedigree->size(); ++i) {
		std::vector<Read*> reads_per_individual;
		for (unsigned int j = 0; j < ploidy; ++j) {
			reads_per_individual.push_back(new Read("superread_"+std::to_string(j)+"_"+std::to_string(i), -1, -1, pedigree->index_to_id(i)));
		}
		superreads.emplace_back(reads_per_individual);
	}

	if (index_backtrace_table.empty()) {
		assert(!input_column_iterator.has_next());
	} else {
		// run through the file again with the input_column_iterator
		unsigned int i = 0; // column index
		while (input_column_iterator.has_next()) {
			const index_and_inheritance_t& v = index_path[i];
			unique_ptr<vector<const Entry *> > column = input_column_iterator.get_next();
			PedigreeColumnCostComputer cost_computer(*column, i, read_sources, pedigree, *pedigree_partitions[v.inheritance_value], distrust_genotypes, allele_counts->at(i));
			cost_computer.set_partitioning(v.index);

			auto population_alleles = cost_computer.get_alleles();
			int n_alleles = allele_counts->at(i);
			
			// TODO: compute proper weights based on likelihoods
			for (unsigned int k = 0; k < pedigree->size(); k++){
				for (unsigned int j = 0; j < ploidy; j++) {
					std::vector<unsigned int> qualities(n_alleles, 10);
					for(int a = 0; a < n_alleles; a++){
						if (a == population_alleles[k].alleles[j]) qualities[a] = 0;
					}
					superreads[k][j]->addVariant(positions->at(i), population_alleles[k].alleles[j], qualities);
				}
			}

			transmission_vector->push_back(v.inheritance_value);
			++i; // next column
		}
	}
	for (unsigned int k = 0; k < pedigree->size(); k++) {
		assert(output_read_set->at(k) != nullptr);
		for (unsigned int j = 0; j < ploidy; j++) {
			output_read_set->at(k)->add(superreads[k][j]);
		}
	}
}

vector<unsigned int>* PedigreeDPTable::get_optimal_partitioning() {
	vector<unsigned int>* partitioning = new vector<unsigned int>(read_set->size(), false);
	for(size_t i = 0; i < index_path.size(); ++i){
		unsigned int index = index_path[i].index;
		for(size_t j = 0; j < indexers[i]->get_read_ids()->size(); ++j){
			unsigned int partition = index % ploidy;
			partitioning->at(indexers[i]->get_read_ids()->at(j)) = partition;
			index /= ploidy;
		}
	}
	return partitioning;
}
