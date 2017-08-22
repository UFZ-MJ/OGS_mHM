/*
 * filterMeshNodes.cpp
 *
 *  Created on: Jan 3, 2011
 *      Author: TF
 */
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

// Base
#include "quicksort.h"
#include "binarySearch.h"

int main (int argc, char *argv[])
{
	if (argc < 7) {
		std::cout << "Usage: " << argv[0] << " --nodes-with-area nodes_with_area_file --nodes-in-polygon nodes_in_polygon --output output_file_name" << std::endl;
		return -1;
	}

	std::string tmp (argv[1]);
	if (tmp.find ("--nodes-with-area") == std::string::npos) {
		std::cout << "could not find switch nodes_with_area" << std::endl;
		return -1;
	}

	tmp = argv[3];
	if (tmp.find ("--nodes-in-polygon") == std::string::npos) {
		std::cout << "could not find switch nodes-in-polygon" << std::endl;
		return -1;
	}

	tmp = argv[5];
	if (tmp.find ("--output") == std::string::npos) {
		std::cout << "missing switch output" << std::endl;
		return -1;
	}

	tmp = argv[2];
	std::ifstream in0 (tmp.c_str());
	tmp = argv[4];
	std::ifstream in1 (tmp.c_str());
	if (!in0 || !in1) {
		std::cout << "could not open one of the input files" << std::endl;
	}

	std::vector<size_t> ids0, ids1;
	std::vector<double> areas;

	size_t tmp_id;
	double tmp_area_val;

	bool all_ok (true);
	while (!in0.eof() && all_ok) {
		in0 >> tmp_id;
		if (in0.fail()) {
			all_ok = false;
		}
		if (all_ok) {
			in0 >> tmp_area_val;
			if (in0.fail ())
				all_ok = false;
			else {
				ids0.push_back (tmp_id);
				areas.push_back (tmp_area_val);
			}
		}
	}
	in0.close ();
	std::cout << "read "<< ids0.size() << " ids and " << areas.size() << " areas from file " << argv[2] << std::endl;

	all_ok = true;
	while (!in1.eof() && all_ok) {
		in1 >> tmp_id;
		if (! in1.fail()) {
			ids1.push_back (tmp_id);
		} else
			all_ok = false;
	}
	in1.close ();
	std::cout << "read "<< ids1.size() << " ids from file " << argv[4] << std::endl;

	std::vector<size_t> perm;
	for (size_t k(0); k<ids0.size(); k++) {
		perm.push_back(k);
	}
	Quicksort<size_t> (ids0, 0, ids0.size(), perm);

	// apply permutation to areas
	std::vector<double> sorted_areas (areas.size());
	for (size_t k(0); k<perm.size(); k++) {
		sorted_areas[k] = areas[perm[k]];
	}

	std::vector<size_t> not_found_ids;
	std::ofstream out (argv[6]);
	if (!out) {
		std::cout << "could not open file " << argv[6] << " for output" << std::endl;
		return -1;
	}
	out.precision (12);
	for (size_t k(0); k<ids1.size(); k++) {
		size_t idx (searchElement (ids1[k], 0, ids0.size(), ids0));
		if (idx != std::numeric_limits<size_t>::max()) {
			out << ids1[k] << " " << sorted_areas[idx] << std::endl;
		} else {
			 not_found_ids.push_back (ids1[k]);
		}
	}
	out.close ();

	if (!not_found_ids.empty()) {
		std::sort (not_found_ids.begin(), not_found_ids.end());
		for (size_t k(0); k<not_found_ids.size(); k++) {
			std::cout << not_found_ids[k] << std::endl;
		}
		std::cout << "number of not found ids: " << not_found_ids.size() << std::endl;
	}
}

