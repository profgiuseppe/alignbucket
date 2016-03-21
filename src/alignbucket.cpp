/**
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation version 2..

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

	Authors:
	- Giuseppe Profiti, <gamma2@users.sourceforge.net>
	- Piero Fariselli, <piero@biocomp.unibo.it>
	*/
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>
#include <list>
#include <utility>
#include <string>
#include <cstdlib>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>
#include <string>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
namespace po = boost::program_options;

//Hypothesis about the number of different lengths in the input
#define HYP_SIZE 35000

/**
 * Reads a file containing the distribution of sequences by length.
 * Each row in the file contains a length, a space and the number of proteins
 * for that length. If a length is not present, the number of
 * proteins is automatically set to 0.
 * @param filename the name of the file to read
 * @param S the vector for the distribution (modified by this function)
 * @param start the minimum length considered (associated to S[0])
 */
void read_distr(const char* filename, vector<unsigned int>* S, int start) {
	// current length
	unsigned int current;
	// last length encountered
	unsigned int lastseed = start - 1;
	// number of sequences associated to current length
	long found = 0;
	string line;

	ifstream myfile(filename, ios_base::in);
	if (myfile.fail()) {
		cerr << "ERROR: Can't open file " << filename << endl;
		exit(-1);
	}
	if (myfile.is_open()) {
		while (getline(myfile, line)) {
			istringstream ss(line);
			ss >> current;
			ss >> found;
			// fills the lengths without elements
			while (lastseed + 1 < current) {
				S->push_back(0);
				lastseed++;
			}
			if (current >= lastseed) {
				S->push_back(found);
				lastseed = current;
			}
		}
		myfile.close();
	}
}

/**
 * Reads a fasta file and populates a vector with the distribution of proteins
 * per length. Also returns a map of protein identifiers associated to lengths.
 * @param filename the file to read
 * @param S vector with length distribution (modified by this function)
 * @param start minimum length to consider (associated to S[0])
 * @return a map of length->vector of identifiers
 */
map<int, vector<string> > read_fasta(const char* filename, vector<unsigned int>* S, int start) {
	map<int, vector<string> > res = map<int, vector<string> >();
	string line;

	ifstream myfile(filename, ios_base::in);
	if (myfile.fail()) {
		cerr << "ERROR: Can't open file " << filename << endl;
		exit(-1);
	}
	if (myfile.is_open()) {
		string lastseen = "";
		int count = 0;
		while (getline(myfile, line)) {
			if (line[0] == '>') {
				if (lastseen.compare("") != 0) {
					res[count].push_back(lastseen);
				}
				lastseen = line.substr(1);
				count = 0;
			}
			else {
				count = count + line.size();
			}
		}
		res[count].push_back(lastseen); //the last one
	}
	myfile.close();

	unsigned int current;
	unsigned int lastseed = start - 1;
	long found = 0;
	for (std::map<int, vector<string> >::iterator it = res.begin(); it != res.end(); ++it) {
		current = it->first;
		found = it->second.size();
		while (lastseed + 1 < current) {
			S->push_back(0);
			lastseed++;
		}
		if (current >= lastseed) {
			S->push_back(found);
			lastseed = current;
		}
	}
	return res;
}

/**
 * Cost function
 * @param lower the lower bound (inclusive)
 * @param upper the upper bound (inclusive)
 * @param sumS vector of partial sums for the distribution
 * @param sumsigma vector of partial sums for length*#elements of that length
 * @return the cost of aligning sequences from lower to upper
 */
mpz_class cost(int lower, int upper, const vector<mpz_class> &sumS, const vector<mpz_class> &sumsigma) {
	mpz_class isigma = sumsigma[upper];
	if (lower > 0)
		isigma = isigma - sumsigma[lower - 1];
	mpz_class total = isigma * (sumS[upper] + 1);
	if (lower > 0)
		total = total - isigma * sumS[lower - 1];
	return total;
}

/**
 * The upper bound for a range
 * @param i the "center" of the range
 * @param maximum the upper bound of the domain
 * @param delta the coverage ratio [0-100]
 * @return the upper bound for the range
 */
int upper(int i, int maximum, int delta) {
	return min(int(floor(i / (delta / 100.0))), maximum);
}


/**
 * Evaluates the intervals identified by the algorithm
 * @param lastpos the last array index to be used
 * @param start the lower bound of lengths
 * @param p the array of interval bounds
 * @param length the length of array p
 * @param delta the coverage ratio
 * @return a vector of pairs representing the buckets
 */
vector<pair<int, int> > generate_intervals(int lastpos, int start, int p[], int length, int delta)
{
	int ll = max(p[lastpos] + start + 1, start);
	vector<pair<int, int> > intervals = vector<pair<int, int> >();
	int end = upper(lastpos + start, length - 1 + start, delta);
	intervals.push_back(pair<int, int>(ll, end));
	int i = p[lastpos];
	while ((ll > start) && (i > 0)) {
		ll = max(p[i] + start + 1, start);
		int end = upper(i + start, length - 1 + start, delta);
		intervals.push_back(pair<int, int>(ll, end));
		i = p[i];
	}

	return intervals;
}


/**
 * Prints the intervals identified by the algorithm
 * @param intervals the pairs representing interval bounds
 */
void print_intervals(vector<pair<int, int> > intervals)
{
	cout << "Intervals: " << endl;
	for (std::vector<pair<int, int> >::iterator it = intervals.begin(); it != intervals.end(); ++it)
	{
		int ll = it->first;
		int end = it->second;
		cout << ll << ", " << end << endl;
	}
}


int main(int argc, const char* argv[]) {
	vector<unsigned int> S;
	vector<mpz_class> sumsigma, sumS;
	S.reserve(HYP_SIZE);
	//length -> identifiers
	map<int, vector<string> > table;
	int length;
	int start = 1;
	string outdir;
	bool verbose = false;
	int delta;

	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("delta", po::value<int>(&delta)->default_value(90), "Coverage percentage needed, please see the paper")
		("start,s", po::value<int>(), "Minimum sequence length to consider")
		("distribution,d", po::value<string>(), "file with length distribution. The content should be <length> <num. of sequences>")
		("fasta,f", po::value<string>(), "file containing the sequences in FASTA format")
		("outdir", po::value<string>(), "output directory")
		("verbose,v", "verbose mode. WARNING: may generate a huge amount of data")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << desc << "\n";
		return 1;
	}

	if (vm.count("start")) {
		start = vm["start"].as<int>();
		cout << "Minimum length set to "
			<< start << ".\n";
	}
	else {
		cout << "Using default minimum length: " << start << endl;
	}

	if (vm.count("distribution")) {
		string fname = vm["distribution"].as<string>();
		cout << "Reading distribution file "
			<< fname << endl;
		read_distr(fname.c_str(), &S, start);
		cout << "Done reading." << endl;
	}
	else if (vm.count("fasta")) {
		string fname = vm["fasta"].as<string>();
		cout << "Reading fasta file "
			<< fname << endl;
		table = read_fasta(fname.c_str(), &S, start);
		cout << "Done reading." << endl;
	}
	else {
		cout << "ERROR: please specify a distribution or a fasta file" << endl;
		cout << desc << endl;
		return 2;
	}

	if (vm.count("outdir")) {
		outdir = vm["outdir"].as<string>();
		if (outdir[outdir.length() - 1] != '/')
			outdir.append("/");
		cout << "Output directory set to "
			<< outdir << endl;
	}
	else {
		outdir = "./";
		cout << "Output directory set to current working directory" << endl;
	}

	if (vm.count("verbose")) {
		verbose = true;
	}


	length = S.size();
	if (verbose)
		cout << "Different lengths considered = " << length << endl;
	sumsigma.resize(length);
	sumS.resize(length);

	sumsigma[0] = start * S[0];
	sumS[0] = S[0];
	if (verbose) {
		cout << "i\tS\tsumS\tsigma\tsumsigma" << endl;
		cout << 0 << "\t" << sumS[0] << "\t" << sumsigma[0] << endl;
	}
	for (int i = 1; i < length; i++) {
		sumS[i] = sumS[i - 1] + S[i];
		sumsigma[i] = sumsigma[i - 1] + (i + start) * S[i];
		if (verbose)
			cout << i << "\t" << sumS[i] << "\t" << sumsigma[i] << endl;
	}

	S.clear();

	vector<mpz_class> B;
	int p[length];
	B.resize(length);
	mpz_class naive = 0;

	ofstream outfile;
	if (verbose) {
		string ff = outdir + "matrix.csv";
		outfile.open(ff.c_str(), ios::out);
		outfile << "i\tB\tp" << endl;
	}


	for (int i = 0; i < length; i++) {
		int a = i + start;
		int upper;
		if (delta == 100)
			upper = min(a, length - 1 + start) - start;
		else
			upper = min(int(floor(a / (delta / 100.0))), length - 1 + start) - start;

		mpz_class M = sumsigma[upper] * (sumS[upper] + 1);
		p[i] = -1;
		B[i] = M;
		for (int k = 1; k <= i; k++) {
			int lower = k - 1;
			M = cost(lower + 1, upper, sumS, sumsigma);
			if (B[i] > M + B[lower]) {
				B[i] = M + B[lower];
				p[i] = lower;
			}
		}
		if (delta == 100)
			naive = naive + cost(max(a, start) - start, upper, sumS, sumsigma);
		else
			naive = naive + cost(max(int(ceil(a * (delta / 100.0))), start) - start, upper, sumS, sumsigma);
		if (verbose)
			outfile << i << "\t" << B[i] << "\t" << p[i] << endl;
	}

	if (verbose)
		outfile.close();

	mpz_class tmptotal = 0;

	cout << "Naive partitioning cost=" << naive << " using " << length << " buckets" << endl;
	cout << "Cost of full interval = " << cost(0, length - 1, sumS, sumsigma) << endl;

	sumsigma.clear();
	sumS.clear();

	cout << "Minimum cost found = " << B[length - 1] << endl;
	vector<pair<int, int> > intervals = generate_intervals(length - 1, start, p, length, delta);
	cout << "For a total of " << intervals.size() << " subsets" << endl;
	if (verbose)
		print_intervals(intervals);

	if (table.size() > 0)
	{ // write the buckets to file
		string fname = outdir + "buckets.list";
		ofstream bucket;
		bucket.open(fname.c_str(), ios::out);
		for (std::vector<pair<int, int> >::iterator it = intervals.begin(); it != intervals.end(); ++it)
		{
			int first = it->first;
			int second = it->second;

			for (int i = first; i <= second; i++)
			{
				for (std::vector<string>::iterator j = table[i].begin(); j != table[i].end(); ++j)
				{
					bucket << *j << "\t" << first << "-" << second << endl;
				}
			}
		}
		bucket.close();
	}

	B.clear();

	return 0;
}
