/*
 *  
 *  CLARK, CLAssifier based on Reduced K-mers.
 *
 *
 */

/*
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2013-2019, Rachid Ounit <clark.ucr.help at gmail.com>
 */

/*
 *  * @author: Rachid Ounit, Ph.D Candidate.
 *  * @project: CLARK, Metagenomic and Genomic Sequences Classification project.
 *  * @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *  *
 *  */


#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include "file.hh"
#include "parameters.hh"

using namespace std;

int main(int argc, char** argv)
{
	if (argc != 4)
	{
		cerr << "<settings file> <k-mer length> <min k-mers frequency (default: 0)>" << endl;
		exit(1);
	}	

	FILE * fd = fopen(argv[1], "r");
	if (fd == NULL)
	{
		cerr << "Failed to open " << argv[1] << endl;
		exit(1);
	}

	vector<string> ele;
	vector<char> sep;
	sep.push_back(' ');
	sep.push_back('\t');
	sep.push_back('\r');

	map<string,uint16_t> db;
	map<string,uint16_t>::iterator it;
	size_t i = 0;
	string line = "";

	getLineFromFile(fd, line);
	ele.clear();
	getElementsFromLine(line, sep, ele);
	string ftarget = ele[1];
	getLineFromFile(fd, line);
	ele.clear();
	getElementsFromLine(line, sep, ele);
	string predbfile = ele[1];

	fclose(fd);
	size_t m_km = atoi(argv[2]);
	fd = fopen(ftarget.c_str(), "r");
	size_t minCt = atoi(argv[3]);
	vector<size_t> countKmers;
	if (m_km < 32 && m_km > 2 && minCt >= 0 && fd != NULL)
	{
		while (getLineFromFile(fd, line))
		{
			ele.clear();
			getElementsFromLine(line, sep, ele);
			string tgt = ele[1];
			it = db.find(tgt);
			if (it == db.end())
			{
				db[tgt] = i;
				i++;
			}
		}
		fclose(fd);
		std::cerr << db.size() << " targets id identified." << endl;
		countKmers.resize(db.size(),0);
		char * flblname = (char*) calloc(10000, sizeof(char));
		sprintf(flblname,"%s/db_central_k%lu_t%lu_s%lu_m%lu.tsk.lb",predbfile.c_str(), m_km, (size_t) db.size(),(size_t) HTSIZE,(size_t) minCt);
		///////////////////////////////////////////////////////////////////////////

		FILE * flabel = fopen(flblname,"r");

		if (flabel == NULL)
		{       cerr << "Failed to open " <<  flabel << endl;
			exit(1);
		}
		size_t _k = m_km, kc = 0;

#define LEN 100000
		uint16_t c[LEN];
		size_t len = fread(c,2,LEN,flabel);
		size_t i =0, fileSize = len;
		while (true)
		{
			while (i < len)
			{
				countKmers[(int) c[i++]]++;
				kc++;
			}
			len = fread(c, 2, LEN, flabel);
			i = 0;
			if (len == 0)
			{       break;  }
		}
		fclose(flabel);
		std::cerr << kc << " k-mers successfully reads from database file:\n" << flblname << endl;
	}
	else
	{
		cerr << "Error with parameters! Please make sure, that: " << endl;
		cerr << " - The database has been set up and created for the k-mers length and min k-mers frequency passed." << endl;
		cerr << " - The parameters passed are correct: the k-mers length is between 2 and 32, and the min frequency is >0." << endl;
		exit(1);
	}
	///////////////////////////////////////////////////////////////////////////
	ofstream fout("./targets.distribution.csv");
	fout << "TaxonomyID,KmersCount," << std::endl;
	size_t t_cpt = 0, t = 0;
	for (it = db.begin(); it != db.end(); it++)
	{
		fout << it->first << "," << countKmers[it->second] << "," << std::endl;
	}
	fout.close();
	cerr << "The k-mers distribution per target is available in the file: targets.distribution.csv" << endl;
	return 0;
	//////////////////////////
}
