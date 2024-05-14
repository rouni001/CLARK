/*
 * CLARK, CLAssifier based on Reduced K-mers.
 */

/*
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Copyright @The Regents of the University of California. All rights reserved.

*/
/*
 * @author: Rachid Ounit, Ph.D.
 * @project: CLARK, Metagenomic and Genomic Sequences Classification project.
 * @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *
 */

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <map>
#include <fstream>
#include <algorithm>
using namespace std;
#include "./file.hh"
#define MXNMLEN 1000

struct node
{
	uint32_t        parent;
	uint8_t         rank;
	node():parent(0),rank(255) {}
};

struct targetAbundance
{
	string 	name;
	string	taxid;
	size_t 	abundance;
	vector<node> lineage;
	targetAbundance():name(""),abundance(0)
	{}	
	bool operator<(const targetAbundance& _e) const
	{
		return name < _e.name;
	}
};

#define NBNODE 8

bool getLineage(const vector<node>& _nodes, const uint32_t& _taxid, vector<node>& _line)
{
	_line.clear();
	_line.resize(NBNODE);
	size_t it = _taxid, tmp = 0;
	if (_nodes[it].parent == 0)
	{
		return false;
	}
	while (true)
	{
		if (_nodes[it].parent == 1)
		{       
			_line[NBNODE-1].parent = 1;
			if (_line[NBNODE-2].parent == 0)
			{ 	
				_line[NBNODE-2].parent = it; 
				_line[NBNODE-2].rank = 0;
			}

			break;
		}
		if (_nodes[it].rank < NBNODE && _line[ _nodes[it].rank ].rank != 0)
		{
			_line[ _nodes[it].rank ].rank = 0;
			_line[ _nodes[it].rank ].parent = it;
		}
		tmp = it;
		it = _nodes[tmp].parent;
	}
	return true;
}

void getNodes(const char* file, vector<node>& nodes)
{
	FILE * fdn = fopen(file, "r");
	if (fdn == NULL)
	{       cerr << "Failed to open " << file << endl;
		exit(-1);
	}
	std::map<std::string,uint8_t> nameTorank;
	std::map<std::string,uint8_t>::iterator it;
	nameTorank["species"] = 0;
	nameTorank["genus"] = 1;
	nameTorank["family"] = 2;
	nameTorank["order"] = 3;
	nameTorank["class"] = 4;
	nameTorank["phylum"] = 5;
	nameTorank["superkingdom"] = 6;
	nameTorank["root"] = 7;
#define MAXNB 200000000
	nodes.clear();
	nodes.resize(MAXNB);

	string line;
	vector<string> ele;
	vector<char> sep;
	sep.push_back(' ');
	sep.push_back('|');
	sep.push_back('\t');
	cerr << "Loading nodes of taxonomy tree... " ;
	int id, idp;
	while (getLineFromFile(fdn, line))
	{
		ele.clear();
		getElementsFromLine(line, sep, ele);
		id = atoi(ele[0].c_str());
		idp = atoi(ele[1].c_str());
		nodes[id].parent = idp;
		it = nameTorank.find(ele[2].c_str());
		if (it != nameTorank.end()  && (ele.size()==3 || ele[3].find("group") == std::string::npos))
		{       nodes[id].rank = it->second; 	}
	}
	cerr<< "done\n";
	fclose(fdn);
}

std::string getmpaFormatted(const std::string& _name)
{
	std::string res = "";
	for(size_t i = 0; i < _name.size(); i++)
	{
		res.push_back(_name[i]==' '?'_':_name[i]);
	}
	return res;
}

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		cerr <<" -c <minConfidenceScore> -g <minGamma> -D <Directory_Path> -F <result1>.csv <result2>.csv ... <result_n>.csv -a <minAbundance> ... "<<endl; 
		cerr <<"\nDefinition of parameters: \n" << endl;
		cerr <<"-c <minConfidenceScore>           \t To filter assignments based on their confidence score (if available) using the\n";
		cerr <<"                                  \t threshold value minConfidenceScore (a value between 0.5 and 1.0). \n";
		cerr <<"                                  \t The abundance estimation for each target will count only\n";
		cerr <<"                                  \t assignments with a confidence score higher than minConfidenceScore.\n";
		cerr <<"                                  \t Assignments that have a confidence score lower than minConfidenceScore\n";
		cerr <<"                                  \t will be marked as unclassified and so counted in the \n";
		cerr <<"                                  \t category UNKNOWN in the abundance estimation report.\n";
		cerr <<"                                  \t The default value is 0.5.\n"<< endl;
		cerr <<"-g <minGamma>                     \t To filter assignments based on their gamma score (if available) using the\n";
		cerr <<"                                  \t threshold value minGamma (a value between 0 and 1.0). \n";
		cerr <<"                                  \t The abundance estimation for each target will count only\n";
		cerr <<"                                  \t assignments with a gamma score higher than minGamma.\n";
		cerr <<"                                  \t Assignments that have a gamma score lower than minGamma\n";
		cerr <<"                                  \t will be marked as unclassified and so counted in the \n";
		cerr <<"                                  \t category UNKNOWN in the abundance estimation report.\n";
		cerr <<"                                  \t The default value is 0.\n"<< endl;
		cerr <<"-D <Directory_Path>               \t The directory path of the database (the same you indicated when calling set_targets.sh)."<< endl;
		cerr <<"                                  \t This parameter is required to load scientific names of all targets ONLY if you pass"<< endl;
		cerr <<"                                  \t results of a metagenomic sample.\n" << endl;
		cerr <<"-F <result1>.csv <result2>.csv ... <result_n>.csv" << endl;
		cerr <<"                                  \t results file or list of results file produced by CLARK.\n";
		cerr <<"                                  \t Important Note: You can pass a results file produced by any mode of execution of CLARK" << endl;
		cerr <<"                                  \t (full, express, spectrum, default), but if you pass several files, make sure they all have\n";
		cerr <<"                                  \t been produced by the same mode. For example, if you pass result1.csv and result2.csv\n";
		cerr <<"                                  \t then result1.csv and result2.csv should be produced through the same mode (e.g., full mode).\n\n";
		cerr <<"-a <minAbundance(%)>              \t To output only estimations that are higher than a certain threshold, minAbundance."<< endl;
		cerr <<"                                  \t minAbundance is a percentage value (between 0 and 100).\n" << endl;
		cerr <<"--highconfidence                  \t To count only high confidence assignments for the abundance estimation. "<< endl;
		cerr <<"                                  \t This option is equivalent to '-c 0.75 -g 0.03', in this version.\n"<< endl;
		cerr <<"--krona				  \t To export results in a 3-column file (i.e., taxon_id,taxon_id,magnitude) for "<<endl;
		cerr <<"                                  \t the Krona tool (Hierarchical data browser, Ondov et al, BMC Bioinformatics, 2011," << endl;
		cerr <<"                                  \t doi: 10.1186/1471-2105-12-385. https://github.com/marbl/Krona/wiki)."<< endl;
		cerr <<"                                  \t With this option on, the program will create the file 'results.krn' containing " << endl;
		cerr <<"                                  \t the unnormalized abundance from CLARK's assignments (with the filtering options if any)." << endl;
		cerr <<"                                  \t Then, to create the Krona diagram, run ktImportTaxonomy with the option '-m 3', for example:" << endl;
		cerr <<"                                  \t $ ktImportTaxonomy-o results.html -m 3 results.krn" << endl;
		cerr <<"--mpa                             \t To export results in the MetaPhlan's mpa format. A two-column file, where the first column," << endl;
		cerr <<"                                  \t is the taxonomy rank and the second column is the total number of reads mapped to that rank or below." << endl;
		cerr << endl;
		exit(1);
	}
	int i_deb = -1, i_end = -1, i_names = -1;
	double minConf = 0.5, minGamma = 0, min = 0;
	bool krona= false;
	bool mpa=false;

	for(int t = 1; t < argc; t++)
	{
		string param(argv[t]);
		if (param == "--highconfidence" || param == "--hc")
		{
			minConf = 0.75;
			minGamma = 0.03;
			continue;
		}
		if (param == "--krona" || param == "--Krona" || param == "--KRONA")
		{
			krona = true;
			continue;
		}
		if (param == "--mpa" || param == "--MPA" || param == "--Mpa")
		{
			mpa = true;
			continue;
		}
		if (param == "-c")
		{
			if (++t >= argc)
			{
				cerr << "Please provide a minimum for the confidence score." << endl;
				exit(1);
			}
			minConf = atof(argv[t]);
			if (minConf <  0.5 || minConf > 1)
			{       cerr << "Please provide a minimum confidence score between 0.5 and 1." << endl;
				exit(1);
			}
			continue;
		}
		if (param == "-g")
		{
			if (++t >= argc)
			{
				cerr << "Please provide a minimum value for the gamma." << endl;
				exit(1);
			}
			minGamma = atof(argv[t]);
			if (minGamma <  0 || minGamma > 1)
			{       cerr << "Please provide a minimum Gamma score between 0 and 1." << endl;
				exit(1);
			}
			continue;
		}
		if (param == "-a")
		{
			if (++t >= argc)
			{
				cerr << "Please provide a minimum value for the abundance." << endl;
				exit(1);
			}
			min = atof(argv[t]);
			if (min <  0 || min > 100)
			{	cerr << "Please provide a minimum abundance between 0 and 100." << endl;
				exit(1);
			}
			continue;
		}
		if (param == "-D")
		{
			if (++t >= argc)
			{
				cerr << "Please provide the directory path containing the database." << endl;
				exit(1);
			}
			i_names = t;
			continue;
		}
		if (param == "-F")
		{
			if (++t >= argc)
			{
				cerr << "Please provide one (or several) CLARK output file(s) (CSV format)." << endl;
				exit(1);
			}
			i_deb = t;
			while (true)
			{
				if (++t == argc || argv[t][0] == '-')
					break;
			}
			i_end = t;
			t--;
			continue;
		}
		cerr << "Failed to recognize option: " << argv[t] << endl;
		exit(1);
	}

	FILE * fd = fopen(argv[i_deb], "r");
	if (fd == NULL)
	{
		cerr << "Failed to open " << argv[i_deb] << endl;
		exit(1);
	}
	string line;
	vector<string> ele;
	vector<char> sep;
	sep.push_back(',');
	sep.push_back('\t');
	sep.push_back('\r');
	std::map<std::string, uint32_t>			idTodDiD;
	std::map<std::string, uint32_t>::iterator 	it;	

	getLineFromFile(fd, line);
	getElementsFromLine(line, sep, ele);
	if (ele.size() < 3)
	{
		cerr << "Failed to extract all data from the file: "<<argv[i_deb]<<". The file does not seem to be a CLARK results file."<< endl;
		exit(1);
	}
	size_t idx = ele.size() == 3 ? 2: ele.size()-5;

	uint32_t i_lbl = 0;
	vector<size_t> abundance;
	vector<std::string> dLabels, dLabelsN;
	vector< vector< node > > lineages;
	map<uint32_t,string> 				idToName;
	map<uint32_t,string>::iterator 			itl;

	size_t total = 0;
	bool admissible = true;
	while (i_deb < i_end)
	{
		cerr << "\rFile: " << argv[i_deb] << "    ";
		while (getLineFromFile(fd, line))
		{
			ele.clear();
			getElementsFromLine(line, sep, ele);
			total++;
			it = idTodDiD.find(ele[idx]);
			admissible = true;
			// check whether the assignment is admissible
			if (idx > 2)
			{
				admissible = atof(ele[idx-1].c_str()) >= minGamma && atof(ele[idx+4].c_str()) >= minConf;
			}
			if (!admissible)
			{
				ele[idx] = "NA";
				it = idTodDiD.find("NA");
			}
			if (it == idTodDiD.end())
			{
				idTodDiD[ele[idx]] = i_lbl++;
				abundance.push_back(1);
				dLabels.push_back(ele[idx]);
				dLabelsN.push_back(ele[idx]);
			}
			else
			{
				abundance[it->second]++;
			}
		}
		fclose(fd);
		i_deb++;
		if (i_deb != i_end)
		{
			fd = fopen(argv[i_deb], "r");	
			getLineFromFile(fd, line);
		}
	}
	cerr <<"\n";
	if (i_names > 0)
	{
		sep.clear();
		sep.push_back('|');
		vector<char> sepc;
		sepc.push_back('\t');
		vector<string> elec;
		char * filename= (char*) calloc(MXNMLEN,sizeof(char));
		sprintf(filename,"%s/taxonomy/nodes.dmp", argv[i_names]);
		vector<node> nodes;
		getNodes(filename,nodes);
		vector<node> lineage;
		lineages.resize(dLabels.size());

		cerr << "Start retrieving lineage for each target identified (" <<  dLabels.size() << ")... " ;
		size_t i = 0;
		while (i < dLabels.size())
		{
			if (dLabels[i] == "NA")
			{       i++;
				continue;
			}
			lineage.clear();
			if (!getLineage(nodes, atoi(dLabels[i].c_str()), lineage))
			{
				cerr << "\nFailed to identify " << dLabels[i] << ": Unknown taxonomy id given the provided taxonomy database."<< endl;		
				dLabels[i]  = "NA";
				dLabelsN[i] = "NA";
				i++;
				continue;
			}                        
			for(size_t t = 0; t < NBNODE-1; t++)
			{
				if (lineage[t].rank == 0 && idToName.find(lineage[t].parent) == idToName.end())
				{	idToName[lineage[t].parent] = "";	}
				lineages[i].push_back(lineage[t]);
			}
			i++;
		}
		cerr << "done." << endl;

		sprintf(filename,"%s/taxonomy/names.dmp", argv[i_names]);
		FILE * fdn = fopen(filename, "r");

		if (fdn == NULL)
		{	cerr << "Failed to open " << filename << endl;
			cerr << "The program will estimates abundance per taxonomy id." << endl;
		}
		else
		{
			size_t cpt = 0;
			cerr << "Retrieving scientific names from taxonomy tree... ";
			while (getLineFromFile(fdn, line) && cpt < dLabels.size() )
			{	
				ele.clear();
				getElementsFromLine(line, sep,ele);
				elec.clear();
				getElementsFromLine(ele[0], sepc, elec);
				it = idTodDiD.find(elec[0]);
				if (it != idTodDiD.end() && ele.size()> 2 && ele[3].find("scientific name") != std::string::npos)
				{
					cpt++;
					elec.clear();
					getElementsFromLine(ele[1], sepc, elec);
					dLabelsN[it->second] = elec[0];
				}
				itl = idToName.find(atoi(elec[0].c_str()));
				if (itl != idToName.end() && ele.size()>2 && ele[3].find("scientific name") != std::string::npos)
				{
					elec.clear();
					getElementsFromLine(ele[1], sepc, elec);
					itl->second = elec[0];
				}	
			}
			cerr << "done." << endl;
			fclose(fdn);
		}
		free(filename);
		filename=NULL;
	}
	vector<targetAbundance> res(dLabels.size());
	for (size_t t= 0; t < dLabels.size(); t++)
	{
		res[t].taxid 		= dLabels[t];
		res[t].name 		= dLabelsN[t];
		res[t].abundance 	= abundance[t];
		if (res[t].name == "NA")
			continue;

		if (lineages.size() > 0)
		{
			for(size_t v = 0; v < lineages[t].size(); v++)
			{	
				res[t].lineage.push_back(lineages[t][v]);	
			}
		}
	}
	std::sort(res.begin(), res.end());	
	if (i_names < 0)
	{ 	cout << "Name,TargetID,";	}
	else
	{	cout << "Name,TaxID,Lineage";	

	}
	cout << ",Count,Proportion_All(%),Proportion_Classified(%)" << endl;
	size_t unk_obj = 0;
	double a = 0, a2 = 0;
	for (size_t t = 0; t < res.size(); t++)
	{
		if (res[t].name == "NA")
		{	unk_obj += res[t].abundance;	
		}
	}
	for (size_t t = 0; t < res.size(); t++)
	{
		if (res[t].name == "NA")
			continue;

		a =  100*( (double) res[t].abundance)/((double) (total)) ;
		a2 = 100*( (double) res[t].abundance)/((double) (total-unk_obj)) ;

		if (a >= min)
		{	cout << res[t].name << "," << res[t].taxid << ",";

			if (res[t].lineage.size() > 0)
			{
				size_t len = res[t].lineage.size();
				cout << idToName[res[t].lineage[len-1].parent] ;
				for(size_t u = len-2 ; u > 0; u--)
					cout << ";" << idToName[res[t].lineage[u].parent];
				cout << ",";
			}			
			cout << res[t].abundance << "," << a << "," << a2 << endl;	
		}
	}
	a =  100*( (double)  unk_obj)/((double) total) ;
	if (a >= min)
	{	if (i_names > 0)
		{	cout << "UNKNOWN,UNKNOWN,UNKNOWN," << unk_obj << "," << a << ",-" << endl; }
		else
		{ 	cout << "UNKNOWN,UNKNOWN," << unk_obj << "," << a << ",-" << endl;}
	}
	if (krona)
	{
		ofstream fout("results.krn", std::ios::binary);
		for(size_t t = 0; t < res.size(); t++)
		{
			if (res[t].name != "NA")
			{
				fout << res[t].taxid << " \t " << res[t].taxid << " \t " << res[t].abundance << endl;
			}
		}
		fout.close();
	}
	if (mpa)
	{
		vector<string> ranks;
		ranks.push_back("s__");		ranks.push_back("g__"); ranks.push_back("f__");         ranks.push_back("o__"); ranks.push_back("c__");         ranks.push_back("p__"); ranks.push_back("d__");
		ofstream fout("results.mpa", std::ios::binary);
		map<uint32_t,int> rankTaken;
		map<uint32_t,int>::iterator it;
		for(size_t t=NBNODE-1; t > 0; t--)
		{
			for(size_t r = 0; r < res.size(); r++)
			{
				if (res[r].lineage.size() <= t || res[r].lineage[t].rank != 0)
					continue;
				uint32_t c_rank = res[r].lineage[t].parent;
				it = rankTaken.find(c_rank);
				if (it != rankTaken.end())
					continue;
				rankTaken[c_rank] = 1;
				int c_count = res[r].abundance;
				size_t len = res[r].lineage.size();
				fout << ranks[len-1] << getmpaFormatted(idToName[res[r].lineage[len-1].parent]);
				for(size_t v = len-2; v >= t; v--)
				{	if (idToName[res[r].lineage[v].parent]  != "")
					{	fout << "|" << ranks[v]  << getmpaFormatted(idToName[res[r].lineage[v].parent]) ;}
				}

				for (size_t s = 0; s < res.size(); s++)
				{
					if (r == s || res[s].lineage.size() <= t)
						continue;

					if (res[s].lineage[t].parent == c_rank)
					{
						c_count += res[s].abundance;
					}
				}
				fout << "\t" << c_count << endl;
			}	
		}
		for(size_t r = 0; r < res.size(); r++)
		{
			if (res[r].name == "NA")
				continue;
			int c_count = res[r].abundance;
			size_t len = res[r].lineage.size();
			fout << ranks[len-1] << getmpaFormatted(idToName[res[r].lineage[len-1].parent]);
			for(size_t v = len-2; v > 0; v--)
			{       if (idToName[res[r].lineage[v].parent]  != "")
				{       fout << "|" << ranks[v]  << getmpaFormatted(idToName[res[r].lineage[v].parent]) ;}
			}
			fout << "|" << ranks[0] << getmpaFormatted(res[r].name);
			fout << "\t" << c_count << endl;
		}
		fout.close();
	}
	return 0;
}
