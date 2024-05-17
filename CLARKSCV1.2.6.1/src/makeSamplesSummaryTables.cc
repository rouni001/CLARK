/*
 *  CLARK, CLAssifier based on Reduced K-mers.
 */

/*
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *    Copyright 2013-2019, Rachid Ounit <clark.ucr.help at gmail.com>
 */

/*
 *  @author: Rachid Ounit, Ph.D.
 *  @project: CLARK, Metagenomic and Genomic Sequences Classification project.
 *  @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *  
 */

#include <cstdlib>
#include <iomanip>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>
#include "file.hh"
#include <string>
#include <stdint.h>

using namespace std;

struct Organism
{
        string 	Name;
        string	Domain;
	size_t	Count;
	size_t	Taxid;
        double 	RProportion;
	double	TProportion;

        Organism(): Name(""), Domain("NA"), Count(-1), Taxid(-1), RProportion(-1.0), TProportion(-1.0)
        {}

        bool operator<(const Organism& _o) const
        {       return (RProportion < _o.RProportion);  }
};

struct FileS
{
	size_t	Idx;
	size_t	Count;
	
	FileS(): Idx(-1), Count(-1)
	{}

	bool operator<(const FileS& _f) const
	{
		return Count < _f.Count;
	}
};


bool CompareVectors(const vector<FileS>& _a, const vector<FileS>& _b, const size_t& _i, const size_t& _size)
{
	if (_i >= _size)
	{	return true;}
	if (_a[_i].Count == _b[_i].Count)
	{	return CompareVectors(_a,_b,_i+1,_size);}
	else
	{	return _a[_i].Count < _b[_i].Count;}
}


struct Target
{
	string 	Name;
	string 	Domain;
	size_t	Taxid;
	size_t	Freq;
	vector<FileS> Counts;
	Target(const size_t& _cap):  Name(""), Domain("NA"), Taxid(-1), Freq(-1)
	{
		Counts.resize(_cap);
	}	
	bool operator<(const Target& _o) const
	{
		if (Freq == _o.Freq)
		{
			return CompareVectors(Counts, _o.Counts, 0, Counts.size());
			//return Name < _o.Name;
		}
		return Freq < _o.Freq;
	}
	void Display() const
	{
		for(size_t t = 0; t < Counts.size(); t++)
		{	cerr << Name << "-" << Counts[t].Idx << "-" <<  Counts[t].Count << "-" << Freq << "\t" ;}
		cerr << endl;
	}
	void Finalize() 
	{
		Freq = 0;
		for(size_t t = 0; t < Counts.size(); t++)
		{
			Freq += Counts[t].Count!=-1?1:0;	
		}
	}
};

struct Sample
{
	string 	Name;
	int	NReads;
	int 	NRAssigned;
	double 	ADiversity;
	vector<Organism> TopOrgs;
	
	Sample(): Name(""), NReads(0), NRAssigned(0), ADiversity(0.0)
	{}
	
	void Add(const Organism& o)
	{
		TopOrgs.push_back(o);
	}	
	
	void Display(const double& abd, ofstream& fout, const bool& isExtended) const
	{
		fout << setprecision(4) <<  Name << "," << NReads << "," << NRAssigned << "," << ADiversity << ",";
		fout << abd << "," << ((double) NRAssigned)/((double) NReads)*100.0 << "%," ;
		for(size_t u = 0; u < TopOrgs.size(); u++)
		{
			if (TopOrgs[u].TProportion >= abd)	
			{	fout << TopOrgs[u].Name;
				if (isExtended)
				{	fout << "," << TopOrgs[u].Taxid;	} 
				fout << "," << TopOrgs[u].RProportion << "%," ;
			}
		}
		fout << std::endl;
	}

	void Update(const double& abd, std::map<size_t, size_t>& taxidToIdx, std::vector< Organism >& orgs, std::vector< std::vector<double> >& records) const
	{
		std::map<size_t, size_t>::iterator it;
		size_t idx = 0;
		for(size_t u = 0; u < TopOrgs.size(); u++)
                {
                        if (TopOrgs[u].TProportion >= abd)
                        {       
				it = taxidToIdx.find(TopOrgs[u].Taxid);
				if (it == taxidToIdx.end())
				{	idx = taxidToIdx.size();
					taxidToIdx[TopOrgs[u].Taxid] = idx;
					std::vector< double > tmp;
					records.push_back(tmp);
					orgs.push_back(TopOrgs[u]);
				}
				else
				{
					idx = it->second;
				}
				records[idx].push_back(TopOrgs[u].TProportion);	
			}
                }
	}
	
	size_t getNumberOrgs() const
	{
		return (TopOrgs.size());
	}
};

size_t hString(const string& _val)
{
	size_t res = 0;
	for(size_t i = 0; i < _val.size() && i < 32 ; i++)
	{	res = res*256+((uint8_t) _val[i]);}
	return res;
}

bool isCLARKReportFile(const char* _filename)
{
	FILE* fd = fopen(_filename,"r");
	if (fd == NULL)
	{	return false;	}
	
	string line = "";
	vector<string> ele;
	vector<char> sep;
	sep.push_back(',');
	getLineFromFile(fd, line);
	getElementsFromLine(line, sep, ele);
	fclose(fd);
	return (ele.size() > 3 && ele[ele.size()-2] == "Proportion_All(%)" && ele[ele.size()-1] == "Proportion_Classified(%)");
}


std::string getFileID(const char* _filename)
{
	string file(_filename);
	vector<char> sep, sep2;
	sep.push_back('/');
	sep2.push_back('.');
	vector<string> ele;
	getElementsFromLine(file, sep, ele);
	if (ele.size() <  1)
	{	std::cerr << "Failed to read the filename: " << _filename << std::endl; exit(1);	}
	string name= ele[ele.size()-1];
	ele.clear();
	getElementsFromLine(name, sep2, ele);
	return ele[0];
}

bool getSummary(const char* filename, const size_t& topk, Sample& S)
{
	string 		line;
	vector<char> 	sep, sepd;
	vector<string> 	ele, eled;
	sep.push_back(',');
	sep.push_back('\r');
	sepd.push_back(';');
	FILE * fd = fopen(filename, "r");
	if (fd == NULL)
	{	std::cerr << "Failed to open " << filename << endl;
		exit(1);
	}
	int nbReads = 0, nbAssigned = 0;
	double diversity = 0;
	getLineFromFile(fd, line);
	
	int incr = (line.find("Lineage") != string::npos)?2:1;
	vector<Organism> Orgs;
	while (getLineFromFile(fd, line))
	{	// Name,TaxID,Count,Proportion_All(%),Proportion_Classified(%)
		// Name,TaxID,Lineage,Count,Proportion_All(%),Proportion_Classified(%)
		ele.clear();
		getElementsFromLine(line, sep, ele);
		if (ele[0] != "UNKNOWN")
		{
			Organism 	e;
			e.Name	 	= ele[0];
			e.Taxid 	= (incr > 0)?(atoi(ele[1].c_str())):(hString(ele[0]));
			e.Count		= atoi(ele[incr+1].c_str());
			e.RProportion 	= atof(ele[incr+3].c_str());
			e.TProportion	= atof(ele[incr+2].c_str());
			if (incr > 1)
			{
				eled.clear();
				getElementsFromLine(ele[2], sepd, eled);
				e.Domain = eled[0];
			}
			nbAssigned 	+= atoi(ele[incr+1].c_str());
			nbReads 	+= atoi(ele[incr+1].c_str());
			diversity	+= (e.RProportion/100.0) * log((e.RProportion/100.0));
			Orgs.push_back(e);
		}
		else
		{
			nbReads 	+= atoi(ele[incr+1].c_str());
		}
	}
	fclose(fd);
	
	sort(Orgs.begin(), Orgs.end());
	S.Name 		= getFileID(filename);
	S.NReads	= nbReads;
	S.NRAssigned 	= nbAssigned;  
	S.ADiversity	= -1 * diversity; 
	for (size_t t = 0; t < topk && t < Orgs.size() ; t++)
	{	
		S.Add(Orgs[Orgs.size()-1-t]);
	}
	return (incr == 2);
}


void DisplaySamples2(const int& topK, const double& abd, const vector<Sample>& Samples, const bool& _isExtended)
{
	vector<Target> 			HitCounts;
	map<string,size_t> 		db;
	map<string,size_t>::iterator 	it;
	
	for(size_t t= 0; t < Samples.size(); t++)
	{
		for(size_t o = 0; o < Samples[t].TopOrgs.size(); o++)
		{
			if (Samples[t].TopOrgs[o].TProportion < abd)
			{	continue;	}
			it = db.find(Samples[t].TopOrgs[o].Name); //TaxID
			if (it == db.end())
			{
				Target e(Samples.size());
				e.Taxid  = Samples[t].TopOrgs[o].Taxid;
				e.Domain = Samples[t].TopOrgs[o].Domain; 
				e.Name   = Samples[t].TopOrgs[o].Name;
				db[e.Name] = HitCounts.size();
				HitCounts.push_back(e);
			}
		}
	}
	for(size_t t= 0; t < Samples.size(); t++)
        {
                for(size_t o = 0; o < Samples[t].TopOrgs.size(); o++)
		{
			if (Samples[t].TopOrgs[o].TProportion < abd)
                        {       continue;       }
                	it = db.find(Samples[t].TopOrgs[o].Name);
			if (it == db.end())
				continue;
                	HitCounts[it->second].Counts[t].Idx = t;
			HitCounts[it->second].Counts[t].Count = Samples[t].TopOrgs[o].Count;
		}
	}
	/* Finalize */
	for(size_t u = 0; u < HitCounts.size(); u++)
        {	HitCounts[u].Finalize();	}	
	for(size_t t = 0; t < Samples.size(); t++)
	{
		size_t c_idx = 0;
		for(size_t u = 0; u < HitCounts.size(); u++)
		{	if (HitCounts[u].Counts[t].Idx != -1)
			{	c_idx = HitCounts[u].Counts[t].Idx;
				break;
			}
		}
		for(size_t u = 0; u < HitCounts.size(); u++)
                {       HitCounts[u].Counts[t].Idx = c_idx;
                }
	}
	/* Sorting */
	sort(HitCounts.begin(), HitCounts.end());
	
	
	/* Display */
	ofstream fout("TableSummary_HitCount.csv");
	string prefix =_isExtended?",,,":"";
	fout << prefix << "Filenames,"  ;
	for(size_t t=0; t < Samples.size(); t++)
	{	
		fout << Samples[HitCounts[0].Counts[t].Idx].Name << "," ;
	}
	fout << endl; 
	fout << prefix <<"#TotalReads," ;
	for(size_t t=0; t < Samples.size(); t++)
        {
                fout << Samples[HitCounts[0].Counts[t].Idx].NReads << "," ;
        }
        fout << endl;
        fout << prefix << "#TotalReadsMapped," ;
	for(size_t t=0; t < Samples.size(); t++)
        {
                fout << Samples[HitCounts[0].Counts[t].Idx].NRAssigned << "," ;
        }
        fout << endl;
       	if (_isExtended)
	{	fout << "Domain,TaxonomyID,Name," ;	}
	fout << "AlphaDiversity,";

	for(size_t t=0; t < Samples.size(); t++)
        {
                fout << Samples[HitCounts[0].Counts[t].Idx].ADiversity << "," ;
        }
        fout << endl;
	for(size_t t=0; t < HitCounts.size(); t++)
        {
		if (_isExtended)
		{	fout << HitCounts[t].Domain<< "," << HitCounts[t].Taxid << "," << HitCounts[t].Name << ",#ReadsMapped," ;}
		else
		{	fout << HitCounts[t].Name << "," ;}
		for(size_t u = 0 ; u < Samples.size(); u++)
		{
			if (HitCounts[t].Counts[u].Count == -1)
			{	fout << ",";}
			else
			{	fout << HitCounts[t].Counts[u].Count << ",";}
		}
		fout << endl;
        }
	fout.close();
}

void DisplaySamples(const int& topK, const double& abd, const vector<Sample>& Samples, const bool& _isExtended)
{
	ofstream fout("./TableSummary_per_Report.csv");
	fout << "SampleName,#Reads,#ClassifiedReads,Alpha-Diversity(S-Index),MinAbundance(%),Proportion(%),";
	
	int maxNOrg = topK, tmpN = 0;
	for (size_t t = 0; t < Samples.size(); t++)
	{
		if (Samples[t].getNumberOrgs() > tmpN)
		{	tmpN = Samples[t].getNumberOrgs();}
	}
	maxNOrg = (topK < tmpN)?topK:tmpN;

	for(size_t i = 0; i < maxNOrg ; i++ )
        {
                fout << "Top" << i+1 << "_Taxon_Name,";
		if (_isExtended)
		{
			fout << "Top" << i+1 << "_Taxon_Taxid,";
		}
		fout << "Top" << i+1 << "_Taxon_RelProp(%),";
	}
	fout << std::endl;
	for(size_t i = 0; i < Samples.size() ; i++ )
	{
		Samples[i].Display(abd, fout, _isExtended);	
	}
	fout.close();
}

void DisplayOrgDistr(const double& abd, const vector<Sample>& Samples, const bool& _isExtended)
{
        std::map<size_t, size_t> TaxidToIdx;
	std::map<size_t, size_t>::iterator it;
	std::vector< Organism > Orgs;
	std::vector< std::vector< double > > Records;

	for(size_t i = 0; i < Samples.size() ; i++ )
        {
                Samples[i].Update(abd, TaxidToIdx, Orgs, Records);
        }
	ofstream fout("./TableSummary_per_Taxon.csv");
	if (_isExtended)
	{
		fout << "Taxon_Name,Taxon_Taxid,Taxon_Domain,MinAbundance(%),#ReportsFound,MinProportion(%),MaxProportion(%),AvgProportion(%),StdDevProportion(%),\n";	
	}
	else
	{
		fout << "Taxon_Name,MinAbundance(%),#ReportsFound,MinProportion(%),MaxProportion(%),AvgProportion(%),StdDevProportion(%),\n";
	}	

	for(it = TaxidToIdx.begin(); it!= TaxidToIdx.end(); it++)
	{
		double avg = 0, stdDev = 0;
		double min = 1000.0, max = -1.0;
		for(size_t v = 0 ; v < Records[it->second].size(); v++)
		{
			if (Records[it->second][v] > max)
			{	max = Records[it->second][v];	}
			if (Records[it->second][v] < min)
                        {       min = Records[it->second][v];   }
			avg += Records[it->second][v];
		}
		avg /= ((double) Records[it->second].size());
		for(size_t v = 0 ; v < Records[it->second].size(); v++)
                {
                        stdDev += (Records[it->second][v] - avg) * (Records[it->second][v] - avg);
                }
		stdDev = sqrt(stdDev);
		stdDev /= ((double) Records[it->second].size());
		fout << Orgs[it->second].Name << "," ;
		if (_isExtended)
		{	fout << Orgs[it->second].Taxid << "," <<  Orgs[it->second].Domain << ","; }
		fout << abd << "," << Records[it->second].size() << "," << min <<"," << max << "," << avg << "," << stdDev << ",\n";
	}
	fout.close();
}

int main(const int argc, const char** argv)
{
	if (argc < 4)
	{	std::cerr << argv[0] << " <top r>  <min abundance in %> <CLARK Report filename1> <CLARK Report filename2> ... "<< std::endl; exit(1);}
	
	vector<Sample> Samples;
	size_t lim = atoi(argv[1]);
	double minAbundance = atof(argv[2]);
	bool isExtended = false;
	for(size_t t = 3; t < argc; t++)
	{
		if (!isCLARKReportFile(argv[t]))
		{
			cerr << "Failed to read " << argv[t] << ". Please check the pathname and that it is a correct CLARK report file." << endl;
			continue;
		}
		Sample S;
		
		isExtended = getSummary(argv[t], lim, S);
		if (S.Name == "" || (S.Name).size() < 1)
		{	continue;}
		Samples.push_back(S);
		
	}
	// orgs
	DisplaySamples(lim, minAbundance, Samples, isExtended);
	DisplaySamples2(lim, minAbundance, Samples, isExtended);
	if (isExtended)
	{
		cerr << "Files: TableSummary_per_Report.csv, TableSummary_HitCount.csv ";
		DisplayOrgDistr(minAbundance, Samples, isExtended);
                cerr << "and TableSummary_per_Taxon.csv were successfully generated." << endl;
                return 0;
	}
	else
	{
		cerr << "Files: TableSummary_per_Report.csv and TableSummary_HitCount.csv were generated." << endl;
		return 0;
	}
}

