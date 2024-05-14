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

#ifndef HASHTABLESTORAGE_HH
#define HASHTABLESTORAGE_HH

#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <string>
#include <stdint.h>
#include <map>

#include "./hashTable_hh.hh"
#include "./kmersConversion.hh"
#include "./dataType.hh"
#include "./spacedKmer.hh"

// CLASS 
// ***************************************************************************************************************************
// Name: Hashtable
// Purpose: The Hashtablestorage class allows a simple data structure for the adding and querying of kmers.
// History:
// Special requirements: The following design is supported and tested for 64-bits operating systems.
// Implementation notes: Hash-table with chaining lists for k-mer storage.
//	 		 This design is efficient and allows fast queries for small, large or big datasets (hundred to billion).
// 

class Hashtable
{
	protected:
		size_t                          m_kmerSize;
		uint8_t				m_k;
		size_t                          m_localIndex;
	public:
		Hashtable();
		Hashtable(const size_t& _kmerSize);
		~Hashtable();
		size_t Size() const;
		size_t kmerSize() const 
		{	return m_kmerSize;} 
};

template <typename HKMERr, typename ELMTr>
class EHashtable: public Hashtable
{
	private:
		hTable<HKMERr, ELMTr>			m_hTable;
		std::vector< spacedKmer >		m_seeds;
		std::vector< std::string >		m_Labels;
		std::map< std::string, ILBL > 		m_mapLbls;
		std::map< std::string, ILBL >::iterator m_it;
		bool					m_checkCentromere;
	
	public:
		EHashtable();	
		EHashtable(const size_t& _kmerSize);
		EHashtable(const size_t& _kmerSize,
                        const std::vector< std::string >&       _labelsA,
                        const std::vector< std::string >&       _labelsC
			);
		EHashtable(const size_t& _kmerSize, 
			const std::vector< std::string >& 	_labelsA, 
			const std::vector< std::string >& 	_labelsC, 
			const std::vector< spacedKmer >& 	_seeds
			);

		~EHashtable();
		bool iskmerLengthValid() const;

		bool addElement(const uint64_t& 		_kmerF, 
			const std::string&			_label,
			const ILBL& 				_rlabel, 
			const size_t& 				_count 
			);

		bool addElement(const uint64_t& 		_kmerF, 
			const std::string& 			_label, 
			const size_t& 				_count
			);
		bool addElement(const uint64_t& 		_kmerF, 
			const uint64_t&	 			_kmerR, 
			const std::string& 			_label, 
			const size_t& 				_count
			);
		bool addElement(const std::string& 		_kmerI, 
			const std::string& 			_label, 
			const size_t& 				_count
			);
		bool addElement(const std::string& 		_kmer, 
			const std::string& 			_label
			);
		bool addElement(const std::string& 		_kmer);

		bool clear();

		bool getTargetID(const std::string&		_label,
			ILBL&					_id
			)
		{
			m_it = m_mapLbls.find(_label);
			if (m_it != m_mapLbls.end())
			{
				_id = m_it->second;
				return true;
			}
			return false;
		}

		void SortAllHashTable(const size_t& 		_iteratorPos = 0)
		{	m_hTable.sortall(_iteratorPos);	}
		
		bool SaveMultiple(const std::vector<std::string>& _filesHT, 
			const std::vector<std::string>& 	_labels, 
			const IOCCR& 				_multiplicity = 1,
			const size_t&				_minCount = 0
			);

		bool SaveIntersectionMultiple(const std::vector<std::string>& _filesHT, 
			const std::vector<std::string>& 	_labels
			);

		bool RemoveCommon(const std::vector<std::string>& _labels_c, 
			const size_t& 				_minCount = 0
			);

		bool queryElement(const IKMER& 			_ikmer, 
			ILBL& 					_iLabel
			) const;

		bool queryElement(const uint64_t& 		_kmerI, 
			ILBL& 					_iLabel
			)
		{	return m_hTable.find(_kmerI, _iLabel);	}		

		bool querySpacedElement(const uint8_t* 		_map, 
			const size_t& 				_i, 
			const size_t& 				_iMax, 
			ILBL& 					_iLabel, 
			const size_t& 				_idHt = 0
			) const
		{
			uint64_t km = 0, kmr = 0, skm = 0;
			if (m_seeds[_idHt].isFwdValid(_map, _i, _iMax, km))
			{
				getSpacedSeed(m_seeds[_idHt].getName(), km, skm);
				if (m_hTable.findFwd(skm, _iLabel, _idHt))
                                {	return true;	}
			}
			km = 0;
			kmr = 0;
			skm = 0;
			if (m_seeds[_idHt].isRvsValid(_map, _i, _iMax, km))
                        {
				getReverseComplement(km, m_seeds[_idHt].getLength(), kmr);
				getSpacedSeed(m_seeds[_idHt].getName(), kmr, skm);
				return (m_hTable.findFwd(skm, _iLabel, _idHt));
			}
			return false; 
		}

		bool queryElement(const uint64_t& 		_kmerIF, 
			const uint64_t& 			_kmerIR, 
			ILBL& 					_iLabel
			) const	
		{	return m_hTable.find(_kmerIF, _kmerIR, _iLabel);	}	

		uint64_t Write(const char * 			_filename, 
			const size_t& 				_iteratorPos, 
			const bool& 				_clearAfter = true
			)
		{	return m_hTable.write(_filename, _iteratorPos, _clearAfter);	}

		bool Read(const char * 				_filename, 
			size_t& 				_sizefile, 
			const size_t& 				_nbCPU, 
			const size_t& 				_modCollision = 1, 
			const bool& 				_mmapLoading = false
			) 
		{	return m_hTable.read(_filename, _sizefile, _nbCPU, _modCollision, _mmapLoading); 	}

		bool Add(const char * 				_filename, 
			size_t& 				_sizefile, 
			const size_t& 				_modCollision = 1, 
			const size_t& 				_idHt = 0
			)
                {       return m_hTable.add(_filename, _sizefile, _modCollision, _idHt);        }

		bool AddDB(const std::vector<std::string>& 	_filesname, 
			size_t& 				_sizefile, 
			const ITYPE& 				_modCollision = 1
			)
                {       return m_hTable.addDB(_filesname, _sizefile, _modCollision);        }

		size_t Size() const
		{	return m_hTable.Load();	}
		
		void Load(const string& 			_fileHT, 
			const std::string& 			_label, 
			const ITYPE& 				_minCount = 0
			);
};

#endif //HASHTABLESTORAGE_HH

// @Author: Rachid Ounit
// @Name: Hashtablestorage class
// @Date: 2013
//
#include <iostream>
#include <iomanip>
#include <cmath>
#include "file.hh"
#include "kmersConversion.hh"
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
Hashtable::Hashtable(): m_kmerSize(0), m_k(0), m_localIndex(0)
{}

Hashtable::Hashtable(const size_t& _kmerSize): m_kmerSize(_kmerSize),  m_k((uint8_t) _kmerSize), m_localIndex(0)
{}

size_t Hashtable::Size() const
{       return m_localIndex;}

Hashtable::~Hashtable(){}

///////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::iskmerLengthValid() const
{
	double beta = log(((double)HTSIZE))/log(4.0);
	if ( sizeof(HKMERr)*4*1 + beta < m_kmerSize )
	{
		cerr << "The k-mer length (" << m_kmerSize  << ") requested is too big against hash-table settings." << endl;
		cerr << "Please choose a k-mer length smaller than "<< (size_t) (( sizeof(HKMERr)*4*1*1.0 + beta ));
		cerr << " or increase the value of HTSIZE and/or change the type of the Cell key (files dataType.hh, parameters.hh)." << endl;
		cerr << "The program must exit now." << endl;
		return false;
	}
	return true;
}

	template <typename HKMERr, typename ELMTr>
EHashtable<HKMERr, ELMTr>::EHashtable(): Hashtable(), m_Labels(1)
{
	m_Labels[0] = "";
	m_mapLbls[""] = 0;
}
	template <typename HKMERr, typename ELMTr>
EHashtable<HKMERr, ELMTr>::EHashtable(const size_t& _kmerSize): Hashtable(_kmerSize), m_Labels(1), m_hTable((uint8_t) _kmerSize)
{	m_Labels[0] = "";
	m_mapLbls[""] = 0;
}

template <typename HKMERr, typename ELMTr>
EHashtable<HKMERr, ELMTr>::EHashtable(const size_t& _kmerSize, const std::vector< std::string >& _labelsA, const std::vector< std::string >& _labelsC): Hashtable(_kmerSize), m_Labels(1+_labelsA.size()+_labelsC.size()), m_hTable((uint8_t) _kmerSize), m_checkCentromere(_labelsC.size()>0)
{
        for(size_t t = 0; t < _labelsA.size(); t++)
        {
                m_Labels[t] = _labelsA[t];
                m_mapLbls[_labelsA[t]] = t;
        }
        for(size_t t = 0; t < _labelsC.size(); t++)
        {
                m_mapLbls[_labelsC[t]] = t + _labelsA.size();
                m_Labels[t + _labelsA.size()] = _labelsC[t];
        }
        m_Labels[_labelsA.size() +  _labelsC.size()] = "";
        m_mapLbls[""] = _labelsA.size() +  _labelsC.size();
}


template <typename HKMERr, typename ELMTr>
EHashtable<HKMERr, ELMTr>::EHashtable(const size_t& _kmerSize, const std::vector< std::string >& _labelsA, const std::vector< std::string >& _labelsC, const std::vector< spacedKmer >& _seeds): Hashtable(_kmerSize), m_Labels(1+_labelsA.size()+_labelsC.size()), m_hTable((uint8_t) _kmerSize)
{
	for(size_t t = 0; t < _labelsA.size(); t++)
	{
		m_Labels[t] = _labelsA[t];
		m_mapLbls[_labelsA[t]] = t;
	}
	for(size_t t = 0; t < _labelsC.size(); t++)
	{
		m_mapLbls[_labelsC[t]] = t + _labelsA.size();
		m_Labels[t + _labelsA.size()] = _labelsC[t];
	}
	m_Labels[_labelsA.size() +  _labelsC.size()] = "";
	m_mapLbls[""] = _labelsA.size() +  _labelsC.size();
	
	for(size_t t = 0; t < _seeds.size(); t++)
		m_seeds.push_back(_seeds[t]);
}

///////////////////////////////////////////////////////////////
	template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::RemoveCommon(const std::vector<std::string>& _labels_c, const size_t& _minCount)
{
	m_hTable.resetIterator();
	ELMTr e;
	uint64_t nbElement = 0;
	uint64_t nbSpec = 0;
	string Lbl, centro, candidate, c_label;
	ILBL i_lbl = 0;
	bool found = false;
	std::map<string,ILBL>::iterator it_Lbl;
	const bool centromereRequested = _labels_c.size() > 0 ;
	while (m_hTable.next())
	{
		m_hTable.elementIterator(e);
		if (e.GetMultiplicity() == 1 && e.GetCount() > _minCount)
		{
			m_hTable.markElementAtIterator();
			nbElement++;
			continue;
		}
		if (centromereRequested && e.GetMultiplicity() == 2  && e.GetCount() > _minCount)
		{
			m_hTable.markElementAtIterator();
			Lbl = m_Labels[e.GetLabel()];
			centro = Lbl.substr(0, Lbl.size()-1);
			found = false;
			for(size_t t = 0; !found && t < _labels_c.size(); t++)
			{
				if (_labels_c[t].size() == Lbl.size())
				{
					candidate = _labels_c[t].substr(0, Lbl.size()-1);
					found = centro.compare(candidate) == 0;
					c_label = _labels_c[t];
				}
			}
			if (found)
			{
				it_Lbl = m_mapLbls.find(c_label);
				m_hTable.updateLabelAtIterator(it_Lbl->second);
				nbSpec++;
			}
		}
	}

	cerr <<"Removal of common k-mers done: "<<nbElement+nbSpec<<" specific "<< m_kmerSize<<"-mers found";
	if (centromereRequested)
	{	cerr << " including "<<nbSpec<<" in centromeres."<<endl;	}
	else 
	{	cerr << "." << endl;	}
	return true;
}

	template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::SaveMultiple(const std::vector<std::string>& _filesHT, const std::vector<std::string>& _labels, const IOCCR& _multiplicity, const size_t& _minCount)
{
	size_t kmerSize = m_kmerSize;
	if (_filesHT.size() < 1 || _labels.size() < 1 || _filesHT.size() != _labels.size() )
	{
		return false;
	}

	// Counting how many kmers satisfies the requirement (label X multiplicity
	vector<FILE*> fds(_filesHT.size());
	std::map<string, size_t> stringToIndex;

	for(size_t t = 0; t < _filesHT.size() ; t++)
	{
		// Writing/Saving data in file
		fds[t] = fopen(_filesHT[t].c_str(), "w+");
		fprintf(fds[t], "#Target specific k-mers labeled %s and appearing strictly more than %lu times.\n", _labels[t].c_str(), _minCount);
		fprintf(fds[t], "#IKMER ICOUNT %lu-MER \n#\n", kmerSize);
		stringToIndex[ _labels[t]] = t;
	}
	m_hTable.resetIterator();
	string kmer;
	while (m_hTable.next())
	{
		ELMTr e;
		uint64_t kmerIndex;
		m_hTable.elementIterator(kmerIndex, e);

		if (e.GetMultiplicity() <= _multiplicity /*&& e.GetCount() > _minCount */)
		{
			if (stringToIndex.find(m_Labels[e.GetLabel()]) != stringToIndex.end())
			{
				size_t i_lbl = stringToIndex.find(m_Labels[e.GetLabel()])->second;
				IndexTovector(kmerIndex, kmerSize, kmer);
				fprintf(fds[i_lbl],"%" PRIu64 "\t%lu\t%s\n", kmerIndex, e.GetCount(), kmer.c_str());
			}
		}
	}

	//Closing all files
	for(size_t t = 0; t < _filesHT.size() ; t++)
	{
		fclose(fds[t]);
	}

	return true;
}

	template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::SaveIntersectionMultiple(const std::vector<std::string>& _filesHT, const std::vector<std::string>& _labels)
{
	size_t kmerSize =  m_kmerSize;
	if (_filesHT.size() < 1 || _labels.size() < 1 || _filesHT.size() != _labels.size() )
	{
		return false;
	}

	//vector<ITYPE> f_indexes(_filesHT.size());
	vector<FILE*> fds(_filesHT.size());
	for(size_t t = 0; t < _filesHT.size() ; t++)
	{
		// Writing/Saving data in file
		fds[t] = fopen(_filesHT[t].c_str(), "w+");
		string name= _labels[t];
		string	centro = name.substr(0, name.size()-1);
		fprintf(fds[t], "#K-mers specific to chromosome-centromere %s\n", centro.c_str());
		fprintf(fds[t], "#IKMER ICOUNT %lu-MER\n#\n", kmerSize);
	}

	ELMTr e;
	uint64_t kmerIndex;
	m_hTable.resetIterator();
	string kmer;
	while (m_hTable.next())
	{
		m_hTable.elementIterator(kmerIndex, e);

		bool label_found = false, sameChr = false;
		size_t i_lbl = 0;
		while (!label_found && i_lbl < _labels.size())
		{
			if ( e.GetMultiplicity() == 2 )
			{
				string label = _labels[i_lbl];
				string Lbl = m_Labels[e.GetLabel()];
				bool sameChr = Lbl.size() == label.size();
				for(size_t t = 0 ; sameChr && t < Lbl.size() - 1 ; t++ )
				{       sameChr = sameChr && Lbl[t] == label[t];}

				if (sameChr)
				{
					label_found = true;
					IndexTovector(kmerIndex, kmerSize, kmer);
					fprintf(fds[i_lbl], "%" PRIu64 "\t%lu\t%s\n",kmerIndex,e.GetCount(), kmer.c_str()); 
				}
				else
				{
					i_lbl++;
				}
			}
			else
			{
				label_found = true;
			}
		}
	}

	//Closing all files
	for(size_t t = 0; t < _filesHT.size() ; t++)
	{
		fclose(fds[t]);
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <typename HKMERr, typename ELMTr>
EHashtable<HKMERr, ELMTr>::~EHashtable()
{
}

	template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::clear()
{
	m_localIndex = 0;
	m_kmerSize = 0;
	m_hTable.clear();
}

	template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::addElement(const std::string&   _kmer)
{
	return addElement(_kmer, "", 1);
}

	template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::addElement(const std::string& _kmer, const string& _label)
{
	return addElement(_kmer,  _label, 1);
}

        template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::addElement(const uint64_t& _kmerF, const std::string& _label, const ILBL& _rlabel, const size_t& _count)
{
	if (m_checkCentromere)
	{	return addElement(_kmerF, _label, _count);	}

        // Checking the forward first
        size_t e_x = 0, e_y = 0;
        ILBL e_l = 0;
        IOCCR mult;
        ICount count;
        if (m_hTable.find(_kmerF, e_x, e_y, e_l, mult, count))
        {
                m_hTable.updateElement(e_x, e_y, _count, _rlabel == e_l);
                return true;
        }
        // Checking the reverse
        uint64_t _kmerR = _kmerF;
        _kmerR = ((_kmerR >> 2)  & 0x3333333333333333UL) | ((_kmerR & 0x3333333333333333UL) << 2);
        _kmerR = ((_kmerR >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_kmerR & 0x0F0F0F0F0F0F0F0FUL) << 4);
        _kmerR = ((_kmerR >> 8)  & 0x00FF00FF00FF00FFUL) | ((_kmerR & 0x00FF00FF00FF00FFUL) << 8);
        _kmerR = ((_kmerR >> 16) & 0x0000FFFF0000FFFFUL) | ((_kmerR & 0x0000FFFF0000FFFFUL) << 16);
        _kmerR = ( _kmerR >> 32                        ) | ( _kmerR                        << 32);
        _kmerR = (((uint64_t)-1) - _kmerR) >> (64 - (m_k << 1));

        if (m_hTable.find(_kmerR, e_x, e_y, e_l, mult, count))
        {
		m_hTable.updateElement(e_x, e_y, _count, _rlabel == e_l);
                return true;
        }

        // Element is not in the table already. Then adding now.
        std::map<string,ILBL>::iterator it_Lbl = m_mapLbls.find(_label);
        m_hTable.insert(_kmerF, it_Lbl->second, _count);
        m_localIndex++;
        return true;
}

	template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::addElement(const uint64_t& _kmerF, const std::string& _label, const size_t& _count)
{
	// Checking the forward first
	size_t e_x = 0, e_y = 0;
	ILBL e_l = 0;
	IOCCR mult;
	ICount count;
	string Lbl;
	bool upLbl, isSameLbl;
	if (m_hTable.find(_kmerF, e_x, e_y, e_l, mult, count))
	{
		Lbl = m_Labels[e_l];
		upLbl = _label[0] == Lbl[0] && Lbl.size() == _label.size() ;
		for(size_t t = 1 ; upLbl && t < Lbl.size() - 1; t++)
		{       upLbl = upLbl && _label[t] == Lbl[t];}
		isSameLbl = upLbl && _label[Lbl.size()-1] == Lbl[Lbl.size()-1];
		m_hTable.updateElement(e_x, e_y, _count, !upLbl, isSameLbl);
		return true;
	}
	// Checking the reverse
	uint64_t _kmerR = _kmerF;
	_kmerR = ((_kmerR >> 2)  & 0x3333333333333333UL) | ((_kmerR & 0x3333333333333333UL) << 2);
	_kmerR = ((_kmerR >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_kmerR & 0x0F0F0F0F0F0F0F0FUL) << 4);
	_kmerR = ((_kmerR >> 8)  & 0x00FF00FF00FF00FFUL) | ((_kmerR & 0x00FF00FF00FF00FFUL) << 8);
	_kmerR = ((_kmerR >> 16) & 0x0000FFFF0000FFFFUL) | ((_kmerR & 0x0000FFFF0000FFFFUL) << 16);
	_kmerR = ( _kmerR >> 32                        ) | ( _kmerR                        << 32);
	_kmerR = (((uint64_t)-1) - _kmerR) >> (64 - (m_k << 1));

	if (m_hTable.find(_kmerR, e_x, e_y, e_l, mult, count))
	{
		Lbl = m_Labels[e_l];
		upLbl = _label[0] == Lbl[0] && Lbl.size() == _label.size();
		for(size_t t = 1 ; upLbl && t < Lbl.size() - 1 ; t++)
		{       upLbl = upLbl && _label[t] == Lbl[t];}
		isSameLbl = upLbl && _label[Lbl.size()-1] == Lbl[Lbl.size()-1];
		m_hTable.updateElement(e_x, e_y, _count, !upLbl, isSameLbl);
		return true;
	}

	// Element is not in the table already. Then adding now.
	std::map<string,ILBL>::iterator it_Lbl = m_mapLbls.find(_label);
	m_hTable.insert(_kmerF, it_Lbl->second, _count);
	m_localIndex++;
	return true;
}

	template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::addElement(const uint64_t& _kmerF, const uint64_t& _kmerR, const std::string& _label, const size_t& _count)
{
	// Checking the forward first
	size_t e_x = 0, e_y = 0;
	ILBL e_l = 0;
	IOCCR mult;
	ICount count;
	string Lbl;
	bool upLbl, isSameLbl;
	if (m_hTable.find(_kmerF, e_x, e_y, e_l, mult, count))
	{
		Lbl = m_Labels[e_l];
		upLbl = _label[0] == Lbl[0] && Lbl.size() == _label.size() ;
		for(size_t t = 1 ; upLbl && t < Lbl.size() - 1; t++)
		{       upLbl = upLbl && _label[t] == Lbl[t];}
		isSameLbl = upLbl && _label[Lbl.size()-1] == Lbl[Lbl.size()-1];
		m_hTable.updateElement(e_x, e_y, _count, !upLbl, isSameLbl);
		return true;
	}
	// Checking the reverse
	if (m_hTable.find(_kmerR, e_x, e_y, e_l, mult, count))
	{
		Lbl = m_Labels[e_l];
		upLbl = _label[0] == Lbl[0] && Lbl.size() == _label.size();
		for(size_t t = 1 ; upLbl && t < Lbl.size() - 1 ; t++)
		{       upLbl = upLbl && _label[t] == Lbl[t];}
		isSameLbl = upLbl && _label[Lbl.size()-1] == Lbl[Lbl.size()-1];
		m_hTable.updateElement(e_x, e_y, _count, !upLbl, isSameLbl);
		return true;
	}

	// Element is not in the table already. Then adding now.
	std::map<string,ILBL>::iterator it_Lbl = m_mapLbls.find(_label);
	m_hTable.insert(_kmerF, it_Lbl->second, _count);
	m_localIndex++;
	return true;
}

	template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::addElement(const std::string& _kmerI, const std::string& _label, const size_t& _count)
{
	string _kmer;
	_kmer = _kmerI;

	// Checking the forward first
	uint64_t kmerIndex = 0;
	vectorToIndex(_kmer, kmerIndex);
	size_t e_x = 0, e_y = 0;
	ILBL e_l = 0;
	IOCCR mult;
	ICount count;

	if (m_hTable.find(kmerIndex, e_x, e_y, e_l, mult, count))
	{
		string Lbl = m_Labels[e_l];
		bool upLbl = _label[0] == Lbl[0] && Lbl.size() == _label.size() ;
		for(size_t t = 1 ; upLbl && t < Lbl.size() - 1; t++)
		{	upLbl = upLbl && _label[t] == Lbl[t];}
		bool isSameLbl = upLbl && _label[Lbl.size()-1] == Lbl[Lbl.size()-1];
		m_hTable.updateElement(e_x, e_y, _count, !upLbl, isSameLbl);
		return true;
	}
	uint64_t rev_kmerIndex = 0;
	// Checking the reverse
	getReverseComplement(_kmer, rev_kmerIndex);

	if (m_hTable.find(rev_kmerIndex, e_x, e_y, e_l, mult, count))
	{
		string Lbl = m_Labels[e_l];
		bool upLbl = _label[0] == Lbl[0] && Lbl.size() == _label.size();
		for(size_t t = 1 ; upLbl && t < Lbl.size() - 1 ; t++)
		{	upLbl = upLbl && _label[t] == Lbl[t];}
		bool isSameLbl = upLbl && _label[Lbl.size()-1] == Lbl[Lbl.size()-1];

		m_hTable.updateElement(e_x, e_y, _count, !upLbl, isSameLbl);
		return true;
	}

	// Element is not in the table already. Then adding now.
	std::map<string,ILBL>::iterator it_Lbl = m_mapLbls.find(_label);
	m_hTable.insert(kmerIndex, it_Lbl->second, _count);
	m_localIndex++;

	return true;
}

template <typename HKMERr, typename ELMTr>
bool EHashtable<HKMERr, ELMTr>::queryElement(const IKMER& _ikmer, ILBL& _iLabel) const
{
	if (m_hTable.find(_ikmer, 0, 1, _iLabel))
	{
		return true;
	}
	return m_hTable.find(_ikmer, 2, 3, _iLabel);
}

	template <typename HKMERr, typename ELMTr>
void EHashtable<HKMERr, ELMTr>::Load(const string& _fileHT, const std::string& _label, const ITYPE& _minCount)
{
	FILE * fd = fopen(_fileHT.c_str(), "r");

	if (fd == NULL)
	{
		cerr << "Failed to open " << _fileHT << endl;
		return;
	}
	string line;

	// Get vector size
	getFirstElementInLineFromFile(fd, line);

	getFirstElementInLineFromFile(fd, line);
	//m_localIndex = 0;

	// Get kmer size
	getFirstElementInLineFromFile(fd, line);
	size_t kSize(atol(line.c_str()));
	m_kmerSize = kSize;

	if (!iskmerLengthValid())
	{       exit(-1);}

	ITYPE count = 0;
	uint64_t kIndex = 0;
	// Populate kmers vector and map
	std::map< string, ILBL >::iterator it_Lbl;
	it_Lbl = m_mapLbls.find(_label);
	while ( getFirstAndSecondElementInLine(fd, kIndex, count) )
	{
		if (count > _minCount)
		{
			m_hTable.insert(kIndex, it_Lbl->second);
			m_localIndex++;
		}
	}
	fclose(fd);
}

