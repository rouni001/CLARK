#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#include "src/file.hh"
#include "src/HashTop.hh"

static int failures = 0;

static void expectTrue(const bool condition, const std::string& message)
{
	if (!condition)
	{
		std::cerr << "FAIL: " << message << std::endl;
		failures++;
	}
}

static void expectString(const std::string& actual, const std::string& expected, const std::string& message)
{
	if (actual != expected)
	{
		std::cerr << "FAIL: " << message << " expected <" << expected << "> got <" << actual << ">" << std::endl;
		failures++;
	}
}

static void expectITYPE(const ITYPE actual, const ITYPE expected, const std::string& message)
{
	if (actual != expected)
	{
		std::cerr << "FAIL: " << message << " expected <" << expected << "> got <" << actual << ">" << std::endl;
		failures++;
	}
}

static std::string joinPath(const std::string& dir, const std::string& name)
{
	if (!dir.empty() && dir[dir.size() - 1] == '/')
	{
		return dir + name;
	}
	return dir + "/" + name;
}

static void writeText(const std::string& path, const std::string& text)
{
	std::ofstream out(path.c_str(), std::ios::binary);
	out << text;
	out.close();
}

static std::string readText(const std::string& path)
{
	std::ifstream in(path.c_str(), std::ios::binary);
	std::ostringstream buffer;
	buffer << in.rdbuf();
	return buffer.str();
}

static void testStringElementParsing()
{
	std::vector<std::string> elements;
	getElementsFromLine(" alpha,beta\tgamma  delta\r\n", 3, elements);
	expectITYPE((ITYPE) elements.size(), 3, "comma/space parser returns max requested elements");
	expectString(elements[0], "alpha", "first comma-separated element");
	expectString(elements[1], "beta", "second comma-separated element");
	expectString(elements[2], "gamma", "third comma-separated element");

	std::vector<char> seps;
	seps.push_back('>');
	seps.push_back('/');
	seps.push_back(' ');
	seps.push_back('\t');
	getElementsFromLine(">read1/2 comment", seps, elements);
	expectITYPE((ITYPE) elements.size(), 3, "custom separator parser extracts all tokens");
	expectString(elements[0], "read1", "custom parser read id");
	expectString(elements[1], "2", "custom parser mate number");
	expectString(elements[2], "comment", "custom parser trailing token");
}

static void testCharElementParsing()
{
	char line[] = "\t123  456\n789";
	char* ptr = line;
	std::vector<std::string> elements;
	getElementsFromLine(ptr, std::strlen(ptr), 2, elements);
	expectITYPE((ITYPE) elements.size(), 2, "char parser honors max elements");
	expectString(elements[0], "123", "char parser first token");
	expectString(elements[1], "456", "char parser second token");
}

static void testFileReaders(const std::string& tmpDir)
{
	const std::string path = joinPath(tmpDir, "readers.txt");
	writeText(path, "first line\n  second\t42\n1234567890123 17\nACGT 29\n");

	FILE* fd = fopen(path.c_str(), "r");
	std::string line;
	expectTrue(getLineFromFile(fd, line), "getLineFromFile reads first line");
	expectString(line, "first line", "getLineFromFile strips newline");
	expectTrue(getFirstElementInLineFromFile(fd, line), "getFirstElementInLineFromFile reads token");
	expectString(line, "second", "getFirstElementInLineFromFile returns first token");

	uint64_t kIndex = 0;
	ITYPE index = 0;
	expectTrue(getFirstAndSecondElementInLine(fd, kIndex, index), "numeric pair reader succeeds");
	expectTrue(kIndex == 1234567890123ULL, "numeric pair reader parses first value");
	expectITYPE(index, 17, "numeric pair reader parses second value");

	ITYPE freq = 0;
	expectTrue(getFirstAndSecondElementInLine(fd, line, freq), "string/frequency reader succeeds");
	expectString(line, "ACGT", "string/frequency reader parses first token");
	expectITYPE(freq, 29, "string/frequency reader parses frequency");
	expectTrue(!getLineFromFile(fd, line), "getLineFromFile reports EOF");
	expectString(line, "", "getLineFromFile clears line at EOF");
	fclose(fd);
}

static void testFileLifecycleAndPairedMerge(const std::string& tmpDir)
{
	const std::string left = joinPath(tmpDir, "left.fastq");
	const std::string right = joinPath(tmpDir, "right.fastq");
	const std::string merged = joinPath(tmpDir, "merged.fa");
	const std::string transient = joinPath(tmpDir, "transient.txt");

	writeText(left, "@read1/1 comment\nACGT\n+\n!!!!\n@read2/1\nGG\n+\n!!\n");
	writeText(right, "@read1/2 comment\nTTAA\n+\n!!!!\n@read2/2\nCC\n+\n!!\n");
	mergePairedFiles(left.c_str(), right.c_str(), merged.c_str());
	expectString(readText(merged), ">read1\nACGTNNNNTTAA\n>read2\nGGNNNNCC\n", "mergePairedFiles concatenates paired FASTQ reads");

	writeText(transient, "temporary\n");
	expectTrue(validFile(transient.c_str()), "validFile accepts readable files");
	deleteFile(transient.c_str());
	expectTrue(!validFile(transient.c_str()), "deleteFile removes readable files");
	deleteFile(NULL);
}

static void testHashTopCounting()
{
	HashTop top;
	ITYPE label = 99;
	ITYPE count = 99;

	top.getBest(label, count);
	expectITYPE(label, 0, "empty HashTop best label is zero");
	expectITYPE(count, 0, "empty HashTop best count is zero");
	top.getSecondBest(label, count);
	expectITYPE(label, 0, "empty HashTop second-best label is zero");
	expectITYPE(count, 0, "empty HashTop second-best count is zero");
	top.getTotal(count);
	expectITYPE(count, 0, "empty HashTop total is zero");

	top.insert(3);
	top.insert(3);
	top.insert(5);
	top.insert(7, 4);
	top.insert(5, 5);

	top.getBest(label, count);
	expectITYPE(label, 5, "HashTop best label follows highest count");
	expectITYPE(count, 6, "HashTop best count combines weighted and unweighted inserts");
	top.getSecondBest(label, count);
	expectITYPE(label, 7, "HashTop second-best label excludes the best label");
	expectITYPE(count, 4, "HashTop second-best count is tracked");
	top.getTotal(count);
	expectITYPE(count, 12, "HashTop total includes weighted inserts");

	std::string scores;
	top.getScoresLine(8, scores);
	expectString(scores, ",0,0,0,2,0,6,0,4", "HashTop scores line reports counts for active token");

	top.next();
	top.getBest(label, count);
	expectITYPE(label, 0, "HashTop next clears best label");
	expectITYPE(count, 0, "HashTop next clears best count");
	top.getTotal(count);
	expectITYPE(count, 0, "HashTop next clears total");
	top.getScoresLine(8, scores);
	expectString(scores, ",0,0,0,0,0,0,0,0", "HashTop next hides stale table values");

	top.insert(2);
	top.insert(4);
	top.getBest(label, count);
	expectITYPE(label, 2, "HashTop keeps first label as best on ties");
	expectITYPE(count, 1, "HashTop tie count remains one");
	top.getSecondBest(label, count);
	expectITYPE(label, 4, "HashTop second-best reports tied non-best label");
	expectITYPE(count, 1, "HashTop second-best tied count remains one");
}

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <temporary-directory>" << std::endl;
		return 2;
	}

	const std::string tmpDir(argv[1]);
	testStringElementParsing();
	testCharElementParsing();
	testFileReaders(tmpDir);
	testFileLifecycleAndPairedMerge(tmpDir);
	testHashTopCounting();

	if (failures != 0)
	{
		std::cerr << failures << " unit test failure(s)" << std::endl;
		return 1;
	}
	std::cout << "PASS: file.cc/file.hh/HashTop.hh unit tests" << std::endl;
	return 0;
}
