#include <Rcpp.h>
using namespace std::placeholders;
using namespace Rcpp;

// The functions below are lifted out from mutscan
// https://github.com/fmicompbio/mutscan/blob/master/src/digestFastqs.cpp

// Translate nucleotide string to amino acid string
std::string translateString(std::string& s) {
    const char* tr = "KQE*TPASRRG*ILVLNHDYTPASSRGCILVFKQE*TPASRRGWMLVLNHDYTPASSRGCILVF";
    std::string aa = "";
    size_t i = 0, j = 0;
    const int pow4[] = {1, 4, 16};
    int codonVal = 0;
    
    for (i = 0; i < s.size(); i++) {
        if (s[i] == '_') {
            j = 0;
            aa.push_back('_');
        } else {
            switch(s[i]) {
            case 'A':
            case 'a':
                // codonVal += 0 * pow4[j];
                break;
                
            case 'C':
            case 'c':
                codonVal += 1 * pow4[j];
                break;
                
            case 'G':
            case 'g':
                codonVal += 2 * pow4[j];
                break;
                
            case 'T':
            case 't':
            case 'u':
            case 'U':
                codonVal += 3 * pow4[j];
                break;
                
            default:
                codonVal += 64;
            break;
            }
            
            if (j == 2) {
                if (codonVal > 63) {
                    aa.push_back('X');
                } else {
                    aa.push_back(tr[codonVal]);
                }
                j = 0;
                codonVal = 0;
            } else {
                j++;
            }
        }
    }
    return aa;
}

// Split string by delimiter (code from https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/)
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Compare positions (to get the right order in the final mutant name)
bool compareCodonPositions(std::string a, std::string b, const char mutNameDelimiter) {
    int posa = std::stoi(split(a, mutNameDelimiter)[1]);
    int posb = std::stoi(split(b, mutNameDelimiter)[1]);
    return (posa < posb);
}

// Convert a nucleotide sequence to a mutant name, in the style of mutscan
// [[Rcpp::export]]
std::string seqToName(const std::string varSeq, const std::string wtSeq,
                      const std::string cType, const std::string codonPrefix) {

    std::set<std::string> mutatedCodons, mutatedBases, mutatedAAs;
    std::set<std::string>::iterator mutatedCodonIt;
    std::string varCodon, wtCodon, varAA, wtAA;
    std::string mutantName = "", mutantNameAA = "";
    std::string mutNameDelimiter = ".";
    
    for (size_t i = 0; i < varSeq.length(); i++) {
        if (varSeq[i] != wtSeq[i]) { // found mismatching base
            // consider codons
            varCodon = varSeq.substr((int)(i / 3) * 3, 3);
            wtCodon = wtSeq.substr((int)(i / 3) * 3, 3);
            mutatedCodons.insert(codonPrefix + mutNameDelimiter +
                std::to_string((int)(i / 3) + 1) + mutNameDelimiter +
                varCodon +
                std::string("_"));
            // consider individual bases
            mutatedBases.insert(codonPrefix + mutNameDelimiter +
                std::to_string((int)(i) + 1) + mutNameDelimiter +
                varSeq.substr((int)(i), 1) +
                std::string("_"));
            varAA = translateString(varCodon);
            wtAA = translateString(wtCodon);
            if (varAA != wtAA) {
                // add to mutatedAA
                mutatedAAs.insert(codonPrefix + mutNameDelimiter + 
                    std::to_string((int)(i / 3) + 1) + mutNameDelimiter +
                    varAA + 
                    std::string("_"));
            }
        }
    }

    // Create name for mutant
    std::vector<std::string> mutatedCodonsOrBasesSorted, mutatedAAsSorted;
    if (cType == "codon") {
        mutatedCodonsOrBasesSorted.assign(mutatedCodons.begin(), mutatedCodons.end());
    } else {
        mutatedCodonsOrBasesSorted.assign(mutatedBases.begin(), mutatedBases.end());
    }
    std::sort(mutatedCodonsOrBasesSorted.begin(), mutatedCodonsOrBasesSorted.end(),
              std::bind(compareCodonPositions, _1, _2, *(mutNameDelimiter.c_str())));
    mutatedAAsSorted.assign(mutatedAAs.begin(), mutatedAAs.end());
    std::sort(mutatedAAsSorted.begin(), mutatedAAsSorted.end(),
              std::bind(compareCodonPositions, _1, _2, *(mutNameDelimiter.c_str())));
    for (size_t i = 0; i < mutatedCodonsOrBasesSorted.size(); i++) {
        mutantName += mutatedCodonsOrBasesSorted[i];
    }
    for (size_t i = 0; i < mutatedAAsSorted.size(); i++) {
        mutantNameAA += mutatedAAsSorted[i];
    }
    // if no mutant codons, name as <codonPrefix>.0.WT
    if (mutatedCodonsOrBasesSorted.size() == 0) {
        mutantName += codonPrefix + mutNameDelimiter + "0" + mutNameDelimiter + "WT_";
    }
    if (mutatedAAsSorted.size() == 0) {
        mutantNameAA += codonPrefix + mutNameDelimiter + "0" + mutNameDelimiter + "WT_";
    }
    
    mutantName.pop_back();
    mutantNameAA.pop_back();
    
    if (cType == "aa") {
        return mutantNameAA;
    } else {
        return mutantName;
    }
}

