#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE seqantest

#include <boost/test/unit_test.hpp>
#include <seqan/graph_msa.h>
#include <seqan/align.h>
using namespace seqan;

using namespace std;

struct Fixture{
	
	Fixture(){
	}
	~Fixture(){
	}
};

struct Zeichen {
	int i;
	Zeichen(const char& c){i =c;}
	Zeichen(){}
	operator unsigned char () const { return '0'+i;}
};
bool operator==(Zeichen const & z, Zeichen const & y){ return z.i == y.i;}

inline double                                                                                                                                    
score( Score<int, Zeichen> const & me,                                                                                                       
		      Zeichen const & left,                                                                                                                            
				Zeichen const & right)                                                                                                                           
{                                                                                                                                                
	return 0.0;
} 

BOOST_FIXTURE_TEST_SUITE( suite, Fixture )

BOOST_AUTO_TEST_CASE( test1 )
{
	 Align<String<Zeichen> > ali;
	 String<Zeichen> x,y;
	 appendValue(x, Zeichen(1));
	 appendValue(x, Zeichen(2));
	 appendValue(x, Zeichen(3));
	 appendValue(x, Zeichen(4));
	 appendValue(y, Zeichen(1));
	 appendValue(y, Zeichen(5));
	 appendValue(y, Zeichen(3));
	 appendValue(y, Zeichen(4));
	 appendValue(rows(ali), x);
	 appendValue(rows(ali), y);
	 int score = localAlignment(ali, Score<int>(3,0,1), SmithWaterman());
	 cout << "Score = " << score << endl;
	 cout << ali; 
	 BOOST_CHECK_EQUAL(11,score);
}




BOOST_AUTO_TEST_CASE( test2 )
{
    typedef String<AminoAcid> TString;
    typedef StringSet<TString, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, unsigned int, Default> > TGraph;
    TString str1 = "GARFIELDTHELASTFATCAT";
    TString str2 = "GARFIELDTHEFASTCAT";
    TString str3 = "GARFIELDTHEVERYFASTCAT";
    TString str4 = "THEFATCAT";
    TStringSet strSet;
    assignValueById(strSet, str1);
    assignValueById(strSet, str2);
    assignValueById(strSet, str3);
    assignValueById(strSet, str4);
    Graph<Alignment<TStringSet, void, WithoutEdgeId> > gOut(strSet);
    globalAlignment(strSet, gOut, MSA_Protein() );
    std::cout << gOut << std::endl;
} 



BOOST_AUTO_TEST_SUITE_END()

