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

struct Zeichen : SimpleType<int, int> {
	int i;
	Zeichen(const char& c){i =c;}
	Zeichen(){}
};
bool operator==(Zeichen const & z, Zeichen const & y){ return z.i == y.i;}

inline double                                                                                                                                    
score( Score<int, Zeichen> const & me,                                                                                                       
		      Zeichen const & left,                                                                                                                            
				Zeichen const & right)                                                                                                                           
{                                                                                                                                                
	return 0.0;
} 

struct MyScore{

};

BOOST_FIXTURE_TEST_SUITE( suite, Fixture )

BOOST_AUTO_TEST_CASE( test1 )
{
	 Align<String<Zeichen> > ali;
	 appendValue(rows(ali), "aphilologicaltheorem");
	 appendValue(rows(ali), "bizarreamphibology");
	 int score = localAlignment(ali, MyScore(), SmithWaterman());
	 cout << "Score = " << score << endl;
	 cout << ali; 
	 BOOST_CHECK_EQUAL(19,score);
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

