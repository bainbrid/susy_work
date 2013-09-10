#include <vector>
#include <string>

using namespace std;

typedef unsigned int uint;
typedef vector<double> vdouble;
typedef vector<string> vstring;
typedef vector<vdouble> vvdouble;

class Pdf {
public:
  double m0;
  double m12;
  std::string name;
  vdouble eff;
  vdouble err;
};

void pdf() {

  Pdf cteq;
  cteq.m0 = 500.;
  cteq.m12 = 500.;
  cteq.name = "CTEQ 6.1";
  cteq.eff.push_back(0.10002);
  cteq.eff.push_back(0.06212);
  cteq.eff.push_back(0.05053);
  cteq.eff.push_back(0.04733);
  cteq.eff.push_back(0.06611);
  cteq.eff.push_back(0.07206);
  cteq.eff.push_back(0.09707);
  cteq.eff.push_back(0.70401);
  cteq.err.push_back(0.000252);
  cteq.err.push_back(0.000144);
  cteq.err.push_back(0.000040);
  cteq.err.push_back(0.000059);
  cteq.err.push_back(0.000176);
  cteq.err.push_back(0.000201);
  cteq.err.push_back(0.000212);
  cteq.err.push_back(0.000891);

  Pdf mtsw;
  mtsw.m0 = 500.;
  mtsw.m12 = 500.;
  mtsw.name = "MSTW 5006 NLO 68%";
  mtsw.eff.push_back(0.011320);
  mtsw.eff.push_back(0.006894);
  mtsw.eff.push_back(0.005171);
  mtsw.eff.push_back(0.005234);
  mtsw.eff.push_back(0.006822);
  mtsw.eff.push_back(0.007392);
  mtsw.eff.push_back(0.008829);
  mtsw.eff.push_back(0.065823);
  mtsw.err.push_back(0.000053);
  mtsw.err.push_back(0.000037);
  mtsw.err.push_back(0.000013);
  mtsw.err.push_back(0.000023);
  mtsw.err.push_back(0.000037);
  mtsw.err.push_back(0.000045);
  mtsw.err.push_back(0.000056);
  mtsw.err.push_back(0.000215);

  Pdf nnpdf;
  nnpdf.m0 = 500.;
  nnpdf.m12 = 500.;
  nnpdf.name = "NNPDF";
  
  
		      
nnpdf.eff.push_back(0.010807  0.000330
nnpdf.eff.push_back(0.006750  0.000190
nnpdf.eff.push_back(0.005031  0.000112
nnpdf.eff.push_back(0.004898  0.000141
nnpdf.eff.push_back(0.007309  0.000283
nnpdf.eff.push_back(0.007990  0.000363
nnpdf.eff.push_back(0.004731  0.001331
nnpdf.eff.push_back(0.066915  0.001995


nnpdf.eff.push_back(0.010807  0.000330
nnpdf.eff.push_back(0.006750  0.000190
nnpdf.eff.push_back(0.005031  0.000112
nnpdf.eff.push_back(0.004898  0.000141
nnpdf.eff.push_back(0.007309  0.000283
nnpdf.eff.push_back(0.007990  0.000363
nnpdf.eff.push_back(0.004731  0.001331
nnpdf.eff.push_back(0.066915  0.001995



}
