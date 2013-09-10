void analyse_closure() {

  //TFile * file = new TFile("RA1_Closure_Spread.root");
  TFile * file[2];
  file[0] = new TFile("RA1_List2_4_6_5_Spread.root");
  file[1] = new TFile("RA1_List2_3_8_8_Spread.root");

  TString histname[8];
  histname[0] = "ht_275";
  histname[1] = "ht_325";
  histname[2] = "ht_375";
  histname[3] = "ht_475";
  histname[4] = "ht_575";
  histname[5] = "ht_675";
  histname[6] = "ht_775";
  histname[7] = "ht_875";

  //TH1D * hists[16];
  TH1D * hists_1D[8];

  TCanvas * canvas[8];


  //std::vector<double> values;
  //std::vector<double> errors;

  double values[8][5];
  double errors[8][5];

  for(unsigned int i=0; i<8; i++) {

    hists_1D[i] = new TH1D(histname[i]+"_1D",";;", 200, -1.005, 0.995);

    //for(unsigned int f=0; f<2; f++) {
    values[i][0] = ((TH1D*)file[0]->Get(histname[i]))->GetBinContent(1);
    values[i][1] = ((TH1D*)file[0]->Get(histname[i]))->GetBinContent(2);
    values[i][2] = ((TH1D*)file[1]->Get(histname[i]))->GetBinContent(1);
    values[i][3] = ((TH1D*)file[1]->Get(histname[i]))->GetBinContent(2);
    values[i][4] = ((TH1D*)file[1]->Get(histname[i]))->GetBinContent(3);

    errors[i][0] = ((TH1D*)file[0]->Get(histname[i]))->GetBinError(1);
    errors[i][1] = ((TH1D*)file[0]->Get(histname[i]))->GetBinError(2);
    errors[i][2] = ((TH1D*)file[1]->Get(histname[i]))->GetBinError(1);
    errors[i][3] = ((TH1D*)file[1]->Get(histname[i]))->GetBinError(2);
    errors[i][4] = ((TH1D*)file[1]->Get(histname[i]))->GetBinError(3);

    for(unsigned int j=0;j<5;j++) {
      if(errors[i][j] > 0) {
	//hists_1D[i]->Fill(values[i][j], 1.0 / (errors[i][j] * errors[i][j])); //with weighting
	hists_1D[i]->Fill(values[i][j]); //without weighting
      }
    }


      //unsigned int index = (f*8)+i;
      //hists[index] = (TH1D*)file->Get(histname[i])->Clone();
      //if(f==0) { //only fill once
      //hists_1D[i] = new TH1D(histname[i]+"_1D",";;", 200, -1.005, 0.995);
      canvas[i] = new TCanvas(histname[i],histname[i],600,600);
      //}
    //}
  }

  double mean[8];
  double sigma[8];
  double numerator[8];
  double denominator[8];

  for(unsigned int i=0; i<8; i++) {

    for(unsigned int j=0; j<5; j++) {
      if(errors[i][j] >0) {
	numerator[i] += values[i][j] / (errors[i][j] * errors[i][j]);
	denominator[i] += 1.0 / (errors[i][j] * errors[i][j]); 
      }
    }
    
    if(denominator[i] >0) {
      mean[i] = numerator[i] / denominator[i];
      sigma[i] = 1.0 / TMath::Sqrt(denominator[i]);
    }
     canvas[i]->cd();
     hists_1D[i]->Draw();
  }

  
//   std::cout << values[2][0] << " +/- " << errors[2][0] << std::endl;
//   std::cout << values[2][1] << " +/- " << errors[2][1] << std::endl;
//   std::cout << values[2][2] << " +/- " << errors[2][2] << std::endl;
//   std::cout << values[2][3] << " +/- " << errors[2][3] << std::endl;
//   std::cout << values[2][4] << " +/- " << errors[2][4] << std::endl;

//   for(unsigned int i=0; i<8; i++) {
//     std::cout << "bin = " << histname[i] << std::endl;
//     for(unsigned int j=1; j<=hists[i]->GetXaxis()->GetNbins(); j++) {
//       if(hists[i]->GetBinContent(j) !=0 && hists[i]->GetBinError(j) !=0) {
// 	hists_1D[i]->Fill(hists[i]->GetBinContent(j), 1.0 / (hists[i]->GetBinError(j) * hists[i]->GetBinError(j)));
// 	std::cout << "entry " << j << ": " << hists[i]->GetBinContent(j) << " +/- " << hists[i]->GetBinError(j) << std::endl;
//       }
//     }

//     canvas[i]->cd();
//     hists_1D[i]->Draw();
//   }

  TH1D * summary = new TH1D("summary",";;",8, -0.5, 7.5);
  TH1D * summary2 = new TH1D("summary2",";;",8, -0.5, 7.5);
  for(unsigned int i=0; i<8; i++) {
    summary->SetBinContent(summary->GetXaxis()->FindBin(i), mean[i]);
    summary->SetBinError(summary->GetXaxis()->FindBin(i), sigma[i]);
    summary2->SetBinContent(summary2->GetXaxis()->FindBin(i), hists_1D[i]->GetMean());
    summary2->SetBinError(summary2->GetXaxis()->FindBin(i), hists_1D[i]->GetRMS());
    std::cout << histname[i] << " : " << hists_1D[i]->GetMean() << " +/- " << hists_1D[i]->GetRMS() << std::endl;
  }
  
  TCanvas * summaryplot = new TCanvas("summaryplot", "summaryplot", 600,600);
  
  summaryplot->cd();
  summary->Draw();
  
  TCanvas * summaryplot2 = new TCanvas("summaryplot2", "summaryplot2", 600,600);
  
  summaryplot2->cd();
  summary2->SetStats(0);
  summary2->GetYaxis()->SetRangeUser(-1.5,1.5);
  summary2->SetLineWidth(3);
  summary2->SetLineColor(1);
  summary2->Draw();

  TBox * box = new TBox(-0.5, -0.1, 7.5, 0.1);
  box->SetFillColor(3);
  box->SetFillStyle(3001);
  box->Draw("same");
  TBox * box2a = new TBox(-0.5, 0.1, 7.5, 0.2);
  box2a->SetFillColor(4);
  box2a->SetFillStyle(3001);
  box2a->Draw("same");
  TBox * box2b = new TBox(-0.5, -0.2, 7.5, -0.1);
  box2b->SetFillColor(4);
  box2b->SetFillStyle(3001);
  box2b->Draw("same");
  TBox * box3a = new TBox(-0.5, 0.2, 7.5, 0.4);
  box3a->SetFillColor(2);
  box3a->SetFillStyle(3002);
  box3a->Draw("same");
  TBox * box3b = new TBox(-0.5, -0.4, 7.5, -0.2);
  box3b->SetFillColor(2);
  box3b->SetFillStyle(3002);
  box3b->Draw("same");  

  //for(unsigned int i=0; i<low_hist->GetXaxis()->GetNbins(); i++) {
  //  low_hist_1D->Fill(low_hist->GetBinContent(i));  
  //}
  
  //for(unsigned int i=0; i<mid_hist->GetXaxis()->GetNbins(); i++) {
  //  mid_hist_1D->Fill(mid_hist->GetBinContent(i));  
  //}

  //for(unsigned int i=0; i<high_hist->GetXaxis()->GetNbins(); i++) {
  //  high_hist_1D->Fill(high_hist->GetBinContent(i));  
  //}
  

  //TCanvas * c1 = new TCanvas("low_hist","low_hist",600,600);
  //c1->cd();
  //low_hist_1D->Draw();
  
  //TCanvas * c2 = new TCanvas("mid_hist","mid_hist",600,600);
  //c2->cd();
  //mid_hist_1D->Draw();
  
  //TCanvas * c3 = new TCanvas("high_hist","high_hist",600,600);
  //c3->cd();
  //high_hist_1D->Draw();
  
  return;
  
}
