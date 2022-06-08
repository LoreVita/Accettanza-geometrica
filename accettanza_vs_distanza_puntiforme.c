//Funzione che resituisce il numero di particelle incidenti ai fini del calcolo
//dell'accettanza geometrica mediante il metodo Montecarlo per una sorgente puntiforme
float accettanza_sorgente_puntiforme(float r_riv, float h, int N)
{
    //Funzione peso per estrarre isotropicamente il valore di theta in un semispazio
    TF1 *f = new TF1("f","TMath::Sin(x)",0, TMath::Pi());
	
	int R = 1;
	float d, theta, phi, x, y, z;
	int Ni = 0; //numero di particelle incidenti sul rivelatore
	for(int i=0; i<N; i++){
	
		//Estrazione isotropa di una direzione (theta, phi) nello spazio a partire dal punto di emissione della particella
		theta = f->GetRandom(0,TMath::Pi());
		phi = gRandom->Uniform(0,2*TMath::Pi());
		
		//Si effettua l'intersezione con l'asse z=d per ottenere la posizione (x,y) della particella sul rivelatore
		d = h/TMath::Cos(theta);
		
		//Si prendono solo le emissioni nel semispazio contenente il rivelatore
		if (d>0.0){
			x = d * TMath::Cos(phi)*TMath::Sin(theta); 
			y = d * TMath::Sin(phi)*TMath::Sin(theta);   
		
			//Verifica dell'incidenza della particella sul rivelatore
			if ((TMath::Power(x,2)+TMath::Power(y,2))<TMath::Power(r_riv,2)) Ni++;
		}
	}
	
	return Ni; 
}

void accettanza_vs_distanza_puntiforme()
{
	//Creazione dei grafici
	TGraphErrors *gr_MC = new TGraphErrors();					
	gr_MC->SetName("Metodo Montecarlo");
	gr_MC->SetLineColor(kRed);
	gr_MC->SetMarkerColor(kRed);
	//gr_MC->SetMarkerStyle(20);
	gr_MC->SetMarkerStyle(1);
	
	TGraph *gr_ap = new TGraph();
	gr_ap->SetName("Metodo approssimato");
	gr_ap->SetLineColor(kGreen+2);
	gr_ap->SetMarkerColor(kGreen+2);
	//gr_ap->SetMarkerStyle(21);
	gr_ap->SetMarkerStyle(1);
	
	TGraph *gr_teo = new TGraph();
	gr_teo->SetName("Metodo analitico");
	gr_teo->SetLineColor(kBlue);
	gr_teo->SetMarkerColor(kBlue);
	//gr_teo->SetMarkerStyle(22);
	gr_teo->SetMarkerStyle(1);
	
	//Inizializzazione variabili che regolano la struttura del ciclo 
	float r_riv = 1; //raggio rivelatore
	float h = 0.0; //distanza iniziale sorgente-rivelatore
	float step = 0.05;
	int n = 1e2; //numero di interazioni del ciclo <-> punti del grafico
	int N = 1e5; //numero di eventi generati per il calcolo dell'accettanza (per un singolo valore)
	float Ni; //numero di particelle incidenti sul rivelatore
	float omega_ap, omega_teo; //angoli solidi
	float acc_MC, acc_ap, acc_teo; //accettanze geometriche
	float err_acc_MC; //errore sull'accettanze geometrica simulata
	
	for(int i=0; i<n; i++){
		
		//Calcolo dell'accetanza mediante tecnica Montecarlo
		Ni = accettanza_sorgente_puntiforme(r_riv,h,N);
		acc_MC = Ni/N;
		gr_MC->SetPoint(i,h,acc_MC);
		//Calcolo errore sull'accetanza simulata
		err_acc_MC = TMath::Sqrt(Ni)/N;
		gr_MC->SetPointError(i,0,err_acc_MC);
		
		//Calcolo dell'accetanza mediante la formula approssimata
		omega_ap = (TMath::Power(r_riv,2)*TMath::Pi())/TMath::Power(h,2);
		acc_ap = omega_ap/(4*TMath::Pi());
		gr_ap->SetPoint(i,h,acc_ap);
		
		//Calcolo dell'accetanza mediante l'espressione analitica
		omega_teo = 2*TMath::Pi()*(1 - h/(TMath::Sqrt(TMath::Power(h,2)+TMath::Power(r_riv,2))));
		acc_teo = omega_teo/(4*TMath::Pi());
		gr_teo->SetPoint(i,h,acc_teo);
		
		h = h+step; //incremento della distanza
		
	}
	
	TCanvas *c1 = new TCanvas("c1","",800,600);
	
	//Creazione Grafico su cui rappresentare i tre set di dati
	TMultiGraph *mg = new TMultiGraph();
	
	//Opzioni grafiche
	mg->SetTitle("""; Distanza sorgente-rivelatore (cm); Accettanza geometrica"); 
	mg->GetXaxis()->CenterTitle(true);	
	mg->GetYaxis()->CenterTitle(true);
	mg->SetMinimum(0.0);
	mg->SetMaximum(0.6);
	mg->GetXaxis()->SetLimits(0.0,5.0);
	
	//Grafici da rappresentare
	mg->Add(gr_ap, "l");
	//mg->Add(gr_ap);
	mg->Add(gr_teo, "l");
	//mg->Add(gr_teo );
	mg->Add(gr_MC); 
	mg->Draw("ap"); // "a" per disegnare gli assi; "p" per disegnare solo i marcatori
	
	//Legenda
	//c1->BuildLegend(0.1,0.7,0.48,0.9);
	c1->BuildLegend(0.9,0.7,0.48,0.9);
    
}