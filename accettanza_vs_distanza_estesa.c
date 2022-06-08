//Funzione per il calcolo dell'angolo solido mediante formula approssimata di Bessel
float omega_bessel(float s, float a, float d)
{
	//(devono essere nelle stesse unità di misura)
	//s=raggio sorgente 
	//a=raggio rivelatore 
	//d=distanza sorgente rivelatore
	
	float A = TMath::Power(s/d,2);
	float B = TMath::Power(a/d,2);
	float F1_1 = (5.0/16.0)*(B/TMath::Power(1+B,(7.0/2.0)));
	float F1_2 = (35.0/64.0)*(TMath::Power(B,2)/TMath::Power(1+B,(9.0/2.0)));
	float F2_1 = (35.0/128.0)*(B/TMath::Power(1+B,(9.0/2.0)));
	float F2_2 = (315.0/256.0)*(TMath::Power(B,2)/TMath::Power(1+B,(11.0/2.0)));
	float F2_3 = (1155.0/1024.0)*(TMath::Power(B,3)/TMath::Power(1+B,(13.0/2.0)));
	float F1 = F1_1 - F1_2;
	float F2 = F2_1 - F2_2 + F2_3;
	
	float omega = 2*TMath::Pi()*(1 - 1/TMath::Power(1+B,(1.0/2.0)) - (3.0/8.0)*(A*B)/TMath::Power(1+B,(5.0/2.0)) + TMath::Power(A,2)*F1 - TMath::Power(A,3)*F2);
	return omega;
}

//Funzione che resituisce il numero di particelle incidenti ai fini del calcolo
//dell'accettanza geometrica mediante il metodo Montecarlo per una sorgente estesa
float accettanza_sorgente_estesa(float r_sorg, float r_riv, float h, int N)
{
	//Funzioni peso
	TF1 *f_s = new TF1("f_s","x",0,r_sorg); //per avere una distribuzione uniforme del raggio nel caso della sorgente estesa
    TF1 *f = new TF1("f","TMath::Sin(x)",0,TMath::Pi()); //per avere una distribuzione uniforme del raggio per la direzione della particella
	
	int Ni = 0; //numero di particelle incidenti sul rivelatore
    double x_s, y_s; //coordinate cartesiane del punto sulla sorgente
    double ro, alfa; //coordinate polari del punto da estrarre sulla sorgente
	int R = 1;
    double d,theta, phi; //angolo zenitale e azimutale
    double x, y, z; //coordinate punto intersezione retta-sfera R=1
	
    for(int i=0; i<N; i++){
		
		//Generazione della posizione di emissione della particella da un generico punto sulla sorgente 
        ro = f_s->GetRandom(0,r_sorg);
        alfa = gRandom->Uniform(0,2*TMath::Pi());
        x_s = ro * TMath::Cos(alfa);
        y_s = ro * TMath::Sin(alfa);
		
		//Estrazione isotropa di una direzione (theta, phi) nello spazio (utile per realizzare il relativo istogramma di controllo)
        theta = f->GetRandom(0,TMath::Pi());
        phi = gRandom->Uniform(0,2*TMath::Pi());
		
		//Si effettua l'intersezione con l'asse z=d per ottenere la posizione (x,y) della particella sul rivelatore
		d = h/TMath::Cos(theta);
		
		//Si prendono solo le emissioni nel semispazio contenente il rivelatore
		if (d>0.0){
	        x = d*TMath::Cos(phi)*TMath::Sin(theta)+ x_s; 
	        y = d*TMath::Sin(phi)*TMath::Sin(theta)+ y_s;
			
			//Verifica dell'incidenza della particella sul rivelatore
			if ((TMath::Power(x,2)+TMath::Power(y,2))<TMath::Power(r_riv,2)) Ni++;
		}
	}
	
	return Ni;
}

void accettanza_vs_distanza_estesa()
{
	//Creazione dei grafici
	TGraphErrors *gr_MC = new TGraphErrors();					
	gr_MC->SetName("Metodo Montecarlo");
	gr_MC->SetLineColor(kRed);
	gr_MC->SetMarkerColor(kRed);
	//gr_MC->SetMarkerStyle(20);
	gr_MC->SetMarkerStyle(1);
	
	TGraph *gr_bes = new TGraph();
	gr_bes->SetName("Metodo numerico di Bessel");
	gr_bes->SetLineColor(kBlue);
	gr_bes->SetMarkerColor(kBlue);
	//gr_bes->SetMarkerStyle(21);
	gr_bes->SetMarkerStyle(1);
	
	//Inizializzazione variabili che regolano la struttura del ciclo 
	float r_sorg = 1;  //raggio sorgente (cm)
	float r_riv = 0.5; //raggio rivelatore (cm)
	float h = 0.0; //distanza iniziale sorgente-rivelatore (cm) (si assume che i centri della sorgnete e del rivelatore stiano sull'asse perpendicolare le due superfici)
	float step = 0.05;
	int n = 1e2; //numero di interazioni del ciclo <-> punti del grafico
	int N = 1e5; //numero di eventi generati per il calcolo dell'accettanza (per un singolo valore)
	float Ni; //numero di particelle incidenti sul rivelatore
	float omega; //angolo solido valutato mediante formula approssimata di Bessel
	float acc_MC, err_acc_MC, acc_bes;
	
	for(int i=0; i<n; i++){
	
		//Calcolo dell'accetanza mediante tecnica Montecarlo
		Ni = accettanza_sorgente_estesa(r_sorg,r_riv,h,N);
		acc_MC = Ni/N;
		gr_MC->SetPoint(i,h,acc_MC);
		//Calcolo errore sull'accetanza simulata
		err_acc_MC = TMath::Sqrt(Ni)/N;
		gr_MC->SetPointError(i,0,err_acc_MC);
	
		//Calcolo dell'accetanza mediante la formula approssimata di Bessel
		omega = omega_bessel(r_sorg,r_riv,h);
		if(h<0.2){
			acc_bes=0.0;
		}
		else{
			acc_bes = omega/(4*TMath::Pi());
		}
		gr_bes->SetPoint(i,h,acc_bes);
	
		h = h+step; //incremento della distanza
	}

	TCanvas *c1 = new TCanvas("c1","c1",800,600);

	//Creazione Grafico su cui rappresentare i tre set di dati
	TMultiGraph *mg = new TMultiGraph();

	//Opzioni grafiche
	mg->SetTitle("""; Distanza sorgente-rivelatore (cm); Accettanza geometrica"); 
	mg->GetXaxis()->CenterTitle(true);	
	mg->GetYaxis()->CenterTitle(true);
	mg->SetMinimum(0.0);
	mg->SetMaximum(0.2);
	mg->GetXaxis()->SetLimits(0.0,2);

    //c1->BuildLegend(0.1,0.7,0.48,0.9);
	c1->BuildLegend(0.9,0.7,0.48,0.9);

	//Grafici da rappresentare
	mg->Add(gr_bes, "l");
	//mg->Add(gr_bes);
	mg->Add(gr_MC); 
	mg->Draw("ap"); // "a" per disegnare gli assi; "p" per disegnare solo i marcatori
	
	//Legenda
	//c1->BuildLegend(0.1,0.7,0.48,0.9);
	c1->BuildLegend(0.9,0.75,0.48,0.9);

}