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

void confronto_errori()
{
	//Creazione dei grafici
	TGraphErrors *gr_punt = new TGraphErrors();		
	gr_punt->SetName("Sorgente puntiforme");
	gr_punt->SetLineColor(kBlue);
	gr_punt->SetMarkerColor(kBlue);
	//gr_punt->SetMarkerStyle(20);
	gr_punt->SetMarkerStyle(1);

	TGraphErrors *gr_estes =new TGraphErrors();		
	gr_estes->SetName("Sorgente estesa");
	gr_estes->SetLineColor(kRed);
	gr_estes->SetMarkerColor(kRed);
	//gr_estes->SetMarkerStyle(21);
	gr_estes->SetMarkerStyle(1);

	//Inizializzazione variabili che regolano la struttura del ciclo
	float r_sorg = 0.5; //raggio sorgente (cm)
	float r_riv = 1; //raggio rivelatore (cm)
	float h = 0.1; //distanza iniziale sorgente-rivelatore (cm) (si assume che i centri della sorgnete e del rivelatore stiano sull'asse perpendicolare le due superfici)
	float step = 0.05;
	int n = 1e2; //numero di interazioni del ciclo <-> punti del grafico
	int N = 1e5; //numero di eventi generati per il calcolo dell'accettanza (per un singolo valore)
	float Ni_punt, Ni_estes; //numero di particelle incidenti sul rivelatore nei due casi
	float acc_punt, acc_estes; //accettanze geometriche
	float err_acc_punt, err_acc_estes; //errori sulle accettanze geometriche
	
	for(int i=0; i<n; i++){
	
		//Calcolo dell'accetanza nel caso di sorgente puntiforme
		Ni_punt = accettanza_sorgente_puntiforme(r_riv,h,N);
		acc_punt = Ni_punt/N;
		gr_punt->SetPoint(i,h,acc_punt);
		//Calcolo errore
		err_acc_punt = TMath::Sqrt(Ni_punt)/N;
		gr_punt->SetPointError(i,0,err_acc_punt);
	
		//Calcolo dell'accetanza nel caso di sorgente estesa
		Ni_estes = accettanza_sorgente_estesa(r_sorg,r_riv,h,N);
		acc_estes = Ni_estes/N;
		gr_estes->SetPoint(i,h,acc_estes);
		//Calcolo errore
		err_acc_estes = TMath::Sqrt(Ni_estes)/N;
		gr_estes->SetPointError(i,0,err_acc_estes);
	
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
	mg->SetMaximum(0.5);
	mg->GetXaxis()->SetLimits(0.0,5.0);

	//Grafici da rappresentare
	mg->Add(gr_punt); 
	mg->Add(gr_estes);
	mg->Draw("alp"); // "a" per disegnare gli assi; "p" per disegnare solo i marcatori
	
	//Legenda
	//c1->BuildLegend(0.1,0.7,0.48,0.9);
	c1->BuildLegend(0.9,0.5,0.48,0.9);
}