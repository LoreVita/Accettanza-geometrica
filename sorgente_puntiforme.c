void sorgente_puntiforme()
{
	float r_riv = 1.0; //raggio rivelatore (cm)
	float h = 5.0; //distanza sorgente-rivelatore (cm)
	
	//Creazione di un istogramma di controllo della particella che incide sul rivelatore
	TH2D *hriv = new TH2D("hriv","",100,-r_riv,r_riv,100,-r_riv,r_riv);
	hriv->SetStats(0);

	//Creazione di un istogramma di controllo della direzione di emissione della particella
	TH3D *hsfera = new TH3D("hsfera","",100,-1,1,100,-1,1,100,-1,1);
	hsfera->SetStats(0);
	
	//Funzione peso per estrarre isotropicamente il valore di theta in tutto lo spazio
	TF1 *f = new TF1("f","TMath::Sin(x)",0, TMath::Pi());
	
	//Inizializzazione variabili che regolano la struttura del ciclo 
	
	int R = 1;
	float d, theta, phi, x, y, z;
	int n = 1E7; //numero di estrazioni <--> numero di particelle emesse
	int Ni = 0; //numero di particelle incidenti sul rivelatore
	
	for(int i=0; i<n; i++){
	
		//Estrazione isotropa di una direzione (theta, phi) nello spazio a partire dal punto di emissione della particella
		theta = f->GetRandom(0,TMath::Pi());
		phi = gRandom->Uniform(0,2*TMath::Pi());
		
		//Transformazione coordiante cartesiane in coordinate polari (utile per realizzare il relativo istogramma di controllo)
		x = R * TMath::Cos(phi)*TMath::Sin(theta); 
		y = R * TMath::Sin(phi)*TMath::Sin(theta);    
		z = R * TMath::Cos(theta);
		hsfera->Fill(x,y,z);
		
		//Si effettua l'intersezione con l'asse z=d per ottenere la posizione (x,y) della particella sul rivelatore
		d = h/TMath::Cos(theta);
		
		//Si prendono solo le emissioni nel semispazio contenente il rivelatore
		if (d>0.0){
			x = d * TMath::Cos(phi)*TMath::Sin(theta); 
			y = d * TMath::Sin(phi)*TMath::Sin(theta);   
		
			//Verifica dell'incidenza della particella sul rivelatore
			if ((TMath::Power(x,2)+TMath::Power(y,2))<TMath::Power(r_riv,2)){
				Ni++;
				hriv->Fill(x,y);
			}
		}
	}
	
	float omega1, omega2, acc,  err_acc;
	
	//Calcolo dell'angolo solido usando la formula approssimata
	omega1 = (TMath::Power(r_riv,2)*TMath::Pi())/TMath::Power(h,2);
	//cout << "L'angolo solido nel caso approssimato è pari a " << omega1 << endl;
	
	//Calcolo dell'angolo solido usando la formula analitica
	omega2 = 2*TMath::Pi()*(1 - h/(TMath::Sqrt(TMath::Power(h,2)+TMath::Power(r_riv,2))));
	//cout << "L'angolo solido nel caso analitico è pari a " << omega2 << endl;
	
	//Accettanza geometrica valutata mediante metodo Montecarlo e relativo errore 
	acc = (float) Ni/n;
	err_acc = (float)(TMath::Sqrt(Ni)/n);
	
	TCanvas *c1 = new TCanvas("c1","c1",500,500);   
	hsfera->Draw();
	
	TCanvas *c2 = new TCanvas("c2","c2",500,500);   
	hriv->Draw("colz");
	
	cout << "Il numero di particelle emesse è: " << n << endl;
	cout << "Il numero di particelle incidenti è: " << Ni << endl;
	cout << "Accettanza geometrica dalla simulazione: " << acc << endl;
	cout << "Errore sull'accettanza geometrica:   +/-  " << err_acc << endl;
	cout << "Accettanza geometrica nell'ipotesi d>>a:  " << omega1/(4*TMath::Pi()) << endl;
	cout << "Accettanza geometrica metodo analitico:   " << omega2/(4*TMath::Pi()) << endl;
	
}