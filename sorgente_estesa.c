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

void sorgente_estesa()
{
	float r_sorg= 0.5; //raggio sorgente (cm)
	float r_riv = 1; //raggio rivelatore (cm)
	float h = 5; //distanza sorgente-rivelatore (cm) (si assume che i centri della sorgnete e del rivelatore stiano sull'asse perpendicolare le due superfici)
	
	//Creazione di un istogramma di controllo della posizione di origine della particella
    TH2D *hsorg = new TH2D("hsorg","",100,-r_sorg,r_sorg,100,-r_sorg,r_sorg);
	hsorg->SetStats(0);
	
	//Creazione di un istogramma di controllo della direzione di emissione della particella
	TH3D *hsfera = new TH3D("hsfera","",100,-1,1,100,-1,1,100,-1,1);
	hsfera->SetStats(0);
	
	//Creazione di un istogramma di controllo della particella che incide sul rivelatore
	TH2D *hriv = new TH2D("hriv","",100,-r_riv,r_riv,100,-r_riv,r_riv);
	hriv->SetStats(0);
   
    //Funzioni peso
	TF1 *f_s = new TF1("f_s","x",0,r_sorg);  //per avere una distribuzione uniforme del raggio nel caso della sorgente estesa
    TF1 *f = new TF1("f","TMath::Sin(x)",0,TMath::Pi()); //per avere una distribuzione uniforme del raggio per la direzione della particella
	
	//Inizializzazione variabili che regolano la struttura del ciclo 
	
	int n = 1E7; //numero di estrazioni <--> numero di particelle emesse
	int Ni = 0; //numero di particelle incidenti sul rivelatore
    double x_s, y_s; //coordinate cartesiane del punto sulla sorgente
    double ro, alfa; //coordinate polari del punto da estrarre sulla sorgente
	int R = 1;
    double d,theta, phi; //angolo zenitale e azimutale
    double x, y, z; //coordinate punto intersezione retta-sfera R=1
	
    for(int i=0; i<n; i++){
		
		//Generazione della posizione di emissione della particella da un generico punto sulla sorgente 
        ro = f_s->GetRandom(0,r_sorg);
        alfa = gRandom->Uniform(0,2*TMath::Pi());
        x_s = ro * TMath::Cos(alfa);
        y_s = ro * TMath::Sin(alfa);
        hsorg->Fill(x_s,y_s);
		
		//Estrazione isotropa di una direzione (theta, phi) nello spazio (utile per realizzare il relativo istogramma di controllo)
        theta = f->GetRandom(0,TMath::Pi());
        phi = gRandom->Uniform(0,2*TMath::Pi());
        
		//Transformazione coordiante cartesiane in coordinate polari
        x = R * TMath::Cos(phi)*TMath::Sin(theta); 
        y = R * TMath::Sin(phi)*TMath::Sin(theta);    
        z = R * TMath::Cos(theta);
        hsfera->Fill(x,y,z);
		
		//Si effettua l'intersezione con l'asse z=d per ottenere la posizione (x,y) della particella sul rivelatore
		d = h/TMath::Cos(theta);
		
		//Si prendono solo le emissioni nel semispazio contenente il rivelatore
		if (d>0.0){
	        x = d*TMath::Cos(phi)*TMath::Sin(theta)+ x_s; 
	        y = d*TMath::Sin(phi)*TMath::Sin(theta)+ y_s;
			
			//Verifica dell'incidenza della particella sul rivelatore
			if ((TMath::Power(x,2)+TMath::Power(y,2))<TMath::Power(r_riv,2)){
				Ni++;
				hriv->Fill(x,y);
			}
		}
	}
	
    TCanvas *c1 = new TCanvas("c1","c1",500,500); 
	hsfera->Draw();

    TCanvas *c2 = new TCanvas("c2","c2",500,500);   
	hriv->Draw("colz");
	
    TCanvas *c3 = new TCanvas("c3","c3",500,500);   
	hsorg->Draw("colz");
	
	cout << "Il numero di particelle emesse è: " << n << endl;
	cout << "Il numero di particelle incidenti è: " << Ni << endl;
	
	//Accettanza geometrica valutata mediante metodo Montecarlo e relativo errore 
	float acc,  err_acc;
	acc = (float) Ni/n;
	err_acc = (float)(TMath::Sqrt(Ni)/n);
	
	cout << "L'accettanza geometrica dalla simulazione risulta essere: " << acc << endl;
	cout << "L'errore sull'accettanza geometrica risulta essere:   +/- " << err_acc << endl;
	
	//Calcolo del valor teorico dell'angolo solido mediante l'opportuna formula approssimata
	float omega = omega_bessel(r_sorg,r_riv,h);
	cout << "L'accettanza geometrica approssimata (Bessel) è pari a:   " << omega/(4*TMath::Pi())<< endl;
	cout << "L'angolo solido teorico è:  " << omega << endl;
	cout << "L'angolo solido simulato è: " << acc*4*TMath::Pi() << endl;
	
}