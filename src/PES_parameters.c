
double LJsigma(int Zi);
double LJsigma(int Zi)
{
	if(Zi==18){return 3.4;}
	if(Zi==10){return 2.782;}
	if(Zi==6){return 3.4;}
	if(Zi==1){return 2.65;}
	if(Zi==14){return 3.7;}
	if(Zi==32){return 3.45;}
	else{return 0.0;}
}

double LJepsilon(int Zi);
double LJepsilon(int Zi)
{
	if(Zi==18){return 0.0102985;}
	if(Zi==10){return 0.003084;}
	if(Zi==6){return 0.00284;}
	if(Zi==1){return 0.00150;}
	if(Zi==14){return 0.01182;}
	if(Zi==32){return 0.0015;}
	else{return 0.0;}
}

double A_coeff(int Zi);
double A_coeff(int Zi)
{
	if(Zi==6){return 1393.6;}
	if(Zi==14){return 1830.8;}
	if(Zi==1){return 86.712;}
	if(Zi==32){return 1769.0;}
	else{return 0.0;}
}

double B_coeff(int Zi);
double B_coeff(int Zi)
{
	if(Zi==6){return 346.74;}
	if(Zi==14){return 471.18;}
	if(Zi==1){return 43.531;}
	if(Zi==32){return 419.23;}
	else{return 0.0;}
}

double lamda_coeff(int Zi);
double lamda_coeff(int Zi)
{
	if(Zi==6){return 3.4879;}
	if(Zi==14){return 2.4799;}
	if(Zi==1){return 3.7879;}
	if(Zi==32){return 2.4451;}
	else{return 0.0;}
}

double mi_coeff(int Zi);
double mi_coeff(int Zi)
{
	if(Zi==6){return 2.2119;}
	if(Zi==14){return 1.7322;}
	if(Zi==1){return 1.98;}
	if(Zi==32){return 1.7047;}
	else{return 0.0;}
} 

double beta_coeff(int Zi);
double beta_coeff(int Zi)
{
	if(Zi==6){return 1.5724e-7;}
	if(Zi==14){return 1.1e-6;}
	if(Zi==1){return 4.0;}
	if(Zi==32){return 9.0166e-7;}
	else{return 0.0;}
}
 
double eta_coeff(int Zi);
double eta_coeff(int Zi)
{
	if(Zi==6){return 0.72751;}
	if(Zi==14){return 0.78734;}
	if(Zi==1){return 1.0;}
	if(Zi==32){return 0.75627;}
	else{return 0.0;}
}

double c_coeff(int Zi);
double c_coeff(int Zi)
{
	if(Zi==6){return 38049.0;}
	if(Zi==14){return 100390.0;}
	if(Zi==1){return 0.0;}
	if(Zi==32){return 106430.0;}
	else{return 0.0;}
}

double d_coeff(int Zi);
double d_coeff(int Zi)
{
	if(Zi==6){return 4.3484;}
	if(Zi==14){return 16.217;}
	if(Zi==1){return 1.0;}
	if(Zi==32){return 15.652;}
	else{return 0.0;}
}

double h_coeff(int Zi);
double h_coeff(int Zi)
{
	if(Zi==6){return -0.57058;}
	if(Zi==14){return -0.59825;}
	if(Zi==1){return 1.0;}
	if(Zi==32){return -0.43884;}
	else{return 0.0;}
}

double R_coeff(int Zi);
double R_coeff(int Zi)
{
	if(Zi==6){return 1.8;}
	if(Zi==14){return 2.7;}
	if(Zi==1){return 0.8;}
	if(Zi==32){return 2.8;}
	else{return 0.0;}
}

double S_coeff(int Zi);
double S_coeff(int Zi)
{
	if(Zi==6){return 2.1;}
	if(Zi==14){return 3.0;}
	if(Zi==1){return 1.0;}
	if(Zi==32){return 3.1;}
	else{return 0.0;}
}

double tersoff_xi(int Zi,int Zj);
double tersoff_xi(int Zi,int Zj)
{
	double xi_Si_C=0.9776;
	double xi_Si_Ge=1.00061;
	double xi_Si_H=0.78;
	double xi_Ge_H=0.76;
	//Si-C interaction
	if (Zi==14 && Zj==6){return xi_Si_C;}
	//C-Si interaction
	else if (Zi==6 && Zj==14){return xi_Si_C;}
	//Si-Ge interaction
	else if (Zi==14 && Zj==32){return xi_Si_Ge;}
	//Ge-Si interaction
	else if (Zi==32 && Zj==14){return xi_Si_Ge;}
	//Si-H interaction
	else if (Zi==14 && Zj==1){return xi_Si_H;}
	//H-Si interaction
	else if (Zi==1 && Zj==14){return xi_Si_H;}
	//Ge-H interaction
	else if (Zi==32 && Zj==1){return xi_Ge_H;}
	//H-Ge interaction
	else if (Zi==1 && Zj==32){return xi_Ge_H;}
	else {return 1.0;}
}
