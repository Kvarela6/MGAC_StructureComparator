#pragma once
#include <set>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <ostream>
#include <algorithm>
#include <iomanip>
#include <map>
#include "SVL.h"

using namespace std;

#define AUtoAng 0.529177249
#define PI 3.1416
#define RydToKCalMol 313.7547345
#define CaltoJoule 4.184
#define atomo1_Angulos "C"  
#define atomo2_Angulos "C"


//elements
typedef enum struct Elem {
	UNK = 0, C, H, N, O, P, Cl, F, S, Br,
}Elem;


Elem getElemType(string in) {
	Elem type;
	//UNK=0,C,H,N,O,P,Cl,F,S,Br,
	if (in == "C")
		type = Elem::C;
	else if (in == "H")
		type = Elem::H;
	else if (in == "N")
		type = Elem::N;
	else if (in == "O")
		type = Elem::O;
	else if (in == "P")
		type = Elem::P;
	else if (in == "Cl")
		type = Elem::Cl;
	else if (in == "F")
		type = Elem::F;
	else if (in == "S")
		type = Elem::S;
	else if (in == "Br")
		type = Elem::Br;
	else
		type = Elem::UNK;

	return type;
}

string getElemName(Elem in) {
	string type;
	//UNK=0,C,H,N,O,P,Cl,F,S,Br,
	if (in == Elem::C)
		type = "C";
	else if (in == Elem::H)
		type = "H";
	else if (in == Elem::N)
		type = "N";
	else if (in == Elem::O)
		type = "O";
	else if (in == Elem::P)
		type = "P";
	else if (in == Elem::Cl)
		type = "Cl";
	else if (in == Elem::F)
		type = "F";
	else if (in == Elem::S)
		type = "S";
	else if (in == Elem::Br)
		type = "Br";
	else
		type = "UNK";

	return type;


}

class ang_dist {
public:



	double R;
	double cos_ang;


	ang_dist() : R(0.0), cos_ang(0.0) {};
	ang_dist(const double R_, const double cos_ang_) {
		R = R_;
		cos_ang = cos_ang_;
	}

	bool operator<(const ang_dist& M) {
		if (cos_ang < M.cos_ang) return true;
		else return false;

	}
	bool operator>(const ang_dist& M) {
		if (cos_ang > M.cos_ang) return true;
		else return false;
	}
	bool operator==(const ang_dist& M) {
		if (cos_ang == M.cos_ang) return true;
		else return false;
	}

	double Rcos() const
	{
		return R * cos_ang;
	}
	double Rsen() const
	{
		return R * sqrt(1 - (cos_ang * cos_ang));
	}
	double dist_Rn(const ang_dist& N) const
	{
		double X_ = this->Rcos() - N.Rcos();
		double Y_ = this->Rsen() - N.Rsen();
		return sqrt((X_ * X_) + (Y_ * Y_));
	}

	// El mapa 
	std::map<Elem, std::vector<ang_dist>> AngulosPorAtomo;


	// Método para mostrar
	void mostrar() const {
		for (const auto& par : AngulosPorAtomo) {
			std::cout << getElemName(par.first) << ": ";
			std::cout << "Size del vector: " << par.second.size();
			///for (double d : par.second) {
			//std::cout << d << " ";
			//}
			//std::cout << std::endl;
		}
	}


};

bool comparacion_angulos_tomando_R(const ang_dist& M, const ang_dist& N)
{
	double tolerancia_R = 1e-3;
	//double tolerancia_angle = 0.01;

	if (abs(M.R - N.R) < tolerancia_R) {
		//if (M.R == N.R) {
		if (M.cos_ang < N.cos_ang) {
			return true;
		}
		else
		{
			return false;
		}
	}
	if (M.R < N.R)
	{
		return true;
	}
	else return false;


	/*	if (M.cos_ang < N.cos_ang) //+ tolerancia_angle)
		{
			if (M.R < N.R)
			//if (-tolerancia_R < M.R - N.R < tolerancia_R)
			{
				//if (M.cos_ang < N.cos_ang+tolerancia_angle) return true;
				return true;

			}
			else return false;
		}
		return false;*/

}


class super_celda {  //Coordenadas en clase vector!
public:
	Mat Celda;
	vector <super_celda> list_coor_frac;
	vector <super_celda> cart_supercelda;
	super_celda() {};    //Constructor por defecto

	super_celda(Mat Coordenadas_) {    //Constructor que tiene parametros
		Celda = Coordenadas_;
	}
};

class Distancias {

public:
	// El mapa 
	std::map<std::string, std::vector<double>> distanciasPorTipo;

	// Constructor vacío
	Distancias() {}

	// Constructor que inicializa con combinaciones si querés
	Distancias(const std::vector<std::string>& combinaciones) {
		for (const auto& tipo : combinaciones) {
			distanciasPorTipo[tipo] = std::vector<double>();  // inicializar vacío
		}
	}

	// Método para agregar distancias
	void agregarDistancia(const std::string& tipo, double valor) {
		distanciasPorTipo[tipo].push_back(valor);
	}

	// Método para mostrar
	void mostrar() const {
		for (const auto& par : distanciasPorTipo) {
			std::cout << par.first << ": ";
			std::cout << "Size del vector: " << par.second.size();
			///for (double d : par.second) {
			//std::cout << d << " ";
			//}
			//std::cout << std::endl;
		}
	}

};


class Parameters { //from output

public:



	string structure_name;
	string grupo_espacial;
	int atomxmol;
	int nat_; //number of atoms per cell
	int sym_; //Simetria
	double energy_ = 0.0; // last energy 
	double energy_kj = 0.0;
	double enthalpy = 0.0;
	double enthalpy_kj = 0.0;
	double alat_;
	double a_, b_, c_; //cell parameters
	double gamma_, alpha_, beta_; //angles
	Vec3 va_, vb_, vc_;
	Vec3 Cart_;
	double diagonal; // diagonal para la supercelda 
	double corte; //
	//vector <string> atomos_cell; //leidos del output, estan en el orden del output

	vector <Elem> atomos_cell;
	vector <string> atoms_combinations;
	vector <Elem> atoms_types;

	vector <string> atomic_type_;  //ordenados para el CIF
	Mat coor_Frac_; // Matriz coordenadas Fraccionarias output
	Mat coor_Cart_; //Matriz coordenadas cartesianas


	vector<ang_dist> ang_C_;
	vector<ang_dist> ang_O_;
	vector<ang_dist> ang_H_;


	std::map<std::string, std::vector<double>> distanciasPorTipo_Param;
	map<Elem, std::vector<ang_dist>> angulosPorAtomos_Param;


	Parameters() : atomxmol(0), nat_(0), sym_(0), alat_(0.0), a_(0.0), b_(0.0), c_(0.0), gamma_(0.0), alpha_(0.0), beta_(0.0), diagonal(0.0), corte(0.0), energy_(0.0) {};
	void cell_parameters(string file_out, super_celda& coor);
	void Armado_Supercelda(super_celda& coor, Parameters& params_);
	void Enlaces(Distancias& dist_, ang_dist& ang_);
	void Distancia_SC(super_celda& coor, Distancias& dist_, const int atomxmol); //Distancias para SUPERCELDA
	//void atomic_labels_(int& Atomxmol_);
	//void get_atoms_for_angles(int& atom_1, int& atom_2, int& Atomxmol_);


	void Angulos(const int atom_1, const int atom_2, const int atomxmol, super_celda& coor, ang_dist ang_);


	friend class Distancias;

};


class estructura {
public:

	vector <Parameters> estructuras;
	//vector <ang_dist> angulos_y_distancias;
	//vector <super_celda> coordenadas;

	void diferencia_Angulos(int i, int j, ofstream& resumen, ofstream& resumen2, double& delta_total_ang, double& Maximo_total_ang);
	ang_dist busqueda_minima_dist(const ang_dist& M, const vector<ang_dist>& angulos_dist, double& distancia_minima_calc, const double corte);
	void diferencia_distancias(int i, int j, ofstream& resumen_angulos, ofstream& resumen2, double& Delta_total, double& Maximo_total);
	void dist_and_angles();

};