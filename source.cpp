
//A crystal structure comparison module developed for the MGAC CSP workflow. 
//to detect equivalent structures and remove duplicates during large - scale evolutionary searches.

//In Construction For details write an email to kristalvarelamuzzati@gmail.com


#include "Header.h"

string get_file_contents(const char* filename);
bool outputVerification(string fileout);
double angle(const Vec3& A, const Vec3& B, const Vec3& C);
Vec3 Frac_2Cart(Parameters& param_, Vec3 Frac);

using namespace std;

int main(int argc, char* argv[])
{
	int Atomxmol_;

	if (argc > 1)
	{
		string sAtomxmol_ = argv[1];
		istringstream ss(sAtomxmol_);
		ss >> Atomxmol_;
		cout << Atomxmol_ << " Numero de atomos por molecula " << endl;
	}
	else
	{
		cout << " falta valor de numero de atomos por molecula como argumento " << endl;
		exit(1);
	}

	//PARA WINDOWS
	string Dir_outs = "D:\\UltimasVersiones_\\Iguales_MOREATOMS\\Bestoutqe\\pru.sys_qe_spcgr*.out";
	string Dir_home = "D:\\UltimasVersiones_\\Iguales_MOREATOMS\\list_out_qe.txt";
	string Dir_best = "D:\\UltimasVersiones_\\Iguales_MOREATOMS\\Bestoutqe\\";
	string preNameOut = "pru.sys_qe_spcgrp";
	stringstream Launcher;
	Launcher << "dir " << Dir_outs << " > " << Dir_home << endl; //Crea lista de los oust del bestoutqe
	//cout << Launcher.str().c_str() << endl;
	system(Launcher.str().c_str());

	string list_out_qe = get_file_contents(Dir_home.c_str()); //Conviernte en string la lista 

	/*
	//PARA LINUX


	system("mkdir done_iguales toCheck");
	//system("mkdir toCheck");

	string preNameOut = "pru.sys_qe_spcgrp";
	string Dir_best = "bestoutqe/";
	system("ls -l bestoutqe/pru.sys_qe_spcgrp*.out > tmp_list_bestout_qe.txt");
	string list_out_qe = get_file_contents("tmp_list_bestout_qe.txt");
	//string dir_best_out = "bestoutqe/";
	*/


	int Atom_1 = 1; //para angulos, deberia ser Carbono 
	int Atom_2 = 13; //para angulos, deberia ser Oxigeno
	bool found = false; //Busca si ya el out fue evaluado previamente. 

	Parameters unit_;
	super_celda coor;
	estructura struc;
	Distancias dist;
	ang_dist ang;

	unit_.atomxmol = Atomxmol_;

	size_t pos = 0;
	while (pos != string::npos)
	{
		pos++;
		pos = list_out_qe.find(preNameOut, pos);
		if (pos == string::npos) break;

		string sfileout = list_out_qe.substr(pos, list_out_qe.find(".out", pos) - pos);
		//cout << sfileout << endl;
		unit_.structure_name = sfileout;

		string _grupo_espacial = sfileout.substr(sfileout.find(preNameOut) + preNameOut.size());
		//cout << prueba_ << endl;
		size_t c_pos = 0;
		_grupo_espacial = _grupo_espacial.substr(0, _grupo_espacial.find("n", c_pos));
		//cout << _grupo_espacial << endl;
		unit_.grupo_espacial = _grupo_espacial;

		sfileout = Dir_best + sfileout + ".out";
		string file_out = get_file_contents(sfileout.c_str());
		//cout << file_out << endl; 

		bool outEstatus = outputVerification(file_out);
		if (outEstatus == false)
		{
			cout << "Out invalido, no termina bien " << unit_.structure_name << endl;
			continue;
		}

		unit_.cell_parameters(file_out, coor);

		unit_.Armado_Supercelda(coor, unit_);

		unit_.Enlaces(dist, ang);

		unit_.Distancia_SC(coor, dist, Atomxmol_);

		unit_.Angulos(Atom_1, Atom_2, Atomxmol_, coor, ang);

		struc.estructuras.push_back(unit_);

	}//END_WHILE. 


	struc.dist_and_angles();
}

void estructura::diferencia_Angulos(int i, int j, ofstream& resumen_angulos, ofstream& resumen2, double& delta_total_ang, double& Maximo_total_ang) //Distancias minimas 
{
	double delta_enlaces_angulos = 0.0;
	double resta_ang = 0.0;
	double resta_cuadrado = 0.0;
	double acumula_ = 0.0;
	int contador = 0;
	int contador2 = 0;
	double acumulacion_total = 0.0;
	double Maximo_resta = 0;
	//double Maximo_total1 = 0.0;

	double angulos_i_max = 0.0;
	double angulos_j_max = 0.0;

	string name_i = estructuras[i].structure_name;
	//cout << name_i << endl;

	string name_j = estructuras[j].structure_name;
	//cout << name_j << endl;
	//cout << struc.estructuras[j].structure_name << endl;

	resumen_angulos << "Estructura i: " << name_i << " Estructura j: " << name_j << endl;

	double corte = 0.0;

	if (estructuras[i].diagonal < estructuras[j].diagonal) corte = estructuras[i].diagonal;
	else corte = estructuras[j].diagonal;
	//double corte2 = corte * 1.3;
	//corte = corte * 1.3;
	resumen_angulos << "corte: " << corte << endl;



	//PARA RECORER MAPAs
	for (const auto& par_i : estructuras[i].angulosPorAtomos_Param) {
		//std::cout << "Clave: " << par_i.first << "\nValores: ";
		const string& clave = getElemName(par_i.first);
		const vector<ang_dist>& vec_i = par_i.second;
		resumen_angulos << "Angulos " << getElemName(par_i.first) << endl;

		auto it_j = estructuras[j].angulosPorAtomos_Param.find(par_i.first);
		if (it_j == estructuras[j].angulosPorAtomos_Param.end()) continue;

		const std::vector<ang_dist>& vec_j = it_j->second;

		for (int z = 0; z < vec_i.size(); z++) {

			if (vec_i[z].R < corte) {
				double distancia_minima = 0.0;
				ang_dist atomo_dist_minima = busqueda_minima_dist(vec_i[z], vec_j, distancia_minima, corte);


				//cout << vec_i[z].cos_ang << " " << atomo_dist_minima.cos_ang << endl; 
				resta_ang = acos(vec_i[z].cos_ang) - acos(atomo_dist_minima.cos_ang);

				resta_cuadrado = resta_ang * resta_ang;
				acumula_ = acumula_ + resta_cuadrado;



				if (Maximo_resta < resta_cuadrado) {
					Maximo_resta = resta_cuadrado;
					if (Maximo_total_ang < Maximo_resta) {

						Maximo_total_ang = Maximo_resta;
						angulos_i_max = vec_i[z].cos_ang;
						angulos_j_max = vec_j[z].cos_ang;

					}
				}

				contador++;
			}

			else continue;

		}
		delta_enlaces_angulos = sqrt(acumula_ / contador);
		resumen_angulos << " RESULTADO delta_Dij PARA " << clave << delta_enlaces_angulos << endl;
		resumen_angulos << " Cantidad de distancias para " << clave << " " << contador << endl;
		resumen_angulos << " Maximo resta cuadrados angulos  " << sqrt(Maximo_resta) << endl;


		acumulacion_total = acumulacion_total + acumula_;
		contador2 = contador2 + contador;
		acumula_ = 0.0;
		contador = 0;
		delta_enlaces_angulos = 0.0;
		Maximo_resta = 0.0;
	}


	delta_total_ang = sqrt((acumulacion_total) / contador2);
	Maximo_total_ang = sqrt(Maximo_total_ang);
	resumen_angulos << " delta_total " << delta_total_ang << endl;
	resumen_angulos << " Cantidad para delta_total para " << " " << contador2 << endl;
	resumen_angulos << " Maximo resta cuadrados angulos  " << Maximo_total_ang << endl;
	resumen_angulos << " angulos correspondientes al maximo:    i: " << angulos_i_max << "    j: " << angulos_j_max << endl << endl;


	resumen2 << delta_total_ang << " " << Maximo_total_ang << endl;



	contador2 = 0;
	acumulacion_total = 0.0;
	angulos_i_max = 0.0;
	angulos_j_max = 0.0;



}





void estructura::diferencia_distancias(int i, int j, ofstream& resumen_dist, ofstream& resumen2, double& Delta_total, double& Maximo_total)
{

	string name_i = estructuras[i].structure_name;
	//cout << name_i << endl;

	string name_j = estructuras[j].structure_name;

	resumen_dist << "Estructura i: " << name_i << " Estructura j: " << name_j << endl;
	resumen_dist << estructuras[i].diagonal << "i " << endl;
	resumen_dist << estructuras[j].diagonal << "j " << endl;
	resumen2 << "Estructura i: " << name_i << " Estructura j: " << name_j << endl;

	double corte = 0.0;
	if (estructuras[i].diagonal < estructuras[j].diagonal) corte = estructuras[i].diagonal;
	else corte = estructuras[j].diagonal;
	//cout << corte << " valor corte tomado " << endl;
	//archivo << corte << " valor corte tomado " << endl;
	resumen_dist << corte << " valor corte tomado " << endl;

	double acumulacion_division = 0.0;
	double acumulacion_divi_total = 0.0;
	double Maximo_divi = 0.0;
	int contador = 0; //cantidad de distancias por tipo de enlace 
	int contador2 = 0;
	double distancia_i_max = 0.0;
	double distancia_j_max = 0.0;
	double delta_Dij = 0.0;


	auto& mapa_i = estructuras[i].distanciasPorTipo_Param;
	auto& mapa_j = estructuras[j].distanciasPorTipo_Param;

	for (const auto& par_i : mapa_i) {
		const string clave = par_i.first;

		auto it_j = mapa_j.find(clave);
		if (it_j != mapa_j.end()) {
			const std::vector<double>& vec_i = par_i.second;
			const std::vector<double>& vec_j = it_j->second;

			for (size_t k = 0; k < vec_i.size() && k < vec_j.size(); ++k) {
				if (vec_i[k] > corte && vec_j[k] > corte) {
					//std::cout << "Break por clave: " << clave << ", índice: " << k << std::endl;
					break;
				}
				else
				{
					double resta = 0.0;
					double suma = 0.0;
					double resta_cuadrado = 0.0;
					double suma_cuadrado = 0.0;
					double division = 0.0;

					resta = vec_i[k] - vec_j[k];
					suma = vec_i[k] + vec_j[k];

					resta_cuadrado = resta * resta;
					suma_cuadrado = suma * suma;
					division = resta_cuadrado / suma_cuadrado;

					acumulacion_division = acumulacion_division + division;
					contador++;

					if (Maximo_divi < division) {
						Maximo_divi = division;
						if (Maximo_total < Maximo_divi) {

							Maximo_total = Maximo_divi;
							distancia_i_max = vec_i[k];
							distancia_j_max = vec_j[k];
						}
					}
				}
			}
		}

		delta_Dij = sqrt(acumulacion_division / contador);
		acumulacion_divi_total = acumulacion_divi_total + acumulacion_division;

		resumen_dist << " RESULTADO delta_Dij PARA " << par_i.first << " " << delta_Dij << endl;
		resumen_dist << " Cantidad de distancias para " << par_i.first << " " << contador << endl;
		resumen_dist << " Maximo Cociente division " << sqrt(Maximo_divi) << endl;

		contador2 = contador2 + contador;
		contador = 0;
		delta_Dij = 0.0;
		acumulacion_division = 0.0;
		Maximo_divi = 0.0;
	}

	Delta_total = sqrt((acumulacion_divi_total) / contador2);
	Maximo_total = sqrt(Maximo_total);
	resumen_dist << " delta_total " << Delta_total << endl;
	resumen_dist << " Maximo total " << Maximo_total << endl;
	resumen_dist << " distancias correspondientes al maximo:    i: " << distancia_i_max << "    j: " << distancia_j_max << endl << endl << endl;
	resumen2 << Delta_total << " " << Maximo_total << " ";



}



void estructura::dist_and_angles()
{
	string done_iguales = "done_iguales";
	string files_toCheck = "toCheck";

	double delta_dist = 0.0;
	double Maximo_total_dist = 0.0;
	double delta_ang = 0.0;
	double Maximo_total_ang = 0.0;

	ofstream resumen_dist, resumen_ang, final_file, iguales, toCheck_;

	resumen_dist.open("valores_distancias.txt", std::ofstream::out | std::ofstream::trunc);
	resumen_ang.open("valores_angulos.txt", std::ofstream::out | std::ofstream::trunc);
	final_file.open("final_file.txt", std::ofstream::out | std::ofstream::trunc);
	iguales.open("iguales.txt", std::ofstream::out | std::ofstream::trunc);
	toCheck_.open("toCheck.txt", std::ofstream::out | std::ofstream::trunc);

	double tolerancia_delta_dist = 1e-2;
	double tolerancia_Maximo_dist = 1e-1;
	double tolerancia_delta_ang = 1e-1; //5e-2
	double tolerancia_maximo_ang = 1e-1;

	double excepcion_max_ang = 0.20;

	for (int i = 0; i < estructuras.size(); i++) {

		for (int j = i + 1; j < estructuras.size(); j++) {

			delta_dist = 0.0;
			Maximo_total_dist = 0.0;
			delta_ang = 0.0;
			Maximo_total_ang = 0.0;

			string name_i = estructuras[i].structure_name;
			string name_j = estructuras[j].structure_name;

			stringstream launcher;

			final_file << estructuras[i].structure_name << "/" << estructuras[j].structure_name << " ";

			diferencia_distancias(i, j, resumen_dist, final_file, delta_dist, Maximo_total_dist);

			//cout << delta_dist << " delta distancias" << endl;
			//cout << Maximo_total_dist << " Maximo total distancias" << endl;
			if ((delta_dist < tolerancia_delta_dist) && (Maximo_total_dist < tolerancia_Maximo_dist))
			{
				diferencia_Angulos(i, j, resumen_ang, final_file, delta_ang, Maximo_total_ang);

				if ((delta_ang < tolerancia_delta_ang) && (Maximo_total_ang < tolerancia_maximo_ang))
				{
					iguales << " i: " << estructuras[i].structure_name << " j: " << estructuras[j].structure_name << endl;
					iguales << " delta_dist " << delta_dist << "; " << "Maximo_dist " << Maximo_total_dist << "; " << "delta_ang " << delta_ang << "; " << "Maximo_ang " << Maximo_total_ang;// << endl << endl;

					//cout << " i: " << estructuras[i].structure_name << " j: " << estructuras[j].structure_name << " son iguales " << endl;

					if (estructuras[j].sym_ < estructuras[i].sym_)
					{
						iguales << " eliminada j, por simetria " << endl << endl;
						launcher << " cp " << "bestoutqe/" << estructuras[j].structure_name << ".out " << done_iguales << endl;
						system(launcher.str().c_str());
						//continue;

					}
					else if (estructuras[j].sym_ == estructuras[i].sym_)
					{
						if (estructuras[i].energy_ < estructuras[j].energy_)
						{
							iguales << " eliminada j, por energia " << endl << endl;
							launcher << " cp " << "bestoutqe/" << estructuras[j].structure_name << ".out " << done_iguales << endl;
							system(launcher.str().c_str());
							//continue;
						}
						else
						{
							iguales << " eliminada i, por energia " << endl << endl;
							launcher << " cp " << "bestoutqe/" << estructuras[i].structure_name << ".out " << done_iguales << endl;
							system(launcher.str().c_str());
							break;
						}
					}
					else
					{
						iguales << " eliminada i, por symetria " << endl << endl;
						launcher << " cp " << "bestoutqe/" << estructuras[i].structure_name << ".out " << done_iguales << endl;
						system(launcher.str().c_str());
						break;
					}

					//cout << delta_ang << " delta_angulos " << endl;
					//cout << Maximo_total_ang << " Maximo total angulos " << endl;
				}
				else if ((delta_ang < tolerancia_delta_ang) && (Maximo_total_ang < excepcion_max_ang))
				{
					toCheck_ << " i: " << estructuras[i].structure_name << " j: " << estructuras[j].structure_name << "  ";
					toCheck_ << " delta_dist " << delta_dist << "; " << "Maximo_dist " << Maximo_total_dist << "; " << "delta_ang " << delta_ang << "; " << "Maximo_ang " << Maximo_total_ang << endl << endl;

					launcher << " cp " << "bestoutqe/" << estructuras[i].structure_name << ".out " << files_toCheck << endl;
					launcher << " cp " << "bestoutqe/" << estructuras[j].structure_name << ".out " << files_toCheck << endl;
					system(launcher.str().c_str());
					//final_file << endl;
					continue;

				}
				else
				{
					continue;
				}
			}
			final_file << endl;
			continue;

		}


	}

	resumen_dist.close();
	resumen_ang.close();
	iguales.close();
	final_file.close();
	toCheck_.close();

}

ang_dist estructura::busqueda_minima_dist(const ang_dist& M, const vector<ang_dist>& angulos_dist, double& distancia_minima_calc, const double corte)
{
	double minimo = 1e3;
	ang_dist punto_minimo;

	for (int i = 0; i < angulos_dist.size(); i++)
	{
		double minimo_temp = M.dist_Rn(angulos_dist[i]);

		if (minimo_temp < minimo)
		{
			minimo = minimo_temp;
			punto_minimo = angulos_dist[i];

		}
		if (angulos_dist[i].R > 2 * corte) break;
	}
	double suma_R = M.R + punto_minimo.R;
	distancia_minima_calc = minimo / fabs(suma_R);
	return punto_minimo;
}


void Parameters::Angulos(const int atom_1, const int atom_2, const int atomxmol, super_celda& coor, ang_dist ang_)
{

	Vec3 V_Atom_1 = Vec3(coor.cart_supercelda[0].Celda[atom_1][0], coor.cart_supercelda[0].Celda[atom_1][1], coor.cart_supercelda[0].Celda[atom_1][2]);
	//Vec3 V_Atom_1_frac = Vec3(coor.list_coor_frac[0].Celda[atom_1][0], coor.list_coor_frac[0].Celda[atom_1][1], coor.list_coor_frac[0].Celda[atom_1][2]);

	Vec3 V_Atom_2 = Vec3(coor.cart_supercelda[0].Celda[atom_2][0], coor.cart_supercelda[0].Celda[atom_2][1], coor.cart_supercelda[0].Celda[atom_2][2]);
	//Vec3 V_Atom_2_frac = Vec3(coor.list_coor_frac[0].Celda[atom_2][0], coor.list_coor_frac[0].Celda[atom_2][1], coor.list_coor_frac[0].Celda[atom_2][2]);


	Vec3 R1 = V_Atom_2 - V_Atom_1;

	for (int j = atomxmol; j < nat_; j++) { //Angulos moleculas internas

		Vec3 V_Atom_n = Vec3(coor.cart_supercelda[0].Celda[j][0], coor.cart_supercelda[0].Celda[j][1], coor.cart_supercelda[0].Celda[j][2]);

		Vec3 Rn = V_Atom_n - V_Atom_1;
		double dist_Rn = len(Rn);

		double cos_ang = dot(Rn, R1) / (len(Rn) * len(R1));

		//string atomo = getElemName(atomos_cell[j]); 

		auto it = ang_.AngulosPorAtomo.find(atomos_cell[j]);

		if (it != ang_.AngulosPorAtomo.end()) {
			// Si existe, agregar la distancia (Valor calculado llamado distancia)
			ang_dist cos(dist_Rn, cos_ang);
			it->second.push_back(cos);
		}

	}

	for (int k = 1; k < coor.cart_supercelda.size(); k++) {

		for (int j = 0; j < nat_; j++) {

			Vec3 V_Atom_n = Vec3(coor.cart_supercelda[k].Celda[j][0], coor.cart_supercelda[k].Celda[j][1], coor.cart_supercelda[k].Celda[j][2]);

			Vec3 Rn = V_Atom_n - V_Atom_1;
			double dist_Rn = len(Rn);

			double cos_ang = dot(Rn, R1) / (len(Rn) * len(R1));

			auto it = ang_.AngulosPorAtomo.find(atomos_cell[j]);

			if (it != ang_.AngulosPorAtomo.end()) {
				// Si existe, agregar la distancia (Valor calculado llamado distancia)
				ang_dist cos(dist_Rn, cos_ang);
				it->second.push_back(cos);
			}

		}
	}

	//ang_.mostrar(); 

	for (auto& par : ang_.AngulosPorAtomo) {
		auto& vec = par.second; // vec es el vector<ang_dist>
		std::sort(vec.begin(), vec.end(), comparacion_angulos_tomando_R);
	}

	angulosPorAtomos_Param = ang_.AngulosPorAtomo;

	//Para imprimir 
	/*
	for (const auto& par : ang_.AngulosPorAtomo) {
		Elem tipo = par.first;
		const auto& vec = par.second;

		// Imprimir el tipo de átomo (ajusta getElemName si tienes otro método)
		std::cout << "Tipo: " << getElemName(tipo) << "\n";

		// Imprimir cada distancia (asumo que ang_dist tiene un campo llamado R)
		for (const auto& ang : vec) {
			std::cout << "  R: " << ang.R << "\n";  // Ajusta esto según tu struct/class ang_dist
		}

		std::cout << "-----------------\n";
	}
	*/


}




void Parameters::Enlaces(Distancias& dist_, ang_dist& ang_)
{

	//Distancias dist; 
	/*
	// Mostrar para comprobar que está bien convertido
	for (Elem e : atomos_cell) {
		cout << "Elemento: " << getElemName(e) << endl;
	}
	*/
	//set<Elem> tipos_unicos;  // Usamos set porque guarda solo únicos
	set<Elem> tipos_unicosSet;  // Guarda solo únicos

	for (size_t i = 0; i < atomxmol; ++i) {  // Solo hasta atomxmol
		tipos_unicosSet.insert(atomos_cell[i]);
	}

	//cout << "Cantidad de tipos atómicos distintos (dentro de la celda asimétrica): " << tipos_unicosSet.size() << endl;

	vector<Elem>tipos_unicos(tipos_unicosSet.begin(), tipos_unicosSet.end());
	atoms_types = tipos_unicos;
	/*
	// Mostrar cuáles son
	//cout << "Tipos encontrados:" << endl;
	for (Elem e : tipos_unicos) {
		cout << getElemName(e) << endl;
	}
	*/

	for (const auto& tipo : tipos_unicos) {
		ang_.AngulosPorAtomo[tipo] = std::vector<ang_dist>();  // vector vacío por ahora
	}
	/*
	// Mostrar lo guardado
	for (const auto& par : ang_.AngulosPorAtomo) {
		std::cout << getElemName(par.first) << ": ";
		//for (ang_dist d : par.second) {
			//std::cout << d << " ";
		//}
		std::cout << std::endl;
	}
	*/

	set <string> DistSet;
	//vector<string> Dist;

	for (int i = 0; i < tipos_unicos.size(); ++i) {
		for (int j = 0; j < tipos_unicos.size(); ++j) {
			string elem1 = getElemName(tipos_unicos[i]);
			string elem2 = getElemName(tipos_unicos[j]);

			string combinacion = (elem1 < elem2) ? elem1 + elem2 : elem2 + elem1;
			DistSet.insert(combinacion);
		}
	}

	vector<string>Dist(DistSet.begin(), DistSet.end());

	for (const auto& combina : Dist) {
		cout << combina << endl;
	}


	atoms_combinations = Dist;

	// El mapa: clave = tipo, valor = vector de distancias
	//std::map<std::string, std::vector<double>> distanciasPorTipo;

	// Inicializar el mapa con claves y vectores vacíos

	for (const auto& tipo : Dist) {
		dist_.distanciasPorTipo[tipo] = std::vector<double>();  // vector vacío por ahora
	}
	/*
	// Mostrar lo guardado
	for (const auto& par : dist_.distanciasPorTipo) {
		std::cout << par.first << ": ";
		for (double d : par.second) {
			std::cout << d << " ";
		}
		std::cout << std::endl;
	}
	*/

	/*
	set<string> enlaces_formados;  // Set para no repetir enlaces
	vector<string> enlaces;     // vector donde guardamos los enlaces


	vector<Elem> tipos_vector(tipos_unicos.begin(), tipos_unicos.end());

	for (int i = 0; i < tipos_unicos.size(); ++i) {
		for (int j = i + 1; j < tipos_unicos.size(); ++j) {
			string elem1 = getElemName(tipos_vector[i]);
			string elem2 = getElemName(tipos_vector[j]);

			// Ordenamos para que "C-H" y "H-C" se consideren iguales (si querés esto)
			string enlace = (elem1 < elem2) ? elem1 + "-" + elem2 : elem2 + "-" + elem1;

			// Solo agregamos si no está ya en el set
			if (enlaces_formados.find(enlace) == enlaces_formados.end()) {
				enlaces_formados.insert(enlace);
				enlaces.push_back(enlace);
				cout << enlace << endl;
			}
		}
	}
	*/
}


void Parameters::Distancia_SC(super_celda& coor, Distancias& dist_, const int atomxmol) //Distancias para SUPERCELDA
{
	double distancia;

	if (nat_ > atomxmol)
	{
		for (int i = 0; i < atomxmol; i++) { //distancias entre moleculas internas

			Vec3 V_1 = Vec3(coor.cart_supercelda[0].Celda[i][0], coor.cart_supercelda[0].Celda[i][1], coor.cart_supercelda[0].Celda[i][2]);

			//Vec3 V_1_ = Vec3(coor.list_coor_frac[0].Celda[i][0], coor.list_coor_frac[0].Celda[i][1], coor.list_coor_frac[0].Celda[i][2]);
			//cout << V_1_ << endl;
			//cout << cl.atomic_species[i] << endl;

			string atom_i = getElemName(atomos_cell[i]);

			for (int z = atomxmol; z < nat_; z++)
			{

				Vec3 V_2 = Vec3(coor.cart_supercelda[0].Celda[z][0], coor.cart_supercelda[0].Celda[z][1], coor.cart_supercelda[0].Celda[z][2]);
				//Vec3 V_2_ = Vec3(coor.list_coor_frac[0].Celda[z][0], coor.list_coor_frac[0].Celda[z][1], coor.list_coor_frac[0].Celda[z][2]);
				//cout << cl.atomic_species[z] << endl;
				//cout << V_2_ << endl;
				distancia = len(V_1 - V_2);

				string atom_z = getElemName(atomos_cell[z]);
				//string combinacion = atom_i + atom_z;

				string combinacion = (atom_i < atom_z) ? atom_i + atom_z : atom_z + atom_i;
				//cout << combinacion << endl; 
				//string combinacion = (elem1 < elem2) ? elem1 + elem2 : elem2 + elem1;
				// Buscar si la combinación existe en el mapa
				//dist.distanciasPorTipo es tu map<string, vector<double>>.
				//find() busca directamente si la combinacion existe.
				//Si la encuentra(!= end()), entonces it->second es el vector<double> asociado, y ahí le haces push_back() para agregar la distancia.
				auto it = dist_.distanciasPorTipo.find(combinacion);

				if (it != dist_.distanciasPorTipo.end()) {
					// Si existe, agregar la distancia (Valor calculado llamado distancia)
					it->second.push_back(distancia);
				}

			}
		}

	}

	for (int j = 1; j < coor.cart_supercelda.size(); j++) {

		for (int i = 0; i < atomxmol; i++) {
			//for (int i = 0; i < cl.nat; i++) {
			Vec3 V_1 = Vec3(coor.cart_supercelda[0].Celda[i][0], coor.cart_supercelda[0].Celda[i][1], coor.cart_supercelda[0].Celda[i][2]);
			//Vec3 V_1_ = Vec3(coor.list_coor_frac[0].Celda[i][0], coor.list_coor_frac[0].Celda[i][1], coor.list_coor_frac[0].Celda[i][2]);
			//cout << V_1_ << endl;
		   // cout << cl.atomic_species[i] << endl;

			string atom_i = getElemName(atomos_cell[i]);

			for (int z = 0; z < nat_; z++) {
				Vec3 V_2 = Vec3(coor.cart_supercelda[j].Celda[z][0], coor.cart_supercelda[j].Celda[z][1], coor.cart_supercelda[j].Celda[z][2]);
				//Vec3 V_2_ = Vec3(coor.list_coor_frac[j].Celda[z][0], coor.list_coor_frac[j].Celda[z][1], coor.list_coor_frac[j].Celda[z][2]);
				//cout << cl.atomic_species[z] << endl;
				//cout << V_2_ << endl;
				distancia = len(V_1 - V_2);

				string atom_z = getElemName(atomos_cell[z]);
				//string combinacion = atom_i + atom_z;
				string combinacion = (atom_i < atom_z) ? atom_i + atom_z : atom_z + atom_i;

				//Buscar si la combinación existe en el mapa
				//dist.distanciasPorTipo es tu map<string, vector<double>>.
				//find() busca directamente si la combinacion existe.
				//Si la encuentra(!= end()), entonces it->second es el vector<double> asociado, y ahí le haces push_back() para agregar la distancia.
				auto it = dist_.distanciasPorTipo.find(combinacion);

				if (it != dist_.distanciasPorTipo.end()) {
					// Si existe, agregar la distancia (Valor calculado llamado distancia)
					it->second.push_back(distancia);
				}

			}
		}
	}


	double dgnal = (sqrt((a_ * a_) + (b_ * b_) + (c_ * c_))) * 1.5; // se multiplica por 1.5 para agrandar el radio de corte
	//cout << dgnal << endl;

	diagonal = dgnal;

	//dist_.mostrar();

	for (auto& par : dist_.distanciasPorTipo) {
		auto& vec = par.second; // El vector que quieres ordenar
		std::sort(vec.begin(), vec.end());
	}

	distanciasPorTipo_Param = dist_.distanciasPorTipo;

}



void Parameters::Armado_Supercelda(super_celda& coor, Parameters& params_)
{
	vector <super_celda> i_list_coor_frac_;
	vector <super_celda> i_cart_supercelda_;

	i_list_coor_frac_.push_back(coor_Frac_);
	const int semi_SC = 7;


	//SUPER CELDA FRAC************ 
	Mat tempCelda(nat_, 3);
	for (int i = -semi_SC; i < semi_SC + 1; i++) {
		for (int j = -semi_SC; j < semi_SC + 1; j++) {
			for (int k = -semi_SC; k < semi_SC + 1; k++) {
				if (i != 0 || j != 0 || k != 0)
				{
					for (int l = 0; l < nat_; l++) {
						tempCelda[l][0] = i_list_coor_frac_[0].Celda[l][0] + i;
						tempCelda[l][1] = i_list_coor_frac_[0].Celda[l][1] + j;
						tempCelda[l][2] = i_list_coor_frac_[0].Celda[l][2] + k;
					}
					i_list_coor_frac_.push_back(tempCelda);
				}
			}
		}
	}

	//SUPERCELDA A CART****** 

	Vec3 temp3;
	Mat tempmat(nat_, 3);
	Vec3 temp3CART;
	Mat3 prueba_CART;
	for (int i = 0; i < i_list_coor_frac_.size(); i++) {
		for (int j = 0; j < nat_; j++) {
			temp3 = Vec3(i_list_coor_frac_[i].Celda[j][0], i_list_coor_frac_[i].Celda[j][1], i_list_coor_frac_[i].Celda[j][2]);
			temp3CART = Frac_2Cart(params_, temp3) * AUtoAng * alat_;
			//temp3CART = Frac_2Cart(va_, vb_, vc_, temp3) * AUtoAng * alat_;
			tempmat[j][0] = temp3CART[0];
			tempmat[j][1] = temp3CART[1];
			tempmat[j][2] = temp3CART[2];
		}
		i_cart_supercelda_.push_back(tempmat);

	}
	coor.list_coor_frac.clear();
	coor.cart_supercelda.clear();

	//cout << coor.list_coor_frac.size() << " 1 deberia ser cero" << endl;
	//cout << coor.cart_supercelda.size() << " 2 deberia ser cero" << endl;

	//coor.list_coor_frac = i_list_coor_frac_;
	coor.cart_supercelda = i_cart_supercelda_;

	//cout << coor.list_coor_frac.size() << " 3 " << endl;
	//cout << coor.cart_supercelda.size() << " 4 " << endl;

	i_list_coor_frac_.clear();
	i_cart_supercelda_.clear();

}
void Parameters::cell_parameters(string file_out, super_celda& coor)
{
	int k = 0;
	double alat = 0.0; // alat
	//Vec3 va, vb, vc; //vectores celda
	Vec3 temp; // Ultimas Coordenadas output (fraccionarias)
	istringstream ss;
	string junk;
	vector <string> atomos;

	const string atoms_cell = "number of atoms/cell      =";
	const string cell_params = "CELL_PARAMETERS";
	const string atom_pos = "ATOMIC_POSITIONS (crystal)";
	const string _energy = "!    total energy              =";
	const string no_symmetry = "No symmetry found";
	const string _symmetry = "Sym. Ops.";
	size_t a_pos = 0;
	size_t c_pos = 0;

	//get nat
	if ((file_out.find(atoms_cell) != string::npos))
	{
		string tempatoms = file_out.substr(file_out.find(atoms_cell) + atoms_cell.size());
		istringstream iss(tempatoms);
		int val;
		iss >> val;
		nat_ = val;
		//cout << "number of atoms/cell= " << nat_ << endl;
	}

	int moleculasxCelda = nat_ / atomxmol;

	cout << fixed << setprecision(5);

	//get energy 
	size_t pos = file_out.rfind(_energy);
	if (pos != string::npos)
	{
		ss.clear();
		pos += _energy.size();
		ss.str(file_out.substr(pos, 160));
		ss >> energy_;
		//cout << endl << "energy:" << energy_ << endl;
	}

	energy_kj = energy_;
	energy_kj *= (RydToKCalMol * CaltoJoule);
	energy_kj = energy_kj / moleculasxCelda;

	//cout << endl << " enegy Kj  " << energy_kj;

	//get Enthalpy
	const string energy_scan = "!    total energy              =";
	double energynotfinal = 0.0;
	double enthalpia = 0.0;
	const string enthalpy_scan = "Final enthalpy           =    ";

	// Buscar la última ocurrencia de energy_scan
	pos = file_out.rfind(energy_scan);
	if (pos != string::npos) {
		size_t pos_second_last;
		pos_second_last = file_out.rfind(energy_scan, pos - 1);
		if (pos_second_last != string::npos) {
			// Extraer la penúltima energía
			ss.clear();
			pos_second_last += energy_scan.size();
			ss.str(file_out.substr(pos_second_last, 160));
			ss >> energynotfinal;
			//cout << endl << "Penultimate energy: " << energynotfinal << endl;
		}
	}

	size_t entalpy_pos = file_out.rfind(enthalpy_scan);
	if (pos != string::npos) {

		ss.clear();
		entalpy_pos += enthalpy_scan.size();
		ss.str(file_out.substr(entalpy_pos, 160));
		ss >> enthalpia;
		//cout << fixed << setprecision(8);
		//cout << endl << " not final enthalpy: " << enthalpia << endl;
	}

	double PV = enthalpia - energynotfinal;

	double FinalEnthalpy = energy_ + PV;
	enthalpy = FinalEnthalpy;
	//cout << endl << " Final enthalpy  " << enthalpy; 

	//energy_KJ = energy/sym;
	enthalpy_kj = FinalEnthalpy;
	enthalpy_kj *= (RydToKCalMol * CaltoJoule);
	enthalpy_kj = enthalpy_kj / moleculasxCelda;

	//cout << endl << " enthalpy Kj  " << enthalpy_kj;

	//get symmetry
	if (file_out.find(no_symmetry) != string::npos) sym_ = 1;
	else
	{
		size_t k_pos = file_out.find(_symmetry);
		//cout << _symmetry.size();
		k_pos = k_pos - 2;
		//k_pos -= _symmetry.size();
		ss.clear();
		ss.str(file_out.substr(k_pos));
		ss >> sym_;
		//cout << sym_ << " simetria" << endl;
	}

	// get cells parameters
	a_pos = file_out.rfind(cell_params);
	if (a_pos != string::npos)
	{
		//get alat
		ss.clear();
		ss.str(file_out.substr(a_pos, 800 * 4));
		//ss >> junk >> junk >> alat >> junk;
		ss >> junk >> junk >> alat_ >> junk;
		//cout << " alat: " << alat << endl;
		//cout << " alat: " << cl.alat1 << endl;

		//get vectors and angles
		ss >> va_[0] >> va_[1] >> va_[2];
		ss >> vb_[0] >> vb_[1] >> vb_[2];
		ss >> vc_[0] >> vc_[1] >> vc_[2];

		//cout << "vectores: " << endl << cl.va[0] << "  " << cl.va[1] << "  " << cl.va[2] << endl;
		//cout << "vectores: " << endl << cl.vb[0] << "  " << cl.vb[1] << "  " << cl.vb[2] << endl;
		//cout << "vectores: " << endl << cl.vc[0] << "  " << cl.vc[1] << "  " << cl.vc[2] << endl;

		a_ = alat_ * AUtoAng * len(va_);
		b_ = alat_ * AUtoAng * len(vb_);
		c_ = alat_ * AUtoAng * len(vc_);
		gamma_ = angle(va_, vl_0, vb_);
		beta_ = angle(va_, vl_0, vc_);
		alpha_ = angle(vb_, vl_0, vc_);
		//cout << a_ << " " << b_ << " " << c_ << endl;
		//cout << alpha_ << " " << beta_ << " " << gamma_ << endl;


		//Mat3 toCart = fracToCart(cl);
		//cout << toCart << endl;
		 // get the coordinates 
		Mat lista_Frac(nat_, 3); // Matriz coordenadas Fraccionarias output

		c_pos = file_out.rfind(atom_pos);
		if (c_pos != string::npos)
		{
			ss.clear();
			ss.str(file_out.substr(c_pos + atom_pos.size(), 800 * nat_));
			ss >> junk >> junk;
			string temporal = ss.str();
			istringstream sOut(temporal);
			int contador = 0;
			while (sOut)
			{
				string final;
				sOut >> final;
				if (contador >= nat_)
					//if (final == "End" || final == "Writing")
					break;
				sOut >> temp[0] >> temp[1] >> temp[2]; //temp coordenadas fraccionarias como Vec3
				//cout << final << temp[0] << " " << temp[1] << " " << temp[2] << endl;
				//cl.atomic_species.push_back(final);
				atomos.push_back(final);

				//  Matriz 1. Matriz Fraccionarias 
				for (int i = 0; i < 3; i++) {
					lista_Frac[k][i] = temp[i];
				}
				k = k + 1;

				//Vec 3 Cart = Frac_2Cart(cl.va, cl.vb, cl.vc, temp);

				//Cart_ = Frac_2Cart(va_, vb_, vc_, temp);
				//cout << cl.Cart << endl;
				//cout << endl;

				// Vec3 Cart_A = Cart * alat * AUtoAng;
				// cout << 2 << endl;
				// cout << Cart_A << endl;
				contador++;
			}
			coor_Frac_ = lista_Frac;
			//cout << cl.coor_Frac_ << endl;

			//vector <super_celda> out_;
			//out_.push_back(lista_Frac);
			//coor.list_coor_frac = out_;


			//coor.list_coor_frac.push_back(cl.coor_Frac_);
			//atomos_cell = atomos;

			for (const string& nombre : atomos) {
				atomos_cell.push_back(getElemType(nombre));
			}

			//cout << endl << endl;
			//cout << "TAMANO" << coor.list_coor_frac.size() << endl;
			//for (int i = 0; i < coor.list_coor_frac.size(); i++) {
			//	cout << coor.list_coor_frac[i].Celda << endl;
			//}
		}
	}
}
string get_file_contents(const char* filename)
{
	std::ifstream in(filename, std::ios::in | std::ios::binary);
	if (in) {
		std::string contents;
		in.seekg(0, std::ios::end);
		contents.resize(in.tellg());
		in.seekg(0, std::ios::beg);
		in.read(&contents[0], contents.size());
		in.close();
		return(contents);
	}
	else// if(in.eof())
	{
		return string("");
	}
}
bool outputVerification(string fileout)
{
	bool outOK = false;

	const string finishOK = "JOB DONE.";
	const string end_coodinates = "End final coordinates";
	const string any_routine_error = "Error in routine";
	const string scf_final = "'scf'";
	const string control_in = "&control"; //para acotar la busqueda de los errores en la ultima salida, y no vaya a buscarlos en los checkpoints previos. 
	const string system_in = "&system"; //asegura que sea el input, buscando &control y &system. y no un error que aparece y coloca &control en el error in routine

	string ultimo_out; //va desde el ultimo input generado del checkpoint e incluye el output final. 

	if (fileout.rfind(control_in) != string::npos && fileout.rfind(system_in) != string::npos)
		ultimo_out = fileout.substr(fileout.rfind(control_in)); //del archivo que sale con los checkpoints, tomamos solo el ultimo input y output. 
	else ultimo_out = fileout;

	if (ultimo_out.size() == 0)
	{
		cout << " outout vacio  " << endl;
		outOK = false;
		return outOK;
	}
	else if (ultimo_out.rfind(any_routine_error) != string::npos)
	{
		cout << " Error in routine " << endl;
		outOK = false;
		return outOK;
	}
	else if (ultimo_out.rfind(end_coodinates) != string::npos && ultimo_out.rfind(finishOK) != string::npos)
	{
		//cout << " job ok " << endl;
		outOK = true;
		return outOK;
	}
	else if (ultimo_out.rfind(scf_final) != string::npos && ultimo_out.rfind(finishOK) != string::npos)
	{
		outOK = true;
		return outOK;
	}
	else
	{
		cout << "ninguno de los anteriores " << endl;
		outOK = false;
		return outOK;
	}

}
//calculate the angle between three points
double angle(const Vec3& A, const Vec3& B, const Vec3& C) {
	double cos_angle = dot(A - B, C - B) / (len(A - B) * len(C - B));

	if (cos_angle > 1.0) cos_angle = 1.0;

	return (acos(cos_angle));
}
Vec3 Frac_2Cart(Parameters& param_, Vec3 Frac) {
	Vec3 Coor_Cart = (param_.va_ * Frac[0]) + (param_.vb_ * Frac[1]) + (param_.vc_ * Frac[2]);
	return Coor_Cart;
}