#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <mpi.h>

#include "model.hpp"
#include "display.hpp"

using namespace std::string_literals;
using namespace std::chrono_literals;

struct ParamsType
{
    double length{1.};
    unsigned discretization{20u};
    std::array<double,2> wind{0.,0.};
    Model::LexicoIndices start{10u,10u};

    void serializable(char* buffer) const {
        memcpy(buffer, this, sizeof(ParamsType));
    }

    void deserializable(const char* buffer) {
        memcpy(this, buffer, sizeof(ParamsType));
    }
};

void analyze_arg( int nargs, char* args[], ParamsType& params )
{
    if (nargs ==0) return;
    std::string key(args[0]);
    if (key == "-l"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour la longueur du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.length = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    auto pos = key.find("--longueur=");
    if (pos < key.size())
    {
        auto subkey = std::string(key,pos+11);
        params.length = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-n"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour le nombre de cases par direction pour la discrétisation du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.discretization = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--number_of_cases=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+18);
        params.discretization = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-w"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la direction du vent !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.wind[0] = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--wind=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+7);
        params.wind[0] = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-s"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la position du foyer initial !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.start.column = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la position du foyer initial" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--start=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+8);
        params.start.column = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }
}

ParamsType parse_arguments( int nargs, char* args[] )
{
    if (nargs == 0) return {};
    if ( (std::string(args[0]) == "--help"s) || (std::string(args[0]) == "-h") )
    {
        std::cout << 
R"RAW(Usage : simulation [option(s)]
  Lance la simulation d'incendie en prenant en compte les [option(s)].
  Les options sont :
    -l, --longueur=LONGUEUR     Définit la taille LONGUEUR (réel en km) du carré représentant la carte de la végétation.
    -n, --number_of_cases=N     Nombre n de cases par direction pour la discrétisation
    -w, --wind=VX,VY            Définit le vecteur vitesse du vent (pas de vent par défaut).
    -s, --start=COL,ROW         Définit les indices I,J de la case où commence l'incendie (milieu de la carte par défaut)
)RAW";
        exit(EXIT_SUCCESS);
    }
    ParamsType params;
    analyze_arg(nargs, args, params);
    return params;
}

bool check_params(ParamsType& params)
{
    bool flag = true;
    if (params.length <= 0)
    {
        std::cerr << "[ERREUR FATALE] La longueur du terrain doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if (params.discretization <= 0)
    {
        std::cerr << "[ERREUR FATALE] Le nombre de cellules par direction doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if ( (params.start.row >= params.discretization) || (params.start.column >= params.discretization) )
    {
        std::cerr << "[ERREUR FATALE] Mauvais indices pour la position initiale du foyer" << std::endl;
        flag = false;
    }
    
    return flag;
}

void display_params(ParamsType const& params)
{
    std::cout << "Parametres définis pour la simulation : \n"
              << "\tTaille du terrain : " << params.length << std::endl 
              << "\tNombre de cellules par direction : " << params.discretization << std::endl 
              << "\tVecteur vitesse : [" << params.wind[0] << ", " << params.wind[1] << "]" << std::endl
              << "\tPosition initiale du foyer (col, ligne) : " << params.start.column << ", " << params.start.row << std::endl;
}

int main(int nargs, char* args[]) {
    MPI_Init(&nargs, &args);
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Création du groupe de calcul
    MPI_Comm compute_comm;
    MPI_Group world_group, compute_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    int exclude[] = {0};
    MPI_Group_excl(world_group, (rank == 0) ? 0 : 1, exclude, &compute_group);
    MPI_Comm_create(MPI_COMM_WORLD, compute_group, &compute_comm);

    ParamsType params;
    char param_buffer[sizeof(ParamsType)];

    if(rank == 0) {
        params = parse_arguments(nargs-1, &args[1]);
        if (!check_params(params)){
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        params.serializable(param_buffer);
            // Préparation pour MPI_Gatherv
            std::vector<int> counts(num_procs-1), displs(num_procs-1);
            int offset = 0;
            for(int i=1; i<num_procs; i++) {
                counts[i-1] = /* calcul de la taille locale */;
                displs[i-1] = offset;
                offset += counts[i-1];
        }
    }
    MPI_Bcast(param_buffer, sizeof(ParamsType), MPI_BYTE, 0, MPI_COMM_WORLD);

    if(rank != 0) {
        params.deserializable(param_buffer);
    }

    if (rank == 0) {
        try {
            // Initialisation de l'affichage
            auto displayer = Displayer::init_instance(params.discretization, params.discretization);
            
            // Préparation pour MPI_Gatherv
            std::vector<int> counts(num_procs-1), displs(num_procs-1);
            int offset = 0;
            for(int i=1; i<num_procs; i++) {
                int local_size = /* calculer la taille locale du processus i */;
                counts[i-1] = local_size;
                displs[i-1] = offset;
                offset += local_size;
            }
        
            std::vector<uint8_t> fire_map(params.discretization * params.discretization);
            std::vector<uint8_t> vegetation_map(params.discretization * params.discretization);
        
            while(running) {
                // Collecte des données
                MPI_Gatherv(nullptr, 0, MPI_DATATYPE_NULL,
                        fire_map.data(), counts.data(), displs.data(), MPI_UNSIGNED_CHAR,
                        0, MPI_COMM_WORLD);
                
                // Mise à jour de l'affichage
                displayer->update(vegetation_map, fire_map);
                
                // Gestion des événements SDL
                if(SDL_PollEvent(&event) && event.type == SDL_QUIT) running = false;
                std::this_thread::sleep_for(10ms);
        }
    }
        catch(const std::exception& e) {
            std::cerr << "Display error: " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }
    else {
        try {
            // Découpage du domaine
            int N = params.discretization;
            int rows_per_proc = N / (num_procs-1);
            int start_row = (rank-1) * rows_per_proc;
            int end_row = (rank == num_procs-1) ? N : start_row + rows_per_proc;

            // Ajout des fantômes
            start_row = std::max(start_row - 1, 0);
            end_row = std::min(end_row + 1, N);

            // Création du modèle local
            Model simu(
                params.length, 
                params.discretization, 
                params.wind, 
                params.start,
                60.0,          // t_max_wind
                start_row,     // t_local_start
                end_row        // t_local_end
            );

            while(simu.update()) {
                simu.exchange_ghosts(); // Échange des cellules fantômes
                
                // Envoi des données au maître via Gatherv
                MPI_Gatherv(simu.fire_map().data(), simu.fire_map().size(), MPI_UNSIGNED_CHAR,
                        nullptr, nullptr, nullptr, MPI_UNSIGNED_CHAR, 
                        0, MPI_COMM_WORLD);
            }

            std::cout << "Temps moyen par iteration: " << total_time.count() / total_steps << " s" << std::endl;

        }
        catch(const std::exception& e) {
            std::cerr << "Compute error: " << e.what() << std::endl;
        }
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}