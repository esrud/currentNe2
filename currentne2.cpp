/*
 * Program: currentNe_v2.0
 *
 * Description: Neural Network for Ne Estimation
 *
 * Authors: Enrique Santiago, Carlos Köpke
 *
 * License: TBD
 */

#include <algorithm>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <random> // std::rand, std::srand
#include <sstream>
#include <stdlib.h>
#include <string>
#include <time.h>

#include "lib/progress.hpp"

#define MAXLOCI 2000000
#define MAXIND 2000
// #define MAXSIBPAIRS 50000
#define MAXCROMO 1000
#define MAXDIST 100000 // desde 0.01cM hasta 1000cM (10M)

// void IntegralTot();
void IntegralUNCROM();
void MixIntegralUNCROM();
void IntegralVARIOSCROM();
// void IntegralCromUMBRALES();
void Ecuacion05();
void MixEcuacion05();
// double Funcion_751_soloLD_mae_ADAM_5();
// double Funcion_desv2logint_manip_log10Lsigmoid_mae_ADAM_4();
double Funcion_desv2logsoloLD_manip_log10Lsigmoid_mae_ADAM_4();
// double CalculaIntervalo_int();
double CalculaIntervalo_soloLD();

char indi[MAXIND][MAXLOCI] = {'\0'}; //** Genotipos dip 0:homo, 1:het, 2:homo
char base[MAXLOCI] = {'\0'};         // Contiene las bases de referencia de cada locus
double posiCM[MAXLOCI] = {0};       // posiciones geneticas en cM del .map
double posiBP[MAXLOCI] = {0};       // posiciones geneticas en cM del .map
double mapdist[MAXDIST] = {0};       // contiene la distribucion de los pares de distancias
char cromo[MAXLOCI];                 // contiene el número del cromosoma de cada locus
int rangocromo[MAXCROMO], ncromos;
double frec[MAXLOCI], homo[MAXLOCI];
double Het[MAXIND], Parent[MAXIND], Omega[MAXIND];
double rangos_05[5],rangos_link[5], rangos_int[5], rangos_soloLD[5];
double chrprop[MAXCROMO] = {0}, chrsize[MAXCROMO] = {0}, cmedchr[MAXCROMO] = {0}; // proporción de SNPs y distancia en M entre extremos
double chrmin[MAXCROMO] = {0}, chrnum[MAXCROMO] = {0};                            // distancia mínima entre parejas consecutivas de marcadores y número de marcadores
bool segrega[MAXLOCI] = {false};
bool validind[MAXIND] = {true};
int ploc[MAXLOCI] = {0};
int pind[MAXIND] = {0};
int *pj, *pi, *prefloc, *prefind, *pk, *pfinloc, *pfinind;
int containd, contaloc, contaloc2, contalocbase, contaseg, eneind, eneloc, enelocsolicitado;
int enecrom, n_sample, n_SNPs, n_threads, ff, gg, ncrom_sample_int, imindist, imaxdist;
double acuD2, acuW, acun, acupq, acuPp2, acur2, acuD205, acuW05, acur205, acun05;
double acuD2link, acuWlink, acur2link, acunlink, effndatalink = 0;
double n, d2s, d2s05, d2slink, d2_poblink, L, fs, fp, backfp, Ch, ele, Ne, Neant, genomesize;
double Ne_integral_crom, Ne_integral_tot, Ne_integral_totcrom, Ne_05, Ne_link,  cM, Ncrom;
double Ne_nohets, Ne_hets, effndata = 0, effndata05 = 0, obsndata = 0, propmiss, n_SNP_pairs, effeneind, effeneind_h;
double acuParent = 0, acuHet = 0, Het_med = 0, Het_esp = 0, Het_var = 0, Het_DT, Het_sesg = 0;
double Parent_med = 0, Parent_var = 0, Parent_DT, Parent_sesg = 0, r_Het_Parent = 0;
double log10_n_sample, log10_nSNPs_solicitado, log10_Ne_obs, ncrom_sample, log10_ncrom_sample, f_pob, d2_pob, d2_pob05;
double DT, lim, posicMacu = 0, posicMant = 0, unomenoschrprop;
double tini, tpas, sumWwithin, sumWbetween;
long int x_containdX[MAXLOCI] = {0};
long int x_contapares[MAXLOCI] = {0};
double xD[MAXLOCI] = {0}, xW[MAXLOCI] = {0}, xr2[MAXLOCI] = {0};
long int x_containdX05[MAXLOCI] = {0};
long int x_contapares05[MAXLOCI] = {0};
long int x_containdXlink[MAXLOCI] = {0};
long int x_contapareslink[MAXLOCI] = {0};
double xD05[MAXLOCI] = {0}, xW05[MAXLOCI] = {0}, xr205[MAXLOCI] = {0};
double xDlink[MAXLOCI] = {0}, xWlink[MAXLOCI] = {0}, xr2link[MAXLOCI] = {0};
bool flagconvergence = true;
int contaelim = 0, enefuerzaajuste = 1;
int nparhermanos = 0, npadrehijo = 0, counthethet = 0, counthethomo = 0, counthomohomo = 0, sumahethomo, countnopadre;
// int hermanos[2][MAXSIBPAIRS];
// int padrehijo[2][MAXIND];
double ratiohets, rationopadre;
double frecmed, nfrecmed;
double d2p = 0, d2p05 = 0, d2pcmed=0;
double Fstmax, increFst;
int fciclo;
bool flagnoestimaks = false, haygenotipos = false;
double Ne05, Nelink, Neinicial;
double cmed = 0;
double maxdistance, increMorgans, sumadist=0;
int maxdistanceindx;
double m, ms, Fst, Dw2, Db2, DbDw,Wt;
double Dw2p, Db2p, DbDwp;
double AA,BB,CC,DD,EE, c, c2, c12,m12,m22;
bool flag_r=false;
bool flag_z=false;
bool flag_Gs=false;
bool flag_chr=false;
bool flag_cM=false;
bool flag_Mb=false;
double Mtot=0, Mbtot=0;

// Variables privadas para paralelizacion:
int *ppi;
int *ppk;
int *ppj;
int ss, _containdX;
double tacuHoHo, tacuHoHetHetHo, tacuHetHet;
double tacui, tacuj;
double D, W, r2, distancia;
int jj2, jj3, id;

struct AppParams
{
    int numThreads;
    int numSample;
    int numSNPs;
    double s;
    double Mchr;
    double umbral;
    double K;
    double r;
    double z;
    double ks;
    double miss;
    bool flagks;
    bool flagnok;
    bool quiet;
    bool printToStdOut;
    bool verbose;
    ProgressStatus progress;
};

struct AppParams params =
    {
        .numThreads = 0,
        .numSample = 0,
        .numSNPs = 0,
        .s = 1,
        .Mchr = 1.0,
        .umbral = 40,
        .K = 0,
        .r = 1,
        .z = 0,
        .ks = 0,
        .miss = 0.2,
        .flagks = false,
        .flagnok = true,
        .quiet = false,
        .printToStdOut = false,
        .verbose = false,
        .progress = ProgressStatus()};

struct PopulationInfo
{
    int numIndividuals;
    int numLoci;
    int numcromo;
};

void readFile_ped(std::string fichinput1, std::string fichinput2, char (&population)[MAXIND][MAXLOCI], PopulationInfo(&popInfo))
{
    /*
     * Takes as input the file name and a pointer to the population
     * matrix and returns the number of individuals in the file
     */
    char base1[1], base2[1];
    int contaLociBase = 0, nline = 0;
    int conta = 0, posi = 0, posi2 = 0, longi = 0, i;
    std::string line;
    std::string cromocod;
    std::string cromocodback = "laksjhbqne";
    std::string str_plink = "AGCTNagctn0123456789";
    std::string str_base = "";
    bool rightletter;

    // READING .ped DATA:
    std::ifstream entrada;
    entrada.open(fichinput1, std::ios::in); // Bucle de lectura del fichero ped
    if (!entrada.good())
    {
        std::cerr << "Could not open \"" << fichinput1 << "\". Does the file exist?" << std::endl;
        exit(EXIT_FAILURE);
    }
    while (std::getline(entrada, line))
    {
        ++nline;
        longi = int(line.length());
        if (longi < 12)
        {
            std::cerr << "Line "<<nline<<" too short in ped file" << std::endl;
            exit(EXIT_FAILURE);
        }
        conta = 0;
        posi = 0;
        while ((posi < longi) && (conta < 6))
        {
            posi2 = posi;
            posi = int(line.find_first_of(" \t", posi2));
            if (posi < 0)
            {
                std::cerr << "Line "<<nline<<" too short in ped file" << std::endl;
                exit(EXIT_FAILURE);
            }
            ++posi;
            ++conta;
        }
        if (conta == 6)
        {
            popInfo.numLoci = 0;
            while (posi < longi)
            { // asigna genot.
                base1[0] = line.at(posi);
                posi2 = posi;
                posi = int(line.find_first_of(" \t", posi2));
                if (posi < 0)
                {
                    break;
                }
                ++posi;
                base2[0] = line.at(posi);
                rightletter = true;
                str_base = base1[0];
                if (str_plink.find(str_base) == str_plink.npos)
                {
                    rightletter = false;
                }
                str_base = base2[0];
                if (str_plink.find(str_base) == str_plink.npos)
                {
                    rightletter = false;
                }
                if (!rightletter)
                {
                    std::cerr << "Wrong allele in line " << nline << " of ped file (maybe also in other lines). Only bases A, G, C, T and N (case-insensitive) and numbers from 0 to 9 are allowed." << std::endl;
                    exit(EXIT_FAILURE);
                }
                if ((base1[0] != '0') && (base2[0] != '0') && (base1[0] != 'N') && (base2[0] != 'N') && (base1[0] != 'n') && (base2[0] != 'n'))
                {
                    if (base[popInfo.numLoci] == '\0')
                    {
                        base[popInfo.numLoci] = base1[0];
                    }
                    if (base1[0] != base[popInfo.numLoci])
                    {
                        base1[0] = 'X';
                    }
                    if (base2[0] != base[popInfo.numLoci])
                    {
                        base2[0] = 'X';
                    }

                    // 0:homo ref, 1:het, 2:homo noref
                    if (base1[0] == base2[0])
                    {
                        if (base1[0] == base[popInfo.numLoci])
                        {
                            population[popInfo.numIndividuals][popInfo.numLoci] = 0;
                        }
                        else
                        {
                            population[popInfo.numIndividuals][popInfo.numLoci] = 2;
                        }
                    }
                    else
                    {
                        population[popInfo.numIndividuals][popInfo.numLoci] = 1;
                    }
                }
                else
                {
                    population[popInfo.numIndividuals][popInfo.numLoci] = 9; // '9' = Genotipo sin asignar
                }
                posi2 = posi;
                posi = int(line.find_first_of(" \t", posi2));
                if (posi < 0)
                {
                    posi = longi;
                }
                ++posi;
                ++popInfo.numLoci;
                if (popInfo.numLoci >= MAXLOCI)
                {
                    std::cerr << "Reached max number of loci (" << MAXLOCI << ")" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }

            if (popInfo.numIndividuals == 0)
            {
                contaLociBase = popInfo.numLoci;
            }

            if (popInfo.numLoci != contaLociBase)
            {
                std::cerr << "Some genomes in the sample are of different sizes" << std::endl;
                exit(EXIT_FAILURE);
            }

            popInfo.numIndividuals++;
            if (popInfo.numIndividuals > MAXIND)
            {
                std::cerr << "Reached limit of sample size (" << MAXIND << ")" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    entrada.close();

    // READING .map DATA:
    flag_chr = false;
    ncromos = 0;
    entrada.open(fichinput2, std::ios::in); // Bucle de lectura del fichero map
    if (entrada.good())
    {
        flag_chr = true;
    }
    if (flag_chr)
    {
        int contalines = 0;
        while (std::getline(entrada, line))
        {
            longi = int(line.length());
            if (longi < 5)
            {
                std::cerr << "Line "<<contalines+1<<" too short in map file" << std::endl;
                exit(EXIT_FAILURE);
            }

            posi = int(line.find_first_of(" \t", 0));
            if (posi <= 0)
            {
                std::cerr << "Empty line in map file" << std::endl;
                exit(EXIT_FAILURE);
            }
            cromocod = line.substr(0, posi);
            if (cromocod != cromocodback)
            {
                cromocodback = cromocod;
                rangocromo[ncromos] = contalines;
                ++ncromos;
            }
            cromo[contalines] = ncromos;

            ++posi;  // Nombre del SNP que no se lee
            posi2 = posi;
            posi=static_cast<int>(line.find_first_of(" \t",posi2));
            if (posi <= 0) {
            std::cerr << "Error in map file (1)" << std::endl;
            exit(EXIT_FAILURE);
            }

            ++posi;  // Localizacion cM
            posi2 = posi;
            posi=static_cast<int>(line.find_first_of(" \t",posi2));
            if (posi <= 0) {
            std::cerr << "Error in map file (2)" << std::endl;
            exit(EXIT_FAILURE);
            }
            posiCM[contalines] = std::stod(line.substr(posi2, posi - posi2));
            if (posiCM[contalines]>0){flag_cM=true;}
            ++posi;

            if (posi > longi) {  // Localizacion bp
            std::cerr << "Error in map file (3)" << std::endl;
            exit(EXIT_FAILURE);
            }
            posiBP[contalines] = std::stoi(line.substr(posi, longi - posi));
            if (posiBP[contalines]>0){flag_Mb=true;}

            ++contalines;
        }

        rangocromo[ncromos] = contalines;

        if (ncromos < 2)
        {
            flag_chr = false;
        }
        entrada.close();
        if (popInfo.numLoci != contalines)
        {
            std::cerr << "Different number of loci in ped and map files" << std::endl;
            exit(EXIT_FAILURE);
        }
        popInfo.numcromo = ncromos;

        if (flag_cM){
            Mtot = 0;
            for (int conta = 0; conta < ncromos; ++conta) {
                Mtot += posiCM[rangocromo[conta+1]-1]-posiCM[rangocromo[conta]];
            }
            Mtot /= 100.0;       // en Morgans
        }

        if (flag_Mb){
            Mbtot = 0;
            for (int conta = 0; conta < ncromos; ++conta) {
                Mbtot += posiBP[rangocromo[conta+1]-1]-posiBP[rangocromo[conta]];
            }
            Mbtot /= 1000000.0;  // en megabases
            if (flag_r){
                for (int conta = 0; conta < popInfo.numLoci; ++conta) {
                    posiCM[conta] = params.r * posiBP[conta] / 1000000;
                }
                Mtot = Mbtot*params.r/100;  // en Morgans
                flag_cM=true;
            }        
        }

        if (flag_cM && params.s>1 && !flag_z){flag_z=true;params.z=5;}

        if (flag_chr)
        {
            unomenoschrprop = 0;
            if(!flag_Gs){genomesize=Mtot;}
            for (i = 0; i < ncromos; ++i)
            {
                chrnum[i] = rangocromo[i + 1] - rangocromo[i]; // número de marcadores en cada cromosoma
                chrprop[i] = chrnum[i] / contalines;           // Proporción de SNPs en cada cromosoma
                if (flag_Gs){
                    chrsize[i] = chrprop[i] * genomesize; // tamaño en Morgans
                }
                else if (flag_cM){
                    chrsize[i] = (posiCM[rangocromo[i+1]-1]-posiCM[rangocromo[i]])/100; // tamaño en Morgans
                }
                unomenoschrprop += chrprop[i];
            }
            unomenoschrprop = 1 - unomenoschrprop;
        }
    }
}

void readFile_tped(std::string fichinput1, char (&population)[MAXIND][MAXLOCI], PopulationInfo(&popInfo))
{
    /*
     * Takes as input the file name and a pointer to the population
     * matrix and returns the number of individuals in the file
     */
    char base1[1], base2[1];
    int contaIndBase = 0, nline = 0;
    int conta = 0, posi = 0, posi2 = 0, longi = 0, i;
    std::string line;
    std::string cromocod;
    std::string cromocodback = "laksjhbqne";
    std::string str_plink = "AGCTNagctn0123456789";
    std::string str_base = "";
    bool rightletter;

    // READING .tped DATA:
    std::ifstream entrada;
    entrada.open(fichinput1, std::ios::in); // Bucle de lectura del fichero tped
    if (!entrada.good())
    {
        std::cerr << "Could not open \"" << fichinput1 << "\". Does the file exist?" << std::endl;
        exit(EXIT_FAILURE);
    }
    flag_chr = false;
    ncromos = 0;
    int contalines = 0;
    popInfo.numLoci = 0;
    while (std::getline(entrada, line))
    {
        ++nline;
        longi = int(line.length());
        if (longi < 12)
        {
            std::cerr << "Line "<<nline<<" too short in tped file" << std::endl;
            exit(EXIT_FAILURE);
        }
        conta = 0;
        posi = 0;
        // chr name
        posi2 = posi;
        posi = int(line.find_first_of(" \t", posi2));
        if (posi < 0)
        {
            std::cerr << "Line "<<nline<<" too short in tped file" << std::endl;
            exit(EXIT_FAILURE);
        }
        cromocod = line.substr(0, posi);
        if (cromocod != cromocodback)
        {
            cromocodback = cromocod;
            rangocromo[ncromos] = contalines;
            ++ncromos;
        }
        cromo[contalines] = ncromos;

        ++contalines; // para crom

        ++posi;
        ++conta;

        // Ignora las primeras columnas
        while ((posi < longi) && (conta < 2))
        {
            posi2 = posi;
            posi = int(line.find_first_of(" \t", posi2));
            if (posi < 0)
            {
                std::cerr << "Line "<<nline<<" too short in tped file" << std::endl;
                exit(EXIT_FAILURE);
            }
            ++posi;
            ++conta;
        }

        // Position in morgans or centimorgans
        posi2 = posi;
        posi = static_cast<int>(line.find_first_of(" \t", posi2));
        posiCM[popInfo.numLoci] =  std::stod(line.substr(posi2, posi-posi2));
        if (posiCM[popInfo.numLoci]>0){flag_cM=true;}
        ++posi;
        ++conta;
        // Base-pair coordinate
        posi2 = posi;
        posi = static_cast<int>(line.find_first_of(" \t", posi2));
        posiBP[popInfo.numLoci] = std::stoi(line.substr(posi2, posi-posi2));
        if (posiBP[popInfo.numLoci]>0){flag_Mb=true;}
        ++posi;
        ++conta;
        if (conta == 4)
        {
            popInfo.numIndividuals = 0;
            while (posi < longi)
            { // asigna genot.
                base1[0] = line.at(posi);
                posi2 = posi;
                posi = int(line.find_first_of(" \t", posi2));
                if (posi < 0)
                {
                    break;
                }
                ++posi;
                base2[0] = line.at(posi);
                rightletter = true;
                str_base = base1[0];
                if (str_plink.find(str_base) == str_plink.npos)
                {
                    rightletter = false;
                }
                str_base = base2[0];
                if (str_plink.find(str_base) == str_plink.npos)
                {
                    rightletter = false;
                }
                if (!rightletter)
                {
                    std::cerr << "Wrong allele in line " << nline << " of tped file (maybe also in other lines). Only bases A, G, C, T and N (case-insensitive) and numbers from 0 to 9 are allowed." << std::endl;
                    exit(EXIT_FAILURE);
                }
                if ((base1[0] != '0') && (base2[0] != '0') && (base1[0] != 'N') && (base2[0] != 'N') && (base1[0] != 'n') && (base2[0] != 'n'))
                {
                    if (base[popInfo.numLoci] == '\0')
                    {
                        base[popInfo.numLoci] = base1[0];
                    }
                    if (base1[0] != base[popInfo.numLoci])
                    {
                        base1[0] = 'X';
                    }
                    if (base2[0] != base[popInfo.numLoci])
                    {
                        base2[0] = 'X';
                    }

                    // 0:homo ref, 1:het, 2:homo noref
                    if (base1[0] == base2[0])
                    {
                        if (base1[0] == base[popInfo.numLoci])
                        {
                            population[popInfo.numIndividuals][popInfo.numLoci] = 0;
                        }
                        else
                        {
                            population[popInfo.numIndividuals][popInfo.numLoci] = 2;
                        }
                    }
                    else
                    {
                        population[popInfo.numIndividuals][popInfo.numLoci] = 1;
                    }
                }
                else
                {
                    population[popInfo.numIndividuals][popInfo.numLoci] = 9; // '9' = Genotipo sin asignar
                }
                posi2 = posi;
                posi = int(line.find_first_of(" \t", posi2));
                if (posi < 0)
                {
                    posi = longi;
                }
                ++posi;
                popInfo.numIndividuals++;
                if (popInfo.numIndividuals > MAXIND)
                {
                    std::cerr << "Reached limit of sample size (" << MAXIND << ")" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }

            if (popInfo.numLoci == 0)
            {
                contaIndBase = popInfo.numIndividuals;
            }

            if (popInfo.numIndividuals != contaIndBase)
            {
                std::cerr << "Some SNP in the sample are not represented in all the individuals" << std::endl;
                exit(EXIT_FAILURE);
            }

            ++popInfo.numLoci;
            if (popInfo.numLoci >= MAXLOCI)
            {
                std::cerr << "Reached max number of loci (" << MAXLOCI << ")" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    if (ncromos > 1)
    {
        flag_chr = true;
    }

    if (flag_cM){
        Mtot = 0;
        for (int conta = 0; conta < ncromos; ++conta) {
            Mtot += posiCM[rangocromo[conta+1]-1]-posiCM[rangocromo[conta]];
        }
        Mtot /= 100.0;       // en Morgans
    }

    if (flag_Mb){
        Mbtot = 0;
        for (int conta = 0; conta < ncromos; ++conta) {
            Mbtot += posiBP[rangocromo[conta+1]-1]-posiBP[rangocromo[conta]];
        }
        Mbtot /= 1000000.0;  // en megabases
        if (flag_r){
            for (int conta = 0; conta < popInfo.numLoci; ++conta) {
                posiCM[conta] = params.r * posiBP[conta] / 1000000;
            }
            Mtot = Mbtot*params.r/100;  // en Morgans
            flag_cM=true;
        }        
    }
        
    if (flag_cM && params.s>1 && !flag_z){flag_z=true;params.z=5;}

//    genomesize = Ncrom * params.Mchr;
    if (flag_chr)
    {
        popInfo.numcromo = ncromos;
        unomenoschrprop = 0;
        if(!flag_Gs){genomesize=Mtot;}
        for (i = 0; i < ncromos; ++i)
        {
            chrnum[i] = rangocromo[i + 1] - rangocromo[i]; // número de marcadores en cada cromosoma
            chrprop[i] = chrnum[i] / contalines;           // Proporción de SNPs en cada cromosoma
            if (flag_Gs){
                chrsize[i] = chrprop[i] * genomesize; // tamaño en Morgans
            }
            else if (flag_cM){
                chrsize[i] = (posiCM[rangocromo[i+1]-1]-posiCM[rangocromo[i]])/100; // tamaño en Morgans
            }
            unomenoschrprop += chrprop[i];
        }
        unomenoschrprop = 1 - unomenoschrprop;
    }

    entrada.close();
}

void readFile_vcf(std::string fichinput1, char (&population)[MAXIND][MAXLOCI], PopulationInfo(&popInfo))
{
    /*
     * Takes as input the file name and a pointer to the population
     * matrix and returns the number of individuals in the file
     */
    char base1[1], base2[1];
    int contaIndBase = 0, nline = 0;
    int conta = 0, posi = 0, posi2 = 0, longi = 0, i, j;
    std::string line;
    std::string cromocod;
    std::string cromocodback = "laksjhbqne";
    std::string str_vcf = "AGCT,agct";
    std::string str_base = "";
    bool rightletter, hayalelos;

    // READING .vcf DATA:
    std::ifstream entrada;
    entrada.open(fichinput1, std::ios::in); // Bucle de lectura del fichero tped
    if (!entrada.good())
    {
        std::cerr << "Could not open \"" << fichinput1 << "\". Does the file exist?" << std::endl;
        exit(EXIT_FAILURE);
    }
    flag_chr = false;
    ncromos = 0;
    int contalines = 0;
    popInfo.numLoci = 0;
    while (std::getline(entrada, line))
    {
        ++nline;
        if (line.at(0) != '#')
        {
            longi = int(line.length());
            if (longi < 12)
            {
                std::cerr << "Line "<<nline<<" too short in vcf file" << std::endl;
                exit(EXIT_FAILURE);
            }
            conta = 0;
            posi = 0;
            // chr name
            posi2 = posi;
            posi = int(line.find_first_of("\t", posi2));
            if (posi < 0)
            {
                std::cerr << "Line "<<nline<<" too short in vcf file" << std::endl;
                exit(EXIT_FAILURE);
            }
            cromocod = line.substr(0, posi);

            ++posi;
            ++conta;


            posi2 = posi;
            posi = static_cast<int>(line.find_first_of("\t", posi2));
            if (posi < 0) {
                std::cerr << "Line too short in vcf file" << std::endl;
                exit(EXIT_FAILURE);
            }
            posiBP[popInfo.numLoci] = std::stoi(line.substr(posi2, posi-posi2));
            if (posiBP[popInfo.numLoci]>0){flag_Mb=true;}
            ++posi;
            ++conta;

            while ((posi < longi) && (conta < 3))
            {
                posi2 = posi;
                posi = int(line.find_first_of("\t", posi2));
                if (posi < 0)
                {
                    std::cerr << "Line "<<nline<<" too short in vcf file" << std::endl;
                    exit(EXIT_FAILURE);
                }
                ++posi;
                ++conta;
            }

            // Mira el alelo REF
            posi2 = posi;
            posi = int(line.find_first_of("\t", posi2));
            j = posi - posi2;
            rightletter = true;
            hayalelos = true;
            for (i = 0; i < j; ++i)
            {
                str_base = line.substr(posi2 + i, 1);
                if (str_base == ".")
                {
                    hayalelos = false;
                }
                else if (str_vcf.find(str_base) == str_vcf.npos)
                {
                    rightletter = false;
                    break;
                }
            }
            if (!rightletter)
            {
                std::cerr << "Wrong reference allele in line " << nline << " of vcf file (maybe also in other lines). Only bases A, G, C and T (case-insensitive) are allowed." << std::endl;
                exit(EXIT_FAILURE);
            }
            ++posi;
            ++conta;

            // Mira el alelo ALT
            posi2 = posi;
            posi = int(line.find_first_of("\t", posi2));
            j = posi - posi2;
            for (i = 0; i < j; ++i)
            {
                str_base = line.substr(posi2 + i, 1);
                if (str_base == ".")
                {
                    hayalelos = false;
                }
                else if (str_vcf.find(str_base) == str_vcf.npos)
                {
                    rightletter = false;
                    break;
                }
            }
        
            if (!rightletter)
            {
                std::cerr << "Wrong alternative allele in line " << nline << " of vcf file. Only bases A, G, C and T (case-insensitive) are allowed." << std::endl;
                exit(EXIT_FAILURE);
            }
            ++posi;
            ++conta;

            while ((posi < longi) && (conta < 8))
            {
                posi2 = posi;
                posi = int(line.find_first_of("\t", posi2));
                if (posi < 0)
                {
                    std::cerr << "Line "<<nline<<" too short in vcf file" << std::endl;
                    exit(EXIT_FAILURE);
                }
                ++posi;
                ++conta;
            }

            // Mira si la posicion 8 es GT (genotipo)
            haygenotipos = false;
            if (line.substr(posi, 2) == "GT")
            {
                haygenotipos = true;
            }

            // Lee los  genotipos:
            if (haygenotipos && hayalelos)
            {
                // contabiliza el cromosoma que habia leido al principio de la linea
                if (cromocod != cromocodback)
                {
                    cromocodback = cromocod;
                    rangocromo[ncromos] = contalines;
                    ++ncromos;
                }
                cromo[contalines] = ncromos;
                ++contalines; // para crom
                // Avanza hasta el primer genotipo:
                posi2 = posi;
                posi = int(line.find_first_of("\t", posi2));
                if (posi < 0)
                {
                    std::cerr << "Line too short in tped file" << std::endl;
                    exit(EXIT_FAILURE);
                }
                ++posi;
                ++conta;
                // empieza a leer genotipos
                popInfo.numIndividuals = 0;
                while (posi < longi)
                { // asigna genot.
                    base1[0] = line.at(posi);
                    posi2 = posi;
                    posi = int(line.find_first_of("/|", posi2));
                    if (posi < 0)
                    {
                        break;
                    }
                    ++posi;
                    base2[0] = line.at(posi);
                    if ((base1[0] != '.') && (base2[0] != '.'))
                    {
                        if (base[popInfo.numLoci] == '\0')
                        {
                            // base[popInfo.numLoci] = '0'; // el de referencia es '0'
                            base[popInfo.numLoci] = base1[0]; // el de referencia es '0'
                        }
                        // if (base1[0] != '0')
                        // {
                        //     base1[0] = '1';
                        // }
                        // if (base2[0] != '0')
                        // {
                        //     base2[0] = '1';
                        // }
                       if (base1[0] != base[popInfo.numLoci])
                        {
                            base1[0] = 'X';
                        }
                        if (base2[0] != base[popInfo.numLoci])
                        {
                            base2[0] = 'X';
                        }

                       // 0:homo ref, 1:het, 2:homo noref
                        if (base1[0] == base2[0])
                        {
                            if (base1[0] == base[popInfo.numLoci])
                            {
                                population[popInfo.numIndividuals][popInfo.numLoci] = 0;
                            }
                            else
                            {
                                population[popInfo.numIndividuals][popInfo.numLoci] = 2;
                            }
                        }
                        else
                        {
                            population[popInfo.numIndividuals][popInfo.numLoci] = 1;
                        }
                    }
                    else
                    {
                        population[popInfo.numIndividuals][popInfo.numLoci] = 9; // '9' = Genotipo sin asignar
                    }
                    posi2 = posi;
                    posi = int(line.find_first_of("\t", posi2));
                    if (posi < 0)
                    {
                        posi = longi;
                    }
                    ++posi;

                    popInfo.numIndividuals++;
                    if (popInfo.numIndividuals > MAXIND)
                    {
                        std::cerr << "Reached limit of sample size (" << MAXIND << ")" << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }

                if (popInfo.numLoci == 0)
                {
                    contaIndBase = popInfo.numIndividuals;
                }

                if (popInfo.numIndividuals != contaIndBase)
                {
                    std::cerr << "Some SNP in the sample are not represented in all the individuals" << std::endl;
                    exit(EXIT_FAILURE);
                }

                ++popInfo.numLoci;
                if (popInfo.numLoci >= MAXLOCI)
                {
                    std::cerr << "Reached max number of loci (" << MAXLOCI << ")" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    if (ncromos > 1)
    {
        flag_chr = true;
    }
//    genomesize = Ncrom * params.Mchr;
    if (flag_Mb){
        Mbtot = 0;
        for (int conta = 0; conta < ncromos; ++conta) {
            Mbtot += posiBP[rangocromo[conta+1]-1]-posiBP[rangocromo[conta]];
        }
        Mbtot /= 1000000.0;  // en megabases
        
        if (flag_r){
            for (int conta = 0; conta < popInfo.numLoci; ++conta) {
                posiCM[conta] = params.r * posiBP[conta] / 1000000;
            }
            Mtot = Mbtot*params.r/100;  // en Morgans
            flag_cM=true;
        }
    }

    if (flag_cM && params.s>1 && !flag_z){flag_z=true;params.z=5;}

    if (flag_chr)
    {
        popInfo.numcromo = ncromos;
        if(!flag_Gs){genomesize=Mtot;}
        unomenoschrprop = 0;
        for (i = 0; i < ncromos; ++i)
        {
            chrnum[i] = rangocromo[i + 1] - rangocromo[i]; // número de marcadores en cada cromosoma
            chrprop[i] = chrnum[i] / contalines;           // Proporción de SNPs en cada cromosoma
            if (flag_Gs){
                chrsize[i] = chrprop[i] * genomesize; // tamaño en Morgans
            }
            else if (flag_cM){
                chrsize[i] = (posiCM[rangocromo[i+1]-1]-posiCM[rangocromo[i]])/100; // tamaño en Morgans
            }
            unomenoschrprop += chrprop[i];
        }
        unomenoschrprop = 1 - unomenoschrprop;
    }

    entrada.close();
}

void printHelp(char *appName)
{
    fprintf(stderr,
            "currentNe2 - Current Ne estimator (v 2.0 - March 2024)\n"
            "Authors: Enrique Santiago - Carlos Köpke\n"
            "\n"
            "USAGE: %s [OPTIONS] <filename_with_extension> [Genome_size_in_Morgans]\n"
            "         where filename is the name of the data file in vcf, ped or tped\n"
            "             format. The filename must include the .vcf, .ped or .tped\n"
            "             extension, depending on the format.\n"
            "         If Genome_size_in_Morgans (~chromosomes number) is specified, \n"
            "         markers are assumed to be evenly distributed in the genome, and:\n"
            "             If chromosome assignments are available in the input \n"
            "               file (or in a accompanying map file when using ped format),\n"
            "               two Ne estimates are calculated: one based on all the\n"
            "               SNP pairs and another based only on SNP pairs located in \n"
            "               different chromosomes. Chromosome sizes are assumed to be\n"
            "               proportional to the number of markers.\n"
            "             If chromosome assignments are not available, chromosomes are \n"
            "               assumed to be one Morgan long and markers are assumed to\n"
            "               be evenly distributed between chromosomes. Only one estimate\n"
            "               based on all SNP pairs is made.\n"
            "         If Genome_size_in_Morgans is not specified, then information on\n"
            "         marker locations within chromosomes is obligatory, and:\n"
            "             If a constant recombination rate in cM/Mb is specified with\n"
            "               the option r, physical locations of markers in bp are \n"
            "               converted to genetic distances in cM and the Ne estimate\n"
            "               using all SNP pairs considers these distances.\n"
            "             If the option r is not used, then the information about marker\n"
            "               locations in a genetic map (i.e., in cM) is used (only  \n"
            "               applicable for ped/map and tped formats).\n\n"
            "OPTIONS:\n"
            "   -h    Print this manual.\n"
            "   -s    Number of SNPs to use in the analysis (all by default).\n"
            "   -k    -If a POSITIVE NUMBER is given, the number of full siblings that\n"
            "         a random individual has IN THE POPULATION (the population is the\n"
            "         set of reproducers). With full lifetime monogamy k=2, with 50%%\n"
            "         of monogamy k=1 and so on. With one litter per multiparous\n"
            "         female k=2, with two litters per female sired by the same father\n"
            "         k=2, but if sired by different fathers k=1, in general, k=2/Le \n"
            "         where Le is the effective number of litters (Santiago et al. 2023).\n"
            "         -If ZERO is specified (i.e., -k 0), each offspring is assumed to\n"
            "         be from a new random pairing.\n"
            "         -If a NEGATIVE NUMBER is specified, the average number of full \n"
            "         siblings observed per individual IN THE SAMPLE. The number k of \n"
            "         full siblings in the population is estimated together with Ne.\n"
            "         -BY DEFAULT, i.e. if the modifier is not used, the average number\n"
            "         of full siblings k is estimated from the input data.\n"
            "   -r    Recombination rate to convert the physical locations in the\n"
            "         input file to genetic locations in cM.\n"
            "   -x    The sample consists of a random set of individuals from\n"
            "         a metapopulation with two subpopulations of equal size.\n"
            "         The k parameter cannot be used when this option is used.\n"
            "         This option is only applicable when the assigments of markers\n"
            "         to chromosomes are available. A single estimate of Ne is generated.\n"
            "   -o    Specifies the output file name. If not specified, the output \n"
            "         filename is built from the name of the input file.\n"
            "   -t    Number of threads (default: %d)\n\n"
            "EXAMPLES: \n"
            "   - Random mating. Full sibs will be identified from the input data.\n"
            "     Genetic distances between markers are obtained from the physical\n"
            "     locations in the input file using a constant rate of 1cM/Mb:\n"
            "         %s -r 1.0 filename\n"
            "   - Random mating and 20 Morgans (equivalent to a genome of 20 chromosomes),\n"
            "     assuming that full siblings are no more common than expected under\n"
            "     random pairing (each offspring from a new random pairing):\n"
            "         %s -k 0 filename 20\n"
            "   - A subsample of 10000 SNPs will be analysed. Individuals are assumed to\n"
            "     come from a metapopulation composed of two subpopulations. The length of\n"
            "     the genome is 30 Morgans and the marker assigments to chromosomes are \n"
            "     available in the input file:\n"
            "         %s -x -s 100000 filename 30\n"
            "   - Full siblings will be identified from the genotyping data. The genetic\n"
            "     locations (in cM) are available in the input file (only for ped/map and\n"
            "     tped formats):\n"
            "         %s filename\n"
            "   - Two full siblings per individual (k = 2) IN THE POPULATION:\n"
            "         %s -k 2 filename 20\n"
            "   - An 80%% of lifetime monogamy in the population. Output filename specified:\n"
            "         %s -k 1.6 -o output_filename filename 20\n"
            "     (with a monogamy rate m = 0.80, the expected number of full\n"
            "     siblings that a random individual has is k = 2*m = 1.6)\n"
            "   - If 0.2 full siblings per individual are OBSERVED IN THE SAMPLE:\n"
            "         %s -k -0.2 filename 20\n"
            "     (NOTE the MINUS SIGN before the number of full sibling 0.2)\n\n",
            appName,
            omp_get_max_threads(),
            appName,
            appName,
            appName,
            appName,
            appName,
            appName,
            appName);
}

int main(int argc, char *argv[])
{
    int i, j, k, j2, j3, conta;
    int containdX, contalocX, mincontaindX;
    double a, b, sx2, sx3, nindmedfs;
    double nDT[5] = {-1.645, -0.6745, 0, 0.6745, 1.645};
    bool superadolimite;
    double aculnc = 0, nlnc = 0, unidades, acusize;

    std::random_device rd;
    std::mt19937 g(rd());
    g.seed(rd());
    std::uniform_real_distribution<> uniforme01(0.0, 1.0);
    for (i = 0; i < MAXIND; ++i)
    {
        validind[i] = true;
    }
    params.K = 0;
    std::string fichspecified = "";
    for (;;)
    {
        switch (getopt(argc, argv, "hxs:k:r:z:m:o:t:qpv"))
        {
        case '?':
        case 'h':
        default:
            printHelp(argv[0]);
            return -1;
        case 'x':
            params.s=2;
            continue;
        case 's':
            params.numSNPs = std::atoi(optarg);
            if (params.numSNPs < 10)
            {
                std::cerr << "Invalid number of SNPs" << std::endl;
                return -1;
            }
            continue;
        case 'k':
            params.K = std::atof(optarg);
            params.flagnok = false;
            if (params.K < 0)
            {
                params.ks = -(params.K);
                params.flagks = true;
            }
            else
            {
                params.flagks = false;
            }
            continue;
        case 'r':
            params.r = std::atof(optarg);
            if (params.r <= 0)
            {
                std::cerr << "Invalid recombination rate" << std::endl;
                return -1;
            }      
            flag_r = true;
            continue;
        case 'z':
            params.z = std::atof(optarg);
            if (params.z < 0)
            {
                std::cerr << "Invalid lowest cut-off genetic distance in cM" << std::endl;
                return -1;
            }      
            flag_z = true;
            continue;
        case 'm':
            params.miss = std::atof(optarg);
            continue;
        case 'o':
            fichspecified = optarg;
            continue;
        case 't':
            params.numThreads = std::atoi(optarg);
            continue;
        case 'q':
            params.quiet = true;
            continue;
        case 'v':
            params.verbose = true;
            continue;
        case 'p':
            params.printToStdOut = true;
            continue;
        case -1:
            break;
        }
        break;
    }
    std::string fich = "";
    std::string extension = "";
    std::string prefijo = "";
    std::string fichinput1 = "";
    std::string fichinput2 = "";
    Ncrom = 0;
    if (optind < argc)
    {
        fich = argv[optind];

        optind++;
        if (optind < argc)
        {
            Ncrom = std::atof(argv[optind]);
            ncromos=int(Ncrom);
            if (Ncrom < 1){
                std::cerr << "Genetic size must be larger than 1 Morgan" << std::endl;
                return -1;
            }
            flag_Gs = true;
            if (flag_r){
                std::cout << "Option r is not compatible with the assignation of a genome size in the command line." << std::endl;
                std::cout << "The genome size given in the command line will be ignored." << std::endl;
                flag_Gs = false;
                Ncrom =ncromos= 0;
            }
        }
    }
    if (flag_z){
        if (flag_Gs){
            std::cerr << "Option z is not compatible with the assignation of a genome size in the command line." << std::endl;
            return -1;
        }
    }

    if (!params.flagnok && params.s>1){
        std::cerr << "Option x is incompatible with option k." << std::endl;
        return -1;
    }

    genomesize = Ncrom * params.Mchr;

    if (fich == "")
    {
        std::cerr << "Missing data file name" << std::endl;
        return -1;
    }

    fichinput1 = fich;
    std::size_t punto = fich.find(".", 0);
    if (punto == std::string::npos)
    {
        std::cerr << "The filename has no extension" << std::endl;
        return -1;
    }
    std::size_t punto2 = 0;
    while (punto2 != std::string::npos)
    {
        punto2 = fich.find(".", punto + 1);
        if (punto2 != std::string::npos)
        {
            punto = punto2;
        }
        else
        {
            break;
        }
    }
    prefijo.assign(fich, 0, punto);
    extension.assign(fich, punto + 1, fich.length() - punto - 1);
    if (prefijo == "")
    {
        std::cerr << "Missing data filename" << std::endl;
        return -1;
    }
    if ((extension != "vcf") && (extension != "ped") && (extension != "tped"))
    {
        std::cerr << "Missing data filename extension. It must be vcf, ped or tped" << std::endl;
        return -1;
    }

    if (extension == "ped")
    {
        fichinput2 = prefijo + ".map";
    }

    std::string fichProgress = fich + "_currentNe2_progress.tmp";
    params.progress.InitTotalTasks(1, fichProgress.c_str());

    n_sample = params.numSample;
    n_SNPs = params.numSNPs;
    n_threads = params.numThreads;

    std::ifstream entrada;

    if (n_threads > 0)
    {
        omp_set_num_threads(n_threads);
    }

    if (!params.quiet)
    {
        std::cout << " A progress report is stored at " << fichProgress << "\n";
        std::cout << " Check it using 'cat " << fichProgress << "'\n";
        std::cout << " Reading file " << fich << std::endl;
    }

    tini = omp_get_wtime();
    // LECTURA DE LOS DATOS:
    struct PopulationInfo popInfo =
        {
            .numIndividuals = 0,
            .numLoci = 0};

    if (extension == "ped")
    {
        fichinput1 = prefijo + ".ped";
        fichinput2 = prefijo + ".map";
        readFile_ped(fichinput1, fichinput2, indi, popInfo);
    }
    else if (extension == "tped")
    {
        fichinput1 = prefijo + ".tped";
        readFile_tped(fichinput1, indi, popInfo);
    }
    else if (extension == "vcf")
    {
        fichinput1 = prefijo + ".vcf";
        readFile_vcf(fichinput1, indi, popInfo);
    }

    double tProcessFile = (omp_get_wtime() - tini);
    if (!params.quiet)
    {
        if ((flag_chr) && (extension == "ped"))
        {
            std::cout << " Reading " << fichinput1 << " and " << fichinput2 << " took " << std::fixed << std::setprecision(2) << tProcessFile << " sec" << std::endl;
        }
        else
        {
            std::cout << " Reading " << fichinput1 << " took " << std::fixed << std::setprecision(2) << tProcessFile << " sec" << std::endl;
        }
    }
    if (!params.quiet)
        std::cout << " Processing file" << std::endl;

    for (j = 0; j < popInfo.numIndividuals; ++j)
    {
        pind[j] = j;
    }
    pfinind = &pind[0] + popInfo.numIndividuals - 1; // último elemento de la matriz
    prefind = pfinind;                               // apunta al último elemento de la matriz para empezar

    for (j = 0; j < popInfo.numLoci; ++j)
    { // Numera los loci en orden consecutivo para los cromosomas elegidos
        ploc[j] = j;
    }

    contaloc2 = popInfo.numLoci;
    std::shuffle(&ploc[0], &ploc[popInfo.numLoci], g); // Aleatoriza los loci
    prefloc = &ploc[0];
    pfinloc = prefloc + popInfo.numLoci;

    //  LOS INDIVIDUOS NO SE ALEATORIZAN
    //    std::shuffle(&pind[0],&pind[popInfo.numIndividuals],g); // Aleatoriza el orden de los individuos

    prefind = &pind[0];
    pfinind = prefind + popInfo.numIndividuals;

    if (n_SNPs > 0)
    {
        enelocsolicitado = eneloc = n_SNPs;
    }
    else
    {
        enelocsolicitado = eneloc = popInfo.numLoci;
    }

    if (n_sample > 0)
    {
        eneind = n_sample;
    }
    else
    {
        eneind = popInfo.numIndividuals;
    }
    if (eneloc > popInfo.numLoci)
    {
        eneloc = popInfo.numLoci;
    }
    if (eneind > popInfo.numIndividuals)
    {
        eneind = popInfo.numIndividuals;
    }

    // Elimina los individuos indicados en el comando
    for (i = eneind - 1; i >= 0; --i)
    {
        if (!validind[i])
        {
            for (j = i; j < eneind; ++j)
            {
                pind[j] = pind[j + 1];
            }
            --eneind;
            ++contaelim;
        }
    }

    // Calculo de frecuencia del alelo noref y de homo noref para todos los loci variables en la muestra
    contaseg = 0;
    containdX = 0;
    pj = prefloc;
    frecmed = 0;
    nfrecmed = 0;
    nindmedfs = 0;
    mincontaindX = int(float(eneind) * 2.0 * (1 - params.miss)); // se admite por defecto un 20% de missing data (params.miss)
    for (j = 0; j < popInfo.numLoci; ++j)
    {
        frec[*pj] = 0;
        segrega[*pj] = false;
        homo[*pj] = 0;
        pk = prefind;
        containdX = 0;
        for (i = 0; i < eneind; i++)
        {
            ff = indi[*pk][*pj];
            if (ff < 9)
            {
                frec[*pj] += ff; // acumulador de frecuencia del alelo noref
                if (ff == 2)
                {
                    ++homo[*pj];
                } // acumulador de homo noref
                ++(++containdX);
            }
            ++pk;
        }
        if (containdX > mincontaindX)
        {
            if ((frec[*pj] > 0) && (frec[*pj] < containdX))
            { // contador de segregantes
                ++contaseg;
                nindmedfs += containdX;
                segrega[*pj] = true;
            }
            frec[*pj] /= (containdX);
            homo[*pj] /= (containdX / 2);
            b = frec[*pj];
            if (b > 0.5)
            {
                b = 1 - b;
            }
            frecmed += b;
            ++nfrecmed;
        }
        ++pj;
        if (contaseg >= eneloc)
        {
            break;
        }
    }

    nindmedfs /= 2 * contaseg;
    frecmed /= nfrecmed;
    if (eneind > popInfo.numIndividuals)
    {
        eneind = popInfo.numIndividuals;
    }
    if (eneloc > contaseg)
    {
        eneloc = contaseg;
    }

    // Calculo de f sample y control del exceso del límite de loci:
    superadolimite = false;
    acupq = 0;
    acuPp2 = 0;
    pj = prefloc;
    conta = 0;
    for (j2 = 0; j2 < eneloc; ++j2)
    {
        for (;;)
        {
            if (segrega[*pj])
            {
                acuPp2 += (homo[*pj] - frec[*pj] * frec[*pj]);
                acupq += (frec[*pj] * (1.0 - frec[*pj]));
                ++pj;
                ++conta;
                break;
            }
            else
            {
                ++pj;
                if (pj > pfinloc)
                {
                    superadolimite = true;
                    break;
                }
            }
        }
        if (superadolimite)
        {
            break;
        }
    }
    eneloc = conta; // Loci segregantes en la muestra
    fs = acuPp2 / acupq;
    Het_esp = acupq / eneloc;

    increMorgans = 0.0001;
    // DISTRIBUTION OF MAP DISTANCES BETWEEN MARKERS ACROSS CHROMOSOMES
    if (flag_chr){
        maxdistance=0;
        maxdistanceindx=0;
        sumadist=0;
        if (flag_Gs){
            for (i = 0; i < ncromos; ++i){
                if (maxdistance<chrsize[i]){
                    maxdistance=chrsize[i];
                }
                acusize = chrsize[i];
                j=int(acusize/increMorgans);
                if (maxdistanceindx<j){
                    maxdistanceindx=j;
                }
                k=1;
                while ((j>=0)){
                    mapdist[j]+=k;
                    sumadist+=k;
                    ++k;
                    --j;
                }
            }
        }
        else if(flag_cM){
            int d=0;
            increMorgans*=100;
            for (i = 0; i < ncromos; ++i){
                for (j=rangocromo[i];j<rangocromo[i+1]-1;++j){
                    for (k = j+1; k < rangocromo[i+1]; ++k){
                        d=int(abs(posiCM[k]-posiCM[j])/increMorgans); // distancia en cM*0.01
                        if (d>maxdistanceindx){maxdistanceindx=d;}
                        ++mapdist[d];
                        ++sumadist;
                    }
                }
                maxdistance = std::max(float(maxdistanceindx)*increMorgans/100.0,maxdistance);
            }
            increMorgans/=100;
        }

        for (j = 0; j <maxdistanceindx; ++j){
            mapdist[j]/=sumadist;
        }
    }
    else{
        sumadist=0;
        // maxdistance = float(Ncrom)/float(ncromos);
        maxdistance = 1;
        maxdistanceindx=int(maxdistance/increMorgans);
        // j=int(((float(eneloc)/Ncrom)/2.0)/increMorgans);
        j=0;
        k=1;
        while ((j<maxdistanceindx)){
            mapdist[j]+=k;
            sumadist+=k;
            ++k;
            ++j;
        }
        for (j = 0; j <maxdistanceindx; ++j){
            mapdist[j]/=sumadist;
        }
    }
    //    std::cout << cmed << std::endl;
    if (maxdistance>10){
            std::cerr << "Chromosome sizes are limited to a maximum of 1000 cM. There is at least" << std::endl;
            std::cerr << "one chromosome larger than this limit." << std::endl;
            return -1;
    }



    fp = (1.0 + fs * (2 * nindmedfs - 1.0)) / (2 * nindmedfs - 1.0 + fs);
    if (std::abs(fp) > 0.08){
        flagnoestimaks = true;
    }
    std::stringstream salidahermanos;
    if (params.flagnok && !flagnoestimaks){
        salidahermanos << std::fixed << std::setprecision(0);
        // En busca hermanos completos usando toda la información disponible
        // SIN ALEATORIZAR:
        a = 4 / (6 - frecmed * (1 - frecmed)) - 0.028;
        // std::cout<<"HERMANOS METODO NUEVO:\n";
        // std::cout<<frecmed<<" ,"<<a<<"\n";
        nparhermanos = 0;
        npadrehijo = 0;
        for (i = 0; i < popInfo.numIndividuals - 1; ++i)
        {
            for (j = i + 1; j < popInfo.numIndividuals; ++j)
            {
                counthethet = 0;
                counthethomo = 0;
                counthomohomo = 0;
                countnopadre = 0;
                for (j2 = 0; j2 < popInfo.numLoci; ++j2)
                {
                    sumahethomo = indi[i][j2] + indi[j][j2];
                    if ((sumahethomo == 1) || (sumahethomo == 3))
                    { // HetHomo
                        ++counthethomo;
                    }
                    else if (sumahethomo == 2)
                    {
                        if (indi[i][j2] == indi[j][j2])
                        {
                            ++counthethet;
                        }
                        else
                        {
                            ++counthomohomo;
                            ++countnopadre;
                        }
                    }
                }
                if (counthethet > 10)
                { // primero mira a ver si hay heterocigotos suficientes para hacer las pruebas
                    rationopadre = float(countnopadre) / (float(counthethet + counthethomo + counthomohomo));
                    b = float(counthethet + counthethomo + counthomohomo);
                    ratiohets = 0;
                    if (b > 0)
                    {
                        ratiohets = (float(counthethet) + float(counthethomo) / 2.0) / b;
                    }
                    if ((rationopadre < 0.001) && (ratiohets > 0.64))
                    { // Es una pareja padre-hijo (error de genotipado del 0.001)
                        // padrehijo[0][npadrehijo] = i;
                        // padrehijo[1][npadrehijo] = j;
                        ++npadrehijo;
                    }
                    else
                    { // si no es una pareja padre-hijo
                        if (ratiohets > a)
                        { // mira a ver si son hermanos completos
                            salidahermanos << i + 1 << ", " << j + 1 << "\n";
                            // hermanos[0][nparhermanos] = i;
                            // hermanos[1][nparhermanos] = j;
                            ++nparhermanos;
                            // std::cout<<i+1<<" ,"<<j+1<<"\n";
                        }
                    }
                }
            }
        }
    }
    if (params.flagnok)
    {
        if (flagnoestimaks)
        {
            params.ks = 0;
        }
        else
        {
            params.ks = float(nparhermanos) * 2.0 / float(popInfo.numIndividuals);
        }
    }

    // Distribucion de heterocigosidades de los individuos (ver Medidas_de_Parent_y_Het.docx)
    pk = prefind;
    for (i = 0; i < eneind; ++i)
    {
        Het[*pk] = 0;
        pj = prefloc;
        contalocX = 0;
        for (j2 = 0; j2 < eneloc; ++j2)
        {
            for (;;)
            {
                if (segrega[*pj])
                {
                    ff = indi[*pk][*pj];
                    if (ff == 1)
                    {
                        ++Het[*pk];
                        ++contalocX;
                    } // acumulador de heterocigosidad de cada indiv
                    else if (ff < 9)
                    {
                        ++contalocX;
                    }
                    ++pj;
                    break;
                }
                else
                {
                    ++pj;
                }
            }
        }
        if (contalocX == 0)
        {
            std::cerr << "There is no genotyping data for at least one individual" << std::endl;
            return -1;
        }
        Het[*pk] /= contalocX; // het de cada individuo
        ++pk;
    }
    acuHet = 0;
    pk = prefind;
    for (i = 0; i < eneind; ++i)
    {
        acuHet += Het[*pk];
        ++pk;
    }
    acuHet /= eneind;
    Het_med = acuHet; // het media de todos los loci y todos los individuos
    pk = prefind;
    sx2 = sx3 = 0;
    pk = prefind;
    for (i = 0; i < eneind; ++i)
    {
        a = (Het[*pk] / acuHet - 1);
        sx2 += a * a;
        sx3 += a * a * a;
        ++pk;
    }
    Het_var = sx2 / (eneind - 1);
    Het_sesg = sx3 / eneind;
    Het_DT = sqrt(Het_var);
    Het_sesg /= (Het_DT * Het_DT * Het_DT);

    // Calculo de D2:
    params.progress.InitCurrentTask(eneloc - 1);
    params.progress.SetCurrentTask(0, "Measuring d²");
    params.progress.SaveProgress();
    acuD2 = acuW = acur2 = acun = 0;
    acuD205 = acuW05 = acur205 = acun05 = 0;
    pj = prefloc;
    int valid_idx[contaloc2] = {0};
    int counter = 0;
    // Inicializamos la tabla de índices válidos (en los que secrega[x] es true)
    for (int idx = 0; idx < contaloc2; idx++)
    {
        if (segrega[ploc[idx]])
        {
            valid_idx[counter] = ploc[idx];
            counter++;
        }
    }

    for (j2 = 0; j2 < eneloc - 1; ++j2)
    {
        ppj = &valid_idx[j2];
#pragma omp parallel for private(tacui,tacuj,tacuHoHo, distancia, id, tacuHoHetHetHo, tacuHetHet, ppk, ppi, D, W, r2, ss, _containdX)
        for (j3 = j2 + 1; j3 < eneloc; ++j3)
        {
            ppi = &valid_idx[j3];
            tacuHoHo = tacuHoHetHetHo = tacuHetHet = 0; // acumuladores de genotipos
            tacui=tacuj=0;
            ppk = prefind;
            _containdX = 0;
            for (int _i = 0; _i < eneind; ++_i)
            {
                ss = indi[*ppk][*ppj] + indi[*ppk][*ppi];
                if (ss < 9)
                {
                    tacui+=indi[*ppk][*ppi];
                    tacuj+=indi[*ppk][*ppj];
                    ++_containdX;
                    if (ss < 2)
                    {
                    }
                    else if (ss == 2)
                    {
                        if (indi[*ppk][*ppj] == indi[*ppk][*ppi])
                        {
                            ++tacuHetHet;
                        }
                    }
                    else if (ss == 3)
                    {
                        ++tacuHoHetHetHo;
                    }
                    else if (ss == 4)
                    {
                        ++tacuHoHo;
                    }
                }
                ++ppk;
            }
            if (_containdX > 1)
            {
                tacui/=(_containdX*2);
                tacuj/=(_containdX*2);
                W = tacui*tacuj;
                D = -2 * W + (2 * tacuHoHo + tacuHoHetHetHo + tacuHetHet / 2) / _containdX;
                D *= D;
                W *= (1 - tacui) * (1 - tacuj);
                if (flag_chr)
                {
                    if ((cromo[*ppi] != cromo[*ppj]))
                    {
                        ++x_contapares05[j3];
                        x_containdX05[j3] += _containdX;
                        xD05[j3] += D;
                        xW05[j3] += W;
                    }
                    else
                    {
                        if (flag_cM){
                            distancia=fabs(posiCM[*ppi] - posiCM[*ppj]);
                            if (distancia>params.z){
                                ++x_contapareslink[j3];
                                x_containdXlink[j3] += _containdX;
                                xDlink[j3] += D;
                                xWlink[j3] += W;
                            }
                        }
                        else{
                            ++x_contapareslink[j3];
                            x_containdXlink[j3] += _containdX;
                            xDlink[j3] += D;
                            xWlink[j3] += W;
                            // xr2link[j3] += D / W; <<<<<XXXXXXX
                        }
                    }
                }
                else
                {
                    ++x_contapares[j3];
                    x_containdX[j3] += _containdX;
                    xD[j3] += D;
                    xW[j3] += W;
                    // xr2[j3] += D / W; <<<<<XXXXXXX
                }
            }


        }
        if (j2 % 1000 == 0)
        {
            params.progress.SetTaskProgress(j2 + 1);
            // params.progress.PrintProgress();
        }
    }

    if (flag_chr)
    {
        for (j3 = 0; j3 < eneloc; ++j3)
        {
            acun05 += x_contapares05[j3];
            effndata05 += x_containdX05[j3];
            acuD205 += xD05[j3];
            acuW05 += xW05[j3];
            // acur205 += xr205[j3]; <<<<<XXXXXXX
        }
        for (j3 = 0; j3 < eneloc; ++j3)
        {
            acunlink += x_contapareslink[j3];
            effndatalink += x_containdXlink[j3];
            acuD2link += xDlink[j3];
            acuWlink += xWlink[j3];
            // acur2link += xr2link[j3]; <<<<<XXXXXXX
        }
        acun = acun05 + acunlink;
        effndata = effndata05 + effndatalink;
        d2s = (acuD205 + acuD2link) / (acuW05 + acuWlink);
        acuD2 = (acuD205 + acuD2link) / acun;
        acuW = (acuW05 + acuWlink) / acun;
        // acur2 = (acur205 + acur2link) / acun; <<<<<XXXXXXX

        d2s05 = acuD205 / acuW05;
        acuD205 /= acun05;
        acuW05 /= acun05;
        // acur205 /= acun05; <<<<<XXXXXXX

        d2slink = acuD2link / acuWlink;
        acuD2link /= acunlink;
        acuWlink /= acunlink;
        // acur2link /= acunlink; <<<<<XXXXXXX
    }
    else
    {
        for (j3 = 0; j3 < eneloc; ++j3)
        {
            acun += x_contapares[j3];
            effndata += x_containdX[j3];
            acuD2 += xD[j3];
            acuW += xW[j3];
            // acur2 += xr2[j3]; <<<<<XXXXXXX
        }
        d2s = acuD2 / acuW;
        acuD2 /= acun;
        acuW /= acun;
        // acur2 /= acun; <<<<<XXXXXXX
    }

    obsndata = double(eneind) * (double(eneloc) * (double(eneloc) - 1.0)) / 2.0;
    // propmiss = 1.0 - effndata / obsndata; <<<<<XXXXXXX
    n_SNP_pairs = (double(eneloc) * double(eneloc - 1)) / 2.0;
    effeneind = effndata / acun;
    propmiss = 1.0 - effeneind / eneind; // <<<<<XXXXXXX
    Wt=acuW/(1-2/(2*effeneind+1));

    std::stringstream salida0;
    std::stringstream salida;
    salida0 << "# (currentNe v 2.0)\n";
    salida0 << "# Command:";
    for (i = 0; i < argc; ++i)
    {
        salida0 << " " << argv[i];
    }
    salida0 << "\n";
    salida0 << "# Running time:";

    salida << "# INPUT PARAMETERS:\n";
    salida << "# Number of chromosomes in the .map file:\n";
    if (flag_chr)
    {
        salida << std::fixed << std::setprecision(2);
        salida << ncromos << "\n";
    }
    else
    {
        salida << "No chromosome information:\n";
    }
    salida << "# Genome size in Morgans:\n";
    salida << std::fixed << std::setprecision(2);
    salida << genomesize << "\n";
    salida << "# Total number of individuals in the input file:\n";
    salida << std::fixed << std::setprecision(0);
    salida << popInfo.numIndividuals << "\n";
    salida << std::fixed << std::setprecision(2);
    salida << "# Effective Number of individuals included in the analysis (excluding missing genotypes):\n";
    salida << effeneind << "\n";
    salida << std::fixed << std::setprecision(0);
    salida << "# Number of SNPs in the input file:\n";
    salida << popInfo.numLoci << "\n";
    salida << "# Number of SNPs included in the analysis (only polymorphic and with less than 20% missing data):\n";
    salida << eneloc << "\n";
    salida << "# Number of SNP pairs included in the analysis:\n";
    salida << n_SNP_pairs << "\n";
    salida << std::fixed << std::setprecision(8);
    salida << "# Proportion of missing data:\n";
    salida << propmiss << "\n";
    salida << std::fixed << std::setprecision(2);
    salida << "# Number of full siblings that a random individual in the population has:\n";
    if (!params.flagks && !params.flagnok){
        salida << params.K << "\n";
    }
    else{
        salida << "Not given\n";
    }
    salida << "# Number of full siblings that a random individuals in the sample has:\n";
    if (params.flagks && !params.flagnok){
        salida << params.ks << "\n#\n";
    }
    else
    {
        salida << "Not given\n#\n";
    }
    salida << std::fixed << std::setprecision(8);
    salida << "# OUTPUT PARAMETERS:\n";
    salida << "# Observed d^2 of the entire sample (weighted correlation of loci pairs):\n";
    salida << d2s << "\n";
    d2_pob = (d2s - (4 * double(effeneind) - 4) / ((2 * double(effeneind) - 1) * (2 * double(effeneind) - 1))) / ((1 - 1 / (2 * double(effeneind))) * 0.25); // APROXIMADO
    if (flag_chr)
    {
        salida << "# Observed d^2 of the sample (only between different chromosomes):\n";
        salida << d2s05 << "\n";
        salida << "# Observed d^2 of the sample (within chromosomes):\n";
        salida << d2slink << "\n";
    }
    salida << "# Expected heterozygosity of the individuals in the sample under H-W eq.:\n";
    salida << 2 * Het_esp << "\n";
    salida << "# Observed heterozygosity of the individuals in the sample:\n";
    salida << Het_med << "\n";
    salida << "# Observed Fis value in the sample (excess of homozygotes):\n";
    salida << fs << "\n";
    salida << "# Estimate of the Fis value in the population (excess of homozygotes):\n";
    salida << fp << "\n";

    params.progress.SetCurrentTask(1, "Analyzing input data");

    if (params.s>1)
    { // Aquí va la población subdividida

        if ((flag_chr) && (d2slink>d2s05))
        {
            // first, search for the maximum possible value of Fst
            double maxFst=1;
            double minFst=0;
            double difmax=0;
            double difmin=0;

            double fp12 = (1 + fp) * (1 + fp);
            double sample_size_h = effeneind * 2;
            double samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
            samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
            samplex += 4 * (sample_size_h - 2.0);
            samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
            double sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));
            double d2= (d2s05 - sampley * fp12) / samplex;
            if (d2>0){
                Fst = std::sqrt(d2);
                Ne05=Nelink=10E9;
                double Ne05ant=10E9;
                double Nelinkant=10E9;
                increFst = Fst;
                for (i=0;i<5;++i){
                    increFst /= 10;
                    while ((Ne05-Nelink)<=0){
                        if (Fst<increFst){break;}
                        Ne05ant=Ne05;
                        Nelinkant=Nelink;
                        Fst-=increFst;
                        MixEcuacion05();
                        Ne05=Ne;
                        MixIntegralUNCROM();
                        Nelink=Ne;
                        if ((Ne05-Nelink)>0){
                            Fst+=increFst;
                            Ne05=Ne05ant;
                            Nelink=Nelinkant;
                            break;
                        }
                    }
                }
                Fst-=increFst;
                MixEcuacion05();
                Ne05=Ne;
                MixIntegralUNCROM();
                Nelink=Ne;
                difmax=Ne05-Nelink;
                maxFst=Fst;
                // then, search for the minimum possible value of Fst
                Fst = 0;
                minFst=Fst;
                MixEcuacion05();
                Ne05=Ne;
                MixIntegralUNCROM();
                Nelink=Ne;
                difmin=Ne05-Nelink;
            }
            else{
                difmin=difmax=1;
            }
            if ((difmin*difmax)<0){
                // Next, search for the confluence of the two estimates of Ne
                int signo=-1;
                increFst =(maxFst-minFst)/10;
                Fst=maxFst;
                MixEcuacion05();
                Ne05=Ne;
                MixIntegralUNCROM();
                Nelink=Ne;
                double dif=Ne05-Nelink;
                int countchanges=0;
                while (countchanges<10){
                    // Check for valid Fst
                    while ((dif*signo<0)){
                        Fst +=increFst*signo; // reset the previous value
                        if (Fst>maxFst){Fst=maxFst;}
                        if (Fst<minFst){Fst=minFst;}
                        MixEcuacion05();
                        Ne05=Ne;
                        MixIntegralUNCROM();
                        Nelink=Ne;
                        dif=Ne05-Nelink;
                    }
                    signo=-signo;
                    ++countchanges;
                    increFst/=10;
                }
                salida << "#\n# Partition of the F-statistics:\n";
                salida << "# Fit (probability of IBD within individuals relatively to the entire metapopulation):\n";
                salida << std::fixed << std::setprecision(6);
                salida << fp << "\n";
                salida << "# Fst (subpopulation differentiation):\n";
                salida << Fst << "\n";
                salida << "# Fis (probability of IBD within individuals in each subpopulation):\n";
                salida << 1-(1-fp)/(1-Fst) << "\n";
                salida << "#\n# Migration rate between subpopulations of equal size:\n";
                salida << std::fixed << std::setprecision(6);
                salida << m << "\n";
                salida << "#\n# Ne of the entire metapopulation:\n";
                salida << std::fixed << std::setprecision(2);
                salida << Nelink * params.s << "\n";
                ncrom_sample_int = int(Ncrom);
                if (ncrom_sample_int < 1)
                {
                    ncrom_sample_int = 1;
                }
                ncrom_sample = Ncrom;
                log10_ncrom_sample = log10(genomesize);
                if (log10_ncrom_sample > log10(60))
                {
                    log10_ncrom_sample = log10(60);
                } // 60 cromosomas es c=0.5
                f_pob = fp;
                log10_n_sample = log10(effeneind);
                log10_nSNPs_solicitado = log10(eneloc);
                Ne_integral_totcrom = Nelink * params.s ;
                log10_Ne_obs = log10(Ne_integral_totcrom);
                for (i = 0; i < 5; ++i)
                {
                    lim = -nDT[i];
                    if (i != 2)
                    {
                        rangos_int[i] = CalculaIntervalo_soloLD();
                    }
                }
                rangos_int[2] = log10_Ne_obs;
                for (i = 0; i < 5; ++i)
                {
                    rangos_int[i] = pow(10, rangos_int[i]);
                }
                salida << "# Lower limit of 50% CI:\n";
                salida << rangos_int[1] << "\n";
                salida << "# Upper limit of 50% CI:\n";
                salida << rangos_int[3] << "\n";
                salida << "# Lower limit of 90% CI:\n";
                salida << rangos_int[0] << "\n";
                salida << "# Upper limit of 90% CI:\n";
                salida << rangos_int[4] << "\n";
                salida << "#\n";
                salida << "# Partition of the expected LD (weighted d² correlation coefficient) of a metapopulation\n";
                salida << std::fixed << std::setprecision(0);
                salida << "# of "<< 2*Nelink<<" individuals divided into two subpopulations of equal size and recirocal\n";
                salida << std::fixed << std::setprecision(8);
                salida << "# migration rate of "<< m<<" :\n";
                salida << "# Rec_rate\td²_wi_subpops\td²_bw_subpops\t2_x_d²_wi_x_bw\td²_entire_metapop:\n";
                // m12=(1-m)*(1-m);
                // m22=(1-2*m)*(1-2*m);
                ms=params.s*m/(params.s-1);
                m12=(1-ms)*(1-ms);
                c=1;
                double increc=0.1;
                for (i=0;i<4;++i){
                    for (j=1;j<10;++j){
                        c-=increc;
                        if (c<0.51){

                            c2 = c * c;
                            c12 = (1 - c) * (1 - c);

                            // AA=(1-m12)*m12*(1-1/(2*Nelink));
                            // BB=2*(1-m12)*(1-2.2/(2*params.s*Nelink))*(1-c);
                            // CC=1-(1-2.2/(2*params.s*Nelink))*c12;
                            // DD=1-m12*(1-1/(2*Nelink))*(1-c);
                            // EE=(1-m12)*(1-m12)*(1-c)*(1-2.2/(2*params.s*Nelink));

                            // Db2=Wt*Fst*Fst;
                            // DbDw=Db2*(AA/DD+(EE/CC+(AA*BB)/(CC*DD))/(2*params.s*Nelink*DD))+Wt*(1-Fst)*(1-Fst)/(4*params.s*params.s*Nelink*Nelink*CC*DD);
                            // Dw2=Db2*(EE/CC+(AA*BB)/(CC*DD))+Wt*(1-Fst)*(1-Fst)*(1+c2)/(2*params.s*Nelink*CC);

                            Db2=(Wt)*Fst*Fst;
                            Dw2=(Wt)*(1-Fst)*(1-Fst)*(1+c2)/(4*Nelink*(1-c12)+2.2*c12);
                            DbDw=(4*(Wt)*Fst*Fst*m12*(m)+(Dw2)/(4*Nelink))/(1-m12*(1-c)*(1-1/2/Nelink));

                            salida << std::fixed << std::setprecision(4);
                            salida << c<<"\t";
                            salida << std::fixed << std::setprecision(8);
                            salida <<Dw2/Wt<<"\t"<<Db2/Wt<<"\t"<<2*DbDw/Wt<<"\t"<<(Dw2+Db2+2*DbDw)/Wt<<"\n";

                        }
                    }
                    increc/=10;
                }
                salida << "# Multiply by the product of genic variances Wt to obtain squared covariances. Wt:\n";
                salida << Wt<<"\n#\n";            
            }
            else{
                if (difmin+difmax == 2){
                    salida << "#\n# The Fst value is too low to be estimated, indicating no signs of population subdivision.\n";
                    // salida << Ne05<<"  "<<Nelink<<"\n";
                    salida << "#\n";  
                }
                else{
                    salida << "#\n# Ne estimate does not converge to a solution within the range of Fst from 0 to 1.\n";
                    // salida << Ne05<<"  "<<Nelink<<"\n";
                    salida << "#\n";  
                }
            }
        }
        else{
            if (!flag_chr){
                salida << "# Ne cannot be estimated because the markers are not assigned to chromosomes.\n";
                salida << "#\n";
            }
            if (d2slink<d2s05){
                salida << "# Ne cannot be estimated because the LD is smaller for linked than for unlinked markers.\n";
                salida << "#\n";
            }
        }
    }
    else
    { // y aquí la población panmíctica
    
        ncrom_sample_int = int(Ncrom);
        if (ncrom_sample_int < 1)
        {
            ncrom_sample_int = 1;
        }
        ncrom_sample = Ncrom;
        log10_ncrom_sample = log10(genomesize);
        if (log10_ncrom_sample > log10(60))
        {
            log10_ncrom_sample = log10(60);
        } // 60 cromosomas es c=0.5
        f_pob = fp;
        log10_n_sample = log10(effeneind);
        log10_nSNPs_solicitado = log10(eneloc);

        backfp = fp;
        for (fciclo = 0; fciclo < 2; ++fciclo)
        { // The second try assumes fp=0
            Neant = 0;
            if (params.flagks || params.flagnok)
            {
                params.K = 0;
                for (i = 0; i < 30; ++i)
                {
                    IntegralVARIOSCROM();
                    if (params.ks == 0)
                    {
                        break;
                    }
                    params.K = Ne / (eneind - 1) * params.ks;
                    if ((std::abs(Ne - Neant) < 0.1) || (std::abs(Ne - Neant) < (Ne / 10000)))
                    {
                        break;
                    }
                    Neant = Ne;
                }
            }
            else
            {
                IntegralVARIOSCROM();
            }
            if (Ne < 100000000)
            {
                break;
            } // No more than 10^9 Ne
            fp = 0;
        }
        fp = backfp;
        if (Ne < 100000000){
            Ne_integral_totcrom = Ne;
            log10_Ne_obs = log10(Ne_integral_totcrom);
            for (i = 0; i < 5; ++i)
            {
                lim = -nDT[i];
                if (i != 2)
                {
                    rangos_int[i] = CalculaIntervalo_soloLD();
                }
            }
            rangos_int[2] = log10_Ne_obs;
            for (i = 0; i < 5; ++i)
            {
                rangos_int[i] = pow(10, rangos_int[i]);
            }

            salida << "#\n"
                << "# Ne estimation by integration over the whole genome (no genetic map available).\n";
            salida << std::fixed << std::setprecision(0);
            salida << "# Based on " << acun << " pairs of SNPs.\n";
            salida << std::fixed << std::setprecision(2);
            salida << "# Ne point estimate:\n";
            salida << rangos_int[2] << "\n";
            salida << "# Lower limit of 50% CI:\n";
            salida << rangos_int[1] << "\n";
            salida << "# Upper bound of 50% CI:\n";
            salida << rangos_int[3] << "\n";
            salida << "# Lower limit of 90% CI:\n";
            salida << rangos_int[0] << "\n";
            salida << "# Upper limit of 90% CI:\n";
            salida << rangos_int[4] << "\n";
            // }
            salida << "# Estimated Number of full siblings that a random individuals has in the population:\n";
            if (params.flagks || params.flagnok)
            {
                if (params.flagnok && flagnoestimaks)
                {
                    salida << "Not estimated because F value is too large\n#\n";
                }
                else
                {
                    salida << params.K << "\n#\n";
                }
            }
            else
            {
                salida << "Not calculated\n#\n";
            }
        }
        else{
            salida << "\nThe solution based on LD between all loci pairs does not converge. \n";
            if ((params.flagks) || params.flagnok || (params.K > 0))
            {
                salida << "\nThe LD data and the number of siblings are incongruent to each other. \n#\n";
            }
        }

        f_pob = fp;
        if (flag_chr)
        {
            effeneind_h = effeneind * 2;
            double samx = (effeneind_h - 2.0) * (effeneind_h - 2.0) * (effeneind_h - 2.0);
            samx += 8.0 / 5.0 * (effeneind_h - 2.0) * (effeneind_h - 2.0);
            samx += 4 * (effeneind_h - 2.0);
            samx /= ((effeneind_h - 1.0) * (effeneind_h - 1.0) * (effeneind_h - 1.0) + (effeneind_h - 1.0) * (effeneind_h - 1.0));
            double samy = (2.0 * effeneind_h - 4.0) / ((effeneind_h - 1.0) * (effeneind_h - 1.0));
            log10_ncrom_sample = log10(60); // 60 cromosomas es c=0.5
            log10_n_sample = log10(effeneind);
            log10_nSNPs_solicitado = log10(eneloc);
            backfp = fp;
            d2_pob05 = (d2s05 - samy) / (samx * 0.25); // APROXIMADO
            for (fciclo = 0; fciclo < 2; ++fciclo)
            { // Two tries in case d2 negative
                Neant = 0;
                if (params.flagks || params.flagnok)
                {
                    params.K = 0;
                    for (i = 0; i < 30; ++i)
                    {
                        Ecuacion05();
                        if (params.ks == 0)
                        {
                            break;
                        }
                        params.K = Ne / (eneind - 1) * params.ks;
                        if ((std::abs(Ne - Neant) < 0.1) || (std::abs(Ne - Neant) < (Ne / 10000)))
                        {
                            break;
                        }
                        Neant = Ne;
                    }
                }
                else
                {
                    Ecuacion05();
                }
                if (Ne < 100000000)
                {
                    break;
                }
                fp = 0;
            }
            fp = backfp;
            Ne_05 = Ne;
            if (Ne < 100000000){
                log10_Ne_obs = log10(Ne_05);
                for (i = 0; i < 5; ++i)
                {
                    lim = -nDT[i];
                    if (i != 2)
                    {
                        rangos_05[i] = CalculaIntervalo_soloLD();
                    }
                }
                rangos_05[2] = log10_Ne_obs;
                for (i = 0; i < 5; ++i)
                {
                    rangos_05[i] = pow(10, rangos_05[i]);
                }

                salida << "#\n"
                        << "# Ne estimation based only on LD between chromosomes:\n";
                // if (rangos_05[2]>100000000){
                //     salida<<"#  It is not possible to calculate Ne: the estimate of d2 in the population is negative!! ("<<d2p05<<")\n";
                // }
                // else{
                salida << std::fixed << std::setprecision(0);
                salida << "# Based on " << acun05 << " pairs of SNPs between chromosomes.\n";
                salida << std::fixed << std::setprecision(2);
                salida << "# Ne point estimate:\n";
                salida << rangos_05[2] << "\n";
                salida << "# Lower limit of 50% CI:\n";
                salida << rangos_05[1] << "\n";
                salida << "# Upper limit of 50% CI:\n";
                salida << rangos_05[3] << "\n";
                salida << "# Lower limit of 90% CI:\n";
                salida << rangos_05[0] << "\n";
                salida << "# Upper limit of 90% CI:\n";
                salida << rangos_05[4] << "\n";
                // }
                salida << "# Estimated Number of full siblings that a random individuals has in the population (c=0.5):\n";
                if (params.flagks || params.flagnok)
                {
                    if (params.flagnok && flagnoestimaks)
                    {
                        salida << "Not estimated because F value is too large\n#\n";
                    }
                    else
                    {
                        salida << params.K << "\n#\n";
                    }
                }
                else
                {
                    salida << "Not calculated\n#\n";
                }
            }
            else{
                salida << "\nThe solution based on the LD between SNPs on different chromosomes does not converge. \n";
                if ((params.flagks) || params.flagnok || (params.K > 0))
                {
                    salida << "\nThe LD data and the number of siblings are incongruent to each other. \n#\n";
                }
            }

            // DENTRO DE CROMOSOMAS:
            log10_ncrom_sample = log10(genomesize/2); // 
            log10_n_sample = log10(effeneind);
            log10_nSNPs_solicitado = log10(eneloc);
            backfp = fp;
            //            d2_poblink = (d2slink - samy) / (samx * (1-cmed)*(1-cmed)); // APROXIMADO
            for (fciclo = 0; fciclo < 2; ++fciclo)
            { // Two tries in case d2 negative
                Neant = 0;
                if (params.flagks || params.flagnok)
                {
                    params.K = 0;
                    for (i = 0; i < 30; ++i)
                    {
                        IntegralUNCROM();
                        if (params.ks == 0)
                        {
                            break;
                        }
                        params.K = Ne / (eneind - 1) * params.ks;
                        if ((std::abs(Ne - Neant) < 0.1) || (std::abs(Ne - Neant) < (Ne / 10000)))
                        {
                            break;
                        }
                        Neant = Ne;
                    }
                }
                else
                {
                    IntegralUNCROM();
                }
                if (Ne < 100000000)
                {
                    break;                
                }
                fp = 0;
            }
            fp = backfp;
            Ne_link = Ne;
            if (Ne < 100000000)
            {
                log10_Ne_obs = log10(Ne_link);
                for (i = 0; i < 5; ++i)
                {
                    lim = -nDT[i];
                    if (i != 2)
                    {
                        rangos_link[i] = CalculaIntervalo_soloLD();
                    }
                }
                rangos_link[2] = log10_Ne_obs;
                for (i = 0; i < 5; ++i)
                {
                    rangos_link[i] = pow(10, rangos_link[i]);
                }

                salida << "#\n"<< "# Ne estimation based on LD within chromosomes:\n";
                salida << std::fixed << std::setprecision(0);
                salida << "# Based on " << acunlink << " pairs of SNPs within chromosomes.\n";
                salida << std::fixed << std::setprecision(2);
                salida << "# Ne point estimate:\n";
                salida << rangos_link[2] << "\n";
                salida << "# Lower limit of 50% CI:\n";
                salida << rangos_link[1] << "\n";
                salida << "# Upper limit of 50% CI:\n";
                salida << rangos_link[3] << "\n";
                salida << "# Lower limit of 90% CI:\n";
                salida << rangos_link[0] << "\n";
                salida << "# Upper limit of 90% CI:\n";
                salida << rangos_link[4] << "\n";
                salida << "# Estimated Number of full siblings that a random individuals has in the population:\n";
                if (params.flagks || params.flagnok)
                {
                    if (params.flagnok && flagnoestimaks)
                    {
                        salida << "Not estimated because F value is too large\n#\n";
                    }
                    else
                    {
                        salida << params.K << "\n#\n";
                    }
                }
                else
                {
                    salida << "Not calculated\n#\n";
                }
            }
            else{
                salida << "\nThe solution based on the LD between SNPs on the same chromosome does not converge. \n";
                if ((params.flagks) || params.flagnok || (params.K > 0))
                {
                    salida << "\nThe LD data and the number of siblings are incongruent to each other. \n#\n";
                }
            }
        }

        if (nparhermanos > 0)
        {
            salida << "# Average number of full siblings per individuals in the sample:\n";
            salida << params.ks<<"\n";
            salida << "# Full sibling pairs (individuals are referenced by their ordinals in the ped file):\n";
            if (flagnoestimaks)
            {
                salida << "Not predicted because F value is too large.\n";
            }
            else
            {
                // for (i = 0; i < nparhermanos; ++i)
                // {
                //     salida << hermanos[0][i] + 1 << ", " << hermanos[1][i] + 1 << "\n";
                // }
                salida<<salidahermanos.str();
            }
            salida << "#\n";
        }
    }

    tpas = (omp_get_wtime() - tini);
    salida0 << (float(tpas)) << "sec\n";
    salida0 << "#\n";

    if (!params.quiet)
    {
        if (params.printToStdOut)
        {
            std::string output = salida0.str()+salida.str();
           std::cout << output << std::endl;
        }
        else
        {
            std::string fichsal = "";
            if (fichspecified == "")
            {
                if (params.s>1){
                    fichsal = prefijo + "_currentNe2_mix_OUTPUT.txt";
                }
                else{
                    fichsal = prefijo + "_currentNe2_OUTPUT.txt";
                }
            }
            else
            {
                fichsal = fichspecified;
            }
            std::ofstream outputFile;
            outputFile.open(fichsal);
            outputFile << salida0.str();
            outputFile << salida.str();
            outputFile.close();
            std::cout << " End of process. Output file " << fichsal << " generated\n";
        }
    }

    std::remove(fichProgress.c_str());
    return 0;
}

//  SOLO ENTRE CROMOSOMAS DISTINTOS:
//  Se calcula cuando hay un fichero .map que tenga al menos la asignación de cada
//  marcador a su cromosoma. No hacen falta que tenga las posiciones físicas o genéticas.
void Ecuacion05()
{
    int MIL, DOSMIL, i, ii, conta;
    double sample_size, sample_size_h, samplex, sampley;
    double d2spred, increNe;
    double fp12, K12;
    bool flagbreak = false;

    K12 = (1.0 + params.K / 4.0);
    MIL = 1000;
    DOSMIL = 2000;
    fp12 = (1 + fp) * (1 + fp);
    sample_size = effeneind;
    sample_size_h = sample_size * 2;
    samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 4 * (sample_size_h - 2.0);
    samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
    sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));
    Ne = 1000;
    for (ii = 0; ii < 2; ++ii)
    {
        increNe = Ne / (4 * (ii + 1));
        for (i = 0; i < 20; ++i)
        {
            for (conta = 0; conta < DOSMIL; ++conta)
            {
                d2p05 = (1 + 0.25 * K12 + 1 / Ne) / (Ne * 1.5 + 0.55);
                d2spred = (d2p05 * 0.25 * samplex + ((params.K / 2 + 2 * 0.25 * K12) / (4 * Ne)) * samplex + sampley) * fp12;
                if (d2spred > d2s05)
                {
                    if (increNe < 0)
                    {
                        increNe = -increNe / 5;
                        break;
                    }
                }
                if (d2spred < d2s05)
                {
                    if (increNe > 0)
                    {
                        increNe = -increNe / 5;
                        break;
                    }
                }
                if ((Ne + increNe) < 3)
                {
                    increNe = increNe / 5;
                    break;
                }
                else
                {
                    Ne += increNe;
                }
                if (abs(increNe) < 0.01)
                {
                    break;
                }
                if (Ne > 100000000)
                {
                    flagbreak = true;
                    break;
                }
                // if (flagbreak){
                //     break;
                // }
            }
            if (abs(increNe) < 0.1)
            {
                break;
            }
            if (flagbreak)
            {
                break;
            }
        }
        if (flagbreak)
        {
            break;
        }
    }

}

//  MIGRACION entre dos subpoblaciones CON c=0.5:
void MixEcuacion05()
{
    int MIL, DOSMIL, i, ii, conta;
    double sample_size, sample_size_h, samplex, sampley;
    double d2spred, increNe;
    double fp12;
    bool flagbreak = false;

    // K12 = (1.0 + params.K / 4.0);
    // MIL = 1000;
    DOSMIL = 2000;
    fp12 = (1 + fp) * (1 + fp);
    sample_size = effeneind;
    sample_size_h = sample_size * 2;
    samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 4 * (sample_size_h - 2.0);
    samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
    sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));
    c=0.5;
    // Ne = Ne05;
    // Ne = std::min(std::max(Ne, 1000.0), 1000000.0);
    Ne = 100000;

    c2 = c * c;
    c12 = (1 - c) * (1 - c);
    m= (1-Fst)/(4*Ne*Fst*pow((params.s/(params.s-1)),2));

    // ALTERNATIVO:
    // m = (1-pow(1-1/(Fst*4*Ne)+1/(4*Ne),0.5))/2;

    if (m>0.5){m=0.5;}
    ms=params.s*m/(params.s-1);

    for (ii = 0; ii < 2; ++ii)
    {
        increNe = Ne / (4 * (ii + 1));
        for (i = 0; i < 20; ++i)
        {
            for (conta = 0; conta < DOSMIL; ++conta)
            {
                m12=(1-ms)*(1-ms);
                
                // AA=(1-m12)*m12*(1-1/(2*Ne));
                // BB=2*(1-m12)*(1-2.2/(2*params.s*Ne))*(1-c);
                // CC=1-(1-2.2/(2*params.s*Ne))*c12;
                // DD=1-m12*(1-1/(2*Ne))*(1-c);
                // EE=(1-m12)*(1-m12)*(1-c)*(1-2.2/(2*params.s*Ne));

                // Db2=Wt*Fst*Fst;
                // DbDw=Db2*(AA/DD+(EE/CC+(AA*BB)/(CC*DD))/(2*params.s*Ne*DD))+Wt*(1-Fst)*(1-Fst)/(4*params.s*params.s*Ne*Ne*CC*DD);
                // Dw2=Db2*(EE/CC+(AA*BB)/(CC*DD))+Wt*(1-Fst)*(1-Fst)*(1+c2)/(2*params.s*Ne*CC);

                Db2=(Wt)*Fst*Fst;
                Dw2=(Wt)*(1-Fst)*(1-Fst)*(1+c2)/(4*Ne*(1-c12)+2.2*c12);
                DbDw=(4*(Wt)*Fst*Fst*m12*(m)+(Dw2)/(4*Ne))/(1-m12*(1-c)*(1-1/2/Ne));

                Dw2p = Dw2*c12 + Db2*(1-m12)*(1-m12) + DbDw*2*(1-m12)*(1-c);
                DbDwp = DbDw*m12*(1-c) + Db2*m12*(1-m12);
                Db2p = Db2*m12*m12;                

                Dw2=Dw2p;
                DbDw=DbDwp;
                Db2=Db2p;

                d2spred=(Dw2+4*Db2+4*DbDw)/Wt*samplex+sampley *fp12;

                if (d2spred > d2s05)
                {
                    if (increNe < 0)
                    {
                        increNe = -increNe / 5;
                        break;
                    }
                }
                if (d2spred < d2s05)
                {
                    if (increNe > 0)
                    {
                        increNe = -increNe / 5;
                        break;
                    }
                }
                if ((Ne + increNe) < 3)
                {
                    increNe = increNe / 5;
                    break;
                }
                else
                {
                    Ne += increNe;
                    m= (1-Fst)/(4*Ne*Fst*pow((params.s/(params.s-1)),2));

                    // ALTERNATIVO:
                    // m = (1-pow(1-1/(Fst*4*Ne)+1/(4*Ne),0.5))/2;

                    if (m>0.5){m=0.5;}
                    ms=params.s*m/(params.s-1);
                }
                if (abs(increNe) < 0.01)
                {
                    break;
                }
                if (Ne > 100000000)
                {
                    flagbreak = true;
                    break;
                }
                // if (flagbreak){
                //     break;
                // }
            }
            if (abs(increNe) < 0.1)
            {
                break;
            }
            if (flagbreak)
            {
                break;
            }
        }
        if (flagbreak)
        {
            break;
        }
    }

}

//  DENTRO Y ENTRE CROMOSOMAS :
void IntegralVARIOSCROM()
{
    int MIL, DOSMIL, i, ii, j, minj, k, conta;
    double sample_size, sample_size_h, samplex, sampley;
    double d2spred, increL, increNe, sumad2spred, sumafrec;
    double fp12, tamacrom, distmin, K12;
    bool flagbreak = false;

    K12 = (1.0 + params.K / 4.0);
    tamacrom = genomesize / float(ncrom_sample_int);
    // distmin = genomesize / eneloc;
    // if (distmin > (tamacrom))
    // {
    //     distmin = tamacrom;
    // }
    // MIL = int(1.0 / distmin);
    // // if (MIL<1000){MIL*=2;}
    // // MIL=int(float(MIL)*tamacrom);
    // if (MIL < 100)
    // {
    //     MIL = 100;
    // }

    increL=increMorgans;
    if (flag_z){
        distmin=increL/2 + params.z/100; // en Morgans
        minj=int(distmin/increMorgans);
    }
    else{
        distmin=increL/2;
        minj=0;
    }
    fp12 = (1 + fp) * (1 + fp);
    DOSMIL = 2000;
    sample_size = effeneind;
    sample_size_h = sample_size * 2;
    samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 4 * (sample_size_h - 2.0);
    samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
    sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));
    Ne = 1000;
    for (ii = 0; ii < 2; ++ii)
    {
        increNe = Ne / (4 * (ii + 1));
        for (i = 0; i < 10; ++i)
        {
            for (conta = 0; conta < DOSMIL; ++conta)
            {
                // increL = tamacrom / MIL;
                sumad2spred = 0;
                distancia = distmin;
                sumafrec = 0;
                j = minj;
                while (distancia < maxdistance)
                {
                    c = (1 - exp(-2 * distancia)) / 2;
                    c2 = c * c;
                    c12 = (1 - c) * (1 - c);
                    d2p = (1 + c2 * K12 + 1 / Ne) / (2 * Ne * (1 - c12) + 2.2 * c12);
                    sumad2spred += mapdist[j] * (d2p * c12 * samplex + ((params.K / 2 + 2 * c2 * K12) / (4 * Ne)) * samplex + sampley) * fp12; // La integral
                    distancia += increL;
                    sumafrec += mapdist[j];
                    ++j;
                }
                d2p = (1 + 0.25 * K12 + 1 / Ne) / (Ne * 1.5 + 0.55);
                d2spred = (d2p * 0.25 * samplex + ((params.K / 2 + 2 * 0.25 * K12) / (4 * Ne)) * samplex + sampley) * fp12;
                if (sumafrec > 0)
                {
                    if (flag_chr){
                        d2spred = float(acun05) / float(acun05+acunlink) * d2spred;
                        d2spred += (sumad2spred / sumafrec) * float(acunlink) / float(acun05+acunlink);
                    }
                    else{
                        d2spred = (Ncrom-1)/Ncrom * d2spred;
                        d2spred += (sumad2spred / sumafrec) / Ncrom;
                    }
                }
                if (d2spred > d2s)
                {
                    if (increNe < 0)
                    {
                        increNe = -increNe / 5;
                        break;
                    }
                }
                if (d2spred < d2s)
                {
                    if (increNe > 0)
                    {
                        increNe = -increNe / 5;
                        break;
                    }
                }
                if ((Ne + increNe) < 3)
                {
                    increNe = increNe / 5;
                    break;
                }
                else
                {
                    Ne += increNe;
                }
                if (abs(increNe) < 0.1)
                {
                    break;
                }
                if (Ne > 100000000)
                {
                    flagbreak = true;
                    break;
                }
                // if (flagbreak){
                //     break;
                // }
            }
            if (abs(increNe) < 0.1)
            {
                break;
            }
            if (flagbreak)
            {
                break;
            }
        }
        if (flagbreak)
        {
            break;
        }
    }
}

//  SOLO DENTRO DE CROMOSOMAS:
//  Se utiliza con la opción -x cuando Fis es >> 0
//  Exige la existencia de mapa
//  Se asume que cada cromosoma es de un morgan aproximadamente.
void IntegralUNCROM()
{
    int MIL, DOSMIL, i, ii, j,minj, k, conta;
    double sample_size, sample_size_h, samplex, sampley;
    double d2spred, increL, increNe, sumad2spred, sumafrec;
    double fp12, tamacrom, distmin, K12;
    bool flagbreak = false;

    K12 = (1.0 + params.K / 4.0);
    // tamacrom = genomesize / float(ncrom_sample_int);
    // distmin = genomesize / eneloc;
    // if (distmin > (tamacrom))
    // {
    //     distmin = tamacrom;
    // }
    // increL=increMorgans;
    // minj=int(distmin/increL);
    // distmin=minj*increL+increL/2;

    increL=increMorgans;
    if (flag_z){
        distmin=increL/2 + params.z/100;
        minj=int(distmin/increMorgans);
    }
    else{
        distmin=increL/2;
        minj=0;
    }
    fp12 = (1 + fp) * (1 + fp);
    DOSMIL = 2000;
    sample_size = effeneind;
    sample_size_h = sample_size * 2;
    samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 4 * (sample_size_h - 2.0);
    samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
    sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));
    Ne = 1000;
    for (ii = 0; ii < 2; ++ii)
    {
        increNe = Ne / (4 * (ii + 1));
        for (i = 0; i < 10; ++i)
        {
            for (conta = 0; conta < DOSMIL; ++conta)
            {
                sumad2spred = 0;
                distancia = distmin;
                sumafrec = 0;
                j = minj;
                while (distancia < maxdistance)
                {
                    c = (1 - exp(-2 * distancia)) / 2;
                    c2 = c * c;
                    c12 = (1 - c) * (1 - c);
                    d2p = (1 + c2 * K12 + 1 / Ne) / (2 * Ne * (1 - c12) + 2.2 * c12);
                    sumad2spred += mapdist[j] * (d2p * c12 * samplex + ((params.K / 2 + 2 * c2 * K12) / (4 * Ne)) * samplex + sampley) * fp12; // La integral
                    distancia += increL;
                    sumafrec += mapdist[j];
                    ++j;
                }
                if (sumafrec > 0)
                {
                    d2spred = (sumad2spred / sumafrec);
                }
                else
                {
                    std::cerr << "Within chromosome integral cannot be computed." << std::endl;
                    exit(EXIT_FAILURE);
                }
                if (d2spred > d2slink)
                {
                    if (increNe < 0)
                    {
                        increNe = -increNe / 5;
                        break;
                    }
                }
                if (d2spred < d2slink)
                {
                    if (increNe > 0)
                    {
                        increNe = -increNe / 5;
                        break;
                    }
                }
                if ((Ne + increNe) < 3)
                {
                    increNe = increNe / 5;
                    break;
                }
                else
                {
                    Ne += increNe;
                }
                if (abs(increNe) < 0.1)
                {
                    break;
                }
                if (Ne > 100000000)
                {
                    flagbreak = true;
                    break;
                }
                // if (flagbreak){
                //     break;
                // }
            }
            if (abs(increNe) < 0.1)
            {
                break;
            }
            if (flagbreak)
            {
                break;
            }
        }
        if (flagbreak)
        {
            break;
        }
    }
}

void MixIntegralUNCROM(){
    int MIL, DOSMIL, i, ii, j,minj, k, conta;
    double sample_size, sample_size_h, samplex, sampley;
    double d2spred, increL, increNe, sumad2spred, sumafrec;
    double fp12, tamacrom, distmin, K12;
    bool flagbreak = false;

    K12 = (1.0 + params.K / 4.0);
    // tamacrom = genomesize / float(ncrom_sample_int);
    // distmin = genomesize / eneloc;
    // if (distmin > (tamacrom))
    // {
    //     distmin = tamacrom;
    // }
    // increL=increMorgans;
    // minj=int(distmin/increL);
    // distmin=minj*increL+increL/2;

    increL=increMorgans;
    if (flag_z){
        distmin=increL/2 + params.z/100;
        minj=int(distmin/increMorgans);
    }
    else{
        distmin=increL/2;
        minj=0;
    }
    fp12 = (1 + fp) * (1 + fp);
    DOSMIL = 2000;
    sample_size = effeneind;
    sample_size_h = sample_size * 2;
    samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 4 * (sample_size_h - 2.0);
    samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
    sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));
    // Ne = Nelink;
    Ne=100000;
    // if (Ne<1000){Ne=1000;}
    // if (Ne>1000000){Ne=1000000;}
    m= (1-Fst)/(4*Ne*Fst*pow((float(params.s)/float(params.s-1)),2));

    // ALTERNATIVO:
    // m = (1-pow(1-1/(Fst*4*Ne)+1/(4*Ne),0.5))/2;

    if (m>0.5){m=0.5;}
    ms=params.s*m/(params.s-1);
    for (ii = 0; ii < 2; ++ii)
    {
        increNe = Ne / (4 * (ii + 1));
        for (i = 0; i < 10; ++i)
        {
            for (conta = 0; conta < DOSMIL; ++conta)
            {
                sumad2spred = 0;
                distancia = distmin;
                sumafrec = 0;
                j = minj;
                while (distancia < maxdistance)
                {
                    c = (1 - exp(-2 * distancia)) / 2;
                    c2 = c * c;
                    c12 = (1 - c) * (1 - c);

                    // m12=(1-m)*(1-m);
                    // m22=(1-2*m)*(1-2*m);
                    m12=(1-ms)*(1-ms);

                    // AA=(1-m12)*m12*(1-1/(2*Ne));
                    // BB=2*(1-m12)*(1-2.2/(2*params.s*Ne))*(1-c);
                    // CC=1-(1-2.2/(2*params.s*Ne))*c12;
                    // DD=1-m12*(1-1/(2*Ne))*(1-c);
                    // EE=(1-m12)*(1-m12)*(1-c)*(1-2.2/(2*params.s*Ne));

                    // Db2=Wt*Fst*Fst;
                    // DbDw=Db2*(AA/DD+(EE/CC+(AA*BB)/(CC*DD))/(2*params.s*Ne*DD))+Wt*(1-Fst)*(1-Fst)/(4*params.s*params.s*Ne*Ne*CC*DD);
                    // Dw2=Db2*(EE/CC+(AA*BB)/(CC*DD))+Wt*(1-Fst)*(1-Fst)*(1+c2)/(2*params.s*Ne*CC);

                    Db2=(Wt)*Fst*Fst;
                    Dw2=(Wt)*(1-Fst)*(1-Fst)*(1+c2)/(4*Ne*(1-c12)+2.2*c12);
                    DbDw=(4*(Wt)*Fst*Fst*m12*(m)+(Dw2)/(4*Ne))/(1-m12*(1-c)*(1-1/2/Ne));

                    Dw2p = Dw2*c12 + Db2*(1-m12)*(1-m12) + DbDw*2*(1-m12)*(1-c);
                    DbDwp = DbDw*m12*(1-c) + Db2*m12*(1-m12);
                    Db2p = Db2*m12*m12;                
                    
                    Dw2=Dw2p;
                    DbDw=DbDwp;
                    Db2=Db2p;

                    // d2spred=(Dw2+4*Db2+4*DbDw)/Wt;
                    d2spred=(Dw2+4*Db2+4*DbDw)/Wt*samplex+sampley * fp12;
                    sumad2spred += mapdist[j] * d2spred; // La integral
                    distancia += increL;
                    sumafrec += mapdist[j];
                    ++j;
                }
                if (sumafrec > 0)
                {
                    // d2spred = (sumad2spred / sumafrec)*samplex+sampley * fp12;
                    d2spred = (sumad2spred / sumafrec);
                }
                else
                {
                    std::cerr << "Within chromosome integral cannot be computed." << std::endl;
                    exit(EXIT_FAILURE);
                }
                if (d2spred > d2slink)
                {
                    if (increNe < 0)
                    {
                        increNe = -increNe / 5;
                        break;
                    }
                }
                if (d2spred < d2slink)
                {
                    if (increNe > 0)
                    {
                        increNe = -increNe / 5;
                        break;
                    }
                }
                if ((Ne + increNe) < 3)
                {
                    increNe = increNe / 5;
                    break;
                }
                else
                {
                    Ne += increNe;
                    m= (1-Fst)/(4*Ne*Fst*pow((float(params.s)/float(params.s-1)),2));

                    // ALTERNATIVO:
                    // m = (1-pow(1-1/(Fst*4*Ne)+1/(4*Ne),0.5))/2;
                    
                    if (m>0.5){m=0.5;}
                    ms=float(params.s)*m/float(params.s-1);
               }
                if (abs(increNe) < 0.1)
                {
                    break;
                }
                if (Ne > 100000000)
                {
                    flagbreak = true;
                    break;
                }
                // if (flagbreak){
                //     break;
                // }
            }
            if (abs(increNe) < 0.1)
            {
                break;
            }
            if (flagbreak)
            {
                break;
            }
        }
        if (flagbreak)
        {
            break;
        }
    }
}

double linear(double xx)
{
    return xx;
}
double sigmoid(double xx)
{
    return 1 / (1 + exp(-xx));
}

double Funcion_desv2logsoloLD_manip_log10Lsigmoid_mae_ADAM_4()
{
    double scaled_log10_n_sample, scaled_log10_nSNPs_sample, scaled_log10_Ne_obs, scaled_log10_ncrom_sample;
    double L0_N0, L0_N1, L0_N2, L0_N3;
    double L1_N0;
    double L2_N0;
    double desv2_log10;

    scaled_log10_n_sample = (log10_n_sample - (1.55102205)) / (0.30404840);
    scaled_log10_nSNPs_sample = (log10_nSNPs_solicitado - (3.66867408)) / (1.379187800);
    scaled_log10_Ne_obs = (log10_Ne_obs - (3.16535400)) / (0.57021052);
    scaled_log10_ncrom_sample = (log10_ncrom_sample - (1.18831731)) / (0.22132427);
    L0_N0 = linear(-2.05769849 + (scaled_log10_n_sample * -1.46369338) + (scaled_log10_nSNPs_sample * -0.34379384) + (scaled_log10_Ne_obs * 0.62112820) + (scaled_log10_ncrom_sample * 0.28521538));
    L0_N1 = linear(-0.00346551 + (scaled_log10_n_sample * 0.04352528) + (scaled_log10_nSNPs_sample * 0.00146236) + (scaled_log10_Ne_obs * 0.04874034) + (scaled_log10_ncrom_sample * 0.02056176));
    L0_N2 = linear(-1.84876776 + (scaled_log10_n_sample * -0.47447950) + (scaled_log10_nSNPs_sample * -0.30161744) + (scaled_log10_Ne_obs * 1.53283715) + (scaled_log10_ncrom_sample * 0.24124029));
    L0_N3 = linear(-0.00046321 + (scaled_log10_n_sample * 0.04229134) + (scaled_log10_nSNPs_sample * 0.00025067) + (scaled_log10_Ne_obs * 0.05397954) + (scaled_log10_ncrom_sample * 0.03016489));
    L1_N0 = sigmoid(1.73418295 + (L0_N0 * -0.70259154) + (L0_N1 * 0.00812468) + (L0_N2 * -0.38448384) + (L0_N3 * 0.01005468));
    L2_N0 = linear(3.26984167 + (L1_N0 * -3.65439820));
    desv2_log10 = L2_N0 * 0.40417856 + (0.15538608);
    return desv2_log10;
}

double CalculaIntervalo_soloLD()
{
    double log10_Nest, desv, incre, nciclos, ciclo, punto;
    int it;

    DT = Funcion_desv2logsoloLD_manip_log10Lsigmoid_mae_ADAM_4();
    if (DT<0){DT=0;}
    DT = pow(DT, 0.5);
    desv = DT * lim;
    incre = DT * lim / 10;
    log10_Nest = log10_Ne_obs - desv;
    if (log10_Nest < 0.7)
    {
        log10_Nest = 0.7;
    }
    nciclos = 4;
    ciclo = 0;
    for (it = 0; it < 100000; ++it)
    {
        log10_Nest += incre;
        if (log10_Nest < 0.7)
        {
            log10_Nest = 0.7;
        }
        DT = Funcion_desv2logsoloLD_manip_log10Lsigmoid_mae_ADAM_4();
        DT = pow(DT, 0.5);
        desv = DT * lim;
        punto = log10_Nest + desv;
        if (((punto > log10_Ne_obs) && (incre > 0)) || ((punto < log10_Ne_obs) && (incre < 0)))
        {
            incre = -incre / 5;
            ciclo += 1;
            if (ciclo > nciclos)
            {
                break;
            }
        }
    }
    return log10_Nest;
}
