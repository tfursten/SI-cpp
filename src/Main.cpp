#include "Main.h"


int main(int ac, char** av)
{
    namespace po = boost::program_options;
    using namespace std;

    static int nGenerations, nMaxX, nMaxY, nPollen, nOvule, nMarkers, nDel, nBurnIn, nSample, nAlleles;
    unsigned int seed;
    static double dSMut, dMMut, dDMut, dSigmaP, dSigmaS, dPdel;
    static float p1, p2;
    static bool f;
    string bound, dist_name, si, infile, outfileName, mut_type;

    ostringstream out;
    try
    {
        po::options_description generic("General Options");
        generic.add_options()
        ("help", "Produce help message")
        ;

        po::options_description config("Configuration");
        config.add_options()
        ("maxX,x", po::value<int>(&nMaxX)->default_value(100),"Set X dimension")
        ("maxY,y", po::value<int>(&nMaxY)->default_value(100),"Set Y dimension")
        ("landscape", po::value<string>(&bound)->default_value("torus"),"Set landscape boundary (torus or rectangle)")
        ("generations,g", po::value<int>(&nGenerations)->default_value(10), "Set number of Generations to run after burn-in")
        ("pollen,p", po::value<int>(&nPollen)->default_value(10), "Set number of pollen produced per individual")
        ("ovule,o", po::value<int>(&nOvule)->default_value(10), "Set number of ovules per individual")
        ("markers,n", po::value<int>(&nMarkers)->default_value(3), "Set number of markers")
        ("del_mark", po::value<int>(&nDel)->default_value(1), "Set number of deleterious markers")
        ("smut,u", po::value<double>(&dSMut)->default_value(0.00001), "Set S locus mutation rate")
        ("mut-type", po::value<string>(&mut_type)->default_value(string("IAM")), "Mutation Model for S-alleles (IAM, KAM), Markers are IAM")
        ("nAllele", po::value<int>(&nAlleles)->default_value(20), "Number of alleles for KAM, Number starting alleles for IAM")
        ("mmut,m", po::value<double>(&dMMut)->default_value(0.00001), "Set marker mutation rate")
        ("pdel", po::value<double>(&dPdel)->default_value(1.0),"Set deleterious selection coefficient")
        ("dmut", po::value<double>(&dDMut)->default_value(0.0001), "Set deleterious mutation rate for unlinked locus")
        ("distribution,d", po::value<string>(&dist_name)->default_value("tri"), "Set Dispersal Distribution")
        ("sigmaP,q", po::value<double>(&dSigmaP)->default_value(2.0), "Set dispersal parameter for pollen")
        ("sigmaS,r", po::value<double>(&dSigmaS)->default_value(2.0), "Set dispersal parameter for seed")
        ("burn,b", po::value<int>(&nBurnIn)->default_value(0),"Set Burn-in Period")
        ("sample,t", po::value<int>(&nSample)->default_value(1),"Sample every n generations after burn-in")
        ("output_file,f", po::value<string>(&outfileName)->default_value(string("data")),"Output File Name")
        ("seed", po::value<unsigned int>(&seed)->default_value(0), "Set PRNG seed")
        ("si,s", po::value<string>(&si)->default_value("nsi"), "Set self-incompatibility system")
        ("pparam", po::value<float>(&p1)->default_value(0),"Extra Parameter for pollen dispersal")
        ("sparam", po::value<float>(&p2)->default_value(0),"Extra Parameter for seed dispersal")
        ("fast", po::value<bool>(&f)->default_value(true),"Use fast dispersal when available")
        ;

        po::options_description hidden("Hidden Options");
        hidden.add_options()
        ("input-file", po::value<string>(&infile), "input file")
        ;

        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);

        po::options_description config_file_options;
        config_file_options.add(config);

        po::options_description visible("Allowed Options");
        visible.add(generic).add(config);

        po::positional_options_description p;
        p.add("input-file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(ac,av).options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);


        if (vm.count("help"))
        {
            cout << visible << "\n";
            return 0;
        }

        if (!infile.empty())
        {
            ifstream ifs(infile.c_str());
            if (!ifs)
            {
                cout << "can not open config file: "<< infile << "\n";
                return 0;
            }
            else
            {
                po::store(parse_config_file(ifs, config_file_options), vm);
                po::notify(vm);
            }
        }

        assert(dSMut>=0);
        assert(dMMut>=0);
        assert(dDMut>=0);
        assert(nMaxX>0);
        assert(nMaxY>0);
        assert(nGenerations>0);
        assert(nSample>0);
        assert(nPollen>0);
        assert(nOvules>0);
        assert(nMarkers>0);
        assert(nDel>0);
        assert(dPdel>=0);


        if(mut_type != "IAM" && mut_type != "KAM"){
            cout << "Not a valid mutation model" << endl;
            throw;
        }

        out << "X dimension set to " << nMaxX << ".\n"
            << "Y dimension set to " << nMaxY << ".\n"
            << "Run for " << nGenerations << " generations.\n"
            << "Burn " << nBurnIn << " generation(s).\n"
            << "Collect data every " << nSample << " Generation(s).\n"
            << "Number of pollen set to " << nPollen << ".\n"
            << "Number of ovules set to " << nOvule << ".\n"
            << "Number of markers set to " << nMarkers << ".\n"
            << "Number of deleterious markers set to " << nDel << ".\n"
            << "S mutation rate set to " << dSMut<< ".\n"
            << "Marker mutation rate set to " << dMMut << ".\n"
            << "Deleterious mutation rate set to " << dDMut << ".\n"
            << "Selection coefficient of deleterious allele is " << dPdel << ".\n"
            << "Mutation model set to " << mut_type << ".\n"
            << "Number of alleles set to " << nAlleles << ".\n";

        assert(dSigmaP>0);
        assert(dSigmaS>0);
        out << "Pollen dispersal parameter set to " << dSigmaP << ".\n"
            << "Extra pollen parameter set to " << p1 << ".\n"
            << "Seed dispersal parameter set to " << dSigmaS << ".\n"
            << "Extra seed paramter set to " << p2 << ".\n";
       

        if (seed)
            out << "User set PRNG seed to: " << seed << ".\n";

    }

    catch(exception& e)
    {
        cout<< e.what() << "\n";
        return 1;
    }

    string datafile = outfileName+".txt";
    string paramfile = outfileName+"_settings.txt";
    cout << "Data saved to: " << datafile << endl;
    cout << "Parameters saved to: " << paramfile << endl;
    ofstream pout;
    ofstream dout;
    pout.open(paramfile);
    dout.open(datafile);
    pout << out.str();
    cout << out.str();
    //Initialize Population
    clock_t start = clock();
    Population pop(pout, dout);
    pop.initialize(nMaxX,nMaxY,bound,nPollen,nOvule, nMarkers, nDel, dSigmaP, dSigmaS, si, dist_name, p1, p2, f, mut_type, nAlleles);
    pop.param(dSMut, dMMut, dDMut,dPdel, seed);
    //Run Simulation
    pop.evolve(nBurnIn, nGenerations, nSample);

    clock_t end = clock();
    float seconds = (float)(end-start)/ CLOCKS_PER_SEC;
    cout << "TIME: " << seconds << endl;

    return 0;

}
