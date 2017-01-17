#include <iostream>
#include <sstream>
#include <vector>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <string>
#include <ctime>
#include <map>

#define MAXA 20
#define MAXN 280

using namespace std;
ofstream flog;                          // file log
string file_name;                       // ten file log va file out
int number_of_round;                    // so round kien chay
int ant_per_round;                      // so kien chay 1 round
int time_limit;                         // gioi han thoi gian chay
int reset_count_down;                   // so lan kq can lap de reset mui
double alpha;                           // alpha
double beta;                            // beta
double rho;                             // toc do bay hoi
int ls_flag;                            // cach local_search
int seed;                               // so khoi tao ham random
int k;                                  // so luong huong duoc luu
string s;                               // xau nhap vao
int n;                                  // do dai xau
double e[MAXA][MAXA];                   // ma tran nang luong MJ
string AA;                              // xau bieu dien cac amino acid
vector<double> T[MAXN];                 // ma tran mui
int x0[12] = { 1,-1,-1, 0, 1,-1, 1, 1, 0, 0,-1, 0};
int y0[12] = { 1,-1, 1, 1, 0, 0,-1, 0, 1,-1, 0,-1};
int z0[12] = { 0, 0, 0,-1, 1,-1, 0,-1, 1,-1, 1, 1};


namespace random_picker
{
    int get_rand (int lim) // get a random number in 0..lim-1
    {
        seed = (12345LL * seed + 67890) % (1LL << 31);
        return seed % lim;
    }

    double get_rand () // get a random number in range [0,1]
    {
        seed = (12345LL * seed + 67890) % (1LL << 31);
        return double(seed) / INT_MAX;
    }

    int pick (int sz, double w[]) // pick a random element from 0..sz-1, w[i]/sum(w) is the chance of i being picked
    {
        double r = get_rand() * accumulate(w, w + sz, 0.0);

        int i = 0;
        while (w[i] < r) r -= w[i++];
        return min(i, sz - 1);
    }
};

struct OutputStream          // synced output stream for standard output and log output
{
    template<typename T>
    OutputStream& operator<< (T v)
    {
        flog << v;
        cout << v; cout.flush();
        return *this;
    }
} fout;

namespace parser
{
    void set_default_parameters ()
    {
        file_name               = "4RXN";
        number_of_round         = 500;
        ant_per_round           = 100;
        time_limit              = INT_MAX;
        reset_count_down        = 10;
        alpha                   = 1.0;
        beta                    = 3.0;
        rho                     = 0.3;
        ls_flag                 = 0;
        seed                    = 0;
        k                       = 1;
    }

    void run (int argc,char *argv[])
    {
        // set parameters to default
        set_default_parameters();
        // check if changes required
        for (int i = 1; i < argc; i += 2)
        {
            if (strcmp(argv[i], "-inst" ) == 0) file_name               = argv[i + 1];
            if (strcmp(argv[i], "-r"    ) == 0) number_of_round         = atoi(argv[i + 1]);
            if (strcmp(argv[i], "-ant"  ) == 0) ant_per_round           = atoi(argv[i + 1]);
            if (strcmp(argv[i], "-tl"   ) == 0) time_limit              = atoi(argv[i + 1]) * CLOCKS_PER_SEC;
            if (strcmp(argv[i], "-reset") == 0) reset_count_down        = atoi(argv[i + 1]);
            if (strcmp(argv[i], "-a"    ) == 0) alpha                   = atof(argv[i + 1]);
            if (strcmp(argv[i], "-b"    ) == 0) beta                    = atof(argv[i + 1]);
            if (strcmp(argv[i], "-rho"  ) == 0) rho                     = atof(argv[i + 1]);
            if (strcmp(argv[i], "-ls"   ) == 0) ls_flag                 = atoi(argv[i + 1]);
            if (strcmp(argv[i], "-seed" ) == 0) seed                    = atoi(argv[i + 1]);
            if (strcmp(argv[i], "-k"    ) == 0) k                       = atoi(argv[i + 1]);
        }

        // input
        freopen((file_name + ".in").c_str(),"r",stdin);
        cin >> s;
        n = s.length();
        // set file_name
        ostringstream seed_str, k_str;
        k_str << k;
        seed_str  << seed;
        cout << "DONE";
        // preparing outputs
        ios::sync_with_stdio(false);
        flog.open((file_name + "_" + k_str.str() +"_"+ seed_str.str() + ".log").c_str());

        // printing current settings
        fout << "ACO for Protein Structure PredictionProblem\n";
        fout << "Test                           : " << file_name                                            << "\n";

        fout << "\nParameter-settings\n";
        fout << "number_of_round                : " << number_of_round                                      << "\n";
        fout << "ant_per_round                  : " << ant_per_round                                        << "\n";
        fout << "time_limit                     : " << time_limit / CLOCKS_PER_SEC                          << " (s)\n";
        fout << "reset_count_down               : " << reset_count_down                                     << "\n";
        fout << "alpha                          : " << alpha                                                << "\n";
        fout << "beta                           : " << beta                                                 << "\n";
        fout << "rho                            : " << rho                                                  << "\n";
        fout << "ls_flag                        : " << ls_flag                                              << "\n";
        fout << "seed                           : " << seed                                                 << "\n";

        // change output stream precision
        flog.precision(2); flog << fixed;
        cout.precision(2); cout << fixed;
        AA = "";
        // gan ma tran MJ
    }
}

namespace convert
{
    int to_int(int x, int y, int z)
    {
        return x * 1000000 + y * 1000 + z;
    }
    void to_xzy(int p, int &x, int &y, int &z)
    {
        z = p % 1000;
        p /= 1000;
        y = p % 1000;
        x = p / 1000;
    }
}

struct Solution
{
    vector<int> X;
    vector<int> Y;
    vector<int> Z;
    double E_MJ;
    map <int,int> visited;

    Solution()
    {
        X = vector<int> (n,-1);
        Y = vector<int> (n,-1);
        Z = vector<int> (n,-1);
        X[0] = Y[0] = Z[0] = 500;
        E_MJ = 0;
    }

    void add(int i, int x, int y, int z, int amino_acid)
    {
        visited[convert::to_int(x,y,z)] = amino_acid;
        X[i] = x;
        Y[i] = y;
        Z[i] = z;
        int p;
        for(int j = 0; j < 12; ++j)
        {
            p = convert::to_int(x + x0[j], y + y0[j], z + z0[j]);
            if(visited.count(p) > 0)
                E_MJ += e[visited[p]][amino_acid];
        }
    }
};

namespace Local_search
{

}

namespace ACO
{
    double Tmin,Tmax;
    int MAXT;
    //Solution Gbest, Ibest;

    void reset_T()
    {

    }

    void update_T()
    {

    }

    void let_ant_run()
    {

    }

    void run()
    {
        Tmax = 1.0;
        Tmin = Tmax / 12 / n;
        MAXT = 1;
        for (int i = 0; i < k; ++i)
            MAXT *= 12;
        for (int i = 0; i < n; ++i)
            T[i] = vector<double> (MAXT,Tmax);
        for (int rnd = 0; rnd < number_of_round; ++rnd)
        {
            for (int ant = 0; ant < ant_per_round; ++ant)
                let_ant_run();
        }
    }
}

int main(int argc,char *argv[])
{
    parser::run(argc, argv);
    ACO::run();
    return 0;
}
