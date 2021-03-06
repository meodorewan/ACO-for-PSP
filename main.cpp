#include <bits/stdc++.h>

#define MAXA 20
#define MAXN 280

using namespace std;
ofstream flog;                                          // file log
string file_name;                                       // ten file log va file out
int number_of_round;                                    // so round kien chay
int ant_per_round;                                      // so kien chay 1 round
int time_limit;                                         // gioi han thoi gian chay
int reset_count_down;                                   // so lan kq can lap de reset mui
double alpha;                                           // alpha
double beta;                                            // beta
double rho;                                             // toc do bay hoi
int ls_flag;                                            // cach local_search
int seed;                                               // so khoi tao ham random
int k;                                                  // so luong huong duoc luu
string s;                                               // xau nhap vao
int n;                                                  // do dai xau
const string AA = "CMFILVWYAGTSQNEDHRKP";               // xau bieu dien cac amino acid
vector<double> T[MAXN];                                 // ma tran mui

// y0 trung ten voi ham trong math.h
namespace directions
{
    int x0[12] = { 1,-1,-1, 0, 1,-1, 1, 1, 0, 0,-1, 0};
    int y0[12] = { 1,-1, 1, 1, 0, 0,-1, 0, 1,-1, 0,-1};
    int z0[12] = { 0, 0, 0,-1, 1,-1, 0,-1, 1,-1, 1, 1};
}

const double MJ_ENERGY[20][20] = {                      //ma tran MJ energy
    {-1.06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.19, 0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-0.23, -0.42, -0.44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.16, -0.28, -0.19, -0.22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-0.08, -0.2, -0.3, -0.41, -0.27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.06, -0.14, -0.22, -0.25, -0.29, -0.29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.08, -0.67, -0.16, 0.02, -0.09, -0.17, -0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.04, -0.13, 0, 0.11, 0.24, 0.02, -0.04, -0.06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0.25, 0.03, -0.22, -0.01, -0.1, -0.09, 0.09, -0.13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-0.08, 0.19, 0.38, 0.25, 0.23, 0.16, 0.18, 0.14, -0.07, -0.38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.19, 0.19, 0.31, 0.14, 0.2, 0.25, 0.22, 0.13, -0.09, -0.26, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-0.02, 0.14, 0.29, 0.21, 0.25, 0.18, 0.34, 0.09, -0.06, -0.16, -0.08, 0.2, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.05, 0.46, 0.49, 0.36, 0.26, 0.24, 0.08, -0.2, 0.08, -0.06, -0.14, -0.14, 0.29, 0, 0, 0, 0, 0, 0, 0},
    {0.13, 0.08, 0.18, 0.53, 0.3, 0.5, 0.06, -0.2, 0.28, -0.14, -0.11, -0.14, -0.25, -0.53, 0, 0, 0, 0, 0, 0},
    {0.69, 0.44, 0.27, 0.35, 0.43, 0.34, 0.29, -0.1, 0.26, 0.25, 0, -0.26, -0.17, -0.32, -0.03, 0, 0, 0, 0, 0},
    {0.03, 0.65, 0.39, 0.59, 0.67, 0.58, 0.24, 0, 0.12, -0.22, -0.29, -0.31, -0.17, -0.3, -0.15, 0.04, 0, 0, 0, 0},
    {-0.19, 0.99, -0.16, 0.49, 0.16, 0.19, -0.12, -0.34, 0.34, 0.2, -0.19, -0.05, -0.02, -0.24, -0.45, -0.39, -0.29, 0, 0, 0},
    {0.24, 0.31, 0.41, 0.42, 0.35, 0.3, -0.16, -0.25, 0.43, -0.04, -0.35, 0.17, -0.52, -0.14, -0.74, -0.72, -0.12, 0.11, 0, 0},
    {0.71, 0, 0.44, 0.36, 0.19, 0.44, 0.22, -0.21, 0.14, 0.11, -0.09, -0.13, -0.38, -0.33, -0.97, -0.76, 0.22, 0.75, 0.25, 0},
    {0, -0.34, 0.2, 0.25, 0.42, 0.09, -0.28, -0.33, 0.1, -0.11, -0.07, 0.01, -0.42, -0.18, -0.1, 0.04, -0.21, -0.38, 0.11, 0.26}
};


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
        //AA = "";
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
            p = convert::to_int(x + directions::x0[j], y + directions::y0[j], z + directions::z0[j]);
            if(visited.count(p) > 0)
                E_MJ += MJ_ENERGY[visited[p]][amino_acid];
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
