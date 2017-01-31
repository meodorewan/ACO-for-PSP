#include <bits/stdc++.h>

#define MAXA 20
#define MAXN 280
#define eps 0.01
#define INF 123456789
#define epsilon 0.0000001

using namespace std;
ofstream flog;                                          // file log
string file_name;                                       // ten file log va file out
string output_file;
int number_of_round;                                    // so round kien chay
int ant_per_round;                                      // so kien chay 1 round
int time_limit;                                         // gioi han thoi gian chay
int reset_count_down;                                   // so lan kq can lap de reset mui
int count_down;                                         // dem nguoc de reset mui
double alpha;                                           // alpha
double beta;                                            // beta
double rho;                                             // toc do bay hoi
int ls_flag;                                            // cach local_search
int seed;                                               // so khoi tao ham random
int k;                                                  // so luong huong duoc luu
string s;                                               // xau nhap vao
int n;                                                  // do dai xau
const string AA = "CMFILVWYAGTSQNEDHRKP";               // xau bieu dien cac amino acid
int a[MAXN];                                            // chuyen tu xau s sang so
int MAXT;                                               // 12 mu k
int algo_flag;                                          // loai thuat toan
vector<double> T[MAXN];                                 // ma tran mui
int x0[12] = { 1,-1,-1, 0, 1,-1, 1, 1, 0, 0,-1, 0};
int y0[12] = { 1,-1, 1, 1, 0, 0,-1, 0, 1,-1, 0,-1};
int z0[12] = { 0, 0, 0,-1, 1,-1, 0,-1, 1,-1, 1, 1};

const double MJ_ENERGY[20][20] = {                      //ma tran MJ energy
    {-1.06, 0.19, -0.23, 0.16, -0.08, 0.06, 0.08, 0.04, 0, -0.08, 0.19, -0.02, 0.05, 0.13, 0.69, 0.03, -0.19, 0.24, 0.71, 0},
    {0.19, 0.04, -0.42, -0.28, -0.2, -0.14, -0.67, -0.13, 0.25, 0.19, 0.19, 0.14, 0.46, 0.08, 0.44, 0.65, 0.99, 0.31, 0, -0.34},
    {-0.23, -0.42, -0.44, -0.19, -0.3, -0.22, -0.16, 0, 0.03, 0.38, 0.31, 0.29, 0.49, 0.18, 0.27, 0.39, -0.16, 0.41, 0.44, 0.2},
    {0.16, -0.28, -0.19, -0.22, -0.41, -0.25, 0.02, 0.11, -0.22, 0.25, 0.14, 0.21, 0.36, 0.53, 0.35, 0.59, 0.49, 0.42, 0.36, 0.25},
    {-0.08, -0.2, -0.3, -0.41, -0.27, -0.29, -0.09, 0.24, -0.01, 0.23, 0.2, 0.25, 0.26, 0.3, 0.43, 0.67, 0.16, 0.35, 0.19, 0.42},
    {0.06, -0.14, -0.22, -0.25, -0.29, -0.29, -0.17, 0.02, -0.1, 0.16, 0.25, 0.18, 0.24, 0.5, 0.34, 0.58, 0.19, 0.3, 0.44, 0.09},
    {0.08, -0.67, -0.16, 0.02, -0.09, -0.17, -0.12, -0.04, -0.09, 0.18, 0.22, 0.34, 0.08, 0.06, 0.29, 0.24, -0.12, -0.16, 0.22, -0.28},
    {0.04, -0.13, 0, 0.11, 0.24, 0.02, -0.04, -0.06, 0.09, 0.14, 0.13, 0.09, -0.2, -0.2, -0.1, 0, -0.34, -0.25, -0.21, -0.33},
    {0, 0.25, 0.03, -0.22, -0.01, -0.1, -0.09, 0.09, -0.13, -0.07, -0.09, -0.06, 0.08, 0.28, 0.26, 0.12, 0.34, 0.43, 0.14, 0.1},
    {-0.08, 0.19, 0.38, 0.25, 0.23, 0.16, 0.18, 0.14, -0.07, -0.38, -0.26, -0.16, -0.06, -0.14, 0.25, -0.22, 0.2, -0.04, 0.11, -0.11},
    {0.19, 0.19, 0.31, 0.14, 0.2, 0.25, 0.22, 0.13, -0.09, -0.26, 0.03, -0.08, -0.14, -0.11, 0, -0.29, -0.19, -0.35, -0.09, -0.07},
    {-0.02, 0.14, 0.29, 0.21, 0.25, 0.18, 0.34, 0.09, -0.06, -0.16, -0.08, 0.2, -0.14, -0.14, -0.26, -0.31, -0.05, 0.17, -0.13, 0.01},
    {0.05, 0.46, 0.49, 0.36, 0.26, 0.24, 0.08, -0.2, 0.08, -0.06, -0.14, -0.14, 0.29, -0.25, -0.17, -0.17, -0.02, -0.52, -0.38, -0.42},
    {0.13, 0.08, 0.18, 0.53, 0.3, 0.5, 0.06, -0.2, 0.28, -0.14, -0.11, -0.14, -0.25, -0.53, -0.32, -0.3, -0.24, -0.14, -0.33, -0.18},
    {0.69, 0.44, 0.27, 0.35, 0.43, 0.34, 0.29, -0.1, 0.26, 0.25, 0, -0.26, -0.17, -0.32, -0.03, -0.15, -0.45, -0.74, -0.97, -0.1},
    {0.03, 0.65, 0.39, 0.59, 0.67, 0.58, 0.24, 0, 0.12, -0.22, -0.29, -0.31, -0.17, -0.3, -0.15, 0.04, -0.39, -0.72, -0.76, 0.04},
    {-0.19, 0.99, -0.16, 0.49, 0.16, 0.19, -0.12, -0.34, 0.34, 0.2, -0.19, -0.05, -0.02, -0.24, -0.45, -0.39, -0.29, -0.12, 0.22, -0.21},
    {0.24, 0.31, 0.41, 0.42, 0.35, 0.3, -0.16, -0.25, 0.43, -0.04, -0.35, 0.17, -0.52, -0.14, -0.74, -0.72, -0.12, 0.11, 0.75, -0.38},
    {0.71, 0, 0.44, 0.36, 0.19, 0.44, 0.22, -0.21, 0.14, 0.11, -0.09, -0.13, -0.38, -0.33, -0.97, -0.76, 0.22, 0.75, 0.25, 0.11},
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
        beta                    = 1.0;
        rho                     = 0.3;
        ls_flag                 = 0;
        algo_flag               = 0;
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
            if (strcmp(argv[i], "-al"   ) == 0) algo_flag               = atoi(argv[i + 1]);
            if (strcmp(argv[i], "-seed" ) == 0) seed                    = atoi(argv[i + 1]);
            if (strcmp(argv[i], "-k"    ) == 0) k                       = atoi(argv[i + 1]);
        }

        // input
        freopen((file_name + ".in").c_str(),"r",stdin);
        cin >> s;
        n = s.length();
        for (int i = 0; i < n; ++i)
            a[i] = AA.find(s[i]);
        // set file_name
        ostringstream seed_str, k_str , algo_flag_str;
        k_str << k;
        seed_str  << seed;
        algo_flag_str << algo_flag;
        // preparing outputs
        ios::sync_with_stdio(false);

        output_file = file_name + "_" + k_str.str() +"_"+ seed_str.str() + "_al = " + algo_flag_str.str();
        flog.open((output_file + ".log").c_str());

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
        fout << "algo_flag                      : " << algo_flag                                            << "\n";
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

    int next_direction(int p, int h)
    {
        p = p % (MAXT/12);
        return p * 12 + h;
    }
    bool compare(double a, double b)
    {
        if (abs(a-b) <= epsilon) return true;
        return false;
    }
}

struct Solution
{
    vector<int> X;
    vector<int> Y;
    vector<int> Z;
    vector<int> H;
    double E_MJ;
    map <int,int> visited;

    Solution()
    {
        X = vector<int> (n,-1);
        Y = vector<int> (n,-1);
        Z = vector<int> (n,-1);
        H = vector<int> (n,-1);
        E_MJ = 0;
        visited.clear();
    }

    void add_first_point(int x, int y, int z)
    {
        X[0] = x;
        Y[0] = y;
        Z[0] = z;
        visited[convert::to_int(x,y,z)] = a[0];
    }

    void add(int i, int j)
    {
        int x,y,z;
        x = X[i-1] + x0[j];
        y = Y[i-1] + y0[j];
        z = Z[i-1] + z0[j];
        visited[convert::to_int(x,y,z)] = a[i];
        X[i] = x;
        Y[i] = y;
        Z[i] = z;
        H[i] = j;
        int p;
        for(int j = 0; j < 12; ++j)
        {
            p = convert::to_int(x + x0[j], y + y0[j], z + z0[j]);
            if(visited.count(p) > 0)
                E_MJ += MJ_ENERGY[visited[p]][a[i]];
        }
    }

    void minus_consecutive_neighbors()
    {
        for (int i = 1; i < n; ++i)
            E_MJ -= MJ_ENERGY[a[i]][a[i-1]];
    }

    int heuristic(int i, int j)
    {
        int x,y,z;
        x = X[i-1] + x0[j];
        y = Y[i-1] + y0[j];
        z = Z[i-1] + z0[j];
        int p;
        p = convert::to_int(x,y,z);
        if (visited.count(p) > 0)
            return INF;
        double res = E_MJ;
        for(int j = 0; j < 12; ++j)
        {
            p = convert::to_int(x + x0[j], y + y0[j], z + z0[j]);
            if(visited.count(p) > 0)
                res += MJ_ENERGY[visited[p]][a[i]];
        }
        return res;
    }
    double recaculate()
    {
        double res = 0;
        int p;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < 12; ++j)
            {
                p = convert::to_int(X[i] + x0[j], Y[i] + y0[j], Z[i] + z0[j]);
                if(visited.count(p) > 0)
                    res += MJ_ENERGY[visited[p]][a[i]];
            }
        }
        res = res / 2;
        for (int i = 1; i < n; ++i)
            res -= MJ_ENERGY[a[i]][a[i-1]];
        return res;
    }

};

namespace Local_search
{

}

namespace ACO
{
    double Tmin,Tmax;
    Solution Gbest, Ibest;
    int unfinished_trip;                                     // so lan kien khong hoan thanh hanh trinh
    double prev_E;                                           // ket qua tim duoc o luot truoc
    void reset_T()
    {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < MAXT; ++j)
                T[i][j] = Tmax;
    }

    void update_T()
    {
        int p = 0;
        for (int i = 1; i < n; ++i)
        {
            for (int j = 0; j < 12; ++j)
                T[i][convert::next_direction(p,j)] = T[i][convert::next_direction(p,j)] * (1-rho) + rho * Tmin;
            p = convert::next_direction(p,Ibest.H[i]);
            T[i][p] += rho * (Tmax - Tmin);
        }
    }

    void let_ant_run()
    {
        Solution sol;
        sol.add_first_point(500,500,500);
        //sol.X[0] = sol.Y[0] = sol.Z[0] = 500;
        double w[12] , max_w , min_w;
        int p = 0, h;
        h = 0;
        sol.add(1,h);
        p = convert::next_direction(p,h);
        for (int i = 2; i < n; ++i)
        {
            for(int j = 0; j < 12; ++j)
                w[j] = sol.heuristic(i,j);
            max_w = -INF;
            min_w = INF;
            for (int j = 0; j < 12; ++j)
                if (w[j] != INF)
                {
                    max_w = max(max_w,w[j]);
                    min_w = min(min_w,w[j]);
                }
            if (max_w == -INF)
            {
                ++unfinished_trip;
                return;
            }
            for (int j = 0; j < 12; ++j)
                if (w[j] == INF)
                    w[j] = 0;
                else w[j] = max_w - w[j] + eps ;//+ (max_w - min_w) / 10;
            for (int j = 0; j < 12; ++j)
                w[j] = pow(w[j],alpha) * pow(T[i][convert::next_direction(p,j)],beta);
            h = random_picker::pick(12,w);
            sol.add(i,h);
            p = convert::next_direction(p,h);
        }
        sol.minus_consecutive_neighbors();
        if (sol.E_MJ < Ibest.E_MJ) Ibest = sol;
    }

    void multi_ants_run()
    {
        Solution sol[2][ant_per_round];
        //sol[0][0].X[0] = sol[0][0].Y[0] = sol[0][0].Z[0] = 500;
        sol[0][0].add_first_point(500,500,500);
        int n_prev_ants = 1, n_cur_ants;
        int p[2][ant_per_round];
        int prev_ant, h;
        double w[ant_per_round][12] , max_w, min_w;
        double sum;
        double r;
        sol[0][0].add(1,0);
        memset(p,0,sizeof(p));
        int cur = 0;
        for(int i = 2; i < n; ++i)
        {
            for (int ant = 0; ant < n_prev_ants; ++ant)
                for (int j = 0; j < 12; ++j)
                    w[ant][j] = sol[cur][ant].heuristic(i,j);
            max_w = -INF;
            min_w = INF;
            for (int ant = 0; ant < n_prev_ants; ++ant)
                for (int j = 0; j < 12; ++j)
                    if (w[ant][j] != INF)
                    {
                        max_w = max(max_w,w[ant][j]);
                        min_w = min(min_w,w[ant][j]);
                    }
            for (int ant = 0; ant < n_prev_ants; ++ant)
                for (int j = 0; j < 12; ++j)
                    if (w[ant][j] != INF)
                        w[ant][j] = max_w - w[ant][j] + eps;
                    else w[ant][j] = 0;
            for (int ant = 0; ant < n_prev_ants; ++ant)
                for (int j = 0; j < 12; ++j)
                    w[ant][j] = pow(w[ant][j],alpha) * pow(T[i][convert::next_direction(p[cur][ant],j)],beta);
            sum = 0;
            for (int ant = 0; ant < n_prev_ants; ++ant)
                for (int j = 0; j < 12; ++j)
                    sum += w[ant][j];
            n_cur_ants = 0;
            for (int ant = 0; ant < ant_per_round; ++ant)
            {
                if (convert::compare(sum,0)) break;
                r = random_picker::get_rand() * sum;
                prev_ant = 0; h = 0;
                while (w[prev_ant][h] < r)
                {
                    r -= w[prev_ant][h];
                    ++h;
                    if (h == 12)
                    {
                        h = 0;
                        ++prev_ant;
                    }
                }
                sol[1-cur][ant] = sol[cur][prev_ant];
                sol[1-cur][ant].add(i,h);
                p[1-cur][ant] = convert::next_direction(p[cur][prev_ant],h);
                sum -= w[prev_ant][h];
                w[prev_ant][h] = 0;
                ++n_cur_ants;
            }
            n_prev_ants = n_cur_ants;
            cur = 1 - cur;
        }
        for (int ant = 0; ant < n_prev_ants; ++ant)
            sol[cur][ant].minus_consecutive_neighbors();
        for (int ant = 0; ant < n_prev_ants; ++ant)
            if (sol[cur][ant].E_MJ < Ibest.E_MJ)
                Ibest = sol[cur][ant];
    }

    void run()
    {
        clock_t st = clock();
        Gbest.E_MJ = INF;
        Tmax = 1.0;
        Tmin = Tmax / 12 / n;
        MAXT = 1;
        for (int i = 0; i < k; ++i)
            MAXT *= 12;
        for (int i = 0; i < n; ++i)
            T[i] = vector<double> (MAXT,Tmax);
        count_down = reset_count_down;
        for (int rnd = 0; rnd < number_of_round; ++rnd)
        {
            fout << "\nRound " << rnd << " :\n";
            Ibest.E_MJ = INF;
            if (algo_flag == 0)
            {
                // tung kien chay
                for (int ant = 0; ant < ant_per_round; ++ant)
                    let_ant_run();
            }
            if (algo_flag == 1)
            {
                // cho tat ca kien chay song song
                multi_ants_run();
            }
            if (Ibest.E_MJ < Gbest.E_MJ) Gbest = Ibest;
            // cap nhat mui
            update_T();
            // khoi tao lai mui khi lien tuc tim duoc ket qua giong nhau
            if (Ibest.E_MJ == prev_E) --count_down;
            else
            {
                count_down = reset_count_down;
                prev_E = Ibest.E_MJ;
            }
            if (count_down == 0)
            {
                reset_T();
                count_down = reset_count_down;
            }
            fout << "Ibest      = " << Ibest.E_MJ << "\n";
            fout << "Gbest      = " << Gbest.E_MJ << "\n";
            fout << "count_down = " << count_down << "\n";
            fout << "unfinished_trip = " << unfinished_trip << "\n";
            fout << "time = " << double(clock() - st) / CLOCKS_PER_SEC << " (s)\n";
            if (double(clock() - st) / CLOCKS_PER_SEC > time_limit)
            {
                fout << "\n - Passed time litmit -\n";
                break;
            }
        }
        flog.close();
        ostringstream seed_str, k_str;
        k_str << k;
        seed_str  << seed;
        flog.open((output_file + ".out").c_str());
        fout << Gbest.E_MJ<< "\n";
        fout << Gbest.recaculate() << "\n";
        for (int i = 0 ; i < n; ++i)
            fout << Gbest.X[i]-500 << " " << Gbest.Y[i]-500 << " " << Gbest.Z[i]-500 << "\n";
    }
}

int main(int argc,char *argv[])
{
    parser::run(argc, argv);
    ACO::run();
    return 0;
}
