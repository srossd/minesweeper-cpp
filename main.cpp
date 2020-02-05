#include <stdio.h>
#include <random>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <queue>
#include <bitset>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <muser2_api.hh>

using namespace std;

typedef struct {
    pair<int, int> coord;
    bool is_mine;
    int core;
    bool success;
} inference;

typedef struct {
    bool won;
    int mines_remaining;
    int max_core;
} game_result;

typedef struct {
    double won;
    double won_std;
    double flagged;
    double flagged_std;
    double max_core;
    double max_core_std;
} sim_data;

vector<vector<int>> buildGrid(int, double, bool = false);
void pretty_print(vector<vector<int>>&);
vector<vector<int>> eff_grid(vector<vector<int>>&);
void flag(vector<vector<int>>&, int, int);
vector<pair<int,int>> frontier_inner(vector<vector<int>>&);
vector<pair<int,int>> frontier_outer(vector<vector<int>>&);
inference infer(vector<vector<int>>&);
game_result inference_play(vector<vector<int>>&, bool=false, bool=false);
vector<pair<int,int>> neighbors(int, int, int);
void reveal(vector<vector<int>>&, vector<vector<int>>&, int, int);
void test_muser();
void get_data(int, int, int);
string exec(const char*);
pair<int, bool> call_muser(vector<pair<int, vector<int>>>, int);
sim_data run_sims(int, double, int);
double mean(vector<double>);
double stdev(vector<double>);
double combine_means(double, double, int, int);
double combine_stds(double, double, double, double, int, int);

int main(int argc, char *argv[]) {
    if(atoi(argv[1]) == 0) {
        int N = atoi(argv[2]);
        double trials = atof(argv[3]);
        int batch = atoi(argv[4]);

        get_data(N, trials, batch);
    }
    else {
        int N = atoi(argv[2]);
        double p = atof(argv[3]);
        int trials = atoi(argv[4]);

        sim_data res = run_sims(N, p, trials);
        printf("%0.5f %0.5f %0.5f %0.5f %0.5f %0.5f\n", res.won, res.won_std, res.flagged, res.flagged_std, res.max_core, res.max_core_std);
    }
}

// pair<int, bool> call_muser(vector<pair<int, vector<int>>> clauses, int nvars) {
//     ofstream gcnf("temp.gcnf");
//     gcnf << "p gcnf " << nvars << " " << clauses.size() << " " << clauses.back().first << endl;
//     for(pair<int,vector<int>> clause : clauses) {
//         gcnf << "{" << clause.first << "} ";
//         for(int var : clause.second)
//             gcnf << var << " ";
//         gcnf << "0" << endl;
//     }
//     gcnf.close();

//     string muser_output = exec("./muser2 -grp -test temp.gcnf");
//     size_t result_idx = muser_output.find("result:");
//     bool sat = muser_output[result_idx+8] == 'S';

//     size_t size_idx = muser_output.find("MUS size:");
//     int size;
//     sscanf(muser_output.substr(size_idx + 10).c_str(), "%d", &size);

//     return pair<int,bool>(size, sat);
// }

string exec(const char* cmd) {
    array<char, 128> buffer;
    string result;
    unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

sim_data run_sims(int N, double p, int trials) {
    vector<double> wins;
    vector<double> completion;
    vector<double> cores;
    
    for(int i = 0; i < trials; i++) {
        vector<vector<int>> grid = buildGrid(N, p, true);
        game_result res = inference_play(grid);
        wins.push_back(res.won);
        completion.push_back(1-res.mines_remaining/(N*N*p));
        cores.push_back(res.max_core);
    }

    return {.won = mean(wins), .won_std = stdev(wins), .flagged = mean(completion), .flagged_std = stdev(completion), .max_core = mean(cores), .max_core_std = stdev(cores)};
}

void get_data(int N, int trials, int batch_size) {
    vector<double> ps(30);
    for(int i = 0; i < ps.size(); i++)
        ps[i] = 0.01 + (0.5-0.01)*i/(ps.size()-1);

    
    vector<double> data(ps.size(), 0);
    vector<double> remaining(ps.size(), 0);
    vector<double> avg_max_core(ps.size(), 0);

    vector<double> data_std(ps.size(), 0);
    vector<double> remaining_std(ps.size(), 0);
    vector<double> avg_max_core_std(ps.size(), 0);

    printf("N = %d\n", N);
    for(int j = 0; j < ps.size(); j++) {
        double p = ps[j];
        printf("p = %.2f (%d/%d)\n",p, (j+1), ps.size());

        int batches = (trials+0.5)/batch_size;
        for(int i = 0; i < batches; i++) {
            printf("Batch %d/%d\n", i+1, batches);
            cout.flush();
            string call = "./main 1 "+to_string(N)+" "+to_string(p)+" "+to_string(batch_size);
            string res = exec(call.c_str());
            double x, xs, y, ys, z, zs;
            sscanf(res.c_str(), "%lf %lf %lf %lf %lf %lf", &x, &xs, &y, &ys, &z, &zs);

            int nseen = i*batch_size;
            double new_data = combine_means(data[j], x, nseen, batch_size);
            double new_remaining = combine_means(remaining[j], y, nseen, batch_size);
            double new_core = combine_means(avg_max_core[j], z, nseen, batch_size);

            data_std[j] = combine_stds(data[j], x, data_std[j], xs, nseen, batch_size);
            remaining_std[j] = combine_stds(remaining[j], y, remaining_std[j], ys, nseen, batch_size);
            avg_max_core_std[j] = combine_stds(avg_max_core[j], z, avg_max_core_std[j], zs, nseen, batch_size);

            data[j] = new_data;
            remaining[j] = new_remaining;
            avg_max_core[j] = new_core;
        }

        data_std[j] /= sqrt(trials);
        remaining_std[j] /= sqrt(trials);
        avg_max_core_std[j] /= sqrt(trials);

        // for(int i = 0; i < trials; i++) {
        //     vector<vector<int>> grid = buildGrid(N, p, true);
        //     game_result res = inference_play(grid);
        //     if(res.won)
        //         data[j]++;
        //     remaining[j] += 1-res.mines_remaining/(N*N*p);
        //     avg_max_core[j] += res.max_core;
        // }
        // remaining[j] /= trials;
        // data[j] /= trials;
        // avg_max_core[j] /= trials;
        printf("Solved %0.0f%% of games, discovered %0.0f%% of mines, average max core = %0.2f\n\n", data[j]*100, remaining[j]*100, avg_max_core[j]);
        cout.flush();
    }

    ofstream data_file("data/solvability_data_"+to_string(N)+".txt");
    for (const auto &e : data) data_file << e << "\n";

    ofstream completion_file("data/completion_data_"+to_string(N)+".txt");
    for (const auto &e : remaining) completion_file << e << "\n";

    ofstream core_file("data/core_data_"+to_string(N)+".txt");
    for (const auto &e : avg_max_core) core_file << e << "\n";

    ofstream data_error_file("data/solvability_data_error_"+to_string(N)+".txt");
    for (const auto &e : data_std) data_error_file << e << "\n";

    ofstream completion_error_file("data/completion_data_error_"+to_string(N)+".txt");
    for (const auto &e : remaining_std) completion_error_file << e << "\n";

    ofstream core_error_file("data/core_data_error_"+to_string(N)+".txt");
    for (const auto &e : avg_max_core_std) core_error_file << e << "\n";
}

int count_mines(vector<vector<int>> grid) {
    int ans = 0;
    for(int i = 0; i < grid.size(); i++)
        for(int j = 0; j < grid.size(); j++)
            ans += grid[i][j] == -1;
    return ans;
}

game_result inference_play(vector<vector<int>>& grid, bool display, bool verbose) {
    int num_mines = count_mines(grid);
    vector<vector<int>> game(grid.size());
    for(int i = 0; i < grid.size(); i++)
        game[i] = vector<int>(grid[0].size(), -2);
    
    int mc = 0;
    reveal(grid, game, 0, 0);
    while(count_mines(game) < num_mines) {
        if(display)
            pretty_print(game);
        
        inference inf = infer(game);
        if(inf.success) {
            if(inf.is_mine) {
                flag(game, inf.coord.first, inf.coord.second);
                if(verbose)
                    printf("A mine is found at (%d, %d)\n", inf.coord.first, inf.coord.second);
            }
            else {
                reveal(grid, game, inf.coord.first, inf.coord.second);
                if(verbose)
                    printf("A safe square is found at (%d, %d)\n", inf.coord.first, inf.coord.second);
            }
            if(inf.core > mc)
                mc = inf.core;
        }
        else {
            int remaining = num_mines - count_mines(game);
            if(verbose)
                printf("No more inferences are possible. %d mines remain hidden\n", remaining);
            return {.won = false, .mines_remaining = remaining, .max_core = mc};
        }
    }
    if(verbose)
        printf("The game is completed by inference\n");
    return {.won = true, .mines_remaining = 0, .max_core = mc};
}

inference infer(vector<vector<int>>& game) {
    vector<vector<int>> eg = eff_grid(game);
    vector<pair<int,int>> f_out = frontier_outer(game);
    vector<pair<int,int>> f_in = frontier_inner(game);

    map<pair<int,int>, int> ind_f_out;
    for(int i = 0; i < f_out.size(); i++)
        ind_f_out[f_out[i]] = i+1;

    muser2 m;
    m.init_all();
    // vector<pair<int,vector<int>>> clauses;
    for(int i = 0; i < f_in.size(); i++) {
        vector<int> bool_vars;
        for(pair<int,int> n : neighbors(f_in[i].first, f_in[i].second, game.size())) {
            if(eg[n.first][n.second] == -2) {
                bool_vars.push_back(ind_f_out[n]);
            }
        }

        if(eg[f_in[i].first][f_in[i].second] == 0) {
            for(pair<int,int> n : neighbors(f_in[i].first, f_in[i].second, game.size())) {
                if(eg[n.first][n.second] == -2) {
                    return {.coord = n, .is_mine = false, .core = 1, .success = true};
                }
            }
        }
        if(eg[f_in[i].first][f_in[i].second] == bool_vars.size()) {
            for(pair<int,int> n : neighbors(f_in[i].first, f_in[i].second, game.size())) {
                if(eg[n.first][n.second] == -2) {
                    return {.coord = n, .is_mine = true, .core = 1, .success = true};
                }
            }
        }

        for(int b = 0; b < 1 << bool_vars.size(); b++) {
            bitset<8> bs(b);
            if(bs.count() != eg[f_in[i].first][f_in[i].second]) {
                vector<muser2::lit> clause;
                for(int i = 0; i < bool_vars.size(); i++)
                    clause.push_back((1-2*bs[i])*bool_vars[i]);

                m.add_clause(clause.data(), clause.data()+clause.size(), i);
                // clauses.push_back(pair<int,vector<int>>(i, clause));
            }
        } 
    }

    muser2::gid extra_group = f_in.size();
    muser2::gid added_grp;
    vector<muser2::lit> is_mine_clause;
    pair<int,bool> res;
    int sat;
    for(int i = 0; i < f_out.size(); i++) {
        is_mine_clause = {i+1};
        added_grp = m.add_clause(is_mine_clause.data(), is_mine_clause.data() + 1, extra_group);
        m.init_run();
        sat = m.test_sat();
        // clauses.push_back(pair<int,vector<int>>(extra_group, is_mine_clause));
        // res = call_muser(clauses, f_out.size());
        if(sat == 20) {
            m.compute_gmus();
            return {.coord = f_out[i], .is_mine = false, .core = (int)m.gmus_gids().size(), .success = true};
        }
        // clauses.pop_back();
        if(added_grp == extra_group)
            m.remove_group(extra_group);
        m.reset_run();
        
        is_mine_clause = {-(i+1)};
        added_grp = m.add_clause(is_mine_clause.data(), is_mine_clause.data() + 1, extra_group);
        m.init_run();
        sat = m.test_sat();
        // clauses.push_back(pair<int,vector<int>>(extra_group, is_mine_clause));
        // res = call_muser(clauses, f_out.size());
        if(sat == 20) {
            m.compute_gmus();
            return {.coord = f_out[i], .is_mine = true, .core = (int)m.gmus_gids().size(), .success = true};
        }
        // clauses.pop_back();
        if(added_grp == extra_group)
            m.remove_group(extra_group);

    }
    m.reset_all();
    
    return {.coord = f_out[0], .is_mine = true, .core = 0, .success = false};
}

void reveal(vector<vector<int>>& grid, vector<vector<int>>& game, int i, int j) {
    queue<pair<int,int>> to_click;
    to_click.push(pair<int,int>(i,j));

    set<pair<int,int>> seen;
    while(!to_click.empty()) {
        pair<int,int> coord = to_click.front();
        to_click.pop();
        if(grid[coord.first][coord.second] != -1) {
            game[coord.first][coord.second] = grid[coord.first][coord.second];

            if(game[coord.first][coord.second] == 0) {
                for(pair<int,int> n : neighbors(coord.first, coord.second, grid.size())) {
                    if(!seen.count(n)) {
                        to_click.push(n);
                        seen.insert(n);
                    }
                }
            }
        }
    }
}

void flag(vector<vector<int>>& game, int i, int j) {
    game[i][j] = -1;
}

vector<pair<int,int>> neighbors(int i, int j, int size) {
    vector<pair<int,int>> n;
    for(int dx = -1; dx <= 1; dx++)
        for(int dy = -1; dy <= 1; dy++)
            if(dx != 0 || dy != 0)
                n.push_back(pair<int,int>((i+dx+size) % size, (j+dy+size) % size));
    return n;
}

vector<pair<int,int>> frontier_inner(vector<vector<int>>& game) {
    vector<pair<int,int>> frontier;
    for(int i = 0; i < game.size(); i++) {
        for(int j = 0; j < game[i].size(); j++) {
            if(game[i][j] >= 0) {
                for(pair<int,int> n : neighbors(i, j, game.size())) {
                    if(game[n.first][n.second] == -2) {
                        frontier.push_back(pair<int,int>(i,j));
                        break;
                    }
                }
            }
        }
    }
    return frontier;
}

vector<pair<int,int>> frontier_outer(vector<vector<int>>& game) {
    vector<pair<int,int>> frontier;
    for(int i = 0; i < game.size(); i++) {
        for(int j = 0; j < game[i].size(); j++) {
            if(game[i][j] == -2) {
                for(pair<int,int> n : neighbors(i, j, game.size())) {
                    if(game[n.first][n.second] != -2) {
                        frontier.push_back(pair<int,int>(i,j));
                        break;
                    }
                }
            }
        }
    }
    return frontier;
}

vector<vector<int>> eff_grid(vector<vector<int>>& grid) {
    vector<vector<int>> eg;
    eg = grid;
    for(int i = 0; i < eg.size(); i++) {
        for(int j = 0; j < eg[0].size(); j++) {
            if(grid[i][j] > 0) {
                for(pair<int,int> n : neighbors(i, j, eg.size())) {
                    if(grid[n.first][n.second] == -1)
                        eg[i][j]--;
                }
            }
        }
    }
    return eg;
}

vector<vector<int>> buildGrid(int size, double p, bool safe_origin) {
    vector<vector<int>> grid(size);
    for(int i = 0; i < size; i++)
        grid[i] = vector<int>(size, 0);
    
    int mines = (int)(p*size*size);

    vector<pair<int,int>> coordinates;
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            if(!safe_origin || (i > 1 && i < (size-1)) || (j > 1 && j < (size-1)))
                coordinates.push_back(pair<int,int>(i, j));
        }
    }

    vector<pair<int,int>> mine_coords;
    std::sample(coordinates.begin(), coordinates.end(), back_inserter(mine_coords), mines, mt19937{random_device{}()});

    for(pair<int,int> coord : mine_coords)
        grid[coord.first][coord.second] = -1;

    vector<vector<int>> new_grid(size);
    for(int i = 0; i < size; i++)
        new_grid[i] = vector<int>(size, 0);

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            if(grid[i][j] != -1) {
                for(pair<int,int> n : neighbors(i, j, size)) {
                    new_grid[i][j] -= grid[n.first][n.second];
                }
            }
            else
                new_grid[i][j] = -1;
        }
    }

    return new_grid;
}

void pretty_print(vector<vector<int>>& game) {
    map<int, char> ch;
    ch[-2] = '?'; ch[-1] = 'X'; ch[0] = ' ';
    for(int i = 1; i <= 8; i++)
        ch[i] = (char)(48+i);
    
    for(int i = 0; i < 6*game.size()+1; i++)
        printf("_");
    printf("\n");
    for(int i = 0; i < game.size(); i++) {
        printf("|");
        for(int j = 0; j < game[i].size(); j++)
            printf("     |");
        printf("\n|");
        for(int j = 0; j < game[i].size(); j++)
            printf("  %c  |", ch[game[i][j]]);
        printf("\n|");
        for(int j = 0; j < game[i].size(); j++)
            printf("_____|");
        printf("\n");
    }
}

double mean(vector<double> v) {
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();

    return mean;
}

double stdev(vector<double> v) {
    double m = mean(v);

    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),
                std::bind2nd(std::minus<double>(), m));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / v.size());

    return stdev;
}

double combine_means(double mu1, double mu2, int n1, int n2) {
    return (n1*mu1 + n2*mu2)/(n1+n2);
}

double combine_stds(double mu1, double mu2, double sigma1, double sigma2, int n1, int n2) {
    double mu = combine_means(mu1, mu2, n1, n2);
    return sqrt((n1*sigma1*sigma1 + n2*sigma2*sigma2 + n1*(mu1-mu)*(mu1-mu) + n2*(mu2-mu)*(mu2-mu))/(n1+n2));
}